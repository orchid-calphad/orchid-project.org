#!/usr/bin/env python3
"""Build the ORCHID bibliography: every paper ORCHID knows about, its roles, evidence and typed relations.

Writes `<orchid-data>/bibliography/` (papers.jsonl/.csv/.parquet, edges.jsonl/.csv/.parquet, systems_index.json,
systems/<SYS>.json, papers/<key>.json, index.json, papers.jsonld, context.jsonld, orchid-bibliography.bib,
orchid-bibliography.csl.json, README.md, SCHEMA.md, CITATION.cff, parquet.json) from local ORCHID sources and a
Crossref metadata cache.

Two steps:

  fetch   (network) fill the Crossref cache: read-only GET https://api.crossref.org/works/<doi> for every DOI not yet
          cached, at most 5 requests/s, User-Agent "ORCHID-bibliography/0.1 (https://orchid-project.org)", no mailto.
          Every response (200 and 404) is stored under build-cache/crossref/ (git-ignored) with its retrieval time.
  build   (offline, default) deterministic: sorted output, no clock reads; the only dates are the cache retrieval
          dates. `--check` rebuilds in memory and exits 1 if any file differs.

Run with an interpreter that has pyarrow (the Parquet files), e.g.
  ~/orchid-py313/bin/python tools/build_bibliography.py fetch
  ~/orchid-py313/bin/python tools/build_bibliography.py build [--check]

Never written: abstracts, full text, local PDF paths. Only Crossref bibliographic fields, ORCHID-derived roles,
counts and edges.
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import threading
import time
import unicodedata
import urllib.error
import urllib.parse
import urllib.request

HOME = os.path.expanduser("~")
SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DEF_DATA = os.environ.get("ORCHID_DATA") or os.path.join(os.path.dirname(SITE), "orchid-data")
DEF_OPEN_TDB = os.path.join(HOME, "projects", "open-tdb")
DEF_CACHE = os.path.join(SITE, "build-cache", "crossref")
DISCOVERY_REF = "discovery/forward-citations"
DISCOVERY_PATH = "catalog/discovery/higher-order-assessment-targets.json"

DATASET = "bibliography"
VERSION = "2026.10"
RELEASE_TAG = "v2026.10"
BIB_VERSION = "0.1"
USER_AGENT = "ORCHID-bibliography/0.1 (https://orchid-project.org)"
CROSSREF = "https://api.crossref.org/works/"
DATA_URL = "https://data.orchid-project.org/"
HALLSTEDT_DOI = "10.1016/j.calphad.2025.102833"
HANDBOOK = "H. Okamoto, Desk Handbook: Phase Diagrams for Binary Alloys, 2nd ed., ASM International (2010)"

ROLES = ["corpus_assessment_read", "cited_original_source", "tdb_source_assessment", "tdb_compilation",
         "handbook_diagram_source", "discovery_lead", "context_source"]
EDGE_TYPES = ["cites", "used_in_assessment", "models_system_in_tdb", "reassesses", "handbook_reproduces",
              "cited_by_crossref"]
EVIDENCE = ["experimental", "model", "dft", "unspecified"]


# ----------------------------------------------------------------------------- helpers
def norm_doi(d):
    if not d:
        return None
    d = str(d).strip()
    d = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:\s*)", "", d, flags=re.I).strip().lower()
    return d if re.match(r"^10\.\d{4,9}/\S+$", d) else None


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def key_of(pid):
    """URL-safe stable file key of a bibliography id: 'p' + first 16 hex of sha256(id)."""
    return "p" + hashlib.sha256(pid.encode("utf-8")).hexdigest()[:16]


def fold(s):
    s = unicodedata.normalize("NFKD", str(s or "").replace("ß", "ss").replace("ø", "o").replace("Ø", "O").replace("ł", "l"))
    s = "".join(c for c in s if not unicodedata.combining(c)).lower()
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def system_id(elements):
    els = sorted({str(e).strip().upper() for e in elements if str(e).strip()})
    return "-".join(els) if els else None


def norm_system(s):
    if not s:
        return None
    if isinstance(s, (list, tuple)):
        return system_id(s)
    return system_id(re.split(r"[-–\s]+", str(s)))


def cache_path(cache, doi):
    h = hashlib.sha256(doi.encode("utf-8")).hexdigest()
    return os.path.join(cache, h[:2], h + ".json")


def jdump(obj):
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def first_author_surname(citation):
    """Folded surname of the first author of a printed reference ('J. E. Kittl, T.B. ...' -> 'kittl',
    'A. Fernandez Guillermet, ...' -> 'fernandez guillermet'). The name must be followed by ',', ';', 'and', 'et al.'
    or '&' (a journal abbreviation such as 'J. Phase Equilib.' is followed by '.'); None otherwise."""
    if not citation:
        return None
    c = citation.strip()
    m = re.match(r"^((?:[A-Z][a-z]?\.\s*-?\s*)+)\s*([^\W\d][^,;()0-9.]{1,40}?)\s*(?:,|;|\band\b|\bet al|&)", c)
    if m:
        return fold(m.group(2)) or None
    m = re.match(r"^([^\W\d][^,;()0-9.]{1,40}?),\s*(?:[A-Z][a-z]?\.\s*-?\s*)+\s*(?:,|;|\band\b|\bet al|&|:)", c)
    if m:
        return fold(m.group(1)) or None
    return None


# ----------------------------------------------------------------------------- sources
class Sources:
    """All local inputs, read once; `hashes` records the sha256 of each input for the release record."""

    def __init__(self, data, open_tdb):
        self.data = data
        self.open_tdb = open_tdb
        self.corpus = os.path.join(open_tdb, "build", "sgte-corpus")
        self.hashes = {}
        lit = os.path.join(data, "literature")
        self.papers = self._json(os.path.join(lit, "papers.json"), "orchid-data/literature/papers.json")
        self.records = []
        p = os.path.join(lit, "orchid-literature.jsonl.gz")
        self.hashes["orchid-data/literature/orchid-literature.jsonl.gz"] = sha256_file(p)
        with gzip.open(p, "rt", encoding="utf-8") as f:
            for line in f:
                self.records.append(json.loads(line))
        self.tdb_links = self._json(os.path.join(self.corpus, "tdb_pdf_links.json"),
                                    "open-tdb/build/sgte-corpus/tdb_pdf_links.json")["links"]
        self.census = self._json(os.path.join(self.corpus, "census.json"), "open-tdb/build/sgte-corpus/census.json")
        self.handbook = self._json(os.path.join(self.corpus, "db-v2.1", "okamoto", "handbook_labels.json"),
                                   "open-tdb/build/sgte-corpus/db-v2.1/okamoto/handbook_labels.json")
        cm = os.path.join(self.corpus, "db-v2.1", "cited-meta")
        self.fulltext_kind = self._json(os.path.join(cm, "fulltext-kind.json"),
                                        "open-tdb/build/sgte-corpus/db-v2.1/cited-meta/fulltext-kind.json")
        self.verified = self._json(os.path.join(open_tdb, "catalog", "verified-assessments.json"),
                                   "open-tdb/catalog/verified-assessments.json")["assessments"]
        self.catalog = self._json(os.path.join(open_tdb, "catalog", "catalog.json"), "open-tdb/catalog/catalog.json")
        # discovery leads: tracked file on the unmerged discovery/forward-citations branch, read through git
        commit = subprocess.run(["git", "-C", open_tdb, "rev-parse", DISCOVERY_REF], capture_output=True, text=True,
                                check=True).stdout.strip()
        raw = subprocess.run(["git", "-C", open_tdb, "show", f"{commit}:{DISCOVERY_PATH}"], capture_output=True,
                             check=True).stdout
        self.hashes[f"open-tdb@{commit[:12]}:{DISCOVERY_PATH}"] = sha256_bytes(raw)
        self.discovery_commit = commit
        self.discovery = json.loads(raw)
        # ESPEI bundle, phase registry, assessed assemblages (optional datasets)
        self.espei = []
        p = os.path.join(data, "espei", "index.csv")
        if os.path.isfile(p):
            self.hashes["orchid-data/espei/index.csv"] = sha256_file(p)
            with open(p, newline="", encoding="utf-8") as f:
                self.espei = list(csv.DictReader(f))
        self.phases = None
        self.phase_systems = {}
        p = os.path.join(data, "phases", "registry.json")
        if os.path.isfile(p):
            self.phases = self._json(p, "orchid-data/phases/registry.json")
            sd = os.path.join(data, "phases", "systems")
            lines = []
            for fn in sorted(os.listdir(sd)) if os.path.isdir(sd) else []:
                if fn.endswith(".json"):
                    fp = os.path.join(sd, fn)
                    lines.append(f"{sha256_file(fp)}  {fn}")
                    with open(fp, encoding="utf-8") as f:
                        self.phase_systems[fn[:-5]] = json.load(f)
            self.hashes["orchid-data/phases/systems/*.json (sha256 of 'sha256  name' lines)"] = sha256_bytes(
                "\n".join(lines).encode())
        self.assemblages = None
        p = os.path.join(data, "assessed-phase-assemblages", "records.csv")
        if os.path.isfile(p):
            self.hashes["orchid-data/assessed-phase-assemblages/records.csv"] = sha256_file(p)
            with open(p, newline="", encoding="utf-8") as f:
                self.assemblages = [{"id": r["id"], "system": r["system"], "original_doi": r["original_doi"]}
                                    for r in csv.DictReader(f)]
        # bundled TDB files (sha256 of the bytes ORCHID distributes)
        self.tdb_dir = os.path.join(SITE, "data")
        if not os.path.isdir(os.path.join(self.tdb_dir, "databases")):
            self.tdb_dir = open_tdb
        self._tdb_sha = {}
        self.catalog_by_name = {c["name"]: c for c in self.catalog}

    def _json(self, path, label):
        with open(path, "rb") as f:
            b = f.read()
        self.hashes[label] = sha256_bytes(b)
        return json.loads(b)

    def tdb_sha(self, rel):
        """sha256 of a bundled TDB (rel = 'databases/sgte-binary/AgAl-04Wit.tdb'); None when not bundled."""
        if rel not in self._tdb_sha:
            p = os.path.join(self.tdb_dir, rel)
            self._tdb_sha[rel] = sha256_file(p) if os.path.isfile(p) else None
        return self._tdb_sha[rel]

    def all_dois(self):
        out = set()
        for pid, p in self.papers.items():
            d = norm_doi(p.get("doi"))
            if d:
                out.add(d)
        for x in self.tdb_links:
            d = norm_doi(x.get("doi"))
            if d:
                out.add(d)
        for sysn, s in self.handbook["systems"].items():
            for src in s.get("sources") or []:
                d = norm_doi(src.get("doi"))
                if d:
                    out.add(d)
        for r in self.discovery["records"]:
            d = norm_doi(r.get("doi"))
            if d:
                out.add(d)
        for a in self.verified:
            d = norm_doi(a.get("doi"))
            if d:
                out.add(d)
        for r in self.espei:
            for d in (r.get("dois") or "").split(";"):
                d = norm_doi(d)
                if d:
                    out.add(d)
        for r in self.records:
            m = r.get("model") or {}
            d = norm_doi(m.get("model_doi"))
            if d:
                out.add(d)
        for r in self.assemblages or []:
            d = norm_doi(r["original_doi"])
            if d:
                out.add(d)
        out.add(HALLSTEDT_DOI)
        return out


# ----------------------------------------------------------------------------- Crossref fetch
class RateLimiter:
    def __init__(self, per_second):
        self.interval = 1.0 / per_second
        self.lock = threading.Lock()
        self.next = 0.0

    def wait(self):
        with self.lock:
            now = time.monotonic()
            t = max(now, self.next)
            self.next = t + self.interval
        if t > now:
            time.sleep(t - now)


def fetch_one(doi, cache, limiter, timeout=30):
    url = CROSSREF + urllib.parse.quote(doi, safe="/")
    limiter.wait()
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
    status, body = None, None
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310 (fixed https host)
            status = resp.status
            body = json.loads(resp.read())
    except urllib.error.HTTPError as exc:
        status = exc.code
        if exc.code == 429:
            time.sleep(10)
        if exc.code not in (404, 400):
            return doi, f"http_{exc.code}"
    except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError) as exc:
        return doi, f"error:{type(exc).__name__}"
    if body is not None and isinstance(body.get("message"), dict):
        body["message"].pop("abstract", None)  # never kept, never published
    rec = {"doi": doi, "url": url, "http_status": status,
           "retrieved_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "body": body}
    p = cache_path(cache, doi)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    tmp = p + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(rec, f, ensure_ascii=False, sort_keys=True)
    os.replace(tmp, p)
    return doi, status


def cmd_fetch(src, cache, rate, workers, limit=None):
    from concurrent.futures import ThreadPoolExecutor
    todo = sorted(d for d in src.all_dois() if not os.path.isfile(cache_path(cache, d)))
    if limit:
        todo = todo[:limit]
    print(f"{len(todo)} DOIs to fetch (rate {rate}/s, {workers} workers)", flush=True)
    limiter = RateLimiter(min(rate, 5.0))
    res = collections.Counter()
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for i, (doi, st) in enumerate(ex.map(lambda d: fetch_one(d, cache, limiter), todo), 1):
            res[str(st)] += 1
            if i % 250 == 0:
                print(f"  {i}/{len(todo)} {dict(res)}", flush=True)
    print(dict(res))
    return 0


# ----------------------------------------------------------------------------- Crossref parsing
CSL_TYPE = {"journal-article": "article-journal", "book-chapter": "chapter", "proceedings-article": "paper-conference",
            "book": "book", "monograph": "book", "edited-book": "book", "reference-book": "book", "report": "report",
            "dissertation": "thesis", "posted-content": "article", "dataset": "dataset", "reference-entry": "entry",
            "book-section": "chapter", "book-part": "chapter", "standard": "standard", "component": "article"}
BIB_TYPE = {"journal-article": "article", "book-chapter": "incollection", "proceedings-article": "inproceedings",
            "book": "book", "monograph": "book", "edited-book": "book", "reference-book": "book", "report": "techreport",
            "dissertation": "phdthesis", "book-section": "incollection", "book-part": "incollection"}


def _date(parts):
    try:
        dp = parts["date-parts"][0]
        return "-".join(f"{int(x):02d}" if i else f"{int(x):04d}" for i, x in enumerate(dp) if x is not None)
    except (KeyError, IndexError, TypeError, ValueError):
        return None


def parse_crossref(rec):
    """Bibliographic fields of a cached Crossref response (no abstract)."""
    if not rec or rec.get("http_status") != 200 or not rec.get("body"):
        return None
    m = rec["body"].get("message") or {}

    def first(k):
        v = m.get(k)
        return v[0] if isinstance(v, list) and v else (v if isinstance(v, str) else None)

    authors = []
    for a in m.get("author") or []:
        orcid = a.get("ORCID")
        if orcid:
            orcid = re.sub(r"^https?://orcid\.org/", "", orcid.strip())
        au = {"given": a.get("given"), "family": a.get("family"), "name": a.get("name"), "orcid": orcid,
              "orcid_authenticated": a.get("authenticated-orcid") if orcid else None}
        authors.append({k: v for k, v in au.items() if v is not None})
    year = None
    for k in ("issued", "published-print", "published-online", "published", "created"):
        d = _date(m.get(k) or {})
        if d:
            year = int(d[:4])
            break
    lic = []
    for x in m.get("license") or []:
        lic.append({k: v for k, v in {"url": x.get("URL"), "content_version": x.get("content-version"),
                                      "start": _date(x.get("start") or {}),
                                      "delay_in_days": x.get("delay-in-days")}.items() if v is not None})
    lic.sort(key=lambda x: (x.get("url") or "", x.get("content_version") or "", x.get("start") or ""))
    refs = []
    for r in m.get("reference") or []:
        d = norm_doi(r.get("DOI"))
        if d:
            refs.append({"doi": d, "asserted_by": r.get("doi-asserted-by")})
    return {
        "doi": norm_doi(m.get("DOI")),
        "title": first("title"),
        "subtitle": first("subtitle"),
        "authors": authors,
        "container_title": first("container-title"),
        "container_title_short": first("short-container-title"),
        "volume": m.get("volume"),
        "issue": m.get("issue"),
        "pages": m.get("page"),
        "article_number": m.get("article-number"),
        "year": year,
        "issued": _date(m.get("issued") or {}),
        "publisher": m.get("publisher"),
        "type": m.get("type"),
        "issn": sorted(set(m.get("ISSN") or [])),
        "isbn": sorted(set(m.get("ISBN") or [])),
        "language": m.get("language"),
        "license": lic,
        "is_referenced_by_count": m.get("is-referenced-by-count"),
        "references_count": m.get("references-count", m.get("reference-count")),
        "references_deposited": len(m.get("reference") or []),
        "references_with_doi": refs,
        "retrieved_at": rec.get("retrieved_at"),
    }


def oa_license(lic):
    """First Creative Commons licence URL in the Crossref `license` field (never guessed from elsewhere)."""
    cc = [x for x in lic if re.search(r"creativecommons\.org/(licenses|publicdomain)/", x.get("url") or "", re.I)]
    if not cc:
        return None
    order = {"vor": 0, "am": 1, "unspecified": 2, "tdm": 3}
    cc.sort(key=lambda x: (order.get(x.get("content_version"), 9), x.get("url")))
    x = cc[0]
    return {"url": x["url"], "content_version": x.get("content_version"), "start": x.get("start"),
            "delay_in_days": x.get("delay_in_days"), "basis": "crossref_license_field"}


def doi_check(local, cr):
    """Compare a DOI's Crossref record with what ORCHID's source printed (first-author surname, year, title)."""
    if not cr:
        return None
    issues = []
    sur = local.get("first_author")
    names = [set(fold(" ".join(x for x in (a.get("given"), a.get("family"), a.get("name")) if x)).split())
             for a in cr["authors"]]
    if sur and names:
        import difflib
        st = [t for t in sur.split() if len(t) > 1]

        def hit(tokens):
            return any(t in tokens or any(difflib.SequenceMatcher(None, t, u).ratio() >= 0.75 for u in tokens) for t in st)
        if st and not any(hit(n) for n in names):
            issues.append("first_author_not_in_crossref_authors")
    y = local.get("year")
    if y and cr.get("year") and abs(int(y) - int(cr["year"])) > 1:
        issues.append("year_differs")
    lt = fold(local.get("title"))
    ct = fold(cr.get("title"))
    if lt and ct and len(lt.split()) >= 4:
        a, b = set(lt.split()), set(ct.split())
        if len(a & b) / max(1, min(len(a), len(b))) < 0.5:
            issues.append("title_differs")
    if not (sur or y or lt):
        return {"status": "unchecked", "issues": []}
    return {"status": "consistent" if not issues else "mismatch_suspect", "issues": issues}


# ----------------------------------------------------------------------------- assembly
YEAR_RE = re.compile(r"(?<!\d)(1[89]\d\d|20[0-3]\d)(?!\d)")


def year_of(text):
    ys = YEAR_RE.findall(text or "")
    return int(ys[-1]) if ys else None


def nodoi_id(first_author, year, title_or_citation):
    basis = f"{fold(first_author)}|{year or ''}|{fold(title_or_citation)}"
    return "nodoi:" + hashlib.sha256(basis.encode("utf-8")).hexdigest()[:16]


class Bib:
    def __init__(self, src, cache):
        self.src = src
        self.cache = cache
        self.P = {}
        self.edges = {}
        self.lit_map = {}
        self.cr = {}
        self.stats = collections.Counter()
        self.tdb_rows = {}

    # -- papers
    def paper(self, pid):
        p = self.P.get(pid)
        if p is None:
            p = self.P[pid] = {"id": pid, "roles": {}, "systems": set(), "local": {}, "aliases": set(),
                               "provenance": set(), "inferred_systems": {}}
        return p

    def doi_paper(self, doi):
        return self.paper("doi:" + doi)

    def local(self, p, source, **kw):
        """Bibliographic fields printed / held locally, first source wins per field."""
        for k, v in kw.items():
            if v not in (None, "", []):
                p["local"].setdefault(k, (v, source))

    def role(self, p, role, **detail):
        d = p["roles"].setdefault(role, {})
        for k, v in detail.items():
            if isinstance(v, set):
                d.setdefault(k, set()).update(v)
            elif isinstance(v, list):
                d.setdefault(k, []).extend(v)
            elif isinstance(v, int) and not isinstance(v, bool) and k.startswith("n_"):
                d[k] = d.get(k, 0) + v
            else:
                d.setdefault(k, v)
        return d

    def edge(self, etype, src, dst, provenance, **attrs):
        k = (etype, src, dst, provenance["source"])
        e = self.edges.get(k)
        if e is None:
            e = self.edges[k] = {"type": etype, "src": src, "dst": dst, "provenance": dict(provenance), "attrs": {}}
        for a, v in attrs.items():
            if isinstance(v, int) and not isinstance(v, bool):
                e["attrs"][a] = e["attrs"].get(a, 0) + v
            elif isinstance(v, collections.Counter):
                e["attrs"].setdefault(a, collections.Counter()).update(v)
            elif isinstance(v, set):
                e["attrs"].setdefault(a, set()).update(v)
            else:
                e["attrs"].setdefault(a, v)
        return e

    # -- steps
    def from_literature_papers(self):
        for lid, rec in sorted(self.src.papers.items()):
            if rec.get("kind") == "unresolved_citation" or lid.startswith("unresolved:"):
                self.stats["literature_unresolved_citation_labels_skipped"] += 1
                continue
            doi = norm_doi(rec.get("doi"))
            authors = rec.get("authors") or []
            fa = None
            if authors:
                fa = fold(str(authors[0]).split()[-1]) if str(authors[0]).split() else None
            elif rec.get("citation"):
                fa = first_author_surname(rec["citation"])
            if doi:
                pid = "doi:" + doi
            else:
                pid = nodoi_id(fa, rec.get("year"), rec.get("title") or rec.get("citation"))
            self.lit_map[lid] = pid
            p = self.paper(pid)
            p["aliases"].add(lid)
            p["provenance"].add("orchid-data/literature/papers.json")
            src_kind = "orchid_literature_" + (rec.get("metadata_source") or "unknown")
            self.local(p, src_kind, title=rec.get("title"), year=rec.get("year"), container=rec.get("container"),
                       citation=rec.get("citation"), authors=authors)
            if rec.get("metadata_source") == "printed_reference" or rec.get("metadata_source") == "tdb_header":
                self.local(p, "printed", check_first_author=fa, check_year=rec.get("year"),
                           check_title=rec.get("title"))
            for s in rec.get("systems") or []:
                p["systems"].add(norm_system(s))
            roles = set(rec.get("roles") or [])
            if rec.get("kind") == "corpus_paper" or "reporting_paper" in roles:
                self.role(p, "corpus_assessment_read", basis="orchid-data/literature (reporting paper)")
            if "cited_source" in roles:
                self.role(p, "cited_original_source", basis="orchid-data/literature (cited source)")

    def from_tdb_links(self):
        census = set(self.src.census)
        for x in sorted(self.src.tdb_links, key=lambda x: (x.get("tdb") or "", x.get("doi") or "")):
            st = x.get("status")
            tdb = x.get("tdb")
            if st == "unparsed_name" or not tdb:
                self.stats["tdb_links_unparsed_name"] += 1
                continue
            sysn = norm_system(x.get("system"))
            doi = norm_doi(x.get("doi"))
            hc = x.get("header_citation")
            if doi:
                p = self.doi_paper(doi)
            elif hc:
                p = self.paper(nodoi_id(first_author_surname(hc), year_of(hc), hc))
                self.stats[f"tdb_links_without_doi_{st}"] += 1
            else:
                continue
            p["provenance"].add("open-tdb/build/sgte-corpus/tdb_pdf_links.json")
            self.local(p, "tdb_header", citation=hc, year=year_of(x.get("year") or hc), container=x.get("journal"))
            self.local(p, "printed", check_first_author=first_author_surname(hc), check_year=year_of(x.get("year") or hc))
            if sysn:
                p["systems"].add(sysn)
            rel = x.get("tdb_path") or f"databases/sgte-binary/{tdb}"
            sha = self.src.tdb_sha(rel)
            row = {"file": tdb, "path": rel, "sha256": sha, "bundled": sha is not None, "system": sysn,
                   "link_status": st, "link_match": x.get("match"), "basis": "sgte_tdb_header_citation"}
            self.role(p, "tdb_source_assessment", tdb_files=[row])
            self.tdb_rows[tdb] = row
            self.edge("models_system_in_tdb", p["id"], "tdb:" + tdb,
                      {"source": "sgte_tdb_header", "dataset": "open-tdb/build/sgte-corpus/tdb_pdf_links.json",
                       "method": f"TDB header citation matched to the Hallstedt 2025 reference list ({st}"
                                 + (f", {x.get('match')}" if x.get("match") else "") + ")"},
                      system=sysn, tdb_sha256=sha, link_status=st)
            pdf = x.get("pdf")
            if pdf and pdf in census:
                self.role(p, "corpus_assessment_read", basis="sgte-corpus census (PDF processed)", in_census=True)
        h = self.doi_paper(HALLSTEDT_DOI)
        n_sgte = sum(1 for c in self.src.catalog if (c.get("local_path") or "").startswith("databases/sgte-binary/"))
        self.role(h, "tdb_compilation", collection="SGTE binary collection", n_tdb_files=n_sgte)
        h["provenance"].add("open-tdb/catalog/catalog.json")

    def from_model_links(self):
        seen = set()
        for r in self.src.records:
            m = r.get("model") or {}
            doi = norm_doi(m.get("model_doi"))
            tdb = m.get("tdb")
            if not doi or not tdb or (doi, tdb) in seen:
                continue
            seen.add((doi, tdb))
            p = self.doi_paper(doi)
            p["provenance"].add("orchid-data/literature (model links)")
            sysn = r.get("system")
            sha = m.get("tdb_sha256")
            row = {"file": tdb, "path": m.get("tdb_path"), "sha256": sha, "bundled": bool(sha), "system": sysn,
                   "link_status": m.get("link_basis"), "link_match": None, "basis": "literature_model_link"}
            if not any(t["file"] == tdb for t in p["roles"].get("tdb_source_assessment", {}).get("tdb_files", [])):
                self.role(p, "tdb_source_assessment", tdb_files=[row])
            self.edge("models_system_in_tdb", p["id"], "tdb:" + tdb,
                      {"source": "orchid_literature_model_link", "dataset": "orchid-data/literature",
                       "method": "record.model.model_doi -> record.model.tdb (" + str(m.get("link_basis")) + ")"},
                      system=sysn, tdb_sha256=sha)

    def from_verified(self):
        for a in self.src.verified:
            doi = norm_doi(a.get("doi"))
            if not doi:
                continue
            p = self.doi_paper(doi)
            p["provenance"].add("open-tdb/catalog/verified-assessments.json")
            sysn = norm_system(a.get("elements"))
            p["systems"].add(sysn)
            self.local(p, "verified_assessments", title=a.get("title"))
            for name in a.get("related_catalog_names") or []:
                c = self.src.catalog_by_name.get(name) or {}
                rel = c.get("local_path") or ""
                sha = self.src.tdb_sha(rel) if rel and c.get("tier") == "bundled" else None
                row = {"file": name, "path": rel or None, "sha256": sha, "bundled": sha is not None, "system": sysn,
                       "link_status": "verified_assessment", "link_match": None, "basis": "verified_assessment",
                       "catalog_tier": c.get("tier"), "license_status": c.get("license_status")}
                if not any(t["file"] == name for t in p["roles"].get("tdb_source_assessment", {}).get("tdb_files", [])):
                    self.role(p, "tdb_source_assessment", tdb_files=[row])
                self.edge("models_system_in_tdb", p["id"], "tdb:" + name,
                          {"source": "verified_assessment", "dataset": "open-tdb/catalog/verified-assessments.json",
                           "method": "hand-verified assessment -> catalog entry"}, system=sysn, tdb_sha256=sha)

    def from_handbook(self):
        for sysn, s in sorted(self.src.handbook["systems"].items()):
            sysn = norm_system(sysn)
            node = f"handbook:okamoto-2010:{sysn}"
            for src in s.get("sources") or []:
                cit = src.get("citation")
                doi = norm_doi(src.get("doi"))
                if doi:
                    p = self.doi_paper(doi)
                elif cit:
                    p = self.paper(nodoi_id(first_author_surname(cit), year_of(cit), cit))
                else:
                    continue
                p["provenance"].add("open-tdb/build/sgte-corpus/db-v2.1/okamoto/handbook_labels.json")
                self.local(p, "handbook_reference", citation=cit, year=year_of(cit))
                if doi:
                    self.local(p, "printed", check_first_author=first_author_surname(cit), check_year=year_of(cit))
                p["systems"].add(sysn)
                self.role(p, "handbook_diagram_source", systems={sysn}, kind=src.get("kind"))
                self.edge("handbook_reproduces", node, p["id"],
                          {"source": "okamoto_handbook_source_line", "dataset": "handbook_labels.json (sources)",
                           "method": "source line printed under the handbook diagram"
                                     + (f"; DOI matched by {src.get('matched_by')}" if doi else "; no DOI")},
                          system=sysn, page=src.get("page"), doi_matched_by=src.get("matched_by"))

    def from_discovery(self):
        for r in sorted(self.src.discovery["records"], key=lambda r: r.get("rank") or 0):
            doi = norm_doi(r.get("doi"))
            if not doi:
                continue
            p = self.doi_paper(doi)
            p["provenance"].add(f"open-tdb@{self.src.discovery_commit[:12]}:{DISCOVERY_PATH}")
            self.local(p, "discovery_snapshot", title=r.get("title"), year=r.get("year"), container=r.get("venue"))
            named = sorted(s for s, b in (r.get("candidate_systems") or {}).items() if b == "named_in_title")
            inferred = {s: b for s, b in sorted((r.get("candidate_systems") or {}).items()) if b != "named_in_title"}
            for s in named:
                p["systems"].add(norm_system(s))
            p["inferred_systems"].update(inferred)
            self.role(p, "discovery_lead", rank=r.get("rank"), score=r.get("score"), max_order=r.get("max_order"),
                      systems_named_in_title=named, candidate_systems_inferred=inferred,
                      discovery_sources=sorted(r.get("sources") or []),
                      crossref_confirmed_seed_count=r.get("crossref_confirmed_seed_count"),
                      license_review=r.get("license_review"), is_seed=bool(r.get("is_seed")),
                      snapshot=f"open-tdb@{self.src.discovery_commit[:12]}:{DISCOVERY_PATH}")

    def from_crossref(self):
        for pid, p in self.P.items():
            if not pid.startswith("doi:"):
                continue
            doi = pid[4:]
            cp = cache_path(self.cache, doi)
            if not os.path.isfile(cp):
                p["crossref_status"] = "not_fetched"
                continue
            with open(cp, encoding="utf-8") as f:
                rec = json.load(f)
            cr = parse_crossref(rec)
            p["crossref_retrieved_at"] = rec.get("retrieved_at")
            if cr is None:
                p["crossref_status"] = "not_found" if rec.get("http_status") == 404 else f"http_{rec.get('http_status')}"
                continue
            p["crossref_status"] = "found"
            self.cr[pid] = cr

    def link_ids(self):
        """Literature ids of every paper id that ORCHID records reference (for edges and evidence)."""
        return self.lit_map

    def from_records(self):
        ev_rep = collections.defaultdict(lambda: {"n": 0, "acc": 0, "ev": collections.Counter(),
                                                  "pg": collections.Counter(), "pc": collections.Counter(),
                                                  "tables": set(), "passages": set(), "tier": collections.Counter(),
                                                  "agree": collections.Counter()})
        ev_cit = collections.defaultdict(lambda: {"n": 0, "acc": 0, "ev": collections.Counter(),
                                                  "pg": collections.Counter(), "pc": collections.Counter(),
                                                  "citing": set()})
        for r in self.src.records:
            pid = self.lit_map.get(r.get("paper"))
            acc = bool((r.get("flags") or {}).get("extraction_accepted"))
            if pid:
                e = ev_rep[pid]
                e["n"] += 1
                e["acc"] += acc
                e["ev"][r.get("evidence_class") or "unspecified"] += 1
                e["pg"][r.get("property_group") or "other property"] += 1
                e["pc"][r.get("property") or "?"] += 1
                loc = r.get("locator") or {}
                where = loc.get("asset") or f"{loc.get('page')}:{loc.get('label')}"
                if loc.get("artifact") == "table":
                    e["tables"].add(where)
                elif loc.get("artifact") == "text_passage":
                    e["passages"].add(r.get("quote") or where)
                e["tier"][r.get("reader_tier") or "unknown"] += 1
                if r.get("reader_agreement"):
                    e["agree"][r["reader_agreement"]] += 1
            seen = set()
            for cs in r.get("cited_sources") or []:
                dst = self.lit_map.get(cs.get("paper"))
                if not dst:
                    if cs.get("paper"):
                        self.stats["cited_source_values_unresolved_label"] += 1
                    continue
                if dst == pid or dst in seen:
                    continue
                seen.add(dst)
                c = ev_cit[dst]
                c["n"] += 1
                c["acc"] += acc
                c["ev"][r.get("evidence_class") or "unspecified"] += 1
                c["pg"][r.get("property_group") or "other property"] += 1
                c["pc"][r.get("property") or "?"] += 1
                if pid:
                    c["citing"].add(pid)
                    prov = {"source": "orchid_extraction", "dataset": "orchid-data/literature/orchid-literature.jsonl.gz",
                            "method": "record.cited_sources (reference printed for the value in the containing paper)"}
                    attrs = dict(n_values=1, n_values_accepted=int(acc),
                                 evidence_class=collections.Counter([r.get("evidence_class") or "unspecified"]),
                                 citation_status=collections.Counter([cs.get("status") or "unknown"]),
                                 property=collections.Counter([r.get("property") or "?"]), systems={r.get("system")})
                    self.edge("cites", pid, dst, prov, **attrs)
                    self.edge("used_in_assessment", dst, pid, prov, **attrs)
        for pid, e in ev_rep.items():
            p = self.P[pid]
            p["evidence_reported"] = e
        for pid, c in ev_cit.items():
            self.P[pid]["evidence_cited"] = c
            self.role(self.P[pid], "cited_original_source", basis="orchid-data/literature (cited source)")

    def from_espei(self):
        for r in self.src.espei:
            ids = set()
            for d in (r.get("dois") or "").split(";"):
                d = norm_doi(d)
                if d:
                    ids.add("doi:" + d)
            for lid in (r.get("papers") or "").split(";"):
                if lid in self.lit_map:
                    ids.add(self.lit_map[lid])
            for pid in ids:
                if pid not in self.P:
                    self.stats["espei_dois_not_in_bibliography"] += 1
                    p = self.paper(pid)
                    p["provenance"].add("orchid-data/espei/index.csv")
                    self.role(p, "espei_reference", basis="cited by an ESPEI dataset of the ORCHID bundle")
                self.P[pid].setdefault("espei", set()).add(r["path"])

    def from_phases(self):
        if not self.src.phases:
            return
        for row in self.src.phases.get("phases") or []:
            for n in row.get("names") or []:
                if n.startswith("paper:"):
                    rest = n[len("paper:"):]
                    m = re.match(r"^(10\.\d{4,9}/.+):([^:]*)$", rest)
                    if m:
                        pid = "doi:" + norm_doi(m.group(1))
                        if pid in self.P:
                            self.P[pid].setdefault("phases", set()).add(row["id"])
        for sysn, doc in self.src.phase_systems.items():
            for ph in doc.get("phases") or []:
                for ref in ph.get("references") or []:
                    d = norm_doi(ref.get("doi"))
                    if d and "doi:" + d in self.P:
                        self.P["doi:" + d].setdefault("phases", set()).add(ph["id"])

    def from_assemblages(self):
        if self.src.assemblages is None:
            return
        for r in self.src.assemblages:
            d = norm_doi(r["original_doi"])
            if d and "doi:" + d in self.P:
                self.P["doi:" + d]["assemblages"] = self.P["doi:" + d].get("assemblages", 0) + 1

    def from_crossref_refs(self):
        for pid, cr in sorted(self.cr.items()):
            for ref in cr["references_with_doi"]:
                dst = "doi:" + ref["doi"]
                if dst == pid:
                    continue
                if dst not in self.P:
                    self.stats["crossref_reference_dois_outside_bibliography"] += 1
                    continue
                self.edge("cited_by_crossref", pid, dst,
                          {"source": "crossref_reference_list", "dataset": "api.crossref.org/works (reference)",
                           "method": "reference list deposited openly by the publisher of the citing work",
                           "retrieved_at": (self.P[pid].get("crossref_retrieved_at") or "")[:10]},
                          doi_asserted_by=ref.get("asserted_by"))

    def reassessments(self):
        """reassesses: A cites B (ORCHID extraction or Crossref reference list) and A, B both model the same exact
        system in a TDB (tdb_source_assessment). Never from dates alone."""
        def tdb_systems(p):
            return {t["system"] for t in p["roles"].get("tdb_source_assessment", {}).get("tdb_files", []) if t["system"]}
        cites = collections.defaultdict(set)
        for (etype, s, d, source), e in self.edges.items():
            if etype in ("cites", "cited_by_crossref"):
                cites[(s, d)].add(source)
        for (s, d), how in sorted(cites.items()):
            ps, pd = self.P.get(s), self.P.get(d)
            if not ps or not pd:
                continue
            common = sorted(tdb_systems(ps) & tdb_systems(pd))
            if not common:
                continue
            ys, yd = self.year(ps), self.year(pd)
            if ys and yd and ys < yd:
                self.stats["reassesses_skipped_citing_older_than_cited"] += 1
                continue
            self.edge("reassesses", s, d,
                      {"source": "citation+same-system", "dataset": "edges cites / cited_by_crossref + TDB systems",
                       "method": "citation+same-system", "citation_basis": sorted(how)},
                      systems=set(common), year_src=ys, year_dst=yd)

    def year(self, p):
        cr = self.cr.get(p["id"])
        if cr and cr.get("year"):
            return cr["year"]
        v = p["local"].get("year")
        return int(v[0]) if v and str(v[0]).isdigit() else None

    def run(self):
        self.from_literature_papers()
        self.from_tdb_links()
        self.from_model_links()
        self.from_verified()
        self.from_handbook()
        self.from_discovery()
        self.from_records()
        self.from_espei()
        self.from_crossref()
        self.from_phases()
        self.from_assemblages()
        self.from_crossref_refs()
        self.reassessments()
        return self


# ----------------------------------------------------------------------------- records
def _counter(c):
    return dict(sorted(c.items(), key=lambda kv: (-kv[1], kv[0])))


def _jsonable(v):
    if isinstance(v, set):
        return sorted(v, key=lambda x: (x is None, str(x)))
    if isinstance(v, collections.Counter):
        return _counter(v)
    if isinstance(v, dict):
        return {k: _jsonable(x) for k, x in sorted(v.items())}
    if isinstance(v, list):
        return [_jsonable(x) for x in v]
    return v


def sys_sort(s):
    return (s.count("-") + 1, s)


def paper_record(b, p):
    pid = p["id"]
    doi = pid[4:] if pid.startswith("doi:") else None
    cr = b.cr.get(pid)
    loc = {k: v[0] for k, v in p["local"].items()}
    loc_src = {k: v[1] for k, v in p["local"].items()}
    rec = {"id": pid, "key": key_of(pid), "doi": doi, "doi_url": f"https://doi.org/{doi}" if doi else None,
           "id_basis": "doi" if doi else "normalized_first_author_year_title",
           "no_doi": doi is None}
    if cr:
        rec.update({k: cr[k] for k in ("title", "subtitle", "authors", "container_title", "container_title_short",
                                       "volume", "issue", "pages", "article_number", "year", "issued", "publisher",
                                       "type", "issn", "isbn", "language")})
        rec["metadata_source"] = "crossref"
    else:
        authors = loc.get("authors") or []
        rec.update({"title": loc.get("title"), "subtitle": None,
                    "authors": [{"name": a} for a in authors] if authors and isinstance(authors[0], str) else authors,
                    "container_title": loc.get("container"), "container_title_short": None, "volume": None,
                    "issue": None, "pages": None, "article_number": None,
                    "year": int(loc["year"]) if str(loc.get("year") or "").isdigit() else None, "issued": None,
                    "publisher": None, "type": None, "issn": [], "isbn": [], "language": None})
        rec["metadata_source"] = loc_src.get("title") or loc_src.get("citation") or loc_src.get("year") or "none"
    cit, csrc = loc.get("citation"), loc_src.get("citation")
    rec["citation_printed"] = cit if cit and (csrc != "orchid_literature_crossref" or not cr) else None
    rec["crossref"] = {"status": p.get("crossref_status", "no_doi" if not doi else "not_fetched"),
                       "retrieved_at": p.get("crossref_retrieved_at"),
                       "references_deposited": cr["references_deposited"] if cr else None,
                       "references_count": cr["references_count"] if cr else None}
    rec["license"] = cr["license"] if cr else []
    rec["is_referenced_by_count"] = cr["is_referenced_by_count"] if cr else None
    rec["is_referenced_by_count_retrieved"] = (p.get("crossref_retrieved_at") or "")[:10] or None if cr else None
    chk = None
    if cr and any(k in loc for k in ("check_first_author", "check_year", "check_title")):
        chk = doi_check({"first_author": loc.get("check_first_author"), "year": loc.get("check_year"),
                         "title": loc.get("check_title")}, cr)
    rec["doi_check"] = chk
    rec["systems"] = sorted({s for s in p["systems"] if s}, key=sys_sort)
    rec["candidate_systems_inferred"] = dict(sorted(p["inferred_systems"].items()))
    rec["roles"] = [r for r in ROLES + ["espei_reference"] if r in p["roles"]]
    rd = {}
    for r in rec["roles"]:
        d = _jsonable(p["roles"][r])
        if r == "tdb_source_assessment":
            d["tdb_files"] = sorted(d["tdb_files"], key=lambda t: t["file"])
        if r == "corpus_assessment_read":
            d["values_released"] = p.get("evidence_reported", {}).get("n", 0)
        rd[r] = d
    rec["role_details"] = rd
    # evidence
    es = {}
    e = p.get("evidence_reported")
    if e:
        es["reported"] = {"n_values": e["n"], "n_values_accepted": e["acc"],
                          "by_evidence_class": _counter(e["ev"]), "by_property_group": _counter(e["pg"]),
                          "by_property": _counter(e["pc"]), "n_tables": len(e["tables"]),
                          "n_text_passages": len(e["passages"]), "n_figures": 0}
        es["reader_agreement"] = {"by_reader_tier": _counter(e["tier"]),
                                  "by_reader_agreement": _counter(e["agree"]) if e["agree"] else None}
    c = p.get("evidence_cited")
    if c:
        es["cited_as_source"] = {"n_values": c["n"], "n_values_accepted": c["acc"],
                                 "by_evidence_class": _counter(c["ev"]), "by_property_group": _counter(c["pg"]),
                                 "by_property": _counter(c["pc"]), "n_citing_papers": len(c["citing"])}
    if p.get("espei"):
        es["espei_datasets"] = sorted(p["espei"])
    if p.get("phases"):
        es["phase_registry_phases"] = sorted(p["phases"])
    if b.src.assemblages is not None and p.get("assemblages"):
        es["assessed_assemblage_records"] = p["assemblages"]
    rec["evidence_summary"] = es
    rd_detail = []
    if "corpus_assessment_read" in p["roles"]:
        rd_detail.append("values_extracted" if (e and e["n"]) else "processed_no_released_values")
    if doi and doi in b.src.fulltext_kind:
        rd_detail.append("full_text_screened_for_origin")
    rec["availability"] = {"open_access_license": oa_license(rec["license"]),
                           "read_by_orchid": "yes" if "corpus_assessment_read" in p["roles"] else "no",
                           "read_detail": rd_detail,
                           "pdf_redistributed": False}
    rec["aliases"] = sorted(p["aliases"])
    rec["provenance"] = sorted(p["provenance"])
    rec["urls"] = {"json": f"{DATA_URL}bibliography/papers/{rec['key']}.json",
                   "html": f"{DATA_URL}bibliography/pages/papers/{rec['key']}.html"}
    return rec


def edge_record(e):
    attrs = _jsonable(e["attrs"])
    for k in ("systems",):
        if k in attrs:
            attrs[k] = sorted((x for x in attrs[k] if x), key=sys_sort)
    eid = "e" + hashlib.sha256(f"{e['type']}|{e['src']}|{e['dst']}|{e['provenance']['source']}".encode()).hexdigest()[:16]
    return {"id": eid, "type": e["type"], "src": e["src"], "dst": e["dst"], "src_kind": node_kind(e["src"]),
            "dst_kind": node_kind(e["dst"]), "count": attrs.get("n_values"), "provenance": e["provenance"],
            "attrs": attrs}


def node_kind(nid):
    return {"doi": "paper", "nodoi": "paper", "tdb": "tdb_file", "handbook": "handbook_diagram"}[nid.split(":", 1)[0]]


# ----------------------------------------------------------------------------- formatting
def author_text(a):
    if a.get("family"):
        return (a.get("given") + " " if a.get("given") else "") + a["family"]
    return a.get("name") or ""


def citation_text(r):
    """One-line citation from the record's metadata (Crossref when present)."""
    if r["metadata_source"] != "crossref" and r.get("citation_printed"):
        return r["citation_printed"] + (f", https://doi.org/{r['doi']}" if r.get("doi") else "")
    au = [author_text(a) for a in r.get("authors") or []]
    au = ", ".join(au[:6]) + (" et al." if len(au) > 6 else "")
    parts = [au] if au else []
    if r.get("title"):
        parts.append(r["title"])
    src = r.get("container_title") or r.get("publisher") or ""
    if r.get("volume"):
        src += f" {r['volume']}"
    if r.get("issue"):
        src += f"({r['issue']})"
    if r.get("year"):
        src += f" ({r['year']})"
    if r.get("pages") or r.get("article_number"):
        src += f" {r.get('pages') or r.get('article_number')}"
    if src.strip():
        parts.append(src.strip())
    if r.get("doi"):
        parts.append(f"https://doi.org/{r['doi']}")
    return ", ".join(parts)


CITE_TITLE = "ORCHID bibliography"
CITE_KEY = "orchid_bibliography_2026"
CITE_ABSTRACT = ("Every paper ORCHID knows about: Crossref bibliographic metadata (CC0), the role of each paper in "
                 "ORCHID, counts of the evidence released from it, and typed relations with provenance (citations of "
                 "values, TDB files, reassessments, handbook diagrams, Crossref reference lists).")
CITE_PAPERS = "the papers themselves when you use their data (doi / citation of each bibliography record)"
CITE_STATUS = ("Version 1, machine-extracted: roles and evidence counts derive from machine-extracted literature values "
               "and link tables; bibliographic fields are Crossref metadata as deposited (DOI mismatches flagged in "
               "doi_check).")


def cite_module():
    """orchid_cite with this dataset registered at run time (same pattern as build_espei_bundle.py), so the
    `cite_as` blocks match the other datasets; see INTEGRATION.md for the permanent registration."""
    import orchid_cite as oc
    if DATASET not in oc.DATASETS:
        oc.THIRD_PARTY.setdefault("papers_bibliography", CITE_PAPERS)
        oc.DATASETS[DATASET] = (CITE_TITLE, CITE_KEY, ["papers_bibliography"])
        oc.EXTRACTED.add(DATASET)
        oc.STATUS_BY.setdefault(DATASET, CITE_STATUS)
        oc.ABSTRACT[DATASET] = CITE_ABSTRACT
        oc.NAMES = tuple(oc.DATASETS)
    return oc


def cite_as(data_dir):
    c = cite_module().cite_as(DATASET, data_dir)
    c["licence_note"] = "Crossref bibliographic metadata: CC0; ORCHID-derived fields (roles, evidence counts, relations): CC BY 4.0"
    return c


def bibtex_escape(s):
    s = str(s)
    for a, b in (("\\", "\\textbackslash{}"), ("{", "\\{"), ("}", "\\}"), ("&", "\\&"), ("%", "\\%"), ("$", "\\$"),
                 ("#", "\\#"), ("_", "\\_"), ("~", "\\textasciitilde{}"), ("^", "\\textasciicircum{}")):
        s = s.replace(a, b)
    return s


def bib_key(r):
    a = (r.get("authors") or [{}])[0]
    fam = fold(a.get("family") or a.get("name") or "anon").split()
    fam = fam[-1] if fam else "anon"
    return f"{fam}{r.get('year') or 'nd'}_{r['key'][1:7]}"


def bibtex_entry(r):
    t = BIB_TYPE.get(r.get("type"), "misc")
    f = []
    au = " and ".join(bibtex_escape(f"{a['family']}, {a['given']}" if a.get("family") and a.get("given")
                                    else a.get("family") or a.get("name") or "") for a in r.get("authors") or [])
    if au:
        f.append(("author", au))
    if r.get("title"):
        f.append(("title", "{" + bibtex_escape(r["title"]) + "}"))
    cont = r.get("container_title")
    if cont:
        f.append(({"article": "journal", "incollection": "booktitle", "inproceedings": "booktitle"}.get(t, "howpublished"),
                  bibtex_escape(cont)))
    for k, bk in (("volume", "volume"), ("issue", "number"), ("pages", "pages"), ("year", "year"),
                  ("publisher", "publisher")):
        if r.get(k):
            v = str(r[k])
            if k == "pages":
                v = v.replace("-", "--")
            f.append((bk, bibtex_escape(v)))
    if r.get("article_number") and not r.get("pages"):
        f.append(("eid", bibtex_escape(r["article_number"])))
    f.append(("doi", r["doi"]))
    lic = (r["availability"].get("open_access_license") or {}).get("url")
    if lic:
        f.append(("license", lic))
    body = ",\n".join(f"  {k:<9} = {{{v}}}" for k, v in f)
    return f"@{t}{{{bib_key(r)},\n{body}\n}}\n"


def csl_item(r):
    it = {"id": r["id"], "type": CSL_TYPE.get(r.get("type"), "article"), "DOI": r["doi"],
          "URL": f"https://doi.org/{r['doi']}"}
    au = []
    for a in r.get("authors") or []:
        x = {}
        if a.get("family"):
            x["family"] = a["family"]
            if a.get("given"):
                x["given"] = a["given"]
        elif a.get("name"):
            x["literal"] = a["name"]
        if x:
            au.append(x)
    if au:
        it["author"] = au
    for k, ck in (("title", "title"), ("container_title", "container-title"),
                  ("container_title_short", "container-title-short"), ("volume", "volume"), ("issue", "issue"),
                  ("pages", "page"), ("publisher", "publisher"), ("article_number", "number"), ("language", "language")):
        if r.get(k):
            it[ck] = r[k]
    if r.get("issued"):
        it["issued"] = {"date-parts": [[int(x) for x in r["issued"].split("-")]]}
    elif r.get("year"):
        it["issued"] = {"date-parts": [[r["year"]]]}
    if r.get("issn"):
        it["ISSN"] = r["issn"][0]
    if r.get("isbn"):
        it["ISBN"] = r["isbn"][0]
    return it


JSONLD_CONTEXT = {
    "@context": {
        "@vocab": "https://schema.org/",
        "schema": "https://schema.org/",
        "orchid": "https://orchid-project.org/vocab/bibliography#",
        "roles": {"@id": "orchid:role", "@container": "@set"},
        "edgeType": "orchid:edgeType",
        "citation": {"@id": "schema:citation", "@type": "@id"},
        "sameAs": {"@id": "schema:sameAs", "@type": "@id"},
        "url": {"@id": "schema:url", "@type": "@id"},
        "license": {"@id": "schema:license", "@type": "@id"},
    }
}


def jsonld_paper(r, cites_out):
    node = {"@type": "ScholarlyArticle" if r.get("type") in (None, "journal-article", "proceedings-article") else
            ("Book" if r.get("type") in ("book", "monograph", "edited-book", "reference-book") else "CreativeWork"),
            "@id": r["urls"]["json"], "identifier": [r["id"]], "name": r.get("title") or r.get("citation_printed"),
            "url": r["urls"]["html"]}
    if r.get("doi"):
        node["identifier"].append({"@type": "PropertyValue", "propertyID": "doi", "value": r["doi"]})
        node["sameAs"] = f"https://doi.org/{r['doi']}"
    au = []
    for a in r.get("authors") or []:
        p = {"@type": "Person", "name": author_text(a)}
        if a.get("family"):
            p["familyName"] = a["family"]
        if a.get("given"):
            p["givenName"] = a["given"]
        if a.get("orcid"):
            p["@id"] = f"https://orcid.org/{a['orcid']}"
        au.append(p)
    if au:
        node["author"] = au
    if r.get("year"):
        node["datePublished"] = r.get("issued") or str(r["year"])
    if r.get("container_title"):
        per = {"@type": "Periodical", "name": r["container_title"]}
        if r.get("issn"):
            per["issn"] = r["issn"]
        part = per
        if r.get("volume"):
            part = {"@type": "PublicationVolume", "volumeNumber": r["volume"], "isPartOf": per}
            if r.get("issue"):
                part = {"@type": "PublicationIssue", "issueNumber": r["issue"], "isPartOf": part}
        node["isPartOf"] = part
    if r.get("pages"):
        node["pagination"] = r["pages"]
    if r.get("publisher"):
        node["publisher"] = {"@type": "Organization", "name": r["publisher"]}
    oa = r["availability"].get("open_access_license")
    if oa:
        node["license"] = oa["url"]
    if cites_out:
        node["citation"] = cites_out
    props = [{"@type": "PropertyValue", "propertyID": "orchid:role", "value": x} for x in r["roles"]]
    props += [{"@type": "PropertyValue", "propertyID": "orchid:system", "value": s} for s in r["systems"]]
    es = r["evidence_summary"]
    for k in ("reported", "cited_as_source"):
        if k in es:
            props.append({"@type": "PropertyValue", "propertyID": f"orchid:values_{k}", "value": es[k]["n_values"]})
            for ev, n in sorted(es[k]["by_evidence_class"].items(), key=lambda kv: (-kv[1], kv[0])):
                props.append({"@type": "PropertyValue", "propertyID": f"orchid:values_{k}:{ev}", "value": n})
    if r.get("is_referenced_by_count") is not None:
        props.append({"@type": "PropertyValue", "propertyID": "crossref:is-referenced-by-count",
                      "value": r["is_referenced_by_count"], "observationDate": r["is_referenced_by_count_retrieved"]})
    props.append({"@type": "PropertyValue", "propertyID": "orchid:read_by_orchid",
                  "value": r["availability"]["read_by_orchid"]})
    node["additionalProperty"] = props
    return node


# ----------------------------------------------------------------------------- tables
CSV_COLS = ["id", "key", "doi", "no_doi", "title", "authors", "orcids", "container_title", "volume", "issue", "pages",
            "year", "publisher", "type", "metadata_source", "citation", "systems", "roles", "tdb_files",
            "handbook_systems", "discovery_rank", "values_reported", "values_reported_accepted",
            "values_reported_experimental", "values_reported_model", "values_reported_dft", "values_reported_unspecified",
            "values_cited", "values_cited_accepted", "property_groups", "properties", "n_tables", "n_text_passages",
            "n_espei_datasets", "n_phase_registry_phases", "n_assessed_assemblage_records", "open_access_license",
            "license_urls", "is_referenced_by_count", "is_referenced_by_count_retrieved", "crossref_status",
            "doi_check", "read_by_orchid", "edges_out", "edges_in", "json_url", "html_url"]


def flat_row(r, ec):
    es = r["evidence_summary"]
    rep = es.get("reported") or {}
    cit = es.get("cited_as_source") or {}
    ev = rep.get("by_evidence_class") or {}
    pg = sorted(set(rep.get("by_property_group") or {}) | set(cit.get("by_property_group") or {}))
    pc = sorted(set(rep.get("by_property") or {}) | set(cit.get("by_property") or {}))
    rd = r["role_details"]
    return {
        "id": r["id"], "key": r["key"], "doi": r["doi"] or "", "no_doi": r["no_doi"], "title": r.get("title") or "",
        "authors": "; ".join(author_text(a) for a in r.get("authors") or []),
        "orcids": "; ".join(a["orcid"] for a in r.get("authors") or [] if a.get("orcid")),
        "container_title": r.get("container_title") or "", "volume": r.get("volume") or "",
        "issue": r.get("issue") or "", "pages": r.get("pages") or r.get("article_number") or "",
        "year": r.get("year"), "publisher": r.get("publisher") or "", "type": r.get("type") or "",
        "metadata_source": r.get("metadata_source") or "", "citation": citation_text(r),
        "systems": ";".join(r["systems"]), "roles": ";".join(r["roles"]),
        "tdb_files": ";".join(t["file"] for t in (rd.get("tdb_source_assessment") or {}).get("tdb_files", [])),
        "handbook_systems": ";".join(sorted((rd.get("handbook_diagram_source") or {}).get("systems", []), key=sys_sort)),
        "discovery_rank": (rd.get("discovery_lead") or {}).get("rank"),
        "values_reported": rep.get("n_values", 0), "values_reported_accepted": rep.get("n_values_accepted", 0),
        "values_reported_experimental": ev.get("experimental", 0), "values_reported_model": ev.get("model", 0),
        "values_reported_dft": ev.get("dft", 0), "values_reported_unspecified": ev.get("unspecified", 0),
        "values_cited": cit.get("n_values", 0), "values_cited_accepted": cit.get("n_values_accepted", 0),
        "property_groups": ";".join(pg), "properties": ";".join(pc),
        "n_tables": rep.get("n_tables", 0), "n_text_passages": rep.get("n_text_passages", 0),
        "n_espei_datasets": len(es.get("espei_datasets") or []),
        "n_phase_registry_phases": len(es.get("phase_registry_phases") or []),
        "n_assessed_assemblage_records": es.get("assessed_assemblage_records", 0),
        "open_access_license": (r["availability"]["open_access_license"] or {}).get("url", ""),
        "license_urls": ";".join(sorted({x.get("url") for x in r["license"] if x.get("url")})),
        "is_referenced_by_count": r.get("is_referenced_by_count"),
        "is_referenced_by_count_retrieved": r.get("is_referenced_by_count_retrieved") or "",
        "crossref_status": r["crossref"]["status"],
        "doi_check": (r.get("doi_check") or {}).get("status", ""),
        "read_by_orchid": r["availability"]["read_by_orchid"],
        "edges_out": ec.get((r["id"], "out"), 0), "edges_in": ec.get((r["id"], "in"), 0),
        "json_url": r["urls"]["json"], "html_url": r["urls"]["html"],
    }


EDGE_COLS = ["id", "type", "src", "dst", "src_kind", "dst_kind", "count", "provenance_source", "provenance_dataset",
             "provenance_method", "systems", "attrs"]


def flat_edge(e):
    a = e["attrs"]
    return {"id": e["id"], "type": e["type"], "src": e["src"], "dst": e["dst"], "src_kind": e["src_kind"],
            "dst_kind": e["dst_kind"], "count": e["count"], "provenance_source": e["provenance"]["source"],
            "provenance_dataset": e["provenance"].get("dataset", ""), "provenance_method": e["provenance"].get("method", ""),
            "systems": ";".join(a.get("systems") or ([a["system"]] if a.get("system") else [])),
            "attrs": jdump({k: v for k, v in a.items() if k not in ("systems", "system")})}


def to_csv(rows, cols):
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({k: ("" if r.get(k) is None else (str(r[k]).lower() if isinstance(r[k], bool) else r[k]))
                    for k in cols})
    return buf.getvalue()


INT_COLS = {"year", "discovery_rank", "values_reported", "values_reported_accepted", "values_reported_experimental",
            "values_reported_model", "values_reported_dft", "values_reported_unspecified", "values_cited",
            "values_cited_accepted", "n_tables", "n_text_passages", "n_espei_datasets", "n_phase_registry_phases",
            "n_assessed_assemblage_records", "is_referenced_by_count", "edges_out", "edges_in", "count"}


def to_parquet(rows, cols):
    import pyarrow as pa
    import pyarrow.parquet as pq
    arrays = []
    fields = []
    for c in cols:
        vals = [r.get(c) for r in rows]
        if c in INT_COLS:
            t = pa.int64()
            vals = [None if v in (None, "") else int(v) for v in vals]
        elif c == "no_doi":
            t = pa.bool_()
        else:
            t = pa.string()
            vals = [None if v is None else str(v) for v in vals]
        arrays.append(pa.array(vals, type=t))
        fields.append(pa.field(c, t))
    table = pa.Table.from_arrays(arrays, schema=pa.schema(fields))
    buf = io.BytesIO()
    pq.write_table(table, buf, compression="zstd", write_statistics=False)
    return buf.getvalue()


# ----------------------------------------------------------------------------- build
def build_outputs(src, cache, with_pages=True):
    b = Bib(src, cache).run()
    lines = "\n".join(f"{v}  {k}" for k, v in sorted(src._tdb_sha.items()) if v)
    src.hashes["TDB files referenced (sha256 of 'sha256  path' lines; orchid-project.org/data)"] = sha256_bytes(lines.encode())
    records = [paper_record(b, p) for p in b.P.values()]
    records.sort(key=lambda r: r["id"])
    by_id = {r["id"]: r for r in records}
    edges = sorted((edge_record(e) for e in b.edges.values()), key=lambda e: (e["type"], e["src"], e["dst"], e["id"]))
    ec = collections.Counter()
    for e in edges:
        ec[(e["src"], "out")] += 1
        ec[(e["dst"], "in")] += 1
    for r in records:
        rc = collections.defaultdict(lambda: {"out": 0, "in": 0})
        r["relation_counts"] = {}
    for e in edges:
        if e["src"] in by_id:
            by_id[e["src"]]["relation_counts"].setdefault(e["type"], {"out": 0, "in": 0})["out"] += 1
        if e["dst"] in by_id:
            by_id[e["dst"]]["relation_counts"].setdefault(e["type"], {"out": 0, "in": 0})["in"] += 1
    for r in records:
        r["relation_counts"] = dict(sorted(r["relation_counts"].items()))
        r["citation"] = citation_text(r)
    out = {}
    cite = cite_as(src.data)
    # counts
    counts = {"papers": len(records), "papers_with_doi": sum(1 for r in records if r["doi"]),
              "papers_without_doi": sum(1 for r in records if not r["doi"]),
              "papers_by_role": {ro: sum(1 for r in records if ro in r["roles"]) for ro in ROLES + ["espei_reference"]},
              "papers_by_role_with_doi": {ro: sum(1 for r in records if ro in r["roles"] and r["doi"])
                                          for ro in ROLES + ["espei_reference"]},
              "edges": len(edges), "edges_by_type": {t: sum(1 for e in edges if e["type"] == t) for t in EDGE_TYPES},
              "crossref": dict(sorted(collections.Counter(r["crossref"]["status"] for r in records if r["doi"]).items())),
              "doi_check": dict(sorted(collections.Counter((r.get("doi_check") or {}).get("status", "not_checked")
                                                           for r in records if r["doi"]).items())),
              "open_access_license": sum(1 for r in records if r["availability"]["open_access_license"]),
              "read_by_orchid": sum(1 for r in records if r["availability"]["read_by_orchid"] == "yes"),
              "authors_with_orcid": sum(1 for r in records for a in r.get("authors") or [] if a.get("orcid")),
              "notes": dict(sorted(b.stats.items()))}
    retrieved = sorted({r["crossref"]["retrieved_at"][:10] for r in records if r["crossref"].get("retrieved_at")})
    counts["crossref_retrieved"] = [retrieved[0], retrieved[-1]] if retrieved else None
    # systems index
    sysidx = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in records:
        for s in r["systems"]:
            for ro in r["roles"] or ["(none)"]:
                sysidx[s][ro].append(r["id"])
    counts["systems"] = len(sysidx)
    counts["systems_by_order"] = dict(sorted(collections.Counter(str(s.count("-") + 1) for s in sysidx).items(),
                                             key=lambda kv: int(kv[0])))
    systems_index = {"dataset": DATASET, "version": VERSION, "note": "system -> role -> paper ids (sorted). Systems "
                     "are sorted element symbols joined by '-'. Per-system detail: systems/<SYS>.json.",
                     "systems": {s: {ro: sorted(v) for ro, v in sorted(d.items())}
                                 for s, d in sorted(sysidx.items(), key=lambda kv: sys_sort(kv[0]))}}
    # papers.jsonl / csv / parquet
    out["papers.jsonl"] = "".join(jdump(r) + "\n" for r in records)
    flat = [flat_row(r, ec) for r in records]
    out["papers.csv"] = to_csv(flat, CSV_COLS)
    out["papers.parquet"] = to_parquet(flat, CSV_COLS)
    out["edges.jsonl"] = "".join(jdump(e) + "\n" for e in edges)
    fe = [flat_edge(e) for e in edges]
    out["edges.csv"] = to_csv(fe, EDGE_COLS)
    out["edges.parquet"] = to_parquet(fe, EDGE_COLS)
    out["parquet.json"] = json.dumps({
        "papers": {"table": "papers.parquet", "rows": len(flat), "id_column": "id",
                   "csv_sha256": sha256_bytes(out["papers.csv"].encode()),
                   "content_sha256": sha256_bytes("\n".join(jdump(x) for x in flat).encode())},
        "edges": {"table": "edges.parquet", "rows": len(fe), "id_column": "id",
                  "csv_sha256": sha256_bytes(out["edges.csv"].encode()),
                  "content_sha256": sha256_bytes("\n".join(jdump(x) for x in fe).encode())},
        "note": "content_sha256 = SHA-256 of the canonical JSON rows (sorted keys); independent of the Parquet encoder",
        "fields": {"papers": [{"name": n, "dataType": t, "description": d} for n, t, d in PAPER_FIELDS],
                   "edges": [{"name": n, "dataType": t, "description": d} for n, t, d in EDGE_FIELDS]}},
        indent=1, ensure_ascii=False) + "\n"
    out["systems_index.json"] = json.dumps(systems_index, ensure_ascii=False, separators=(",", ":")) + "\n"
    # per-paper relations
    rel = collections.defaultdict(list)
    for e in edges:
        rel[e["src"]].append(("out", e))
        rel[e["dst"]].append(("in", e))

    def brief(nid):
        r = by_id.get(nid)
        if r is None:
            return {"id": nid, "kind": node_kind(nid)}
        return {"id": nid, "key": r["key"], "doi": r["doi"], "title": r.get("title") or r.get("citation_printed"),
                "year": r.get("year")}

    for r in records:
        rels = []
        for direction, e in rel.get(r["id"], []):
            other = e["dst"] if direction == "out" else e["src"]
            x = {"direction": direction, "type": e["type"], "edge_id": e["id"], "other": brief(other),
                 "count": e["count"], "provenance": {"source": e["provenance"]["source"]}}
            if e["type"] != "cited_by_crossref":  # full provenance of every edge: edges.jsonl (by edge_id)
                x["provenance"] = e["provenance"]
                x["attrs"] = e["attrs"]
            rels.append(x)
        rels.sort(key=lambda x: (x["type"], x["direction"], x["other"]["id"]))
        doc = {"paper": r, "relations": rels, "cite_as": cite,
               "licence": "Crossref bibliographic metadata: CC0. ORCHID-derived fields (roles, evidence, relations): "
                          "CC BY 4.0."}
        out[f"papers/{r['key']}.json"] = json.dumps(doc, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n"
    # per-system files (assessments with the data papers they used)
    for s, d in systems_index["systems"].items():
        ids = sorted({i for v in d.values() for i in v})
        assessments = []
        for i in ids:
            r = by_id[i]
            if "tdb_source_assessment" not in r["roles"] and "corpus_assessment_read" not in r["roles"]:
                continue
            used = []
            for direction, e in rel.get(i, []):
                if direction == "out" and e["type"] == "cites" and s in (e["attrs"].get("systems") or []):
                    used.append({**brief(e["dst"]), "n_values": e["attrs"].get("n_values"),
                                 "n_values_accepted": e["attrs"].get("n_values_accepted"),
                                 "evidence_class": e["attrs"].get("evidence_class")})
            used.sort(key=lambda x: (-(x["n_values"] or 0), x["id"]))
            reass = sorted({e["dst"] for dr, e in rel.get(i, []) if dr == "out" and e["type"] == "reassesses"
                            and s in (e["attrs"].get("systems") or [])})
            tdbs = [t["file"] for t in (r["role_details"].get("tdb_source_assessment") or {}).get("tdb_files", [])
                    if t.get("system") == s]
            rep = (r["evidence_summary"].get("reported") or {})
            assessments.append({**brief(i), "tdb_files": tdbs, "values_reported": rep.get("n_values", 0),
                                "data_papers_used": used, "n_data_papers_used": len(used),
                                "n_values_from_cited_papers": sum(u["n_values"] or 0 for u in used),
                                "reassesses": reass})
        assessments.sort(key=lambda a: (a.get("year") or 0, a["id"]))
        doc = {"system": s, "elements": s.split("-"), "papers_by_role": d,
               "papers": [brief(i) for i in ids], "assessments": assessments, "cite_as": cite}
        out[f"systems/{s}.json"] = json.dumps(doc, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n"
    # search index (compact rows)
    cols = ["key", "id", "doi", "title", "year", "container", "roles", "systems", "properties", "evidence_classes",
            "values_reported", "values_cited", "authors"]
    rows = []
    for r in records:
        es = r["evidence_summary"]
        props = sorted(set((es.get("reported") or {}).get("by_property") or {}) |
                       set((es.get("cited_as_source") or {}).get("by_property") or {}))
        evs = sorted(set((es.get("reported") or {}).get("by_evidence_class") or {}) |
                     set((es.get("cited_as_source") or {}).get("by_evidence_class") or {}))
        rows.append([r["key"], r["id"], r["doi"], r.get("title") or r.get("citation_printed"), r.get("year"),
                     r.get("container_title"), r["roles"], r["systems"], props, evs,
                     (es.get("reported") or {}).get("n_values", 0), (es.get("cited_as_source") or {}).get("n_values", 0),
                     [author_text(a) for a in (r.get("authors") or [])[:3]]])
    out["index.json"] = json.dumps({"dataset": DATASET, "version": VERSION, "release": RELEASE_TAG,
                                    "bibliography_version": BIB_VERSION, "counts": counts,
                                    "inputs": dict(sorted(src.hashes.items())),
                                    "id_rule": "id = 'doi:<lowercase DOI>' or 'nodoi:<first 16 hex of sha256(folded "
                                               "first-author surname|year|folded title or printed citation)>'; key = "
                                               "'p' + first 16 hex of sha256(id); one record per id at papers/<key>.json",
                                    "urls": {"paper": DATA_URL + "bibliography/papers/{key}.json",
                                             "paper_page": DATA_URL + "bibliography/pages/papers/{key}.html",
                                             "system": DATA_URL + "bibliography/systems/{system}.json",
                                             "system_page": DATA_URL + "bibliography/pages/systems/{system}.html"},
                                    "columns": cols, "rows": rows, "cite_as": cite},
                                   ensure_ascii=False, separators=(",", ":")) + "\n"
    # BibTeX / CSL-JSON (Crossref-backed records only)
    crs = [r for r in records if r["metadata_source"] == "crossref"]
    out["orchid-bibliography.bib"] = (
        "% ORCHID bibliography " + VERSION + ": entries generated from Crossref metadata (CC0); "
        f"{len(crs)} records with a DOI and a Crossref record.\n\n" + "\n".join(bibtex_entry(r) for r in crs))
    out["orchid-bibliography.csl.json"] = json.dumps([csl_item(r) for r in crs], ensure_ascii=False, indent=1) + "\n"
    # JSON-LD
    cites_out = collections.defaultdict(list)
    for e in edges:
        if e["type"] in ("cites", "cited_by_crossref") and e["dst"] in by_id:
            u = by_id[e["dst"]]["urls"]["json"]
            if u not in cites_out[e["src"]]:
                cites_out[e["src"]].append(u)
    out["context.jsonld"] = json.dumps(JSONLD_CONTEXT, indent=1) + "\n"
    graph = [jsonld_paper(r, sorted(cites_out.get(r["id"], []))) for r in records]
    out["papers.jsonld"] = json.dumps({"@context": JSONLD_CONTEXT["@context"], "@type": "Dataset",
                                       "@id": DATA_URL + "bibliography/papers.jsonld", "name": "ORCHID bibliography",
                                       "license": "https://creativecommons.org/licenses/by/4.0/",
                                       "hasPart": graph}, ensure_ascii=False, separators=(",", ":")) + "\n"
    out.update(docs(counts, cite))
    if with_pages:
        import bibliography_pages
        bibliography_pages._SYSDOCS.clear()
        bibliography_pages._SYSDOCS.update({s: out[f"systems/{s}.json"] for s in systems_index["systems"]})
        out.update(bibliography_pages.render(records, rel, by_id, systems_index["systems"], cite, jsonld_paper,
                                             cites_out, citation_text))
    return out, counts


def cmd_build(src, args):
    out, counts = build_outputs(src, args.cache, with_pages=not args.no_pages)
    root = os.path.join(args.data, "bibliography")
    managed_dirs = ("papers", "systems", "pages")
    if args.check:
        bad = []
        for rel, content in sorted(out.items()):
            p = os.path.join(root, rel)
            data = content if isinstance(content, bytes) else content.encode("utf-8")
            if not os.path.isfile(p):
                bad.append(f"missing {rel}")
            else:
                with open(p, "rb") as f:
                    if f.read() != data:
                        bad.append(f"differs {rel}")
        for d in managed_dirs:
            for dp, _, fns in os.walk(os.path.join(root, d)):
                for fn in fns:
                    rel = os.path.relpath(os.path.join(dp, fn), root)
                    if rel not in out:
                        bad.append(f"extra {rel}")
        for x in bad[:30]:
            print(x)
        print(f"check: {len(out)} files, {len(bad)} problems")
        return 1 if bad else 0
    os.makedirs(root, exist_ok=True)
    for d in managed_dirs:  # remove stale generated files
        for dp, _, fns in os.walk(os.path.join(root, d)):
            for fn in fns:
                rel = os.path.relpath(os.path.join(dp, fn), root)
                if rel not in out:
                    os.remove(os.path.join(dp, fn))
    for rel, content in out.items():
        p = os.path.join(root, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        data = content if isinstance(content, bytes) else content.encode("utf-8")
        if os.path.isfile(p):
            with open(p, "rb") as f:
                if f.read() == data:
                    continue
        with open(p, "wb") as f:
            f.write(data)
    print(json.dumps(counts, indent=1))
    return 0


# ----------------------------------------------------------------------------- documentation
PAPER_FIELDS = [  # (column, Croissant dataType, description) for papers.csv / papers.parquet
    ("id", "sc:Text", "Stable bibliography id: 'doi:<lowercase DOI>' or 'nodoi:<16 hex>' (sha256 of folded first-author surname|year|folded title or printed reference). Primary key."),
    ("key", "sc:Text", "URL-safe file key: 'p' + first 16 hex of sha256(id). Record at papers/<key>.json, page at pages/papers/<key>.html."),
    ("doi", "sc:Text", "DOI in lower case (empty when the work has none)."),
    ("no_doi", "sc:Boolean", "True when the work has no DOI and the id is derived from first author, year and title."),
    ("title", "sc:Text", "Title (Crossref; for works without a Crossref record: the title held by ORCHID's sources, may be empty)."),
    ("authors", "sc:Text", "Authors as 'Given Family', '; '-joined, in Crossref order."),
    ("orcids", "sc:Text", "ORCID iDs deposited with Crossref for the authors, '; '-joined."),
    ("container_title", "sc:Text", "Journal, book or proceedings title."),
    ("volume", "sc:Text", "Volume."),
    ("issue", "sc:Text", "Issue."),
    ("pages", "sc:Text", "Page range or article number."),
    ("year", "sc:Integer", "Publication year (Crossref 'issued'; otherwise as printed)."),
    ("publisher", "sc:Text", "Publisher (Crossref)."),
    ("type", "sc:Text", "Crossref work type (journal-article, book-chapter, ...)."),
    ("metadata_source", "sc:Text", "Where the bibliographic fields come from: 'crossref', or the ORCHID source that printed the reference (orchid_literature_printed_reference, tdb_header, handbook_reference, discovery_snapshot, ...)."),
    ("citation", "sc:Text", "One-line citation built from the metadata (printed reference for works without a Crossref record)."),
    ("systems", "sc:Text", "Chemical systems the work covers in ORCHID (sorted element symbols joined by '-'), ';'-joined. Discovery leads contribute only systems named in their title."),
    ("roles", "sc:Text", "Roles in ORCHID, ';'-joined: corpus_assessment_read, cited_original_source, tdb_source_assessment, tdb_compilation, handbook_diagram_source, discovery_lead, context_source, espei_reference."),
    ("tdb_files", "sc:Text", "TDB files whose header cites this work as the assessment (tdb_source_assessment), ';'-joined."),
    ("handbook_systems", "sc:Text", "Systems whose assessed (handbook) diagram names this work as its source, ';'-joined."),
    ("discovery_rank", "sc:Integer", "Rank in the Crossref higher-order assessment discovery (discovery_lead)."),
    ("values_reported", "sc:Integer", "Literature values ORCHID read from this paper (all released records, accepted or not)."),
    ("values_reported_accepted", "sc:Integer", "Of those, values with flags.extraction_accepted."),
    ("values_reported_experimental", "sc:Integer", "Values read from this paper with evidence class experimental."),
    ("values_reported_model", "sc:Integer", "... evidence class model."),
    ("values_reported_dft", "sc:Integer", "... evidence class dft."),
    ("values_reported_unspecified", "sc:Integer", "... evidence class unspecified."),
    ("values_cited", "sc:Integer", "Values in other papers whose printed source (cited_sources) is this work."),
    ("values_cited_accepted", "sc:Integer", "Of those, accepted values."),
    ("property_groups", "sc:Text", "Property groups of the values reported in or citing this work, ';'-joined."),
    ("properties", "sc:Text", "Property codes (HM_FORM, INVARIANT_TEMPERATURE, ...) of those values, ';'-joined."),
    ("n_tables", "sc:Integer", "Distinct tables of this paper the released values come from."),
    ("n_text_passages", "sc:Integer", "Distinct sentences (running-text quotes) of this paper the released values come from."),
    ("n_espei_datasets", "sc:Integer", "ESPEI datasets of the ORCHID bundle derived from values of this work (paths in papers.jsonl)."),
    ("n_phase_registry_phases", "sc:Integer", "Phase-registry phases that cite this work (paper label or reference)."),
    ("n_assessed_assemblage_records", "sc:Integer", "Assessed phase-assemblage records whose original evaluation is this work."),
    ("open_access_license", "sc:URL", "First Creative Commons licence URL in the Crossref license field (empty: none deposited). Never inferred from elsewhere."),
    ("license_urls", "sc:Text", "Every licence URL in the Crossref license field (TDM, publisher, CC), ';'-joined."),
    ("is_referenced_by_count", "sc:Integer", "Crossref is-referenced-by-count on the retrieval date."),
    ("is_referenced_by_count_retrieved", "sc:Date", "Date the Crossref record was retrieved."),
    ("crossref_status", "sc:Text", "found | not_found (404 at Crossref) | no_doi | not_fetched."),
    ("doi_check", "sc:Text", "consistent | mismatch_suspect: Crossref first author / year / title compared with the reference printed in ORCHID's source (empty: nothing printed to compare)."),
    ("read_by_orchid", "sc:Text", "yes: ORCHID processed the paper's full text in its extraction corpus; no otherwise."),
    ("edges_out", "sc:Integer", "Edges with this paper as src (edges.*)."),
    ("edges_in", "sc:Integer", "Edges with this paper as dst."),
    ("json_url", "sc:URL", "Resolvable JSON record of this paper."),
    ("html_url", "sc:URL", "Static HTML page of this paper."),
]
EDGE_FIELDS = [
    ("id", "sc:Text", "Stable edge id: 'e' + first 16 hex of sha256(type|src|dst|provenance source). Primary key."),
    ("type", "sc:Text", "cites | used_in_assessment | models_system_in_tdb | reassesses | handbook_reproduces | cited_by_crossref."),
    ("src", "sc:Text", "Source node id (paper id, or handbook:okamoto-2010:<SYS> for handbook_reproduces)."),
    ("dst", "sc:Text", "Target node id (paper id, or tdb:<file> for models_system_in_tdb)."),
    ("src_kind", "sc:Text", "paper | handbook_diagram."),
    ("dst_kind", "sc:Text", "paper | tdb_file."),
    ("count", "sc:Integer", "cites / used_in_assessment: number of values carrying that citation."),
    ("provenance_source", "sc:Text", "orchid_extraction | crossref_reference_list | sgte_tdb_header | orchid_literature_model_link | verified_assessment | okamoto_handbook_source_line | citation+same-system."),
    ("provenance_dataset", "sc:Text", "Input the edge was read from."),
    ("provenance_method", "sc:Text", "How the edge was established."),
    ("systems", "sc:Text", "Chemical system(s) of the edge, ';'-joined."),
    ("attrs", "sc:Text", "Remaining attributes as compact JSON (evidence-class and citation-status counts, TDB sha256, page, ...)."),
]


def _field_table(fields):
    rows = ["| column | type | description |", "|---|---|---|"]
    rows += [f"| `{c}` | {t} | {d} |" for c, t, d in fields]
    return "\n".join(rows)


def docs(counts, cite):
    c = counts
    role_rows = "\n".join(f"| `{r}` | {c['papers_by_role'][r]} | {c['papers_by_role_with_doi'][r]} |"
                          for r in ROLES + ["espei_reference"])
    edge_rows = "\n".join(f"| `{t}` | {c['edges_by_type'][t]} |" for t in EDGE_TYPES)
    cr = c["crossref"]
    readme = f"""# ORCHID bibliography

Every paper ORCHID knows about — the assessments and data papers it read, the original sources its values cite, the
assessments behind its TDB files, the original evaluations behind the assessed (handbook) diagrams, and discovery
leads — with Crossref bibliographic metadata, the role each paper plays in ORCHID, counts of the evidence ORCHID
released from it, and typed relations between papers, TDB files and handbook diagrams.

Release `{RELEASE_TAG}` (bibliography {BIB_VERSION}). Built by `tools/build_bibliography.py` (orchid-project.org;
deterministic, offline from a Crossref response cache; `--check` verifies these files).

## Counts

| | papers | with DOI |
|---|---:|---:|
| all | {c['papers']} | {c['papers_with_doi']} |
{role_rows}

Papers without a DOI: {c['papers_without_doi']} (ids derived from first author, year and title, flagged `no_doi`).
Crossref records: {cr.get('found', 0)} found, {cr.get('not_found', 0)} DOIs unknown to Crossref (404),
{cr.get('not_fetched', 0)} not fetched; retrieved {' to '.join(c['crossref_retrieved'] or ['-'])}.
DOI consistency check against the printed reference: {c['doi_check'].get('consistent', 0)} consistent,
{c['doi_check'].get('mismatch_suspect', 0)} mismatch suspects (flagged, not removed).
Open-access licence deposited with Crossref: {c['open_access_license']} papers. Chemical systems: {c['systems']}.

| edge type | edges |
|---|---:|
{edge_rows}

## Files

| file | content |
|---|---|
| `papers.jsonl` | one full record per paper (fields in SCHEMA.md) |
| `papers.csv`, `papers.parquet` | flat table, one row per paper (stable `id`) |
| `edges.jsonl`, `edges.csv`, `edges.parquet` | typed relations with provenance (stable `id`) |
| `papers/<key>.json` | one small JSON per paper: record + its relations (resolvable URL) |
| `index.json` | release info, counts, input hashes, id rule, URL templates and a compact search table (`columns`, `rows`) |
| `systems_index.json` | system -> role -> paper ids |
| `systems/<SYS>.json` | one system: papers by role, assessments with the data papers they used (and counts) |
| `orchid-bibliography.bib`, `orchid-bibliography.csl.json` | BibTeX and CSL-JSON generated from Crossref metadata (Crossref-backed records only) |
| `papers.jsonld`, `context.jsonld` | schema.org ScholarlyArticle graph (ORCID author ids, isPartOf journal, schema:citation edges, roles / evidence as additionalProperty) |
| `pages/` | static HTML pages per paper and per system (no JavaScript; JSON-LD and COinS embedded) |
| `parquet.json` | row counts, content hashes and Croissant field descriptions of the Parquet tables |
| `SCHEMA.md`, `CITATION.cff` | field reference; how to cite |

## How to use

```python
import pandas as pd
papers = pd.read_parquet("https://data.orchid-project.org/bibliography/papers.parquet")
edges = pd.read_parquet("https://data.orchid-project.org/bibliography/edges.parquet")
# original data papers used by the assessments of Al-Ni
import json, urllib.request
doc = json.load(urllib.request.urlopen("https://data.orchid-project.org/bibliography/systems/AL-NI.json"))
for a in doc["assessments"]:
    print(a["year"], a["title"], a["n_data_papers_used"])
```

One record: `https://data.orchid-project.org/bibliography/papers/<key>.json` with `key = "p" + sha256(id)[:16]`
(look up `key` by DOI in `index.json`). The ORCHID MCP server (`orchid-mcp`) exposes `get_paper`,
`search_papers`, `get_paper_relations` and `get_assessment_evidence` over these files.

## What is and is not here

- Bibliographic metadata only (no abstracts, no full text, no PDFs, no local paths). Crossref metadata is used as
  deposited; `doi_check` flags DOIs whose Crossref record disagrees with the reference ORCHID's source printed.
- Roles and evidence counts are derived from the ORCHID release files (literature values, ESPEI bundle, phase
  registry, assessed phase assemblages) and from the SGTE-corpus link tables; see `provenance` in each record.
- Edges are never inferred: `cites` / `used_in_assessment` come from the citation printed for each value ORCHID
  extracted; `cited_by_crossref` only from reference lists publishers deposited openly with Crossref, kept
  separate from ORCHID's own extraction and restricted to works in this bibliography; `reassesses` only when the
  later assessment cites the earlier one and both model the same exact system in a TDB file (method
  `citation+same-system`, never from dates alone). `supersedes` is not emitted (no source states supersession);
  `same_system` is not emitted as edges (use `systems_index.json`). No `context_source` papers yet: no same-system
  context manifest exists in the ORCHID sources.
- `discovery_lead` papers were found by the Crossref higher-order assessment discovery and were not read; their
  `candidate_systems_inferred` (citation closures) are hypotheses, only title-named systems enter `systems`.
- Open-access status is taken only from the Crossref `license` field.

## Licence

Crossref bibliographic metadata: CC0 (as distributed by Crossref). ORCHID-derived fields (roles, evidence counts,
relations, systems, checks): CC BY 4.0 — cite this bibliography: {cite['text']}.
"""
    schema = f"""# ORCHID bibliography: schema

Release `{RELEASE_TAG}`, bibliography {BIB_VERSION}. Generated by `tools/build_bibliography.py` (deterministic;
`--check`). Croissant-ready field descriptions: the column tables below (also in `parquet.json`).

## Ids

- Paper: `doi:<lowercase DOI>`; works without DOI `nodoi:<16 hex>` = sha256 of `folded first-author surname|year|folded
  title (or printed reference)`; `no_doi: true`. Literature ids of the same work (`ref:...`, `doi:...` in
  `literature/papers.json`) are listed in `aliases`.
- File key: `p` + first 16 hex of sha256(id); `papers/<key>.json`, `pages/papers/<key>.html`.
- Other nodes (edges only): `tdb:<file>` (a TDB file of the ORCHID catalog), `handbook:okamoto-2010:<SYS>` (the
  assessed diagram of a binary system in {HANDBOOK}).
- Edge: `e` + first 16 hex of sha256(`type|src|dst|provenance source`).

## Paper record (`papers.jsonl`, `papers/<key>.json` -> `paper`)

| field | meaning |
|---|---|
| `id`, `key`, `doi`, `doi_url`, `id_basis`, `no_doi` | identifiers (above) |
| `title`, `subtitle`, `authors[{{given, family, name, orcid, orcid_authenticated}}]`, `container_title`, `container_title_short`, `volume`, `issue`, `pages`, `article_number`, `year`, `issued`, `publisher`, `type`, `issn[]`, `isbn[]`, `language` | Crossref metadata (`metadata_source: crossref`); otherwise the fields ORCHID's sources hold |
| `metadata_source` | `crossref` or the ORCHID source of the fields |
| `citation` | one-line citation; `citation_printed`: the reference as printed in an ORCHID source when it is not Crossref-derived |
| `crossref` | `{{status, retrieved_at, references_deposited, references_count}}` |
| `license[]` | Crossref `license` entries `{{url, content_version, start, delay_in_days}}` |
| `is_referenced_by_count`, `is_referenced_by_count_retrieved` | Crossref citation count and its retrieval date |
| `doi_check` | `{{status: consistent|mismatch_suspect|unchecked, issues[]}}`: Crossref first author / year / title vs the reference printed in ORCHID's source; null when nothing printed to compare |
| `systems[]` | chemical systems (sorted element symbols joined by `-`) |
| `candidate_systems_inferred` | discovery leads: candidate systems from citation closures, with basis (hypotheses; not in `systems`) |
| `roles[]`, `role_details{{}}` | roles (below) and their details |
| `evidence_summary` | `reported` (values read from this paper: `n_values`, `n_values_accepted`, `by_evidence_class`, `by_property_group`, `by_property`, `n_tables`, `n_text_passages`, `n_figures` (0: no released value comes from a figure)), `reader_agreement` (`by_reader_tier`, `by_reader_agreement` when the release carries it), `cited_as_source` (values in other papers citing it: counts as above and `n_citing_papers`), `espei_datasets[]` (paths in the ESPEI bundle), `phase_registry_phases[]`, `assessed_assemblage_records` |
| `availability` | `open_access_license` (`{{url, content_version, start, delay_in_days, basis: crossref_license_field}}` or null), `read_by_orchid` (`yes`/`no`), `read_detail[]` (`values_extracted`, `processed_no_released_values`, `full_text_screened_for_origin`), `pdf_redistributed: false` |
| `aliases[]`, `provenance[]` | literature ids of the same work; inputs that contributed to the record |
| `relation_counts` | `{{edge type: {{out, in}}}}` |
| `urls` | `json`, `html` |

`papers/<key>.json` adds `relations[]`: `{{direction: out|in, type, edge_id, other{{id, key, doi, title, year, roles}}, count, provenance, attrs}}`.

## Roles

| role | meaning, source |
|---|---|
| `corpus_assessment_read` | ORCHID processed the paper's full text (SGTE-corpus census) and/or released values read from it (`literature/papers.json` reporting paper); `values_released` |
| `cited_original_source` | cited by a paper ORCHID read as the source of a value (`cited_sources` of the literature records) |
| `tdb_source_assessment` | the assessment a TDB file's header cites (SGTE collection: header citation matched to the Hallstedt 2025 reference list, `tdb_pdf_links.json`); also literature model links and hand-verified assessments. `tdb_files[{{file, path, sha256, bundled, system, link_status, link_match, basis}}]` |
| `tdb_compilation` | the publication of a TDB collection (Hallstedt 2025, SGTE binary collection) |
| `handbook_diagram_source` | original evaluation printed under an assessed diagram of {HANDBOOK}; `systems` (citations only; no handbook text) |
| `discovery_lead` | found by the Crossref higher-order assessment discovery (`open-tdb` branch `{DISCOVERY_REF}`, `{DISCOVERY_PATH}`), not read: `rank`, `score`, `max_order`, `systems_named_in_title`, `candidate_systems_inferred`, `discovery_sources`, `license_review` |
| `context_source` | same-system context source (none in this release) |
| `espei_reference` | cited by an ESPEI dataset and not otherwise known |

## Edges (`edges.jsonl`)

`{{id, type, src, dst, src_kind, dst_kind, count, provenance{{source, dataset, method, ...}}, attrs{{}}}}`

| type | meaning |
|---|---|
| `cites` | paper ORCHID read -> work printed as the source of its values; `count` = values carrying that citation; `attrs`: `n_values_accepted`, `evidence_class`, `citation_status` (resolved / author_mismatch / surname_mismatch), `property`, `systems` |
| `used_in_assessment` | inverse of `cites`: original data paper -> the paper (assessment) that used it, same attributes |
| `models_system_in_tdb` | paper -> `tdb:<file>`; `system`, `tdb_sha256`, link status |
| `reassesses` | later assessment -> earlier assessment of the same exact system; only when the later cites the earlier (`citation_basis`: orchid_extraction and/or crossref_reference_list) and both have a TDB of that system; method `citation+same-system` |
| `handbook_reproduces` | `handbook:okamoto-2010:<SYS>` -> original evaluation; `page`, `doi_matched_by` |
| `cited_by_crossref` | src's publisher-deposited Crossref reference list contains dst (both in this bibliography); `doi_asserted_by`; kept apart from ORCHID's own extraction |

## Flat tables

### papers.csv / papers.parquet

{_field_table(PAPER_FIELDS)}

### edges.csv / edges.parquet

{_field_table(EDGE_FIELDS)}

## Licence

Crossref metadata CC0; ORCHID-derived fields CC BY 4.0.
"""
    cff = f"""cff-version: 1.2.0
message: "If you use these data, please cite this Zenodo record and the papers whose data you use."
title: "ORCHID bibliography"
type: dataset
authors:
  - family-names: "Arróyave"
    given-names: "Raymundo"
    orcid: "https://orcid.org/0000-0001-7548-8686"
    affiliation: "Texas A&M University"
  - family-names: "Saal"
    given-names: "James"
    orcid: "https://orcid.org/0000-0003-0935-158X"
    affiliation: "Citrine Informatics"
  - family-names: "Shin"
    given-names: "Dongwon"
    orcid: "https://orcid.org/0000-0002-5797-3423"
    affiliation: "Oak Ridge National Laboratory"
version: "{VERSION}"
doi: "{cite['doi']}"
license: CC-BY-4.0
url: "https://data.orchid-project.org/bibliography/"
repository-code: "https://github.com/orchid-calphad/orchid-data"
abstract: "Every paper the ORCHID CALPHAD data project knows about ({c['papers']} works, {c['papers_with_doi']} with DOI): Crossref bibliographic metadata (CC0), the role each paper plays in ORCHID (read, cited original source, assessment behind a TDB file, original evaluation behind an assessed diagram, discovery lead), counts of the evidence ORCHID released from it, and {c['edges']} typed relations with provenance (citations of values, TDB files, reassessments, handbook diagrams, Crossref reference lists). Version 1, machine-extracted."
identifiers:
  - type: doi
    value: "{cite['part_of']['doi']}"
    description: "Umbrella record ORCHID data (all datasets of release {RELEASE_TAG})"
keywords: [CALPHAD, thermodynamics, bibliography, citation graph, phase diagrams, AI-ready data]
"""
    return {"README.md": readme, "SCHEMA.md": schema, "CITATION.cff": cff}


# ----------------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("command", nargs="?", default="build", choices=["build", "fetch", "dois"])
    ap.add_argument("--data", default=DEF_DATA, help="orchid-data checkout (writes bibliography/)")
    ap.add_argument("--open-tdb", default=DEF_OPEN_TDB, help="open-tdb checkout (read only)")
    ap.add_argument("--cache", default=DEF_CACHE, help="Crossref response cache (git-ignored)")
    ap.add_argument("--rate", type=float, default=4.0, help="fetch: requests per second (capped at 5)")
    ap.add_argument("--workers", type=int, default=4, help="fetch: concurrent requests")
    ap.add_argument("--limit", type=int, default=None, help="fetch: at most this many DOIs")
    ap.add_argument("--check", action="store_true", help="build: exit 1 if bibliography/ differs from a fresh build")
    ap.add_argument("--no-pages", action="store_true", help="build: skip the static HTML pages (tools/bibliography_pages.py)")
    args = ap.parse_args(argv)
    src = Sources(args.data, args.open_tdb)
    if args.command == "dois":
        ds = sorted(src.all_dois())
        cached = sum(os.path.isfile(cache_path(args.cache, d)) for d in ds)
        print(f"{len(ds)} DOIs, {cached} cached")
        return 0
    if args.command == "fetch":
        return cmd_fetch(src, args.cache, args.rate, args.workers, args.limit)
    return cmd_build(src, args)


if __name__ == "__main__":
    raise SystemExit(main())
