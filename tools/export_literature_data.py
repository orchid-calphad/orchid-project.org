#!/usr/bin/env python3
"""Export the public ORCHID literature-data release from the SGTE-corpus record store.

Reads (never writes) the db-v2.1 record store produced by the open-tdb SGTE-corpus campaign:

  invariants/master/*.jsonl   invariant-reaction tables (source, artifact, reading, observation records)
  general/records/*.jsonl     general thermodynamic / phase-equilibrium / structure tables
  text/*.jsonl                running-text values, each linked to one verbatim quote and page
  context/*.jsonl             paper-level experimental context (samples, treatments, measurements)

and writes derived data only (values, units, conditions, phases, compositions, evidence class, method,
citations, one verbatim quote per value with its page) under data/literature/:

  index.json                  systems -> counts by evidence class / property group, file list
  systems/<SYS>.json          records of one chemical system (+ the papers and contexts they reference)
  papers.json                 paper id -> citation, DOI, year, title, kind, systems
  contexts.json               paper id -> experimental context items (quote + page each)
  orchid-literature.jsonl.gz  every record, one per line
  orchid-literature.csv.gz    every record, flattened
  SCHEMA.md                   field documentation (written by this script)

Never exported: page images, crops, renders, full text, full table text, PDFs, reader notes, the content of
dropped (unquoted) fields. Okamoto handbook material is out of scope.

Deterministic: same inputs give byte-identical outputs (sorted keys, gzip mtime 0). Standard library only.

  python3 tools/export_literature_data.py            # write
  python3 tools/export_literature_data.py --check    # verify the files on disk are current (exit 1 if not)
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import hashlib
import importlib.util
import io
import json
import os
import re
import sys
from collections import Counter, defaultdict

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import doi_corrections  # noqa: E402
import literature_quality as quality  # noqa: E402
import literature_gates as gates  # noqa: E402
import phase_rule  # noqa: E402
import phase_label_check  # noqa: E402
import orchid_cite  # noqa: E402
import uncertainty_budget as ub  # noqa: E402
import ml_folds  # noqa: E402

CORR = doi_corrections.load()   # shared DOI correction table (tools/doi_corrections.json)

DEFAULT_CORPUS = os.path.expanduser("~/projects/open-tdb/build/sgte-corpus")
DEFAULT_PDFS = os.path.expanduser("~/data/PD-extraction")
DEFAULT_OPEN_TDB = os.path.expanduser("~/projects/open-tdb")
SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_REL = os.path.join("data", "literature")
RELEASE = "orchid-literature 2026-10-23 (record store db-v2.1: invariant/general v2.1.19r, text+context fixup-20261022)"
MAX_SYSTEM_BYTES = 1_900_000

EVIDENCE = {"experimental": "experimental", "model_value": "model", "model": "model",
            "computed": "dft", "dft": "dft", "unspecified": "unspecified"}

TEXT_QUANTITY = {"formation_enthalpy": "HM_FORM", "mixing_enthalpy": "HM_MIX", "gibbs_energy": "GIBBS_ENERGY",
                 "entropy": "ENTROPY", "partial_enthalpy": "PARTIAL_ENTHALPY", "heat_capacity": "CPM",
                 "heat_content": "HEAT_CONTENT", "emf": "EMF", "activity": "ACTIVITY",
                 "activity_coefficient": "ACTIVITY_COEFFICIENT", "vapor_pressure": "VAPOR_PRESSURE",
                 "chemical_potential": "CHEMICAL_POTENTIAL", "lattice_parameter": "LATTICE_UNSPECIFIED_AXIS",
                 "other": "OTHER"}

GENERIC_METHODS = {"as_reported", "as_reported_in_text", "as_reported_in_assessment_table", ""}


def property_group(obs_type: str, q: str | None) -> str:
    if obs_type == "invariant_reaction":
        return "invariant reaction"
    if obs_type == "phase_transition":
        return "phase transition"
    if obs_type in ("phase_region_vertex_set", "phase_assemblage_node"):
        return "phase equilibrium"
    if obs_type == "structural_property":
        return "crystal structure"
    q = (q or "").upper()
    if obs_type == "property_expression":
        return "model expression"
    if q.startswith(("HM_FORM", "SM_FORM", "GM_FORM")) or "FORM" in q:
        return "formation"
    if q.startswith(("HM_MIX", "GM_MIX", "SM_MIX", "PARTIAL", "EXCESS")) or "MIX" in q:
        return "mixing / partial"
    if q.startswith(("CPM", "CP", "HEAT_CONTENT", "SM", "ENTROPY", "HM")):
        return "heat capacity / entropy / enthalpy"
    if q.startswith(("ACR", "ACTIVITY", "EMF", "CHEMICAL_POTENTIAL", "MU", "VAPOR", "GIBBS", "GM")):
        return "activity / EMF / Gibbs energy"
    if q.startswith(("LATTICE", "MOLAR_VOLUME", "SITE", "VOLUME", "DENSITY")):
        return "crystal structure"
    return "other property"


# ----------------------------------------------------------------------------------------------- helpers

def jl(path):
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                yield json.loads(line)


def num(v):
    if isinstance(v, bool) or v is None:
        return None
    if isinstance(v, (int, float)):
        if v != v or v in (float("inf"), float("-inf")):
            return None
        return round(float(v), 10) if isinstance(v, float) else v
    return None


def clean(obj):
    """Drop None / empty containers recursively, round floats (stable output)."""
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            v = clean(v)
            if v is None or v == {} or v == [] or v == "":
                continue
            out[k] = v
        return out
    if isinstance(obj, list):
        return [clean(v) for v in obj if clean(v) not in (None, {}, [], "")]
    if isinstance(obj, float):
        if obj != obj:
            return None
        return float(f"{obj:.10g}")
    return obj


def doi_norm(d):
    if not d or not isinstance(d, str):
        return None
    d = d.strip()
    d = re.sub(r"^(https?://(dx\.)?doi\.org/|doi:)", "", d, flags=re.I)
    return d.lower() or None


SHA2DOI = {}  # pdf sha256 prefix (16 hex) -> DOI, for corpus sources whose record lacks the DOI
SHA2CORR = {}  # pdf sha256 prefix -> DOI correction entry (PDF filed under a DOI that names no article)
KEY2CORR = {}  # paper key -> DOI correction entry whose DOI the paper lost


def lit_doi(d):
    """(DOI, correction entry): the DOI is None when tools/doi_corrections.json drops it for literature papers."""
    d = doi_norm(d)
    e = CORR.literature(d) if d else None
    return (None, e) if e else (d, None)


def paper_key(src: dict | None):
    if not src:
        return None
    d, corr = lit_doi(src.get("doi"))
    if d:
        return "doi:" + d
    sha = src.get("file_sha256")
    if sha:
        if sha[:16] in SHA2DOI:
            return "doi:" + SHA2DOI[sha[:16]]
        key = "sha256:" + sha[:16]
        corr = SHA2CORR.get(sha[:16]) or corr
    elif corr:
        key = corr["replacement_id"] if corr["context"] == "literature_cited_source" else \
            "nodoi:" + hashlib.sha256(json.dumps(src.get("citation"), sort_keys=True).encode()).hexdigest()[:16]
    else:
        return src.get("pid") or src.get("id")
    if corr:
        KEY2CORR.setdefault(key, corr)
    return key


def system_name(elems):
    return "-".join(sorted({e.upper() for e in elems}))


def reader_of(reading_id: str) -> str:
    r = reading_id.rsplit(":", 1)[-1]
    r = re.sub(r"-(\d+)-(\d+)(?=-|$)", r"-\1.\2", r)
    # v2.1.19r: the blind second Claude read 'claude-opus-5.5#2' (reading-id slug claude-opus-5-5-2) keeps its '#2'
    r = re.sub(r"^(claude-[a-z]+-\d+\.\d+)-(\d+)$", r"\1#\2", r)
    return "text_layer" if r == "text-layer" else r


# _export.reader_agreement (record builder, annotate_reader_agreement.py): how many independent visual readers carry
# the value. Published per value with a human label; the data explorer shows a "single reader" badge.
READER_AGREEMENT_LABELS = {
    "multi_reader": "confirmed (two or more readers)",
    "single_reader": "unconfirmed (single reader)",
    "single_reader_text_layer": "single reader + PDF text layer",
}


def reader_agreement(exp):
    """(code, label) from an observation's _export, or (None, None) when the record builder did not set it."""
    code = (exp or {}).get("reader_agreement")
    if code not in READER_AGREEMENT_LABELS:
        return None, None
    return code, READER_AGREEMENT_LABELS[code]


# ----------------------------------------------------------------------------------------------- bibliography

class Biblio:
    def __init__(self, corpus, pdfs):
        self.cache_dir = os.path.join(corpus, "crossref-cache")
        self.header = {}          # doi -> header citation (TDB header)
        self.header_by_sha = {}   # file sha16 -> header citation
        links = os.path.join(corpus, "tdb_pdf_links.json")
        if os.path.exists(links):
            data = json.load(open(links, encoding="utf-8"))
            for ln in data.get("links", []):
                d = doi_norm(ln.get("doi"))
                cit = ln.get("header_citation")
                # a header citation whose matched DOI is a mismatch (tools/doi_corrections.json) does not describe it
                tdb_corr = CORR.tdb(ln.get("tdb"), d) if d else None
                if d and cit and not tdb_corr:
                    self.header.setdefault(d, (cit, ln.get("year")))
                if ln.get("pdf") and ln.get("status") == "linked_pdf_local":
                    p = os.path.join(pdfs, ln["pdf"])
                    if os.path.exists(p):
                        h = hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
                        # the PDF was fetched under the DOI, so it is that DOI's work, unless the DOI names no article
                        pdf_doi = None if (d and CORR.anywhere(d)) else d
                        if pdf_doi:
                            SHA2DOI.setdefault(h, pdf_doi)
                        elif cit:
                            self.header_by_sha.setdefault(h, (cit, ln.get("year")))
                            if d:
                                SHA2CORR.setdefault(h, tdb_corr or CORR.anywhere(d))
        self.meta = {}
        mp = os.path.join(corpus, "db-v2.1", "cited-meta", "meta.json")
        if os.path.exists(mp):
            for v in json.load(open(mp, encoding="utf-8")).values():
                d = doi_norm(v.get("doi"))
                if d and v.get("status") == "found":
                    self.meta.setdefault(d, v)

    def crossref(self, doi):
        key = hashlib.sha256(doi.lower().encode()).hexdigest()[:20]
        p = os.path.join(self.cache_dir, key + ".json")
        if not os.path.exists(p):
            return None
        try:
            d = json.load(open(p, encoding="utf-8"))
        except Exception:
            return None
        if doi_norm(d.get("DOI")) != doi:
            return None
        return d

    def describe(self, doi):
        """-> dict(title, authors, year, container, volume, page) from local caches only."""
        if not doi:
            return {}
        cr = self.crossref(doi)
        if cr:
            auth = []
            for a in cr.get("author", []) or []:
                nm = " ".join(x for x in (a.get("given"), a.get("family")) if x) or a.get("name")
                if nm:
                    auth.append(nm)
            year = None
            for k in ("published-print", "published-online", "issued", "created"):
                dp = (cr.get(k) or {}).get("date-parts")
                if dp and dp[0] and dp[0][0]:
                    year = dp[0][0]
                    break
            return {"title": (cr.get("title") or [None])[0], "authors": auth, "year": year,
                    "container": (cr.get("short-container-title") or cr.get("container-title") or [None])[0],
                    "volume": cr.get("volume"), "page": cr.get("page"), "metadata_source": "crossref"}
        m = self.meta.get(doi)
        if m:
            return {"title": m.get("title"), "authors": [m["first_author"]] if m.get("first_author") else [],
                    "year": m.get("year"), "container": m.get("venue"), "metadata_source": "openalex/crossref"}
        return {}


def repair_replacement_chars(obj, names):
    """Source records carry U+FFFD where an accented letter was lost upstream (e.g. 'Fern\ufffdndez').
    Restore each damaged word from a known name (Crossref author names in the local cache) when exactly one
    candidate matches with the damaged positions as wildcards; otherwise drop the U+FFFD (never invent letters)."""
    def fix(s):
        if "\ufffd" not in s:
            return s
        def word(m):
            w = m.group(0)
            pat = re.compile("^" + "".join("." if ch == "\ufffd" else re.escape(ch) for ch in w) + "$")
            idx = [i for i, ch in enumerate(w) if ch == "\ufffd"]
            hits = {n for n in names if pat.match(n) and all(ord(n[i]) > 127 for i in idx)}  # the lost letter was non-ASCII
            return hits.pop() if len(hits) == 1 else w.replace("\ufffd", "")
        return re.sub(r"[^\s,;.()\[\]\"]*\ufffd[^\s,;.()\[\]\"]*", word, s)
    if isinstance(obj, str):
        return fix(obj)
    if isinstance(obj, list):
        return [repair_replacement_chars(x, names) for x in obj]
    if isinstance(obj, dict):
        return {k: repair_replacement_chars(v, names) for k, v in obj.items()}
    return obj


def crossref_name_words(cache_dir):
    words = set()
    for f in os.listdir(cache_dir):
        if not f.endswith(".json"):
            continue
        try:
            d = json.load(open(os.path.join(cache_dir, f), encoding="utf-8"))
        except Exception:
            continue
        for a in d.get("author", []) or []:
            for part in (a.get("given"), a.get("family")):
                for w in re.split(r"[\s,]+", part or ""):
                    if w and "\ufffd" not in w:
                        words.add(w)
    return words


# Energy-unit normalization (formation / mixing / partial, kJ/mol-atom only when the per-atom basis is explicit)
# lives in tools/literature_quality.py (normalize_energy_unit).


def format_citation(meta, fallback=None):
    if meta.get("title") or meta.get("authors"):
        a = meta.get("authors") or []
        auth = ", ".join(a[:6]) + (" et al." if len(a) > 6 else "")
        parts = [p for p in (auth, meta.get("title")) if p]
        venue = " ".join(str(x) for x in (meta.get("container"), meta.get("volume")) if x)
        if meta.get("page"):
            venue += (", " if venue else "") + str(meta["page"])
        if meta.get("year"):
            venue += f" ({meta['year']})"
        if venue.strip():
            parts.append(venue.strip())
        return ", ".join(parts) + "."
    return fallback


# ----------------------------------------------------------------------------------------------- loading

class Store:
    def __init__(self, corpus, pdfs):
        self.corpus = corpus
        self.db = os.path.join(corpus, "db-v2.1")
        self.bib = Biblio(corpus, pdfs)
        self.sources = {}         # source id/alias -> source record
        self.papers = {}          # paper key -> paper dict
        self.excluded = Counter()
        self.notes = Counter()

    # -- sources / papers
    def add_source(self, s):
        d, sha = doi_norm(s.get("doi")), s.get("file_sha256")
        if d and CORR.anywhere(d):
            d = None   # a DOI that names no article does not identify the PDF
        if d and sha and s.get("kind") in ("journal_article", "book_chapter"):
            SHA2DOI.setdefault(sha[:16], d)
        for k in [s["id"]] + list(s.get("aliases") or []):
            self.sources.setdefault(k, s)

    def source_in_file(self, sid, local):
        """Source record for id `sid`. Ids derived from a DOI that names no article (tools/doi_corrections.json,
        scope everywhere: e.g. source:doi:10-1557-jmr) collide across different PDFs; for those the file's own source
        record is used."""
        s = self.sources.get(sid)
        if s is not None and sid in local and CORR.anywhere(doi_norm(s.get("doi"))):
            return local[sid]
        return s

    def paper_for(self, src, role):
        key = paper_key(src)
        if not key:
            return None
        p = self.papers.get(key)
        if p is None:
            p = self.papers[key] = {"id": key, "kind": set(), "systems": set(), "role": set(),
                                    "_src": src, "doi": key[4:] if key.startswith("doi:") else None}
        p["kind"].add(src.get("kind") or "unknown")
        p["role"].add(role)
        return key

    def paper_from_cited(self, c):
        """A cited work printed in a table/text reference list (dict from _export.cited / attribution.cited)."""
        d, corr = lit_doi(c.get("doi"))
        if corr:   # DOI dropped as a mismatch: the printed reference is kept under a no-DOI key
            key = corr["replacement_id"] if corr["context"] == "literature_cited_source" else \
                "nodoi:" + hashlib.sha256((c.get("text") or "").encode()).hexdigest()[:16]
            KEY2CORR.setdefault(key, corr)
        else:
            key = ("doi:" + d) if d else None
        if key is None:
            return None
        p = self.papers.get(key)
        if p is None:
            p = self.papers[key] = {"id": key, "kind": {"cited_work"}, "systems": set(), "role": set(),
                                    "_src": {"doi": d, "citation": {"printed_reference": c.get("text")},
                                             "kind": "cited_work"}, "doi": d}
        p["role"].add("cited_source")
        if c.get("text") and not p["_src"].get("citation", {}).get("printed_reference"):
            p["_src"].setdefault("citation", {})["printed_reference"] = c.get("text")
        return key

    def finalize_papers(self):
        out = {}
        for key, p in self.papers.items():
            src = p["_src"]
            doi = p["doi"]
            meta = self.bib.describe(doi) if doi else {}
            cit = src.get("citation") or {}
            printed = cit.get("printed_reference") if isinstance(cit, dict) else None
            if not printed and src.get("kind") == "unresolved_citation":
                ca = (src.get("cited_as") or [{}])[0]
                printed = f"reference [{ca.get('label')}] of {ca.get('in', '').replace('source:doi:', 'doi:')} (no printed entry could be read)"
            header = None
            if doi and doi in self.bib.header:
                header = self.bib.header[doi]
            elif key.startswith("sha256:") and key[7:] in self.bib.header_by_sha:
                header = self.bib.header_by_sha[key[7:]]
            title_field = cit.get("title") if isinstance(cit, dict) else None
            fallback = printed or (header[0] if header else None) or (title_field or None)
            year = meta.get("year") or (header[1] if header else None) or (cit.get("year") if isinstance(cit, dict) else None)
            try:
                year = int(str(year)[:4]) if year else None
            except ValueError:
                year = None
            kinds = sorted(p["kind"])
            kind = "corpus_paper" if any(k in ("journal_article", "book_chapter") for k in kinds) else (
                "cited_work" if "cited_work" in kinds else kinds[0])
            corr = KEY2CORR.get(key) if not doi else None
            out[key] = clean({
                "id": key,
                "doi": doi,
                "doi_status": doi_corrections.STATUS if corr else None,
                "doi_correction": corr["id"] if corr else None,
                "doi_url": f"https://doi.org/{doi}" if doi else None,
                "citation": format_citation(meta, fallback) or (f"doi:{doi}" if doi else "(citation not resolved)"),
                "title": meta.get("title"),
                "authors": meta.get("authors"),
                "year": year,
                "container": meta.get("container"),
                "kind": kind,
                "source_kinds": kinds,
                "roles": sorted(p["role"]),
                "systems": sorted(p["systems"]),
                "metadata_source": meta.get("metadata_source") or ("tdb_header" if header else ("printed_reference" if printed else None)),
                "pdf_sha256_prefix": (src.get("file_sha256") or "")[:16] or None,
            })
        return dict(sorted(out.items()))


# ----------------------------------------------------------------------------------------------- records

def phases_of(obs):
    out = []
    statephases = {}
    for sp in ((obs.get("state") or {}).get("phases") or []):
        if sp.get("printed_label") is not None:
            statephases.setdefault(sp["printed_label"], sp)
    mp = (obs.get("_export") or {}).get("model_phases") or {}
    for ph in obs.get("phases") or []:
        lab = ph.get("printed_label")
        sp = statephases.get(lab) or {}
        out.append(clean({"label": lab, "model_phase": sp.get("model_phase") or mp.get(lab),
                          "mapping_confidence": ph.get("mapping_confidence") or sp.get("mapping_confidence"),
                          "pearson_symbol": sp.get("pearson_symbol"), "prototype": sp.get("prototype")}))
    return out


def comp_basis(conversions):
    """Composition basis of stored phase compositions, from the builder's conversion log."""
    c = {x.split(":")[0] for x in conversions or []}
    if c & {"composition_basis_disputed", "composition_basis_mixed", "composition_basis_other",
            "composition_range_basis_disputed", "composition_range_basis_other",
            "composition_range_basis_weight_percent"}:
        return "unresolved"
    if "mass_percent_to_mole_fraction" in c:
        return "mole_fraction (converted from mass percent)"
    if "atomic_percent_to_fraction" in c:
        return "mole_fraction (converted from atomic percent)"
    return "mole_fraction"


def phase_comps(obs, basis):
    out = []
    for pc in obs.get("phase_compositions") or []:
        comp = {k.upper(): num(v) for k, v in (pc.get("composition") or {}).items() if num(v) is not None}
        if not comp:
            continue
        b = basis
        if any(v > 1.0001 or v < -1e-9 for v in comp.values()):
            b = "unresolved"
        out.append(clean({"phase": (pc.get("phase") or {}).get("printed_label"), "composition": comp,
                          "basis": b, "composition_source": pc.get("composition_basis")}))
    return out


def model_link(exp):
    me = (exp or {}).get("model_evidence")
    if not me:
        return None
    tdb = me.get("tdb") or {}
    mc = me.get("model_citation") or {}
    mdoi = doi_norm(mc.get("doi"))
    corr = CORR.tdb(tdb.get("name"), mdoi) if mdoi else None   # TDB header citation -> DOI mismatch
    return clean({"tdb": tdb.get("name"), "tdb_path": tdb.get("path"), "tdb_sha256": tdb.get("sha256"),
                  "catalog_entry": tdb.get("catalog_entry"), "link_basis": me.get("link_basis"),
                  "model_citation": ((mc.get("citation") or {}).get("title") or None),
                  "model_doi": None if corr else mdoi,
                  "model_doi_status": doi_corrections.STATUS if corr else None,
                  "model_doi_correction": corr["id"] if corr else None,
                  "reproduced_by_calculation": me.get("model_reproduction")})


def cited_list(cited, store):
    out = []
    for c in cited or []:
        pk = store.paper_from_cited(c)
        d, corr = lit_doi(c.get("doi"))
        out.append(clean({"label": c.get("label"), "reference": c.get("text"), "doi": d,
                          "doi_status": doi_corrections.STATUS if corr else None,
                          "doi_correction": corr["id"] if corr else None,
                          "first_author": c.get("first_author"), "year": c.get("year"),
                          "kind": c.get("work_kind") or c.get("kind"), "status": c.get("status"),
                          "paper": pk}))
    return out


# tolerance of the record builder's CALPHAD cross-check (master_store.PLAUSIBLE_DT_K / PLAUSIBLE_DX), written when the
# record does not carry its own
PLAUSIBILITY_TOLERANCE = {"dT_K": 25.0, "dx": 0.05}
CALPHAD_CHECK = {"pass": "agrees", "fail": "disagrees"}


def espei_fields(records, data):
    """espei_ready (record builder's readiness AND the physical phase rule, tools/phase_rule.py), espei_certified
    (espei_ready AND agreement with the CALPHAD prediction) and calphad_check on every record. Values never change."""
    rule = phase_rule.PhaseRule(data)
    counts = Counter()
    log = []
    for r in records:
        fl = r["flags"]
        blocked = [f for f in phase_rule.BLOCKING_FLAGS if fl.get(f)]
        if fl.get("espei_ready") and blocked:
            fl["espei_ready"] = False
            fl["espei_ready_upstream"] = True
            fl["espei_blocked_by"] = blocked
            counts["blocked_by_flag"] += 1
            for f in blocked:
                counts["blocked_by_flag_" + f] += 1
            log.append({"id": r["id"], "system": r["system"], "property": r.get("property"),
                        "T_K": r.get("temperature_K"), "labels": [p.get("label") for p in r.get("phases") or []],
                        "reason": "blocked_by_flag", "detail": ", ".join(blocked)})
        elif fl.get("espei_ready"):
            c = rule.classify(r)
            if c["result"] == "fail":
                fl["espei_ready"] = False
                fl["espei_ready_upstream"] = True
                fl["espei_phase_rule"] = clean({"result": "fail", "reason": c["reason"], "detail": c.get("detail")})
                counts["phase_rule_fail_" + c["reason"]] += 1
                log.append({"id": r["id"], "system": r["system"], "property": r.get("property"),
                            "T_K": r.get("temperature_K"), "labels": [p.get("label") for p in r.get("phases") or []],
                            "reason": c["reason"], "detail": c.get("detail")})
            else:
                fl["espei_phase_rule"] = clean({"result": "pass", "basis": c.get("basis")})
                counts["phase_rule_pass"] += 1
        pl = (fl.get("plausibility") or {}).get("result")
        fl["calphad_check"] = CALPHAD_CHECK.get(pl, "not_available")
        fl["espei_certified"] = bool(fl.get("espei_ready")) and pl == "pass" and not any(
            g.get("effect") == "not_certified" for g in r.get("gate_failures") or [])
        if fl["espei_certified"]:
            fl["espei_certified_basis"] = "pycalphad_invariant_agreement"
        counts["espei_ready"] += bool(fl.get("espei_ready"))
        counts["espei_certified"] += fl["espei_certified"]
        counts["espei_ready_calphad_disagrees"] += bool(fl.get("espei_ready")) and pl == "fail"
        counts["espei_ready_no_prediction"] += bool(fl.get("espei_ready")) and pl not in ("pass", "fail")
    return dict(sorted(counts.items())), log


def quality_flags(obs):
    q = obs.get("quality") or {}
    st = q.get("status") or {}
    pl = st.get("plausibility") or {}
    return clean({
        "extraction_accepted": bool(st.get("extraction_accepted")),
        "row_state": st.get("row_state") or st.get("consensus_state"),
        "review_state": q.get("review_state"),
        "espei_ready": st.get("espei_ready"),
        "state_fully_specified": st.get("state_fully_specified"),
        "specification_gaps": st.get("specification_gaps"),
        "plausibility": clean({"result": pl.get("result"), "model": pl.get("model"),
                               "dT_K": num(pl.get("dT_K")),
                               "tolerance": pl.get("tolerance") or PLAUSIBILITY_TOLERANCE})
        if pl.get("result") not in (None, "not_evaluated") else None,
        "source_disputed": bool(obs.get("source_disputed")) or None,
        "judge_adjudicated": bool(st.get("adjudicated_by")) or None,
        "magnitude_suspect": bool(st.get("magnitude_suspect")) or None,
        "system_mismatch": bool(st.get("system_mismatch")) or None,
    })


def base_record(obs, kind, store, paper, art, readers):
    cond = obs.get("conditions") or {}
    exp = obs.get("_export") or {}
    ev = obs.get("evidence_class") or EVIDENCE.get(obs.get("origin"), "unspecified")
    ev = EVIDENCE.get(ev, ev)
    sysn = system_name(obs.get("system") or [])
    rec = {
        "_oid": obs.get("id"),
        "id": obs.get("pid_compact"),
        "pid": obs.get("pid"),
        "system": sysn,
        "elements": sorted({e.upper() for e in obs.get("system") or []}),
        "record_kind": kind,
        "observation_type": obs.get("observation_type"),
        "evidence_class": ev,
        "temperature_K": num(cond.get("T_K")),
        "temperature_printed": cond.get("T_printed"),
        "temperature_basis": cond.get("T_basis"),
        "_T_measurement_K": num(cond.get("T_measurement_K")),   # row temperature (audit round 2; popped by the gates)
        "temperature_max_K": num(cond.get("T_max_K")),
        "temperature_range_K": cond.get("T_range_K"),
        "temperature_bound": cond.get("T_bound"),
        "temperature_approximate": cond.get("T_approximate"),
        "temperature_uncertainty": clean({"value": num((obs.get("temperature_uncertainty") or {}).get("value")),
                                          "basis": (obs.get("temperature_uncertainty") or {}).get("basis")}),
        "pressure_Pa": num(cond.get("P_Pa")),
        "pressure_basis": cond.get("P_basis"),
        "overall_composition": {k.upper(): num(v) for k, v in (cond.get("overall_composition") or {}).items()
                                if num(v) is not None} or None,
        "overall_composition_basis": cond.get("composition_basis"),
        "phases": phases_of(obs),
        "reaction": clean({"printed": (obs.get("reaction") or {}).get("printed"),
                           "type": (obs.get("reaction") or {}).get("type")}),
        "paper": paper,
        "locator": art,
        "readers": readers,
        "reader_tier": (obs.get("quality") or {}).get("status", {}).get("reader_tier") or obs.get("reader_tier"),
        "tool_release": obs.get("tool_release") or obs.get("_tool_release"),
    }
    return rec, exp, ev


def origin_block(obs, store, paper, measured_src):
    measured_key = paper_key(measured_src) if measured_src else None
    own = bool(measured_key and paper and measured_key == paper)
    return own, clean({
        "class": obs.get("origin"),
        "basis": obs.get("origin_basis"),
        "measured_in": measured_key,
        "measured_in_basis": obs.get("measured_in_basis"),
        "own_work_of_paper": own,
    })


def conditions_block(own, ev, paper, contexts):
    if own and ev == "experimental":
        if paper in contexts:
            return {"kind": "paper_context", "paper": paper}
        return {"kind": "paper_context_not_extracted", "paper": paper}
    if own:
        return {"kind": "not_applicable_own_calculation" if ev in ("model", "dft") else "own_work_origin_unspecified"}
    return {"kind": "cited_source", "note": "conditions in the cited source"}


def finish(rec, exp, obs, store, own, ev, contexts, origin):
    rec["origin"] = origin
    rec["conditions"] = conditions_block(own, ev, rec["paper"], contexts)
    rec["model"] = model_link(exp) if ev == "model" else None
    if ev == "model" and rec["model"] is None and exp.get("model_evidence_unlinked"):
        rec["model_unlinked_reason"] = exp.get("model_evidence_unlinked")
    return rec


def load_contexts(store):
    out = {}
    for f in sorted(glob.glob(os.path.join(store.db, "context", "*.jsonl"))):
        for r in jl(f):
            if r.get("record_type") != "context":
                continue
            sp = r.get("source_pid") or ""
            if sp.startswith("doi:") and lit_doi(sp[4:])[0]:
                key = "doi:" + doi_norm(sp[4:])
            else:
                key = paper_key(store.sources.get(r.get("source"))) or sp
            items = []
            for group, fields in (("sample", ("label", "nominal_composition", "measured_composition",
                                               "composition_basis", "starting_materials", "purity", "synthesis")),
                                  ("treatment", ("purpose", "applies_to", "temperature_printed", "temperature",
                                                 "temperature_status", "duration", "duration_parsed",
                                                 "atmosphere_or_container", "cooling")),
                                  ("measurement", ("method", "instrument", "applies_to", "calibration",
                                                   "heating_or_cooling_rate", "atmosphere", "crucible",
                                                   "stated_uncertainty"))):
                for it in r.get(group + "s") or []:
                    d = {"type": group}
                    for k in fields:
                        if k in it:
                            d[k] = it[k]
                    d.update({"quote": it.get("quote"), "page": it.get("page"), "id": it.get("pid_compact"),
                              "extraction_accepted": bool(it.get("extraction_accepted")),
                              "second_reader_agrees": bool(it.get("second_reader_agrees")),
                              "has_dropped_fields": bool(it.get("fields_not_in_quote") or it.get("second_reader_dropped_fields")) or None,
                              "review_flags": sorted({x.get("flag") for x in it.get("review_flags") or [] if x.get("flag")}) or None})
                    items.append(clean(d))
            if not items:
                continue
            uniq = {}
            for x in items:
                uniq.setdefault(dumps(x), x)
            store.notes["context_items_duplicate_dropped"] += len(items) - len(uniq)
            items = list(uniq.values())
            items.sort(key=lambda x: (x["type"], x.get("page") or 0, x.get("id") or ""))
            out[key] = {"paper": key, "context_id": r.get("id"), "section_basis": r.get("section_basis"),
                        "reader_tier": r.get("reader_tier"),
                        "readers": [x for x in ((r.get("reader") or {}).get("identity"), r.get("second_reader_identity")) if x],
                        "items": items}
            store.notes["context_items"] += len(items)
            store.notes["context_items_accepted"] += sum(1 for x in items if x.get("extraction_accepted"))
    return out


def excluded_file(path):
    b = os.path.basename(path).lower()
    if "rejected" in b and "wrong-match" in b:
        return "pdf_wrong_match"
    if "okamoto" in b or "handbook" in b:
        return "okamoto_out_of_scope"
    return None


def restricted(src):
    r = (src or {}).get("rights") or {}
    return r.get("status") == "restricted" or r.get("redistribute_extracted_values") in (False, "no", "denied")


def load_tables(store, sub, kind, contexts, records):
    files = sorted(glob.glob(os.path.join(store.db, sub, "*.jsonl")))
    for f in files:
        why = excluded_file(f)
        recs = list(jl(f))
        obs_list = [r for r in recs if r.get("record_type") == "observation"]
        if why:
            store.excluded[why] += len(obs_list)
            continue
        arts = {}
        corpus_src = None
        local = {}
        for r in recs:
            if r.get("record_type") == "source":
                store.add_source(r)
                for k in [r["id"]] + list(r.get("aliases") or []):
                    local.setdefault(k, r)
                if r.get("kind") in ("journal_article", "book_chapter") and corpus_src is None:
                    corpus_src = r
            elif r.get("record_type") == "artifact":
                arts[r["id"]] = r
        for obs in obs_list:
            m = re.match(r"obs:(artifact:.+?:p\d+:[^:]+)(?::row-(\d+))?", obs["id"])
            art = arts.get(m.group(1)) if m else None
            if art is None and len(arts) == 1:
                art = next(iter(arts.values()))
            src = store.source_in_file((art or {}).get("source"), local) or corpus_src
            if restricted(src):
                store.excluded["restricted_rights"] += 1
                continue
            if not obs.get("system"):
                store.excluded["no_chemical_system"] += 1
                continue
            paper = store.paper_for(src, "reporting_paper") if src else None
            row = int(m.group(2)) if m and m.group(2) else None
            loc = clean({"artifact": "table", "label": (art or {}).get("label"), "caption": (art or {}).get("caption"),
                         "page": (art or {}).get("page"), "row": row + 1 if row is not None else None,
                         "asset": obs.get("asset")})
            readers = sorted({reader_of(x) for x in obs.get("readings") or []})
            rec, exp, ev = base_record(obs, kind, store, paper, loc, readers)
            ot = obs.get("observation_type")
            basis = comp_basis(exp.get("conversions"))
            rec["phase_compositions"] = phase_comps(obs, basis)
            if rec["overall_composition"] and any(v > 1.0001 for v in rec["overall_composition"].values()):
                rec["overall_composition_basis"] = "unresolved"
            # value
            if ot == "invariant_reaction":
                rec["property"], rec["value"], rec["unit"] = "INVARIANT_TEMPERATURE", rec["temperature_K"], "K"
            elif ot == "phase_transition" and obs.get("value") is None:
                rec["property"], rec["value"], rec["unit"] = "TRANSITION_TEMPERATURE", rec["temperature_K"], "K"
            elif ot in ("phase_region_vertex_set", "phase_assemblage_node") and obs.get("value") is None:
                rec["property"] = "PHASE_COMPOSITION" if ot == "phase_region_vertex_set" else "PHASE_ASSEMBLAGE"
            else:
                rec["property"] = obs.get("quantity") or "OTHER"
                rec["value"] = num(obs.get("value"))
                rec["unit"] = obs.get("unit")
                pr = exp.get("printed") or {}
                rec["property_printed"] = obs.get("printed_quantity") or pr.get("quantity")
                rec["value_printed"] = pr.get("value") if not isinstance(pr.get("value"), (dict, list)) else None
                rec["unit_printed"] = pr.get("unit")
                rec["unit_basis"] = quality.unit_basis_from_conversions(exp.get("conversions"))
                div = [c.split(":/")[1] for c in exp.get("conversions") or [] if c.startswith("per_formula_unit_to_per_atom:/")]
                if div:
                    try:
                        rec["per_formula_unit_divisor"] = float(div[0])
                    except ValueError:
                        pass
                u = obs.get("uncertainty") or {}
                rec["uncertainty"] = clean({"value": num(u.get("value")), "basis": u.get("basis"), "note": u.get("note")})
                if obs.get("value_range"):
                    vr = obs["value_range"]
                    rec["value_range"] = clean({"min": num(vr.get("min")), "max": num(vr.get("max"))})
                if obs.get("value_expression"):
                    rec["value_expression"] = (obs["value_expression"] or {}).get("printed")
                rs = obs.get("reference_state")
                if rs:
                    rec["reference_state"] = rs
                if obs.get("component"):
                    rec["component"] = obs["component"]
                if obs.get("structure"):
                    rec["structure"] = obs["structure"]
                if obs.get("calculation_settings"):
                    rec["calculation_settings"] = obs["calculation_settings"]
            rec["property_group"] = property_group(ot, rec.get("property"))
            meth = obs.get("method") or ""
            rec["method"] = None if meth in GENERIC_METHODS else meth
            rec["method_basis"] = "printed in the table row/column" if rec["method"] else None
            measured_src = store.source_in_file(obs.get("measured_in"), local)
            if measured_src is None and obs.get("measured_in") and src and obs["measured_in"] in ([src["id"]] + list(src.get("aliases") or [])):
                measured_src = src
            own, origin = origin_block(obs, store, paper, measured_src)
            rec["cited_sources"] = cited_list(exp.get("cited"), store)
            if not own and measured_src is not None:
                store.paper_for(measured_src, "cited_source")
            rec["quote"] = None
            rec["quote_page"] = None
            fl = quality_flags(obs)
            fl["unsigned_in_print"] = None
            fl["sign_pdf_page_verified"] = None
            rec["reader_agreement"], rec["reader_agreement_label"] = reader_agreement(exp)
            rec["flags"] = clean(fl)
            rec["flags"]["extraction_accepted"] = bool((obs.get("quality") or {}).get("status", {}).get("extraction_accepted"))
            rec = finish(rec, exp, obs, store, own, ev, contexts, origin)
            records.append(rec)


def load_text(store, contexts, records):
    for f in sorted(glob.glob(os.path.join(store.db, "text", "*.jsonl"))):
        why = excluded_file(f)
        recs = list(jl(f))
        if why:
            store.excluded[why] += sum(1 for r in recs if r.get("record_type") == "observation")
            continue
        arts, reads, local = {}, {}, {}
        for r in recs:
            t = r.get("record_type")
            if t == "source":
                store.add_source(r)
                for k in [r["id"]] + list(r.get("aliases") or []):
                    local.setdefault(k, r)
            elif t == "artifact":
                arts[r["id"]] = r
            elif t == "reading":
                reads[r["id"]] = r
        for obs in (r for r in recs if r.get("record_type") == "observation"):
            att = obs.get("attribution") or {}
            src = store.source_in_file(att.get("reported_in"), local) or \
                store.source_in_file((arts.get(obs.get("artifact")) or {}).get("source"), local)
            if restricted(src):
                store.excluded["restricted_rights"] += 1
                continue
            if not obs.get("system"):
                store.excluded["no_chemical_system"] += 1
                continue
            paper = store.paper_for(src, "reporting_paper") if src else None
            art = arts.get(obs.get("artifact")) or {}
            rd = reads.get(obs.get("locator_reading")) or {}
            pr = rd.get("printed") or {}
            quote = pr.get("quote")
            if not quote:
                store.excluded["text_without_quote"] += 1
                continue
            loc = clean({"artifact": "text_passage", "label": art.get("label"), "page": art.get("page"),
                         "asset": obs.get("asset")})
            readers = sorted({(reads.get(x) or {}).get("reader", {}).get("identity", reader_of(x)).split(" / ")[0]
                              for x in obs.get("readings") or []})
            rec, exp, ev = base_record(obs, "running_text", store, paper, loc, readers)
            rt = obs.get("reported_temperature") or {}
            if rt.get("printed") or rt.get("value_printed") is not None:
                rec["temperature_printed"] = rt.get("printed") or (f"{rt.get('value_printed')} {rt.get('unit_printed') or ''}".strip())
            rec["temperature_bound"] = rt.get("bound")
            rec["temperature_max_K"] = num(rt.get("range_max"))
            if obs.get("composition_basis"):
                rec["overall_composition_basis"] = "mole_fraction (printed basis: " + ", ".join(obs["composition_basis"]) + ")"
            rec["phase_compositions"] = []
            ot = obs.get("observation_type")
            if ot == "phase_transition":
                rec["property"], rec["value"], rec["unit"] = "TRANSITION_TEMPERATURE", rec["temperature_K"], "K"
                if (obs.get("reaction") or {}).get("type") in ("eutectic", "peritectic", "monotectic", "eutectoid",
                                                              "peritectoid", "syntectic", "metatectic", "catatectic"):
                    rec["property"] = "INVARIANT_TEMPERATURE"
            elif ot == "phase_assemblage_node":
                rec["property"] = "PHASE_ASSEMBLAGE"
            else:
                q = obs.get("quantity") or "other"
                rec["property"] = TEXT_QUANTITY.get(q, q.upper())
                rec["property_printed"] = obs.get("quantity_printed")
                rec["value"] = num(obs.get("value"))
                rec["unit"] = obs.get("unit")
                rec["reference_state"] = obs.get("reference_state_printed")
                u = obs.get("uncertainty") or {}
                rec["uncertainty"] = clean({"value": num(u.get("value")), "basis": u.get("basis"), "note": u.get("note")})
            rec["property_group"] = property_group(ot if ot != "equilibrium_property" else "x", rec.get("property")) \
                if ot != "phase_transition" else ("invariant reaction" if rec["property"] == "INVARIANT_TEMPERATURE" else "phase transition")
            meth = obs.get("method") or ""
            rec["method"] = None if meth in GENERIC_METHODS else meth
            rec["method_basis"] = "printed in the quoted sentence (fix-up 2026-10-22 gate)" if rec["method"] else None
            measured_src = store.source_in_file(obs.get("measured_in"), local)
            if measured_src is None and src and obs.get("measured_in") == att.get("reported_in"):
                measured_src = src
            own, origin = origin_block(obs, store, paper, measured_src)
            if not own and obs.get("measured_in") and obs.get("measured_in") == att.get("reported_in"):
                own = True
                origin["own_work_of_paper"] = True
                origin["measured_in"] = paper
            if att.get("note"):
                origin["attribution_note"] = att["note"]
            if not own and measured_src is not None:
                store.paper_for(measured_src, "cited_source")
            rec["link_rule"] = obs.get("link_rule")
            rec["cited_sources"] = cited_list(att.get("cited"), store)
            rec["citations_printed"] = obs.get("citations_printed") or None
            rec["quote"] = quote
            rec["quote_page"] = pr.get("page") or art.get("page")
            fl = quality_flags(obs)
            fl["extraction_accepted"] = bool((obs.get("quality") or {}).get("status", {}).get("extraction_accepted"))
            fl["value_in_quote"] = obs.get("value_in_quote")
            fl["unsigned_in_print"] = obs.get("sign_as_printed") == "unsigned_in_print" or None
            fl["sign_pdf_page_verified"] = obs.get("value_source") == "pdf_page_verified" or None
            if fl["sign_pdf_page_verified"]:
                fl["value_first_reader"] = num(obs.get("value_first_reader"))
            fl["has_dropped_fields"] = bool(obs.get("dropped_fields")) or None
            fl["dropped_field_count"] = len(obs.get("dropped_fields") or []) or None
            fl["review_flags"] = sorted({x.get("flag") for x in obs.get("review_flags") or [] if x.get("flag")}) or None
            fl["second_reader"] = obs.get("second_reader")
            fl["unit_conflict_between_readers"] = bool(obs.get("unit_conflict_reader")) or None
            fl["composition_conflict_between_readers"] = bool(obs.get("composition_conflict_reader")) or None
            rec["flags"] = clean(fl)
            rec["flags"]["extraction_accepted"] = fl["extraction_accepted"]
            rec = finish(rec, exp, obs, store, own, ev, contexts, origin)
            records.append(rec)


# ----------------------------------------------------------------------------------------------- output

CSV_COLUMNS = ["id", "system", "record_kind", "observation_type", "property_group", "property", "property_code_raw",
               "property_printed", "value", "unit", "unit_basis", "per_formula_unit_divisor", "value_before_quality_fix",
               "basis_correction", "sign_heading_printed", "value_stored", "unit_stored", "value_printed", "unit_printed", "uncertainty", "value_range", "value_expression",
               "temperature_K", "temperature_printed", "temperature_basis", "temperature_bound", "temperature_max_K",
               "temperature_uncertainty", "uncertainty_budget", "sigma_effective", "sigma_method_class", "pressure_Pa", "pressure_basis", "phases", "phase_compositions",
               "overall_composition", "overall_composition_basis", "reference_state", "reaction", "transition_kind",
               "transition_kind_printed", "evidence_class",
               "origin", "model", "model_unlinked_reason", "method", "method_basis", "conditions",
               "conditions_summary", "cited_sources", "citations_printed", "paper", "paper_citation", "paper_doi", "quote", "quote_page",
               "locator", "link_rule", "flags", "extraction_accepted", "espei_ready", "espei_certified", "calphad_check", "readers", "reader_tier", "reader_agreement",
               "reader_agreement_label", "tool_release", "pid",
               "merged_ids", "merged_pids", "corrections", "gate_failures", "reader_agreement_basis",
               "temperature_K_stored", "system_stored", "property_code_before_gates", "evidence_class_before_gates",
               "phases_stored", "phases_dropped", "unit_stored", "pressure_basis_stored", "phase_compositions_dropped",
               "transition_kind_stored", "uncertainty_stored", "reference_state_stored", "related_records",
               "temperature_basis_stored", "chemical_system", "paper_id", "fold_lso", "fold_lpo", "source_system",
               "overall_composition_stored"]


def summarize_context(ctx):
    if not ctx:
        return None
    parts = []
    for it in ctx["items"]:
        if not it.get("extraction_accepted"):
            continue
        if it["type"] == "sample":
            s = "; ".join(f"{k} {it[k]}" for k in ("purity", "starting_materials", "synthesis") if it.get(k))
        elif it["type"] == "treatment":
            s = "; ".join(str(it[k]) for k in ("temperature_printed", "duration", "atmosphere_or_container", "cooling") if it.get(k))
        else:
            s = "; ".join(str(it[k]) for k in ("method", "instrument", "heating_or_cooling_rate", "atmosphere") if it.get(k))
        if s:
            parts.append(f"{it['type']} (p{it.get('page')}): {s}")
    return " | ".join(parts) or None


def csv_value(v):
    if v is None:
        return ""
    if isinstance(v, (dict, list)):
        return json.dumps(v, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    if isinstance(v, bool):
        return "true" if v else "false"
    return str(v)


def dumps(obj, compact=True):
    if compact:
        return json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return json.dumps(obj, ensure_ascii=False, sort_keys=True, indent=1)


def gz(data: bytes) -> bytes:
    buf = io.BytesIO()
    with gzip.GzipFile(filename="", mode="wb", fileobj=buf, mtime=0, compresslevel=9) as g:
        g.write(data)
    return buf.getvalue()


SCHEMA_MD = """# ORCHID literature data: schema

Release: `{release}`. Generated by `tools/export_literature_data.py` (deterministic; `--check` verifies the files).

Derived data extracted from published papers (tables and running text) by independent machine readers
(Gemini Flash / Pro, Claude Opus; text-layer reading for tables) and accepted only on reader agreement plus
deterministic checks. **Licence:** derived data CC BY 4.0 with attribution to ORCHID and to the cited source
papers; each `quote` is a short excerpt kept to make the value verifiable, which remains under the publisher's
copyright and is not relicensed. The source
papers, page images, crops and full text are **not** redistributed.

**Status.** {status}

## Files

| file | content |
|---|---|
| `index.json` | release info, totals, per-system counts by evidence class and property group, file list |
| `systems/<SYS>.json` | `{{system, elements, records[], artifacts{{}}, models{{}}, papers{{}}, contexts{{}}}}` for one chemical system. To keep files small, a record's `locator` is `{{ref, row}}` with `ref` a key of `artifacts` (table/passage label, caption, page) and its `model` is a key of `models`; `elements` is stored once per file. The `.jsonl.gz` dump carries the full records; `<SYS>` = elements sorted alphabetically, uppercase, joined by `-` (as in `data/phase-diagrams`). Large systems are split into `<SYS>.json`, `<SYS>.part2.json`, ... (listed in `index.json`) |
| `papers.json` | paper id -> bibliographic record |
| `contexts.json` | paper id -> experimental-context items (samples, treatments, measurements), each with its quote and page |
| `orchid-literature.jsonl.gz` | all records, one JSON object per line (same fields as `systems/*.json`) |
| `orchid-literature.csv.gz` | all records flattened; nested fields are compact JSON strings |
| `CITATION.txt`, `CITATION.cff` | how to cite this dataset (the CSV carries no metadata; a leading comment line would break CSV readers). `index.json` and every `systems/*.json` also carry a top-level `cite_as` object: text citation, BibTeX and Zenodo DOI of this dataset and of the umbrella record, and the sources to cite as well |
| `QUALITY_FIXES.md`, `quality_fixes.json` | data-hardening report: duplicates merged, vocabulary maps, garbled fields, DFT-parity outlier classes (every changed record with before / after / evidence) |

## Record fields

| field | meaning |
|---|---|
| `id` | stable compact ORCHID id (`ORC-...`) |
| `pid` | persistent locator id: paper / page / table or passage / row |
| `chemical_system`, `paper_id`, `fold_lso`, `fold_lpo` | machine-learning folds: `chemical_system` (sorted elements, = `system`), `paper_id` (stable source id), `fold_lso` = leave-chemical-system-out fold int(sha256(chemical_system), 16) mod 5, `fold_lpo` = leave-paper-out fold int(sha256(paper_id), 16) mod 5 (tools/ml_folds.py). Subsystems are their own systems (a held-out ternary may have its binaries in training); drop subsystems of held-out systems for strict chemistry extrapolation |
| `system`, `elements` | chemical system of the record's own state (sorted uppercase elements joined by `-`): a binary compound's formation enthalpy or lattice parameter, or a pure element's transition, printed in a ternary / multi-component paper is filed under its binary / unary subsystem (audit round 2, A1 C4: `corrections` rule `system_of_the_state`, flag `system_is_state_subsystem`) |
| `source_system` | the chemical system of the source table / paper the record was read from (equal to `system` unless the state is a lower-order subsystem) |
| `record_kind` | `invariant_table` (invariant-reaction table), `general_table` (thermochemical / phase-equilibrium / structure table), `running_text` (value stated in a sentence) |
| `observation_type` | `invariant_reaction`, `phase_transition`, `equilibrium_property`, `structural_property`, `phase_region_vertex_set`, `phase_assemblage_node`, `property_expression` |
| `property` | quantity kind: `INVARIANT_TEMPERATURE`, `TRANSITION_TEMPERATURE`, `PHASE_COMPOSITION`, `PHASE_ASSEMBLAGE`, or a thermodynamic / structural quantity (`HM_FORM`, `HM_MIX`, `CPM`, `SM`, `GIBBS_ENERGY`, `EMF`, `ACR_<EL>`, `LATTICE_A`, ...); a controlled code (see Controlled vocabularies). `property_code_raw`: the code as the record builder wrote it, when it was malformed and cleaned |
| `property_group` | coarse filter group (invariant reaction, phase transition, phase equilibrium, formation, mixing / partial, heat capacity / entropy / enthalpy, activity / EMF / Gibbs energy, crystal structure, model expression, other property) |
| `property_printed` | the quantity symbol or name as printed |
| `value`, `unit` | value in the normalized unit (temperatures: K; formation enthalpies and Gibbs energies of formation, mixing enthalpies, Gibbs energies of mixing and partial molar enthalpies: kJ/mol-atom when the per-atom basis is explicit, otherwise kJ/mol (basis unstated); lattice lengths a, b, c (`LATTICE_A`, `LATTICE_B`, `LATTICE_C`, `LATTICE_A0`, `LATTICE_C0`, `LATTICE_UNSPECIFIED_AXIS`, `LATTICE_A_R`): angstrom (Å) when the stored unit is nm, pm or Å, unconverted when no unit was printed (`not_printed`, `unknown`), flag `lattice_length_outside_typical_range` outside 1-100 Å; other properties (e.g. Gibbs energies, which include the coefficients B of printed expressions A + BT in J/mol/K) keep the unit as read). `value_stored`, `unit_stored`: the value before that conversion. `value_printed`, `unit_printed`: as printed (tables) |
| `unit_basis` | how the per-atom / per-formula basis of an energy was established: `printed_unit` (in the printed unit), `column_header`, `caption`, `reviewer`, `state_completion_<where>` (a quoted header / caption / footnote / sentence), `eV_per_atom_printed`, `magnitude_heuristic` (inferred by the record builder from the magnitude; formation values flagged `basis_inferred_from_magnitude`), `assumed_per_mole_of_atoms` (solution phase, not printed), `not_printed`, `unstated`. For mixing / partial quantities only an explicit basis gives kJ/mol-atom. `per_formula_unit_divisor`: the record builder divided a per-formula-unit value by this number of atoms |
| `value_before_quality_fix`, `basis_correction`, `sign_heading_printed` | set when a DFT-parity outlier was corrected (QUALITY_FIXES.md): the value before the fix; `{{formula, atoms_per_formula_unit, operation}}` for a per-formula value that was labelled per atom; the printed heading that carries a minus sign |
| `corrections` | literature-level gate corrections (QUALITY_FIXES.md section 7): `[{{gate, field, before, after, rule, evidence, visual_check}}]`; `visual_check` = `{{page, pdf_sha256_prefix, image_sha256, verdict, seen, checked}}` (the PDF page was rendered and inspected; no image or page text is published). The value before a correction is kept in `value_stored`, `temperature_K_stored`, `unit_stored`, `system_stored`, `property_code_before_gates`, `evidence_class_before_gates`, `phases_stored` (with `phases_dropped` [{{label, reason}}]), `pressure_basis_stored`, `overall_composition_stored`, and per phase composition `composition_stored` / `basis_stored` (A1 audit round 2: fractions under an at.% header, header naming the wrong element, compositions assigned to the wrong phase; gate `A1`) |
| `gate_failures` | `[{{gate, code, effect, detail, visual_check}}]`: a gate demoted (`effect` `not_accepted`: `flags.extraction_accepted` false) or flagged (`not_espei_ready`: `espei_ready` and `espei_certified` false; `flag_only`) the value; codes also in `flags.gate_failures`; see QUALITY_FIXES.md section 7 |
| `phase_compositions_dropped` | compositions printed under a label that the page check found not to be a phase of the row (a column header); kept, not used |
| `transition_kind_stored` | the reaction type before a page-checked correction to the printed type |
| `uncertainty_stored`, `reference_state_stored`, `temperature_basis_stored` | the stored uncertainty / reference state before a page-checked correction (e.g. a scatter-aware sigma for compounds whose independent measurements disagree, `uncertainty.basis = scatter_between_sources`; a printed sigma combined with the composition offset of a value measured at a phase's historical formula, `uncertainty.basis = printed_with_composition_offset`; the reference state the page states; a value converted to the standard reference state carries the conversion in `corrections[].evidence.conversion`; primary works of `cited_sources` resolved from the paper's reference list keep the unresolved entries in `corrections[].before`) |
| `property_expression` | (audit round 2) `{{form, term, printed, temperature_range_K}}` when the value is one coefficient of a printed expression in T (`A + B T`, a Cp polynomial); the property code then ends in `_EXPRESSION_TERM` and the record has no single temperature (`temperature_K_stored` keeps the one stored before) |
| `calorimeter_temperature_K`, `reference_temperature_K` | (audit round 2) the calorimeter temperature of a standard enthalpy of formation stored at 298.15 K (direct-synthesis / solute-solvent drop calorimetry); the reference temperature of a heat content H(T) - H(298) stored at its row temperature T. Property codes `HM_FROM_COMPOUNDS` / `GM_FROM_COMPOUNDS` / `SM_FROM_COMPOUNDS`: formation quantities whose printed reference is compounds (e.g. oxides), not the elements; `flags.non_crystalline_state`: liquid / amorphous / glassy state printed at room temperature |
| `related_records` | ids of records of the same measurement listed in another paper (a compilation), or of the independent measurements behind a scatter-aware sigma |
| `merged_ids`, `merged_pids` | ids (and persistent locators) of records that were exact duplicates of this one and were merged into it |
| `value_range` | `{{min,max}}` when a range was printed; `value_expression`: a printed expression (model parameters) |
| `uncertainty` | `{{value, basis, note}}`; `basis: printed` is the source's own; `stated_default` / `reading_resolution` are ORCHID placeholders (e.g. "5% default", "10% placeholder; none printed") and are **not** measured uncertainties |
| `uncertainty_budget` | uncertainty budget of the record's measured quantities, components kept separate (PI decision 2026-10-24; tools/uncertainty_budget.py; calibration in UNCERTAINTY_CALIBRATION.json / .md): `{{calibration (version and sha256 of UNCERTAINTY_CALIBRATION.json), sigma_rule, quantities{{value, composition, temperature_K}}}}`; `value` = the record's value, `composition` = its printed phase compositions (mole fraction), `temperature_K` = the equilibration temperature of a phase-composition / assemblage observation. Per quantity, in its `unit`: `sigma_printed` `{{value, kind, taken_as}}` as printed (null if none; ORCHID placeholders are not printed sigmas), `sigma_reading` `{{value, basis}}` half of the last printed digit (`last_printed_digit`; compositions: `last_stored_digit`), `sigma_method` `{{class, value, floor, floor_unit, basis (calibrated / prior / calibrated vs experiment), n_pairs, ci95, entry}}` the method-class floor (robust between-study scatter of independent same-method measurements in ORCHID, or a documented prior; DFT: calibrated sigma_DFT; model values: model-vs-experiment deviation), `sigma_composition` (value measured at a historical formula: `{{value, dx, slope}}`), `sigma_scatter` (page-checked mutually inconsistent group: `{{value, n_sources, group}}`), `sigma_conversion` (reference-state conversion), then `sigma_effective` = max(sqrt(max(printed, reading, method)^2 + composition^2 + conversion^2), scatter) (`sigma_effective_reason` when none applies). Null = the component does not apply. CSV: the budget as JSON plus `sigma_effective` and `sigma_method_class` of the value |
| `temperature_K` | equilibrium / measurement temperature in K; `temperature_printed` as printed (unit as printed, often °C); `temperature_basis` how it was obtained (`printed`, `quantity_label`, `caption`, `assumed_room_temperature`, `not_printed`, ...); `temperature_bound` `{{op, value_K, printed}}` for a printed bound (`op` `<` upper bound, `>` lower bound), `temperature_max_K` / `temperature_range_K` for ranges; `temperature_uncertainty` `{{value,basis}}` |
| `pressure_Pa`, `pressure_basis` | pressure and how it was set: `printed`, or `assumed_default_1_atm_not_printed` (1 atm assumed, no pressure printed; the record builder's label kept in `pressure_basis_stored`) |
| `phases` | `[{{label, model_phase, mapping_confidence, pearson_symbol, prototype, label_variants, label_source}}]`: printed phase label and, where mapped, the CALPHAD model phase; `label_variants`: other readers' spellings of the same printed label (kept when duplicate readings were merged, e.g. 'NiAl (B2)' beside 'B2'); `label_source: page_check` for a phase added from the printed reaction on a page check; `structure_basis: mapped_model_phase`: the Pearson symbol / prototype are those of the mapped model phase, not printed; `model_phase_rejected` {{model_phase, reason, pearson_symbol, prototype}}: a corpus mapping the label cannot carry (tools/phase_label_check.py); `model_phase_basis: system_alias:<id>`: the model phase comes from (or is confirmed by) the per-system alias table of conventional phase names (tools/phase_aliases.json, e.g. Cu-Sn ε -> CU3SN_L, Cu-Zn β -> BCC_A2, Cu-Si κ -> HCP_A3), applied when the printed composition / temperature agree (`model_phase_alias_check`: composition, temperature, composition+temperature or label_only); `structure_basis: system_alias`: Pearson symbol / prototype are the handbook structure of that alias |
| `phase_compositions` | `[{{phase, composition{{EL: x}}, basis, composition_source}}]`. `basis` is `mole_fraction` (possibly "converted from atomic/mass percent") or `unresolved` (printed basis disputed/mixed, or values outside 0..1: do not use as mole fractions). `composition_source`: `printed`, `stoichiometry`, `formula_nominal`, `congruent_same_composition`, `printed_range`, `printed_bound`, ... |
| `overall_composition`, `overall_composition_basis` | alloy / sample composition (mole fraction) where given |
| `reference_state` | reference state of a thermodynamic quantity, always an object: per element `{{EL: {{phase, fixed_state_variables}}}}`, or `{{printed: "..."}}` when only free text was read |
| `reaction` | `{{printed, type}}`: reaction as printed and its class (`type` = `transition_kind`) |
| `transition_kind`, `transition_kind_printed` | controlled reaction / transition class (see Controlled vocabularies) and the class as printed or read |
| `evidence_class` | `experimental`, `model` (value from a CALPHAD/thermodynamic model), `dft` (first-principles), `unspecified` (origin could not be settled). Never mixed or averaged |
| `origin` | `{{class, basis, measured_in, measured_in_basis, own_work_of_paper, attribution_note}}`: how the evidence class was decided (e.g. `printed_source_experimental`, `cited_work_before_1980_assumed_experimental`), the paper the value was measured/calculated in, and whether that is the containing paper itself |
| `model` | for `model` values linked unambiguously to a TDB: `{{tdb, tdb_path, tdb_sha256, catalog_entry, link_basis, model_citation, model_doi, reproduced_by_calculation}}`; `model_unlinked_reason` otherwise. `model_doi` is null with `model_doi_status: dropped_mismatch` and `model_doi_correction` when the TDB header citation's DOI is a mismatch (`tools/doi_corrections.json`) |
| `method` | row-level method only when printed (table cell/column, or inside the value's quoted sentence); `method_basis` says which |
| `conditions` | `{{kind: paper_context, paper}}`: the containing paper's own experimental context applies (see `contexts`); `paper_context_not_extracted`; `cited_source` ("conditions in the cited source"): the value is cited from another work and the containing paper's context does **not** apply; `not_applicable_own_calculation` |
| `conditions_summary` | (CSV only) accepted context items of that paper, condensed |
| `cited_sources` | the value's own source(s) as printed in the paper's reference list: `{{label, reference, doi, first_author, year, kind, status, paper}}`; `year` is re-derived from the printed `reference` (a parenthesised year first, never a page number; tools/literature_quality.py `cited_year`), and where this changed the reader's year that year is kept in `year_read` with the rule in `year_basis` (`printed_parenthesised_year`, `printed_year_not_page`, `page_number_not_year`, `ambiguous_several_years`, `after_citing_paper`; a null year means no year could be read safely); a DOI dropped as a mismatch is null with `doi_status: dropped_mismatch` and `doi_correction` (the printed `reference` is kept) |
| `citations_printed` | (text) citation labels printed in the sentence |
| `paper` | id of the containing paper (`doi:<doi>` or `sha256:<pdf prefix>`) -> `papers.json`; CSV adds `paper_citation`, `paper_doi` |
| `quote`, `quote_page` | (running text) the short verbatim excerpt backing the value (usually the sentence that states it; some excerpts span several sentences) and its page (`null` with `flags.page_unreadable` when the page read was not a page number); for tables `locator` gives table label, caption, page and row instead |
| `locator` | `{{artifact: table|text_passage, label, caption, page, row, asset}}` |
| `flags` | `extraction_accepted` (reader agreement + checks passed; the explorer hides `false` by default), `row_state` (`settled`, `awaiting_reader`, `disputed_in_review_queue`), `review_state`, `espei_ready` (the value CAN be used in ESPEI: complete state, phases mapped to a model, phase-rule / unit checks, proper format, and the physical phase rule of QUALITY_FIXES.md section 6; `espei_ready_upstream` true and `espei_phase_rule` {{result, reason, detail}} when only the phase rule turned it false), `espei_certified` (`espei_ready` AND agreement with a CALPHAD prediction: `plausibility.result` = pass), `espei_certified_basis` (`pycalphad_invariant_agreement`; later bases may be added, e.g. `espei_driving_force`), `calphad_check` (`agrees`, `disagrees`, `not_available`), `state_fully_specified`, `specification_gaps`, `plausibility` (CALPHAD cross-check: result, model, dT_K, dx {{phase: mole-fraction difference}}, calculation, tolerance {{dT_K, dx}}, `fail_reason` for fails with |dT| < 5 K), `plausibility_before_correction` (the cross-check of a value whose temperature a gate corrected; computed on the old temperature), `source_disputed`, `judge_adjudicated`, `magnitude_suspect`, `system_mismatch`, `value_in_quote`, `unsigned_in_print` (minus sign absent in print, value as printed), `sign_pdf_page_verified` (sign corrected from the PDF page image; `value_first_reader` kept), `has_dropped_fields` / `dropped_field_count` (fields removed because they were not printed in the quote; content not published), `review_flags` (codes only), `second_reader`, reader conflicts; data-hardening flags (QUALITY_FIXES.md): `page_unreadable`, `reference_state_dropped_not_in_quote`, `basis_inferred_from_magnitude`, `basis_corrected_from_formula`, `sign_from_heading`, `dft_outlier_unexplained` |
| `link_rule` | (text) the rule that linked the value to its quote: the quote is verbatim on one page/passage and the printed value is inside it |
| `reader_agreement`, `reader_agreement_label` | (tables) how many independent visual readers carry the value: `multi_reader` ("confirmed (two or more readers)"), `single_reader` ("unconfirmed (single reader)": one visual reader read the row and no second reader has confirmed it; the explorer shows a "single reader" badge), `single_reader_text_layer` ("single reader + PDF text layer": one visual reader, every value also printed in the PDF text layer). Two blind, independent Claude reads (`claude-opus-5.5` and `claude-opus-5.5#2` in `readers`) count as two readers. Where the record builder did not classify the value (invariant tables, running text) it is derived from `readers` (`reader_agreement_basis: derived_from_readers`, settled rows only) |
| `readers`, `reader_tier`, `tool_release` | extraction provenance: the readers whose readings agreed (`text_layer` = the PDF text layer), the reader tier, a pipeline status of the table or passage (which reader passes ran), not an agreement claim (`preliminary_flash+claude`: two independent model readers agree; `confirmed_pro`: confirmed by a third reader; `+fable_judge`: a judge chose between disagreeing readers; `preliminary_flash` for text: first reader, with the second reader's agreement in `flags.second_reader`), the record-builder release |

## Controlled vocabularies

`transition_kind`: {kinds}. Mapping (deterministic, `tools/literature_quality.py`): typos repaired (Peritetic,
Peritectiod, ...); keywords (eutectoid before eutectic, magnetic before order-disorder, metatectic and catatectic
together); Rhines letters (E/e eutectic, P/p peritectic, U and class II quasi-peritectic; ternary class I eutectic);
"C.M.", "melting point", "congruent melting" congruent melting; a bare "congruent" (point, reaction) is congruent
melting when a liquid phase is among the record's phases and congruent transformation otherwise; "polymorphic",
"allotropic", "solid state transition" polymorphic; printed alternatives ("eutectic (or peritectic)"), liquidus /
solidus points, extrema (max, min, saddle), "C" and anything else are `other`. Every printed variant, its term and
the rule are listed in QUALITY_FIXES.md.

Property codes: malformed builder codes are cleaned ({propmap}; `LATTICE__` -> `LATTICE_UNSPECIFIED_ANGLE` for
angles, `LATTICE_UNSPECIFIED_AXIS` otherwise; other `LATTICE_LATTICE_PARAMETER_*` and trailing / doubled
underscores stripped). The raw code is kept in `property_code_raw`.

## Data hardening (`QUALITY_FIXES.md`, `quality_fixes.json`)

Duplicates merged (`merged_ids`), vocabularies, garbled pages / bounds / reference states, DFT-parity outlier
classes and energy-unit conversions, with every changed record (id, before, after, evidence) for pages,
reference states and outliers, and the complete variant tables for the vocabularies.

## Paper fields (`papers.json`)

Ids are `doi:<doi>` (corpus papers and cited works with a DOI), `sha256:<16 hex>` (corpus papers without DOI),
`ref:<hash>` (a cited work without DOI, citation = the printed reference), `nodoi:<hash>` (a cited work whose matched
DOI was dropped as a mismatch, same id as in the bibliography) or `unresolved:...` (a reference label the
citing paper's list could not resolve). Fields: `id`, `doi`, `doi_status` / `doi_correction` (`dropped_mismatch` and
the entry id of `tools/doi_corrections.json` when the DOI an upstream matcher assigned was dropped because its Crossref
record describes a different work; `doi` is then null and the citation is the printed one), `doi_url`, `citation` (formatted from local Crossref metadata, the TDB header, or the printed reference),
`title`, `authors`, `year`, `container`, `kind` (`corpus_paper`: a paper the values were read from; `cited_work`:
a work cited for a value), `roles` (`reporting_paper`, `cited_source`), `systems`, `metadata_source`,
`pdf_sha256_prefix` (first 16 hex of the source PDF's SHA-256, an identifier only; the PDF is not distributed).

## Context items (`contexts.json`, `systems/*.json` -> `contexts`)

`{{paper, context_id, section_basis, reader_tier, readers, items[]}}`; each item has `type` (`sample`, `treatment`,
`measurement`), its fields (`purity`, `starting_materials`, `synthesis`, `nominal_composition`; `temperature_printed`,
`duration`, `atmosphere_or_container`, `cooling`, `purpose`; `method`, `instrument`, `heating_or_cooling_rate`,
`atmosphere`, `crucible`, `calibration`, `stated_uncertainty`, `applies_to`), its `quote` and `page`, and
`extraction_accepted`, `second_reader_agrees`, `has_dropped_fields`, `review_flags`. A context is attached only to values
whose origin is the paper's own experimental work.
"""


def apply_contributions(records, data):
    """Accepted contributed corrections (orchid-data contributions/corrections/*.jsonl; contributions/SCHEMA.md).

    A correction is applied only when the maintainers' page check `contributed_correction:<record id>:<field>` in
    tools/literature_visual_checks.json has verdict `confirmed` and its `old` value still matches; the value before
    is kept in the record's `corrections[]` with the contributor, who is added to the record's `contributors`.
    Returns the applied check keys (empty: the release is unchanged).
    """
    tool = os.path.join(data, "contributions", "tools", "orchid_contrib.py")
    files = sorted(glob.glob(os.path.join(data, "contributions", "corrections", "*.jsonl")))
    if not os.path.exists(tool) or not files:
        return []
    spec = importlib.util.spec_from_file_location("orchid_contrib", tool)
    oc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(oc)
    corrections = [json.loads(line) for f in files for line in open(f, encoding="utf-8") if line.strip()]
    by_id = {r["id"]: r for r in records}
    return oc.apply_corrections(by_id, corrections, gates.load_visual())


def build(corpus, pdfs, open_tdb=DEFAULT_OPEN_TDB, data=None):
    store = Store(corpus, pdfs)
    # pass 1: all sources (cited works resolve across stores)
    for sub in ("invariants/master", "general/records", "text"):
        for f in sorted(glob.glob(os.path.join(store.db, sub, "*.jsonl"))):
            for r in jl(f):
                if r.get("record_type") == "source":
                    store.add_source(r)
    contexts = load_contexts(store)
    records = []
    load_tables(store, "invariants/master", "invariant_table", contexts, records)
    load_tables(store, "general/records", "general_table", contexts, records)
    load_text(store, contexts, records)

    # ids unique: the compact id names a table row; several values of one row share it, so a short hash of
    # the internal observation id is appended (stable, independent of record order)
    seen = Counter(r["id"] for r in records)
    for r in records:
        oid = r.pop("_oid")
        if not r["id"] or seen[r["id"]] > 1:
            r["id"] = (r["id"] or "ORC") + "." + hashlib.sha256(oid.encode()).hexdigest()[:6]
            store.notes["row_ids_with_value_suffix"] += 1
    assert len({r["id"] for r in records}) == len(records)
    records = [clean(r) for r in records]
    records, qreport = quality.apply(records, open_tdb)
    for r in records:
        r["flags"]["extraction_accepted"] = bool(r["flags"].get("extraction_accepted"))
    # literature-level consistency gates (tools/literature_gates.py): corrections only when evidenced and confirmed on
    # the page (tools/literature_visual_checks.json), otherwise demotions / flags; never a printed value or quote
    units, plaus = gates.load_store_evidence(corpus)
    import phase_rule as _pr
    records, qreport["gates"] = gates.apply(
        records, _pr.PhaseRule(data or orchid_cite.DEFAULT_DATA), store_units=units, plaus_store=plaus,
        describe=lambda pk: store.bib.describe(pk[4:]) if pk and pk.startswith("doi:") else {})
    for r in records:
        r["flags"]["extraction_accepted"] = bool(r["flags"].get("extraction_accepted"))
    # model-phase sanity of the corpus label mapping (tools/phase_label_check.py; audit round 2, A2 class 1): a Greek
    # letter or a compound mapped to an element structure, a liquid label on a solid model, a printed composition
    # outside the model phase's field -> mapping removed (kept in model_phase_rejected)
    qreport["phase_labels"] = phase_label_check.apply(records, _pr.PhaseRule(data or orchid_cite.DEFAULT_DATA))
    applied = apply_contributions(records, data or orchid_cite.DEFAULT_DATA)
    if applied:   # accepted contributed corrections (page-checked); none: the release is byte-identical
        qreport["contributed_corrections"] = applied
    qreport["espei"], qreport["espei_phase_rule_failures"] = espei_fields(records, data or orchid_cite.DEFAULT_DATA)
    qreport["espei_phase_rule"] = phase_rule.RULE_TEXT
    records.sort(key=lambda r: (r["system"], r["record_kind"], r.get("pid") or "", r["id"]))
    # uncertainty budget of every value (tools/uncertainty_budget.py; PI decision 2026-10-24): components kept
    # separate, method floors calibrated from the release's own independent measurements; printed values and printed
    # uncertainties are not altered
    store.uncertainty_calibration = ub.calibrate(records, ub.load_dft_calibration(data or orchid_cite.DEFAULT_DATA))
    ub.annotate(records, store.uncertainty_calibration)
    for r in records:   # machine-learning folds (tools/ml_folds.py): leave-system-out and leave-paper-out, 5 folds
        ml_folds.annotate(r)

    for r in records:
        for pk in [r.get("paper")] + [c.get("paper") for c in r.get("cited_sources") or []] + [(r.get("origin") or {}).get("measured_in")]:
            if pk and pk in store.papers:
                store.papers[pk]["systems"].add(r["system"])
    papers = store.finalize_papers()
    # year of each cited source re-derived from its printed reference (page numbers read as years; audit A5 F2)
    cy = Counter()
    quality.fix_cited_years(records, lambda pk: (papers.get(pk) or {}).get("year") if pk else None, cy)
    qreport.setdefault("counts", {}).update({k: v for k, v in sorted(cy.items())})
    used_ctx = {r["paper"] for r in records if (r.get("conditions") or {}).get("kind") == "paper_context"}
    contexts_out = {k: contexts[k] for k in sorted(contexts) if k in used_ctx}
    store.notes["contexts_unused_no_own_experimental_values"] = len(contexts) - len(contexts_out)
    store.quality = qreport
    return store, records, papers, contexts_out


PD_SYSTEMS = set()


def outputs(store, records, papers, contexts):
    pdi = os.path.join(SITE, "data", "phase-diagrams", "index.json")
    if os.path.exists(pdi):
        PD_SYSTEMS.update(json.load(open(pdi, encoding="utf-8")).get("systems", {}))
    files = {}
    cite = orchid_cite.cite_as("literature")
    by_sys = defaultdict(list)
    for r in records:
        by_sys[r["system"]].append(r)
    index_systems = {}
    for sysn in sorted(by_sys):
        recs = by_sys[sysn]
        chunks, cur, size = [], [], 0
        for r in recs:
            b = len(dumps(r).encode()) + 1
            if cur and size + b > MAX_SYSTEM_BYTES:
                chunks.append(cur)
                cur, size = [], 0
            cur.append(r)
            size += b
        chunks.append(cur)
        names = []
        for i, ch in enumerate(chunks):
            name = f"{sysn}.json" if i == 0 else f"{sysn}.part{i + 1}.json"
            pk = set()
            for r in ch:
                pk.add(r.get("paper"))
                pk.add((r.get("origin") or {}).get("measured_in"))
                for c in r.get("cited_sources") or []:
                    pk.add(c.get("paper"))
            ck = {r["paper"] for r in ch if (r.get("conditions") or {}).get("kind") == "paper_context"}
            arts, models, compact = {}, {}, []
            for r in ch:
                r2 = dict(r)
                loc = dict(r2.pop("locator", {}) or {})
                row = loc.pop("row", None)
                akey = loc.pop("asset", None) or r2.get("pid", "").split("#")[0]
                arts.setdefault(akey, loc)
                r2["locator"] = clean({"ref": akey, "row": row})
                if r2.get("model"):
                    mkey = "M" + hashlib.sha256(dumps(r2["model"]).encode()).hexdigest()[:8]
                    models.setdefault(mkey, r2["model"])
                    r2["model"] = mkey
                r2.pop("elements", None)
                compact.append(r2)
            doc = {"system": sysn, "elements": recs[0]["elements"], "part": i + 1, "parts": len(chunks),
                   "cite_as": cite,
                   "release": RELEASE, "licence": "CC BY 4.0 (ORCHID data; quotes are short excerpts kept for verification, the publishers' copyright, not relicensed)",
                   "note": "records[].locator.ref -> artifacts; records[].model -> models (full records: orchid-literature.jsonl.gz)",
                   "records": compact, "artifacts": dict(sorted(arts.items())), "models": dict(sorted(models.items())),
                   "papers": {k: {x: y for x, y in papers[k].items() if x != "systems"} for k in sorted(p for p in pk if p and p in papers)},
                   "contexts": {k: contexts[k] for k in sorted(ck) if k in contexts}}
            files[f"systems/{name}"] = dumps(doc).encode()
            names.append(f"systems/{name}")
        ev = Counter(r["evidence_class"] for r in recs)
        pg = Counter(r["property_group"] for r in recs)
        index_systems[sysn] = {"elements": recs[0]["elements"], "n": len(recs),
                               "phase_diagram": sysn in PD_SYSTEMS,
                               "accepted": sum(1 for r in recs if r["flags"]["extraction_accepted"]),
                               "by_evidence_class": dict(sorted(ev.items())),
                               "by_property_group": dict(sorted(pg.items())),
                               "by_record_kind": dict(sorted(Counter(r["record_kind"] for r in recs).items())),
                               "files": names}
    tot = {
        "records": len(records),
        "accepted": sum(1 for r in records if r["flags"]["extraction_accepted"]),
        "systems": len(by_sys),
        "papers_reporting": sum(1 for p in papers.values() if "reporting_paper" in (p.get("roles") or [])),
        "papers_total": len(papers),
        "by_evidence_class": dict(sorted(Counter(r["evidence_class"] for r in records).items())),
        "by_evidence_class_accepted": dict(sorted(Counter(r["evidence_class"] for r in records if r["flags"]["extraction_accepted"]).items())),
        "by_property_group": dict(sorted(Counter(r["property_group"] for r in records).items())),
        "by_property": dict(sorted(Counter(r.get("property") or "?" for r in records).most_common(60))),
        "by_record_kind": dict(sorted(Counter(r["record_kind"] for r in records).items())),
        "with_quote": sum(1 for r in records if r.get("quote")),
        "with_paper_context": sum(1 for r in records if (r.get("conditions") or {}).get("kind") == "paper_context"),
        "cited_source_conditions": sum(1 for r in records if (r.get("conditions") or {}).get("kind") == "cited_source"),
        "model_linked_to_tdb": sum(1 for r in records if r.get("model")),
        "contexts": len(contexts),
        "context_items": sum(len(c["items"]) for c in contexts.values()),
        "excluded": dict(sorted(store.excluded.items())),
        "espei": store.quality["espei"],
        "uncertainty_budget": ub.totals(records, store.uncertainty_calibration),
        "gates": {k: store.quality["gates"][k] for k in ("release", "records_in", "records_out", "accepted_before",
                                                          "accepted_after", "espei_ready_before", "espei_ready_after",
                                                          "espei_ready_after_note", "records_corrected",
                                                          "records_with_gate_failures", "visual_checks",
                                                          "visual_check_entries", "visual_checks_note")},
        "quality": dict(store.quality["counts"], duplicates_merged=store.quality["duplicates"].get("records_removed", 0),
                        duplicate_groups=store.quality["duplicates"].get("groups", 0),
                        formation_basis_inferred_from_magnitude=sum(
                            1 for r in records if (r.get("flags") or {}).get("basis_inferred_from_magnitude")),
                        by_unit_basis=dict(sorted(Counter(r["unit_basis"] for r in records if r.get("unit_basis")
                                                          and r.get("property") in quality.ENERGY_KJ).items()))),
    }
    index = {"release": RELEASE, "licence": {"derived_data": "CC BY 4.0", "quotes": "short excerpts kept so that each value can be verified; publishers' copyright, not relicensed (LICENSES.md)",
                                             "source_papers": "not redistributed"},
             "schema": "SCHEMA.md", "totals": tot, "cite_as": cite,
             "doi_corrections": dict(CORR.provenance(), papers_affected=sum(1 for p in papers.values() if p.get("doi_status")),
                                     values_affected=sum(1 for r in records if any(
                                         x.get("doi_status") for x in r.get("cited_sources") or []) or (r.get("model") or {}).get("model_doi_status")
                                         or (papers.get(r.get("paper")) or {}).get("doi_status"))),
             "files": {"papers": "papers.json", "contexts": "contexts.json", "jsonl": "orchid-literature.jsonl.gz",
                       "csv": "orchid-literature.csv.gz", "quality_fixes": "QUALITY_FIXES.md",
                       "quality_fixes_json": "quality_fixes.json", "citation_txt": "CITATION.txt",
                       "citation_cff": "CITATION.cff", "uncertainty_calibration": "UNCERTAINTY_CALIBRATION.json",
                       "uncertainty_calibration_md": "UNCERTAINTY_CALIBRATION.md"},
             "systems": index_systems}
    files["index.json"] = dumps(index).encode()
    papers = repair_replacement_chars(papers, crossref_name_words(store.bib.cache_dir))
    files["papers.json"] = dumps(papers).encode()
    files["contexts.json"] = dumps(contexts).encode()
    files["orchid-literature.jsonl.gz"] = gz(("\n".join(dumps(r) for r in records) + "\n").encode())
    sbuf = io.StringIO()
    w = csv.writer(sbuf, lineterminator="\n")
    w.writerow(CSV_COLUMNS)
    for r in records:
        row = dict(r)
        p = papers.get(r.get("paper")) or {}
        row["paper_citation"] = p.get("citation")
        row["paper_doi"] = p.get("doi")
        row["extraction_accepted"] = r["flags"]["extraction_accepted"]
        for k in ("espei_ready", "espei_certified", "calphad_check"):
            row[k] = r["flags"].get(k)
        qv = ((r.get("uncertainty_budget") or {}).get("quantities") or {}).get("value") or {}
        row["sigma_effective"] = qv.get("sigma_effective")
        row["sigma_method_class"] = (qv.get("sigma_method") or {}).get("class")
        row["conditions_summary"] = summarize_context(contexts.get(r["paper"])) if (r.get("conditions") or {}).get("kind") == "paper_context" else (
            "conditions in the cited source" if (r.get("conditions") or {}).get("kind") == "cited_source" else None)
        w.writerow([csv_value(row.get(c)) for c in CSV_COLUMNS])
    files["orchid-literature.csv.gz"] = gz(sbuf.getvalue().encode())
    files["SCHEMA.md"] = SCHEMA_MD.format(release=RELEASE, status=orchid_cite.STATUS, kinds=", ".join("`%s`" % k for k in quality.TRANSITION_KINDS),
                                          propmap=", ".join("`%s` -> `%s`" % (k, v) for k, v in sorted(quality.PROPERTY_FIX.items()))).encode()
    files["QUALITY_FIXES.md"] = (quality.quality_md(store.quality, RELEASE) + "\n" + gates.gates_md(store.quality["gates"])
                                + (phase_label_check.report_md(store.quality["phase_labels"]) if store.quality.get("phase_labels") else "")).encode()
    files["quality_fixes.json"] = dumps(dict(store.quality, release=RELEASE)).encode()
    # the CSV cannot carry metadata (a leading comment line breaks CSV readers): citation in a file beside it
    files["CITATION.txt"] = orchid_cite.citation_txt("literature", intro=(
        "ORCHID literature values: orchid-literature.csv.gz, orchid-literature.jsonl.gz, index.json, "
        "systems/*.json, papers.json, contexts.json (release %s)." % RELEASE)).encode()
    files["CITATION.cff"] = orchid_cite.cff("literature").encode()
    files["UNCERTAINTY_CALIBRATION.json"] = ub.calibration_bytes(store.uncertainty_calibration)
    files["UNCERTAINTY_CALIBRATION.md"] = ub.calibration_md(store.uncertainty_calibration, records).encode()
    return files, tot


LIT_FIELDS = [("id", "sc:Text"), ("system", "sc:Text"), ("record_kind", "sc:Text"), ("observation_type", "sc:Text"),
              ("property_group", "sc:Text"), ("property", "sc:Text"), ("value", "sc:Float"), ("unit", "sc:Text"),
              ("uncertainty", "sc:Text"), ("uncertainty_budget", "sc:Text"), ("sigma_effective", "sc:Float"),
              ("temperature_K", "sc:Float"), ("temperature_printed", "sc:Text"),
              ("pressure_Pa", "sc:Float"), ("phases", "sc:Text"), ("phase_compositions", "sc:Text"),
              ("reaction", "sc:Text"), ("evidence_class", "sc:Text"), ("origin", "sc:Text"), ("model", "sc:Text"),
              ("method", "sc:Text"), ("conditions", "sc:Text"), ("cited_sources", "sc:Text"), ("paper", "sc:Text"),
              ("paper_citation", "sc:Text"), ("paper_doi", "sc:Text"), ("quote", "sc:Text"), ("quote_page", "sc:Integer"),
              ("flags", "sc:Text"), ("extraction_accepted", "sc:Boolean"), ("chemical_system", "sc:Text"),
              ("paper_id", "sc:Text"), ("fold_lso", "sc:Integer"), ("fold_lpo", "sc:Integer")]


def croissant_update(path, files):
    """Add / refresh the literature entries of data/croissant.json (other entries untouched)."""
    d = json.load(open(path, encoding="utf-8"))
    dist = [x for x in d.get("distribution", []) if not str(x.get("@id", "")).startswith("literature")]
    sha = lambda b: hashlib.sha256(b).hexdigest()
    base = "https://orchid-project.org/data/literature/"
    dist += [
        {"@type": "cr:FileObject", "@id": "literature-csv-gz", "name": "orchid-literature.csv.gz",
         "contentUrl": base + "orchid-literature.csv.gz", "encodingFormat": "application/gzip",
         "sha256": sha(files["orchid-literature.csv.gz"])},
        {"@type": "cr:FileObject", "@id": "literature-csv", "name": "orchid-literature.csv",
         "containedIn": {"@id": "literature-csv-gz"}, "encodingFormat": "text/csv"},
        {"@type": "cr:FileObject", "@id": "literature-jsonl-gz", "name": "orchid-literature.jsonl.gz",
         "contentUrl": base + "orchid-literature.jsonl.gz", "encodingFormat": "application/gzip",
         "sha256": sha(files["orchid-literature.jsonl.gz"])},
        {"@type": "cr:FileObject", "@id": "literature-index", "name": "literature-index.json",
         "contentUrl": base + "index.json", "encodingFormat": "application/json", "sha256": sha(files["index.json"])},
        {"@type": "cr:FileObject", "@id": "literature-papers", "name": "literature-papers.json",
         "contentUrl": base + "papers.json", "encodingFormat": "application/json", "sha256": sha(files["papers.json"])},
        {"@type": "cr:FileObject", "@id": "literature-contexts", "name": "literature-contexts.json",
         "contentUrl": base + "contexts.json", "encodingFormat": "application/json", "sha256": sha(files["contexts.json"])},
        {"@type": "cr:FileSet", "@id": "literature-systems", "name": "literature-systems",
         "encodingFormat": "application/json", "includes": "data/literature/systems/*.json"},
    ]
    d["distribution"] = dist
    rs = [x for x in d.get("recordSet", []) if x.get("@id") != "literature"]
    rs.append({"@type": "cr:RecordSet", "@id": "literature", "name": "literature",
               "description": "Values extracted from published papers (CC BY 4.0 data; quotes are short excerpts kept for "
                              "verification and remain the publishers' copyright). "
                              "Nested fields are compact JSON strings; see data/literature/SCHEMA.md.",
               "field": [{"@type": "cr:Field", "@id": f"literature/{n}", "name": n, "dataType": t,
                          "source": {"fileObject": {"@id": "literature-csv"}, "extract": {"column": n}}}
                         for n, t in LIT_FIELDS]})
    d["recordSet"] = rs
    d["version"] = orchid_cite.VERSION
    d["citeAs"] = orchid_cite.bibtex_block([])
    if "literature" not in d.get("description", ""):
        d["description"] = d["description"].rstrip(".") + ", and literature data (invariant reactions, transitions, phase compositions, thermochemical and structural values with evidence class, citation and quote)."
    return json.dumps(d, indent=1, ensure_ascii=False).encode()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--corpus", default=DEFAULT_CORPUS, help="SGTE-corpus build directory (read only)")
    ap.add_argument("--pdfs", default=DEFAULT_PDFS, help="PDF corpus root (read only; hashed for no-DOI citations)")
    ap.add_argument("--out", default=os.path.join(SITE, OUT_REL))
    ap.add_argument("--open-tdb", default=DEFAULT_OPEN_TDB, help="open-tdb checkout holding hulls/snapshot*/ (read only)")
    ap.add_argument("--data", default=orchid_cite.DEFAULT_DATA,
                    help="orchid-data checkout (read only: phases/ registry and phase-diagrams/ for the ESPEI phase rule)")
    ap.add_argument("--check", action="store_true", help="verify the files on disk match a fresh export")
    a = ap.parse_args(argv)
    store, records, papers, contexts = build(a.corpus, a.pdfs, a.open_tdb, a.data)
    files, tot = outputs(store, records, papers, contexts)
    cpath = os.path.join(SITE, "data", "croissant.json")
    croissant = croissant_update(cpath, files) if os.path.exists(cpath) else None
    if a.check:
        bad = []
        if croissant is not None and open(cpath, "rb").read() != croissant:
            bad.append("../croissant.json")
        for rel, data in files.items():
            p = os.path.join(a.out, rel)
            if not os.path.exists(p) or open(p, "rb").read() != data:
                bad.append(rel)
        existing = {os.path.relpath(p, a.out) for p in glob.glob(os.path.join(a.out, "systems", "*.json"))}
        stale = sorted(existing - set(files))
        if bad or stale:
            print(f"out of date: {len(bad)} changed/missing, {len(stale)} stale", file=sys.stderr)
            for x in (bad + stale)[:20]:
                print("  " + x, file=sys.stderr)
            return 1
        print(f"ok: {len(files)} files current")
        return 0
    os.makedirs(os.path.join(a.out, "systems"), exist_ok=True)
    for p in glob.glob(os.path.join(a.out, "systems", "*.json")):
        if os.path.relpath(p, a.out) not in files:
            os.remove(p)
    for rel, data in files.items():
        p = os.path.join(a.out, rel)
        tmp = p + ".tmp"
        with open(tmp, "wb") as fh:
            fh.write(data)
        os.replace(tmp, p)
    if croissant is not None:
        with open(cpath + ".tmp", "wb") as fh:
            fh.write(croissant)
        os.replace(cpath + ".tmp", cpath)
    sizes = {k: len(v) for k, v in files.items()}
    big = max((v, k) for k, v in sizes.items() if k.startswith("systems/"))
    print(json.dumps({"totals": tot, "notes": dict(store.notes),
                      "sizes": {k: v for k, v in sorted(sizes.items()) if not k.startswith("systems/")},
                      "system_files": sum(1 for k in sizes if k.startswith("systems/")),
                      "system_bytes_total": sum(v for k, v in sizes.items() if k.startswith("systems/")),
                      "largest_system_file": big}, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
