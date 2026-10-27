#!/usr/bin/env python3
"""Prepare the ORCHID data release for Zenodo: eight citable datasets plus an umbrella record.

Records (names used in DOI placeholders "10.5281/zenodo.PENDING-<name>" and in <orchid-data>/zenodo/):

  orchid-data      umbrella "ORCHID data": documentation, licences, Parquet manifest, Zenodo metadata (hasPart)
  literature       ORCHID literature values                       literature/ (minus the DFT benchmark files) + Parquet
  phase-diagrams   ORCHID calculated binary phase diagrams        phase-diagrams/ + Parquet
  hulls            ORCHID binary convex hulls (MP, OQMD, CALPHAD) hulls/ + Parquet
  dft-benchmark    ORCHID DFT-vs-experiment formation-enthalpy benchmark   dft_parity, Kim 2017 cross-check + Parquet
  phase-registry   ORCHID phase registry (binary systems)         phases/ (registry JSON/JSON-LD/CSV/Parquet, structure types)
  assessed-phase-assemblages  ORCHID assessed phase assemblages   assessed-phase-assemblages/ (points in assessed
                   phase diagrams: JSONL/CSV/Parquet, ESPEI ZPF datasets; tools/build_assessed_assemblages.py)
  espei            ORCHID ESPEI datasets from the literature      espei/ (ESPEI JSON per evidence class, index, calibration;
                   tools/build_espei_bundle.py)
  bibliography     ORCHID bibliography                            bibliography/ (papers, typed edges, Crossref metadata,
                   DOI checks; tools/build_bibliography.py; bibliography/pages/ are web views, not released)

A tenth record, orchid-examples (software), archives the runnable examples repository (orchid-calphad/orchid-examples):
its files are not in the data repository but in a pinned commit of that repository (`--examples`, `--examples-ref`);
`manifest` writes zenodo/orchid-examples.files.json from `git ls-tree`, `bundle --record orchid-examples` builds the zip
with `git archive`. The ORCHID benchmark (orchid-bench) is a code repository with its own Zenodo GitHub-integration
record, not part of this data release.

Per record the data repository holds zenodo/<name>.metadata.json (Zenodo deposition metadata) and
zenodo/<name>.files.json (the files of the record with sizes and SHA-256, written by `manifest`).

  python3 tools/zenodo_release.py manifest [--check]   # (re)write zenodo/<name>.files.json from the RECORDS rules
  python3 tools/zenodo_release.py card [--check]       # Parquet sizes in DATASET_CARD.md; copy to hf/README.md
  python3 tools/zenodo_release.py metadata [--check]   # zenodo/<name>.metadata.json + metadata_numbers.json (build_zenodo_metadata.py)
  python3 tools/zenodo_release.py scan                 # restricted types/paths, secret patterns, local paths in every record
  python3 tools/zenodo_release.py plan                 # records, files, sizes; every file present and unchanged, < 50 GB
  python3 tools/zenodo_release.py bundle [--record N]  # dist/orchid-<name>-v2026.10.zip (deterministic, MANIFEST.sha256)
  python3 tools/zenodo_release.py writeback --doi literature=10.5281/zenodo.123 [--doi ...] [--no-rebuild]
  python3 tools/zenodo_release.py credit [--check]     # accepted, opted-in contributors -> CONTRIBUTORS.md, Zenodo
                                                       # metadata (and CITATION.cff if contributors.json says so),
                                                       # then refresh the manifests (run before `bundle`)

`writeback` records the reserved DOIs in zenodo/dois.json, replaces the placeholders in every text file of both
repositories (data repository, website pages, CITATION.cff files, llms.txt, croissant.json, DATASET_CARD.md,
Zenodo metadata, ...), then reruns the builders (which read the DOIs through tools/orchid_cite.py) so derived
hashes stay right, refreshes the file manifests and runs every --check. TDB files, gzip and Parquet bytes are never
edited as text. Nothing here talks to Zenodo: the PI uploads the bundles (zenodo/README.md in the data repository).

Standard library only (the Parquet rebuild inside writeback runs tools/build_parquet.py with --parquet-python).
"""
from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orchid_cite  # noqa: E402

SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAX_BYTES = 50 * 10**9  # Zenodo per-record limit
ZIP_DATE = (1980, 1, 1, 0, 0, 0)
STORED = (".gz", ".parquet", ".zip", ".png", ".jpg")

DOCS = ["README.md", "CHANGELOG.md", "LICENSE", "LICENSES.md", "NOTICE.md", "CITATIONS.md", "CITATIONS.json"]
RECORDS = {  # name -> (include patterns, exclude patterns); "dir/**" = everything below dir
    "orchid-data": (DOCS + ["CITATION.cff", "DATASET_CARD.md", "parquet/MANIFEST.json", "zenodo/README.md",
                            "zenodo/dois.json", "zenodo/*.metadata.json", "zenodo/metadata_numbers.json"], []),
    "literature": (DOCS + ["DATASET_CARD.md", "zenodo/literature.metadata.json", "literature/**",
                           "parquet/literature_values.parquet", "parquet/papers.parquet",
                           "parquet/experimental_context.parquet", "parquet/MANIFEST.json"],
                   ["literature/dft_parity.json", "literature/kim2017_crosscheck.*"]),
    "phase-diagrams": (DOCS + ["zenodo/phase-diagrams.metadata.json", "phase-diagrams/**",
                               "parquet/phase_diagram_index.parquet", "parquet/phase_diagram_invariants.parquet",
                               "parquet/MANIFEST.json"], []),
    "hulls": (DOCS + ["zenodo/hulls.metadata.json", "hulls/**", "parquet/hull_vertices.parquet",
                      "parquet/MANIFEST.json"], []),
    "dft-benchmark": (DOCS + ["zenodo/dft-benchmark.metadata.json", "literature/dft_parity.json",
                              "literature/kim2017_crosscheck.json", "literature/kim2017_crosscheck.md",
                              "literature/QUALITY_FIXES.md", "parquet/dft_parity.parquet", "parquet/MANIFEST.json"], []),
    "phase-registry": (DOCS + ["zenodo/phase-registry.metadata.json", "phases/**"], ["phases/pages/**"]),
    "assessed-phase-assemblages": (DOCS + ["zenodo/assessed-phase-assemblages.metadata.json",
                                           "assessed-phase-assemblages/**"], []),
    "espei": (DOCS + ["zenodo/espei.metadata.json", "espei/**"], []),
    "bibliography": (DOCS + ["zenodo/bibliography.metadata.json", "bibliography/**"], ["bibliography/pages/**"]),
}
assert set(RECORDS) == set(orchid_cite.NAMES)
DATA_DIRS = ("literature", "phase-diagrams", "hulls", "parquet", "phases", "assessed-phase-assemblages", "espei",
             "bibliography")  # every file here must be in some record
WEB_VIEWS = ("phases/pages/", "bibliography/pages/")  # generated static pages inside a data dir (tools/build_static_pages.py): not released

TEXT_EXT = (".json", ".md", ".html", ".txt", ".cff", ".xml", ".csv", ".yml", ".yaml")
SKIP_DIRS = {".git", "dist", "__pycache__", "node_modules"}
SKIP_PREFIX = ("data/databases/", "review/", "tools/")  # TDB bytes, frozen review kit, code (reads dois.json)


# ------------------------------------------------------------------------------------------- helpers --
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def all_files(root):
    out = []
    for d, dirs, files in os.walk(root):
        dirs[:] = sorted(x for x in dirs if x not in SKIP_DIRS and not x.startswith("."))
        for f in files:
            if (f.startswith(".") and f != ".zenodo.json") or f.endswith(".tmp"):
                continue
            out.append(os.path.relpath(os.path.join(d, f), root).replace(os.sep, "/"))
    return sorted(out)


def matches(path, pat):
    if pat.endswith("/**"):
        return path.startswith(pat[:-2])
    return fnmatch.fnmatchcase(path, pat)


def select(files, name):
    inc, exc = RECORDS[name]
    sel = [f for f in files if any(matches(f, p) for p in inc) and not any(matches(f, p) for p in exc)]
    missing = [p for p in inc if "*" not in p and p not in sel]
    return sel, missing


# ------------------------------------------------------------------------------- examples (software) record --
EXAMPLES = "orchid-examples"
DEFAULT_EXAMPLES = os.environ.get("ORCHID_EXAMPLES") or os.path.join(os.path.dirname(SITE), "orchid-examples")


def git(repo, *args, binary=False):
    r = subprocess.run(["git", "-C", repo] + list(args), capture_output=True)
    if r.returncode:
        raise SystemExit("git %s failed in %s: %s" % (" ".join(args), repo, r.stderr.decode().strip()))
    return r.stdout if binary else r.stdout.decode()


def examples_manifest(repo, ref):
    """Files of the examples record = the tree of one commit (what `git archive` puts into the zip)."""
    commit = git(repo, "rev-parse", "--verify", ref + "^{commit}").strip()
    entries = []
    for line in git(repo, "ls-tree", "-r", "-l", "-z", commit).split("\0"):
        if not line:
            continue
        meta, path = line.split("\t", 1)
        mode, kind, blob, size = meta.split()
        if kind != "blob":
            raise SystemExit("examples tree holds a %s (submodule?): %s" % (kind, path))
        data = git(repo, "cat-file", "blob", blob, binary=True)
        entries.append({"path": path, "bytes": int(size), "sha256": hashlib.sha256(data).hexdigest()})
    entries.sort(key=lambda e: e["path"])
    return {"record": EXAMPLES, "title": orchid_cite.SOFTWARE[EXAMPLES][0], "version": orchid_cite.VERSION,
            "bundle": bundle_name(EXAMPLES), "source": "git archive of github.com/orchid-calphad/orchid-examples",
            "commit": commit, "n_files": len(entries), "total_bytes": sum(e["bytes"] for e in entries), "files": entries}


def manifest_path(data, name):
    return os.path.join(data, "zenodo", "%s.files.json" % name)


def bundle_name(name):
    return "%s-%s.zip" % (name if name.startswith("orchid-") else "orchid-" + name, orchid_cite.RELEASE_TAG)


def build_manifest(data, name, files):
    sel, missing = select(files, name)
    if missing:
        raise SystemExit("%s: listed file(s) missing: %s" % (name, ", ".join(missing)))
    entries = [{"path": p, "bytes": os.path.getsize(os.path.join(data, p)), "sha256": sha256_file(os.path.join(data, p))}
               for p in sel]
    return {"record": name, "title": orchid_cite.entry(name, data)["title"], "version": orchid_cite.VERSION,
            "bundle": bundle_name(name), "include": RECORDS[name][0], "exclude": RECORDS[name][1],
            "n_files": len(entries), "total_bytes": sum(e["bytes"] for e in entries), "files": entries}


def dump(obj):
    return json.dumps(obj, indent=1, ensure_ascii=False, sort_keys=False) + "\n"


def human(n):
    for u in ("B", "kB", "MB", "GB"):
        if n < 1000 or u == "GB":
            return "%.1f %s" % (n, u) if u != "B" else "%d B" % n
        n /= 1000.0


# -------------------------------------------------------------------------------------------- manifest --
def cmd_manifest_examples(a):
    """zenodo/orchid-examples.files.json, pinned to a commit of the examples repository (--examples-ref, default HEAD)."""
    p = manifest_path(a.data, EXAMPLES)
    if not os.path.exists(os.path.join(a.examples, ".git")):
        print("examples repository not found at %s: skipped" % a.examples)
        return 0
    cur = json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None
    ref = a.examples_ref or (cur["commit"] if cur else "HEAD")  # keep the pin; re-pin only with --examples-ref
    text = dump(examples_manifest(a.examples, ref))
    if (open(p, encoding="utf-8").read() if os.path.exists(p) else None) != text:
        if a.check:
            print("stale manifest: " + os.path.relpath(p, a.data), file=sys.stderr)
            return 1
        with open(p + ".tmp", "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(p + ".tmp", p)
        print("wrote " + os.path.relpath(p, a.data))
    return 0


def cmd_manifest(a):
    files = all_files(a.data)
    bad = []
    for name in RECORDS:
        text = dump(build_manifest(a.data, name, files))
        p = manifest_path(a.data, name)
        cur = open(p, encoding="utf-8").read() if os.path.exists(p) else None
        if cur != text:
            if a.check:
                bad.append(p)
            else:
                with open(p + ".tmp", "w", encoding="utf-8") as fh:
                    fh.write(text)
                os.replace(p + ".tmp", p)
                print("wrote " + os.path.relpath(p, a.data))
    rc = cmd_manifest_examples(a)
    if a.check:
        if bad:
            print("stale manifest(s): " + ", ".join(os.path.relpath(p, a.data) for p in bad), file=sys.stderr)
            return 1
        print("zenodo manifests current (%d data records + examples)" % len(RECORDS))
    return rc


def cmd_metadata(a):
    import build_zenodo_metadata
    return build_zenodo_metadata.run(a.data, check=a.check, export=getattr(a, "export", None))


def cmd_credit(a):
    """Release hook: credit accepted contributors (contributions/tools/orchid_contrib.py credit in the data repository)."""
    tool = os.path.join(a.data, "contributions", "tools", "orchid_contrib.py")
    if not os.path.exists(tool):
        print("no contributions/tools/orchid_contrib.py in %s: nothing to credit" % a.data)
        return 0
    rc = subprocess.run([sys.executable, tool, "credit", "--data", a.data] + (["--check"] if a.check else [])).returncode
    if rc or a.check:
        return rc
    a.check = False
    return cmd_manifest(a)


def load_manifest(data, name):
    p = manifest_path(data, name)
    if not os.path.exists(p):
        raise SystemExit("missing %s: run `zenodo_release.py manifest`" % p)
    return json.load(open(p, encoding="utf-8"))


def verify(data, name, files=None):
    """Problems of one record: manifest vs rules, and every listed file present with its size and SHA-256."""
    m = load_manifest(data, name)
    probs = []
    if m.get("include") != RECORDS[name][0] or m.get("exclude") != RECORDS[name][1]:
        probs.append("manifest rules differ from RECORDS (rerun `manifest`)")
    if files is not None:
        sel, missing = select(files, name)
        probs += ["listed but missing: " + p for p in missing]
        listed = [e["path"] for e in m["files"]]
        if sel != listed:
            extra, gone = sorted(set(sel) - set(listed)), sorted(set(listed) - set(sel))
            probs += ["new file not in manifest: " + p for p in extra[:10]]
            probs += ["manifest file no longer selected: " + p for p in gone[:10]]
    for e in m["files"]:
        p = os.path.join(data, e["path"])
        if not os.path.exists(p):
            probs.append("missing: " + e["path"])
        elif os.path.getsize(p) != e["bytes"] or sha256_file(p) != e["sha256"]:
            probs.append("changed since manifest: " + e["path"])
    if m["total_bytes"] >= MAX_BYTES:
        probs.append("total %s >= 50 GB Zenodo limit" % human(m["total_bytes"]))
    return m, probs


# ------------------------------------------------------------------------------------------------ card --
def card_text(data):
    """DATASET_CARD.md with num_examples / download_size taken from parquet/MANIFEST.json and phases/parquet.json."""
    man = json.load(open(os.path.join(data, "parquet", "MANIFEST.json"), encoding="utf-8"))["tables"]
    pj = os.path.join(data, "phases", "parquet.json")  # phase registry tables (tools/build_phase_registry.py)
    if os.path.exists(pj):
        pm = json.load(open(pj, encoding="utf-8"))
        for cfg, ent in (("phases", pm), ("structure_types", pm.get("structure_types"))):
            if ent and os.path.exists(os.path.join(data, "phases", ent["table"])):
                man[cfg] = {"rows": ent["rows"], "bytes": os.path.getsize(os.path.join(data, "phases", ent["table"]))}
    ap_ = os.path.join(data, "assessed-phase-assemblages", "parquet.json")  # tools/build_assessed_assemblages.py
    if os.path.exists(ap_):
        ent = json.load(open(ap_, encoding="utf-8"))
        tp = os.path.join(data, "assessed-phase-assemblages", ent["table"])
        if os.path.exists(tp):
            man["assessed_phase_assemblages"] = {"rows": ent["rows"], "bytes": os.path.getsize(tp)}
    bp = os.path.join(data, "bibliography", "parquet.json")  # tools/build_bibliography.py: papers and edges tables
    if os.path.exists(bp):
        bj = json.load(open(bp, encoding="utf-8"))
        for cfg, key in (("bibliography_papers", "papers"), ("bibliography_edges", "edges")):
            ent = bj.get(key)
            tp = os.path.join(data, "bibliography", ent["table"]) if ent else None
            if tp and os.path.exists(tp):
                man[cfg] = {"rows": ent["rows"], "bytes": os.path.getsize(tp)}
    out, cur, info = [], None, False
    for line in open(os.path.join(data, "DATASET_CARD.md"), encoding="utf-8").read().split("\n"):
        if line.startswith("dataset_info:"):
            info = True
        elif info and line.startswith("---"):
            info = False
        tm = re.match(r"^\| `(\w+)` \| [\d,]+ \|", line)   # prose table of the card: rows of each config
        if not info and tm and tm.group(1) in man:
            line = re.sub(r"^(\| `\w+` \| )[\d,]+( \|)", lambda m_: m_.group(1) + "{:,}".format(man[tm.group(1)]["rows"]) + m_.group(2), line)
        m = re.match(r"^- config_name: (\S+)$", line)
        if info and m:
            cur = m.group(1)
        if info and cur in man:
            if re.match(r"^    num_examples: \d+$", line):
                line = "    num_examples: %d" % man[cur]["rows"]
            elif re.match(r"^  download_size: \d+$", line):
                line = "  download_size: %d" % man[cur]["bytes"]
        out.append(line)
    return "\n".join(out)


def cmd_card(a):
    """Refresh the Parquet sizes in DATASET_CARD.md and copy it to hf/README.md (the mirror card)."""
    text = card_text(a.data)
    bad = []
    for rel in ("DATASET_CARD.md", "hf/README.md"):
        p = os.path.join(a.data, rel)
        cur = open(p, encoding="utf-8").read() if os.path.exists(p) else None
        if cur != text:
            if a.check:
                bad.append(rel)
            else:
                with open(p, "w", encoding="utf-8") as fh:
                    fh.write(text)
                print("wrote " + rel)
    if a.check:
        if bad:
            print("stale: " + ", ".join(bad), file=sys.stderr)
            return 1
        print("dataset card current (sizes match parquet/MANIFEST.json; hf/README.md identical)")
    return 0


# ------------------------------------------------------------------------------------------------ plan --
def plan_examples(a, dois):
    """The software record: manifest present, pinned commit exists and still matches, examples HEAD not ahead of the pin."""
    probs = []
    p = manifest_path(a.data, EXAMPLES)
    print("%-15s %-72s" % (EXAMPLES, orchid_cite.SOFTWARE[EXAMPLES][0][:72]))
    doi = dois[EXAMPLES]
    print("    DOI %s%s" % (doi, "  (placeholder: reserve on Zenodo, then writeback)" if orchid_cite.is_placeholder(doi) else ""))
    if not os.path.exists(os.path.join(a.data, "zenodo", EXAMPLES + ".metadata.json")):
        probs.append("missing zenodo/%s.metadata.json" % EXAMPLES)
    if not os.path.exists(p):
        probs.append("missing %s: run `zenodo_release.py manifest`" % p)
    elif not os.path.exists(os.path.join(a.examples, ".git")):
        probs.append("examples repository not found at %s (--examples)" % a.examples)
    else:
        m = json.load(open(p, encoding="utf-8"))
        print("    %d files, %s -> dist/%s (git archive of commit %s)" % (m["n_files"], human(m["total_bytes"]), m["bundle"], m["commit"][:10]))
        if dump(examples_manifest(a.examples, m["commit"])) != open(p, encoding="utf-8").read():
            probs.append("examples manifest does not match its pinned commit (rerun `manifest`)")
        head = git(a.examples, "rev-parse", "HEAD").strip()
        if head != m["commit"]:
            probs.append("examples HEAD %s differs from the pinned commit %s (rerun `manifest --examples-ref <tag>` to re-pin)"
                         % (head[:10], m["commit"][:10]))
        if git(a.examples, "status", "--porcelain").strip():
            probs.append("examples working tree has uncommitted changes (the zip is built from the commit, not the tree)")
        if m["n_files"] >= 100000 or m["total_bytes"] >= MAX_BYTES:
            probs.append("examples record too large")
    for q in probs:
        print("    PROBLEM: " + q)
    print()
    return len(probs)


def cmd_plan(a):
    files = all_files(a.data)
    dois = orchid_cite.dois(a.data)
    bad = 0
    covered = set()
    print("ORCHID data release %s (version %s), data repository %s\n" % (orchid_cite.RELEASE_TAG, orchid_cite.VERSION, a.data))
    for name in RECORDS:
        m, probs = verify(a.data, name, files)
        covered.update(e["path"] for e in m["files"])
        meta_p = os.path.join(a.data, "zenodo", "%s.metadata.json" % name)
        if not os.path.exists(meta_p):
            probs.append("missing zenodo/%s.metadata.json" % name)
        doi = dois[name]
        print("%-15s %-72s" % (name, m["title"]))
        print("    DOI %s%s" % (doi, "  (placeholder: reserve on Zenodo, then writeback)" if orchid_cite.is_placeholder(doi) else ""))
        print("    %d files, %s -> %s" % (m["n_files"], human(m["total_bytes"]), os.path.join("dist", m["bundle"])))
        tops = {}
        for e in m["files"]:
            k = e["path"].split("/")[0] if "/" in e["path"] else "(top level)"
            n, b = tops.get(k, (0, 0))
            tops[k] = (n + 1, b + e["bytes"])
        for k in sorted(tops):
            print("      %-14s %6d files %10s" % (k, tops[k][0], human(tops[k][1])))
        for p in probs[:20]:
            print("    PROBLEM: " + p)
        bad += len(probs)
        print()
    bad += plan_examples(a, dois)
    loose = [f for f in files if f.split("/")[0] in DATA_DIRS and f not in covered and not f.startswith(WEB_VIEWS)]
    for f in loose[:20]:
        print("PROBLEM: data file in no record: " + f)
    bad += len(loose)
    not_released = sorted({f.split("/")[0] for f in files if f not in covered and "/" in f} - set(DATA_DIRS))
    print("Not in any record (static web views, the hf mirror card, the zenodo file manifests): %s" % ", ".join(not_released))
    print("Limit check: every record < 50 GB; upload one zip per record (Zenodo also caps a record at 100 files).")
    print("OK" if not bad else "FAIL: %d problem(s)" % bad)
    return 1 if bad else 0


# ------------------------------------------------------------------------------------------------- scan --
RESTRICTED_EXT = (".pdf", ".png", ".jpg", ".jpeg", ".tif", ".tiff", ".gif", ".bmp", ".webp", ".djvu", ".epub", ".docx", ".doc",
                  ".pptx", ".xlsx", ".key", ".pem", ".p12", ".pfx", ".sqlite", ".db", ".pkl", ".npy", ".npz", ".h5", ".ipynb_checkpoints")
RESTRICTED_PATH = re.compile(r"(^|/)(private|secret|secrets|\.env[^/]*|id_rsa[^/]*|id_ed25519[^/]*|credentials[^/]*|"
                             r"[^/]*handbook[-_ ]?(page|image|crop)[^/]*|[^/]*okamoto[-_ ]?(page|crop|image)[^/]*|pages_png|page_images?)(/|$)", re.I)
SECRET_RE = re.compile(r"AKIA[0-9A-Z]{16}|gh[pousr]_[A-Za-z0-9]{36,}|github_pat_[A-Za-z0-9_]{20,}|sk-[A-Za-z0-9_-]{24,}|hf_[A-Za-z0-9]{30,}|"
                       r"pypi-AgE[A-Za-z0-9_-]{20,}|xox[baprs]-[A-Za-z0-9-]{10,}|AIza[0-9A-Za-z_-]{35}|"
                       r"-----BEGIN [A-Z ]*PRIVATE KEY-----|X-ELS-APIKey\s*[:=]\s*\S{8,}|Bearer\s+[A-Za-z0-9._-]{30,}|"
                       r"(?i:api[_-]?key|secret|token|passwd|password)[\"']?\s*[:=]\s*[\"'][A-Za-z0-9/+_-]{20,}[\"']")
LOCAL_RE = re.compile(r"/Users/[A-Za-z0-9._-]+/|/home/[a-z][a-z0-9_-]{2,}/|C:\\Users\\")
SCAN_MAX = 200 * 2**20  # bytes of one file scanned as text


def _iter_text(path):
    """Text of a file for scanning: plain, or gunzipped (.gz); None when binary/too large."""
    import gzip
    try:
        if path.endswith(".gz"):
            with gzip.open(path, "rb") as fh:
                raw = fh.read(SCAN_MAX)
        elif path.endswith(TEXT_EXT + (".jsonl", ".tdb", ".bib", ".tsv", ".py", ".sh", ".toml", ".cfg", ".ini", ".tex", ".yml", ".ipynb")):
            if os.path.getsize(path) > SCAN_MAX:
                return None
            with open(path, "rb") as fh:
                raw = fh.read()
        else:
            return None
    except OSError:
        return None
    return raw.decode("utf-8", "replace")


def scan_files(paths, root):
    """-> (restricted, secrets, local_paths): findings over repository-relative paths."""
    restricted, secrets, local = [], [], []
    for rel in paths:
        low = rel.lower()
        if low.endswith(RESTRICTED_EXT) or RESTRICTED_PATH.search(rel):
            restricted.append(rel)
            continue
        t = _iter_text(os.path.join(root, rel))
        if t is None:
            continue
        m = SECRET_RE.search(t)
        if m:
            secrets.append("%s: %s..." % (rel, m.group(0)[:12]))
        m = LOCAL_RE.search(t)
        if m:
            local.append(rel)
    return restricted, secrets, local


def cmd_scan(a):
    """Pre-publication scan of every record's file list: restricted file types / paths, secret patterns, local paths."""
    files = all_files(a.data)
    bad = 0
    for name in RECORDS:
        m, _ = verify(a.data, name, files)
        r, sec, loc = scan_files([e["path"] for e in m["files"]], a.data)
        exts = {}
        for e in m["files"]:
            x = os.path.splitext(e["path"])[1].lower() or "(none)"
            if e["path"].endswith(".gz"):
                x = ".%s.gz" % os.path.splitext(e["path"][:-3])[1].lstrip(".")
            exts[x] = exts.get(x, 0) + 1
        print("%-27s %6d files; types: %s" % (name, m["n_files"], ", ".join("%s %d" % kv for kv in sorted(exts.items(), key=lambda kv: -kv[1])[:8])))
        for x in r[:10]:
            print("    RESTRICTED type/path: " + x)
        for x in sec[:10]:
            print("    SECRET pattern: " + x)
        if loc:
            print("    local absolute path in %d file(s), e.g. %s" % (len(loc), ", ".join(loc[:3])))
        bad += len(r) + len(sec)
    if os.path.exists(os.path.join(a.examples, ".git")):  # examples: the files of the pinned commit
        m = json.load(open(manifest_path(a.data, EXAMPLES), encoding="utf-8"))
        tmp = None
        import tempfile
        tmp = tempfile.mkdtemp(prefix="orchid-examples-scan-")
        try:
            zp = os.path.join(tmp, "x.zip")
            git(a.examples, "archive", "--format=zip", "-o", zp, m["commit"])
            with zipfile.ZipFile(zp) as zf:
                zf.extractall(os.path.join(tmp, "t"))
            r, sec, loc = scan_files([e["path"] for e in m["files"]], os.path.join(tmp, "t"))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        print("%-27s %6d files (commit %s)" % (EXAMPLES, m["n_files"], m["commit"][:10]))
        own = [x for x in r if (x.split("/")[1:2] == ["figures"] and x.endswith((".pdf", ".png"))) or x.endswith(".npz")]
        r = [x for x in r if x not in own]   # the examples' own figures and numpy input tables: reviewed once, not a finding
        if own:
            print("    review once (the examples' own figures and .npz input tables, not publisher material): %d files, e.g. %s"
                  % (len(own), ", ".join(own[:3])))
        for x in r[:10]:
            print("    RESTRICTED type/path: " + x)
        for x in sec[:10]:
            print("    SECRET pattern: " + x)
        if loc:
            print("    local absolute path in %d file(s), e.g. %s" % (len(loc), ", ".join(loc[:3])))
        bad += len(r) + len(sec)
    print("SCAN OK (no restricted file type or path, no secret pattern)" if not bad else "SCAN FAIL: %d finding(s)" % bad)
    print("Not checked here: verbatim quotes are capped by a PI decision (see the runbook); licences per file are in LICENSES.md / NOTICE.md.")
    return 1 if bad else 0


# ---------------------------------------------------------------------------------------------- bundle --
def zip_record(data, name, out_dir):
    m, probs = verify(data, name)
    if probs:
        raise SystemExit("%s: %s" % (name, "; ".join(probs[:5])))
    top = bundle_name(name)[:-4] + "/"
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, bundle_name(name))
    sums = "".join("%s  %s\n" % (e["sha256"], e["path"]) for e in m["files"])
    tmp = path + ".tmp"

    def info(arc, stored):
        zi = zipfile.ZipInfo(top + arc, date_time=ZIP_DATE)
        zi.create_system = 3
        zi.external_attr = 0o100644 << 16
        zi.compress_type = zipfile.ZIP_STORED if stored else zipfile.ZIP_DEFLATED
        return zi

    with zipfile.ZipFile(tmp, "w", allowZip64=True) as zf:
        zf.writestr(info("MANIFEST.sha256", False), sums.encode(), compresslevel=6)
        for e in m["files"]:
            zi = info(e["path"], e["path"].endswith(STORED))
            h = hashlib.sha256()
            with open(os.path.join(data, e["path"]), "rb") as src, zf.open(zi, "w", force_zip64=e["bytes"] > 2**31) as dst:
                for b in iter(lambda: src.read(1 << 20), b""):
                    h.update(b)
                    dst.write(b)
            if h.hexdigest() != e["sha256"]:
                raise SystemExit("%s changed while bundling" % e["path"])
    os.replace(tmp, path)
    return path


def zip_examples(a, out_dir):
    """Deterministic zip of the pinned examples commit (git archive: file times = commit time)."""
    m = json.load(open(manifest_path(a.data, EXAMPLES), encoding="utf-8"))
    if dump(examples_manifest(a.examples, m["commit"])) != open(manifest_path(a.data, EXAMPLES), encoding="utf-8").read():
        raise SystemExit("examples manifest does not match its pinned commit: rerun `manifest`")
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, bundle_name(EXAMPLES))
    git(a.examples, "archive", "--format=zip", "--prefix=%s/" % bundle_name(EXAMPLES)[:-4], "-o", path, m["commit"])
    return path


def cmd_bundle(a):
    out = a.out or os.path.join(SITE, "dist")
    names = a.record or list(RECORDS) + [EXAMPLES]
    lines, md5s = [], []
    for name in names:
        if name not in RECORDS and name != EXAMPLES:
            raise SystemExit("unknown record %r (one of %s, %s)" % (name, ", ".join(RECORDS), EXAMPLES))
        p = zip_examples(a, out) if name == EXAMPLES else zip_record(a.data, name, out)
        d = sha256_file(p)
        lines.append("%s  %s\n" % (d, os.path.basename(p)))
        md5s.append("%s  %s\n" % (md5_file(p), os.path.basename(p)))   # Zenodo shows the MD5 of an uploaded file
        print("%-15s %s  %10s  sha256 %s" % (name, os.path.relpath(p, SITE), human(os.path.getsize(p)), d))
    if not a.record:
        with open(os.path.join(out, "SHA256SUMS"), "w") as fh:
            fh.writelines(lines)
        with open(os.path.join(out, "MD5SUMS"), "w") as fh:
            fh.writelines(md5s)
    return 0


# ------------------------------------------------------------------------------------------- writeback --
DOI_RE = re.compile(r"^10\.\d{4,9}/[^\s\"'<>]+$")


def replace_in_tree(root, mapping, dry):
    pat = re.compile("(%s)(?![\\w-])" % "|".join(re.escape(orchid_cite.placeholder(n)) for n in mapping))
    changed = []
    for rel in all_files(root):
        if not rel.endswith(TEXT_EXT) or rel.startswith(SKIP_PREFIX):
            continue
        p = os.path.join(root, rel)
        with open(p, "rb") as fh:
            raw = fh.read()
        if orchid_cite.PLACEHOLDER_PREFIX.encode() not in raw:
            continue
        text = raw.decode("utf-8")
        new = pat.sub(lambda m: mapping[m.group(1)[len(orchid_cite.PLACEHOLDER_PREFIX):]], text)
        if new != text:
            changed.append(rel)
            if not dry:
                with open(p + ".tmp", "wb") as fh:
                    fh.write(new.encode("utf-8"))
                os.replace(p + ".tmp", p)
    return changed


def sync_literature(data):
    """Mirror <site>/data/literature (where the literature builders write) into <orchid-data>/literature."""
    src, dst = os.path.join(SITE, "data", "literature"), os.path.join(data, "literature")
    sf, df = set(all_files(src)), set(all_files(dst))
    n = 0
    for rel in sorted(sf):
        s, d = os.path.join(src, rel), os.path.join(dst, rel)
        if not os.path.exists(d) or open(s, "rb").read() != open(d, "rb").read():
            os.makedirs(os.path.dirname(d), exist_ok=True)
            shutil.copyfile(s, d)
            n += 1
    for rel in sorted(df - sf):
        os.remove(os.path.join(dst, rel))
        n += 1
    print("synced literature/: %d file(s) updated" % n)


def run(cmd):
    print("$ " + " ".join(cmd), flush=True)
    r = subprocess.run(cmd, cwd=SITE)
    return r.returncode


def cmd_writeback(a):
    mapping = {}
    for item in a.doi:
        if "=" not in item:
            raise SystemExit("--doi takes name=DOI, got %r" % item)
        name, doi = item.split("=", 1)
        doi = doi.strip().removeprefix("https://doi.org/").removeprefix("doi:")
        if name not in RECORDS and name != EXAMPLES:
            raise SystemExit("unknown record %r (one of %s, %s)" % (name, ", ".join(RECORDS), EXAMPLES))
        if not DOI_RE.match(doi) or orchid_cite.is_placeholder(doi):
            raise SystemExit("not a DOI: %r" % doi)
        mapping[name] = doi
    if len(set(mapping.values())) != len(mapping):
        raise SystemExit("the same DOI was given to two records")
    dp = os.path.join(a.data, "zenodo", "dois.json")
    dj = json.load(open(dp, encoding="utf-8"))
    for name, doi in mapping.items():
        old = dj["dois"].get(name)
        if old and not orchid_cite.is_placeholder(old) and old != doi:
            raise SystemExit("%s already has DOI %s; edit zenodo/dois.json by hand to change it" % (name, old))
    roots = [a.data, SITE] + ([a.examples] if os.path.exists(os.path.join(a.examples, ".git")) else [])
    for root in roots:  # (the examples repository: README and CITATION.cff carry its own DOI; commit there, then re-pin)
        ch = replace_in_tree(root, mapping, a.dry_run)
        print("%s %d file(s) in %s" % ("would change" if a.dry_run else "replaced placeholders in", len(ch), root))
        for rel in ch[:15]:
            print("    " + rel)
        if len(ch) > 15:
            print("    ... %d more" % (len(ch) - 15))
    if a.dry_run:
        return 0
    dj = json.load(open(dp, encoding="utf-8"))  # (the text pass already rewrote it; make sure)
    dj["dois"].update(mapping)
    with open(dp, "w", encoding="utf-8") as fh:
        fh.write(json.dumps(dj, indent=1, ensure_ascii=False, sort_keys=True) + "\n")
    if a.no_rebuild:
        print("--no-rebuild: run the builders and `zenodo_release.py manifest` before bundling")
        return 0
    py, ppy = sys.executable, a.parquet_python
    steps = [[py, "tools/export_literature_data.py"], [py, "tools/build_dft_parity.py"], "sync",
             [py, "tools/add_cite_as.py", "--data", a.data], [ppy, "tools/build_parquet.py", "--data", a.data],
             [py, "tools/build_static_pages.py", "--data", a.data],
             [py, "tools/build_phase_registry.py", "--data", a.data, "--parquet-python", ppy],
             [ppy, "tools/build_assessed_assemblages.py", "--data", a.data],   # before espei: espei/assessed reads it
             [py, "tools/build_espei_bundle.py", "--data", a.data],
             [ppy, "tools/build_bibliography.py", "build", "--data", a.data], "card", "metadata", "manifest"]
    for st in steps:
        if st == "sync":
            sync_literature(a.data)
        elif st == "card":
            cmd_card(argparse.Namespace(data=a.data, check=False))
        elif st == "metadata":
            cmd_metadata(argparse.Namespace(data=a.data, check=False))
        elif st == "manifest":
            cmd_manifest(argparse.Namespace(data=a.data, check=False, examples=a.examples, examples_ref=a.examples_ref))
        elif run(st):
            raise SystemExit("step failed: " + " ".join(st))
    checks = [[py, "tools/export_literature_data.py", "--check"], [py, "tools/build_dft_parity.py", "--check"],
              [py, "tools/add_cite_as.py", "--check", "--data", a.data],
              [ppy, "tools/build_parquet.py", "--data", a.data, "--check"],
              [py, "tools/build_static_pages.py", "--data", a.data, "--check"],
              [py, "tools/build_phase_registry.py", "--data", a.data, "--parquet-python", ppy, "--check"],
              [py, "tools/build_espei_bundle.py", "--data", a.data, "--check"],
              [ppy, "tools/build_assessed_assemblages.py", "--data", a.data, "--check"],
              [ppy, "tools/build_bibliography.py", "build", "--data", a.data, "--check"],
              [py, "tools/build_precision_sample.py", "--check"],   # verifies the frozen sample; never regenerates it
              [py, "tools/zenodo_release.py", "--data", a.data, "card", "--check"],
              [py, "tools/zenodo_release.py", "--data", a.data, "metadata", "--check"],
              [py, "tools/zenodo_release.py", "--data", a.data, "--examples", a.examples, "manifest", "--check"]]
    fails = [" ".join(c) for c in checks if run(c)]
    left = []
    for root in roots:
        for name in mapping:
            for rel in replace_in_tree(root, {name: "X"}, True):
                left.append("%s: %s" % (os.path.basename(root), rel))
    for x in left:
        print("placeholder still present: " + x)
    if fails or left:
        print("FAIL: %s" % "; ".join(fails + (["placeholders left"] if left else [])))
        return 1
    print("writeback done: %s. Review both repositories (git diff), commit, then bundle and publish." %
          ", ".join("%s=%s" % kv for kv in sorted(mapping.items())))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=orchid_cite.DEFAULT_DATA, help="orchid-data repository (default: %(default)s)")
    ap.add_argument("--examples", default=DEFAULT_EXAMPLES, help="orchid-examples repository (default: %(default)s)")
    ap.add_argument("--examples-ref", default=None, help="manifest: commit or tag of the examples repository to pin (default HEAD)")
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("manifest", help="write zenodo/<name>.files.json (sizes, SHA-256) from the record rules")
    s.add_argument("--check", action="store_true", help="exit 1 if a manifest is stale")
    s = sub.add_parser("metadata", help="(re)write zenodo/<name>.metadata.json and metadata_numbers.json (build_zenodo_metadata.py)")
    s.add_argument("--check", action="store_true", help="exit 1 if a metadata file is stale")
    s.add_argument("--export", metavar="DIR", help="also write one HTML form sheet per record into DIR (for filling the web form)")
    sub.add_parser("scan", help="pre-publication scan of every record's files: restricted types/paths, secret patterns, local paths")
    s = sub.add_parser("card", help="refresh Parquet sizes in DATASET_CARD.md and copy it to hf/README.md")
    s.add_argument("--check", action="store_true", help="exit 1 if the card or its mirror is stale")
    sub.add_parser("plan", help="print records, files and sizes; check presence, hashes and the 50 GB limit")
    s = sub.add_parser("bundle", help="build one deterministic zip per record into dist/")
    s.add_argument("--record", action="append", help="only this record (repeatable)")
    s.add_argument("--out", help="output directory (default: <site>/dist, git-ignored)")
    s = sub.add_parser("writeback", help="replace the DOI placeholders everywhere, rebuild and check")
    s.add_argument("--doi", action="append", required=True, metavar="NAME=DOI", help="e.g. literature=10.5281/zenodo.123")
    s.add_argument("--dry-run", action="store_true", help="list the files that would change; write nothing")
    s.add_argument("--no-rebuild", action="store_true", help="only replace text and update zenodo/dois.json")
    default_ppy = os.path.expanduser("~/orchid-py313/bin/python")
    s.add_argument("--parquet-python", default=default_ppy if os.path.exists(default_ppy) else sys.executable,
                   help="interpreter with pyarrow for tools/build_parquet.py (default: %(default)s)")
    s = sub.add_parser("credit", help="credit accepted contributors (CONTRIBUTORS.md, Zenodo metadata), then refresh manifests")
    s.add_argument("--check", action="store_true", help="exit 1 if the credit files are stale")
    a = ap.parse_args(argv)
    return {"manifest": cmd_manifest, "metadata": cmd_metadata, "scan": cmd_scan, "card": cmd_card, "plan": cmd_plan, "bundle": cmd_bundle, "writeback": cmd_writeback,
            "credit": cmd_credit}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
