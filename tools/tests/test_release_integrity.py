#!/usr/bin/env python3
"""Permanent tests for the error classes of audits A5 (figures) and A6 (paper numbers and claims), 2026-10-24.

Unit tests (no release needed) and, when the release is present, invariants of the release:
  1. calculated invariants: phase_diagram_invariants has one row per distinct invariant (tools/calc_invariants.py);
     n_mapper_hits sums to the hits in phase-diagrams/sgte/*.json; index n_invariants agrees
  2. cited-source years: never the trailing page number of the printed reference, never after the citing paper
  3. lattice lengths a, b, c: in angstrom wherever a unit was printed (no nm / pm left)
  4. pressure: no value keeps the assumed 1 atm when the page prints another pressure
  5. running text: every value has quote and quote_page
  6. signs read on the page by gate G3 carry flag_sign_pdf_page_verified
  7. quality_fixes.json: espei.espei_ready equals the released flag count; the gate-level count is labelled
  8. licence statements: no "fair use" / "single sentence" claims in released text files
  9. bibliography counts: zenodo/bibliography.metadata.json and DATASET_CARD.md agree with bibliography/index.json
 10. hull catalogs: the energy corrections applied upstream are stated (MP2020 for the MP GGA_GGA+U snapshot)
 11. quotes: no released quote (literature records, context items, Parquet, any JSON `quote` field of the release)
     is longer than 500 characters; a clipped quote carries flags.quote_clipped, still contains the printed value and
     keeps its quote_page (PI decision 2026-10-24; clipping happens at export, tools/quote_clip.py)
usage: ~/orchid-py313/bin/python tools/tests/test_release_integrity.py [--data ~/projects/orchid-data]"""
import argparse
import glob
import gzip
import json
import os
import re
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import calc_invariants as CI  # noqa: E402
import literature_quality as Q  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--data", default=os.path.expanduser("~/projects/orchid-data"))
A = ap.parse_args()
fails = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (("  " + str(detail)) if detail and not cond else ""))
    if not cond:
        fails.append(name)


# ------------------------------------------------------------------ unit: distinct invariants
P = lambda ph, x: {"phase": ph, "x": x}  # noqa: E731
invs = [{"T": 1000.0, "points": [P("LIQUID", 0.3), P("FCC_A1", 0.1), P("AB", 0.5)]},
        {"T": 1000.04, "points": [P("FCC_A1", 0.1003), P("LIQUID", 0.3002), P("AB", 0.5)]},     # repeat hit
        {"T": 900.0, "points": [P("FCC_A1", 0.2), P("FCC_A1#2", 0.8), P("LIQUID", 0.5)]},
        {"T": 900.0, "points": [P("FCC_A1#2", 0.2), P("FCC_A1", 0.8), P("LIQUID", 0.5)]},       # #2 label swapped
        {"T": 1500.0, "points": [P("AB_B2", 0.5), P("AB_B33", 0.5), P("LIQUID", 0.6)]},
        {"T": 1500.0, "points": [P("LIQUID", 0.4), P("AB_B2", 0.5), P("AB_B33", 0.5)]}]         # other side: distinct
d = CI.distinct(invs)
check("distinct: repeated hit merged, label swap merged, two-sided polymorph kept apart",
      d == [(0, 2), (2, 2), (4, 1), (5, 1)], d)

# ------------------------------------------------------------------ unit: cited years
CASES = [("M. Fitzner, O.J. Kleppa, Metall. Trans. A 24 (1993) 1827-1834.", "1834", 2000, 1993),
         ("Geist et al., Acta Cryst. B 26 (1970) 1800–1802", "1802", 2005, 1970),
         ("Y.K. Rao and B. V. Patil: Metall. Trans., 1971, vol. 2, pp. 1829-35.", "1829", 1987, 1971),
         ("P.S. Herle, and J.D. Corbett, Inorg. Chem., 2001, 40(8), p 1858–1864.", "1864", 2023, 2001),
         ("V. Regnault, Ann. Chim. Phys. 1 (1841) 129.", "1841", 2000, 1841),
         ("E.A. Owen, D.P. Morris, J. Inst. Met. 76 (1949–1950) 145–168.", "1950", 2006, 1950),
         ("Borzone G, Raggio R and Delsante S 2003 Intermetallics 11 1217–22", "2003", 2019, 2003),
         ("K.A. Gschneidner, Metall. Trans. B 1 (1970) 1961.", "1961", 2009, 1970)]
for ref, y, py, want in CASES:
    got, basis = Q.cited_year(ref, y, py)
    check("cited_year %s -> %s" % (ref[:40], want), got == want, (got, basis))
check("cited_year: a database year range (1991-2006) is kept as a year",
      Q.cited_year("PURE 4.4 SGTE Pure Elements (Unary) Database (1991-2006).", "2006", 2010)[0] in (1991, 2006))

# ------------------------------------------------------------------ unit: lattice lengths
r = {"property": "LATTICE_A", "value": 0.4094, "unit": "nm", "uncertainty": {"value": 0.0002, "basis": "printed"}}
Q.normalize_length_unit(r, __import__("collections").Counter())
check("lattice nm -> angstrom, stored value kept", r["unit"] == "Å" and abs(r["value"] - 4.094) < 1e-9
      and r["value_stored"] == 0.4094 and r["unit_stored"] == "nm" and abs(r["uncertainty"]["value"] - 0.002) < 1e-12, r)
r = {"property": "LATTICE_C", "value": 562.6, "unit": "pm"}
Q.normalize_length_unit(r, __import__("collections").Counter())
check("lattice pm -> angstrom", r["unit"] == "Å" and abs(r["value"] - 5.626) < 1e-9, r)
r = {"property": "LATTICE_A", "value": 4.1, "unit": "not_printed"}
Q.normalize_length_unit(r, __import__("collections").Counter())
check("lattice without a printed unit is not converted", r["unit"] == "not_printed" and r["value"] == 4.1, r)

# ------------------------------------------------------------------ release
D = A.data
if not os.path.isdir(os.path.join(D, "literature")):
    print("release not found: unit tests only")
    sys.exit(1 if fails else 0)
try:
    import pandas as pd
except ImportError:
    pd = None

if pd is not None:
    inv = pd.read_parquet(os.path.join(D, "parquet", "phase_diagram_invariants.parquet"))
    hits = 0
    per_file = {}
    for fp in glob.glob(os.path.join(D, "phase-diagrams", "sgte", "*.json")):
        j = json.load(open(fp))
        n = len(j.get("invariants") or [])
        hits += n
        per_file[os.path.basename(fp)] = len(CI.distinct(j.get("invariants") or []))
    check("invariants: n_mapper_hits sums to the hits in the JSON files", int(inv.n_mapper_hits.sum()) == hits,
          (int(inv.n_mapper_hits.sum()), hits))
    check("invariants: one row per distinct invariant of each file",
          inv.groupby("file").size().to_dict() == {k: v for k, v in per_file.items() if v},)
    inv["key"] = [(f, round(t, 1), tuple(sorted(re.sub(r"#\d+$", "", p) for p in ph)),
                   tuple(round(c["x"], 2) for c in sorted(cs, key=lambda c: c["x"])))
                  for f, t, ph, cs in zip(inv.file, inv.T_K, inv.phases, inv.compositions)]
    check("invariants: no two rows share file, T (0.1 K), phases and compositions (0.01)", not inv.key.duplicated().any(),
          int(inv.key.duplicated().sum()))
    pi = pd.read_parquet(os.path.join(D, "parquet", "phase_diagram_index.parquet"))
    check("invariants: phase_diagram_index.n_invariants counts distinct invariants",
          int(pi.n_invariants.sum()) == len(inv), (int(pi.n_invariants.sum()), len(inv)))

papers = json.load(open(os.path.join(D, "literature", "papers.json")))
bad_year, after, n_cs = [], [], 0
for line in gzip.open(os.path.join(D, "literature", "orchid-literature.jsonl.gz"), "rt", encoding="utf-8"):
    rec = json.loads(line)
    py = (papers.get(rec.get("paper")) or {}).get("year")
    for c in rec.get("cited_sources") or []:
        y = c.get("year")
        if y is None:
            continue
        n_cs += 1
        t = c.get("reference") or ""
        ps = Q._page_starts(t)
        allm = [(int(m.group(1)), m.start(1)) for m in re.finditer(r"(?<!\d)(1[6-9]\d\d|20[0-3]\d)(?!\d)", t)]
        years = {v for v, s in allm if s not in ps}
        pages = {v for v, s in allm if s in ps}
        if int(y) in pages and int(y) not in years:
            bad_year.append((rec["id"], y, t[:80]))
        if py and int(y) > int(py):
            after.append((rec["id"], y, py))
    if rec.get("record_kind") == "running_text" and (not rec.get("quote") or rec.get("quote_page") is None):
        fails.append("running text without quote/quote_page: " + rec["id"])
    if rec.get("property") in ("LATTICE_A", "LATTICE_B", "LATTICE_C") and re.match(r"^(nm|pm)\b", rec.get("unit") or ""):
        fails.append("lattice length not in angstrom: " + rec["id"])
    if any(x.get("code") == "pressure_not_1_atm_on_page" for x in rec.get("gate_failures") or []) \
            and rec.get("pressure_Pa") == 101325:
        fails.append("assumed 1 atm kept although the page prints another pressure: " + rec["id"])
    if any((c.get("gate"), c.get("rule")) == ("G3", "sign_read_from_page") for c in rec.get("corrections") or []) \
            and not (rec.get("flags") or {}).get("sign_pdf_page_verified"):
        fails.append("page-read sign without flag_sign_pdf_page_verified: " + rec["id"])
check("cited years: none is a page number of the printed reference (%d cited sources)" % n_cs, not bad_year,
      bad_year[:3])
check("cited years: none after the citing paper", not after, after[:3])
check("running text, lattice units, pressure, page-read signs (record-level checks above)",
      not any(f.startswith(("running text", "lattice", "assumed 1 atm", "page-read")) for f in fails))

q = json.load(open(os.path.join(D, "literature", "quality_fixes.json")))
if pd is not None:
    v = pd.read_parquet(os.path.join(D, "parquet", "literature_values.parquet"), columns=["flag_espei_ready"])
    check("espei_ready: quality_fixes espei.espei_ready = released flag count",
          q["espei"]["espei_ready"] == int(v.flag_espei_ready.sum()), (q["espei"]["espei_ready"], int(v.flag_espei_ready.sum())))
check("espei_ready: the gate-level count is labelled as before the phase rule",
      "before the ESPEI phase rule" in (q["gates"].get("espei_ready_after_note") or ""))
led = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                  "literature_visual_checks.json")))
check("page-check ledger: entries counted as such", q["gates"].get("visual_check_entries") == len(led["checks"]),
      (q["gates"].get("visual_check_entries"), len(led["checks"])))
check("page-check ledger: says who inspected the pages", "AI agent" in led["description"] and "not by a person" in led["description"])

BAN = re.compile(r"fair use|single sentence|one sentence per value|one verbatim sentence", re.I)
for rel in ("LICENSES.md", "NOTICE.md", "llms.txt", "index.md", "README.md", "DATASET_CARD.md", "literature/CITATION.txt",
            "literature/index.json", "literature/SCHEMA.md", "hf/README.md"):
    p = os.path.join(D, rel)
    if os.path.exists(p):
        txt = open(p, encoding="utf-8").read()
        check("licence wording: no fair-use / single-sentence claim in %s" % rel, not BAN.search(txt),
              (BAN.search(txt) or [""])[0])

# review pages committed to the public site repository carry the same licence wording
SITE = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
for p in sorted(glob.glob(os.path.join(SITE, "review", "**", "*.html"), recursive=True)):
    txt = open(p, encoding="utf-8", errors="replace").read()
    check("licence wording: no fair-use claim in %s" % os.path.relpath(p, SITE), not re.search(r"fair use", txt, re.I))

bi = json.load(open(os.path.join(D, "bibliography", "index.json")))["counts"]
nb = sum(1 for line in open(os.path.join(D, "bibliography", "papers.jsonl"), encoding="utf-8") if line.strip())
check("bibliography: index papers = papers.jsonl lines", bi.get("papers", nb) == nb, (bi.get("papers"), nb))
meta = json.dumps(json.load(open(os.path.join(D, "zenodo", "bibliography.metadata.json"))))
card = open(os.path.join(D, "DATASET_CARD.md"), encoding="utf-8").read()
for label, n in (("papers", nb), ("edges", bi["edges"]), ("with DOI", bi["papers_with_doi"])):
    s = "{:,}".format(n)
    check("bibliography %s (%s) stated in zenodo/bibliography.metadata.json" % (label, s), s in meta)
check("DATASET_CARD bibliography_papers row = papers.jsonl", ("| `bibliography_papers` | %s |" % "{:,}".format(nb)) in card)
check("DATASET_CARD bibliography_edges row = edges", ("| `bibliography_edges` | %s |" % "{:,}".format(bi["edges"])) in card)
if pd is not None:
    check("DATASET_CARD phase_diagram_invariants row = Parquet rows",
          ("| `phase_diagram_invariants` | %s |" % "{:,}".format(len(inv))) in card)

for layer in ("dft-mp", "dft-oqmd"):
    pv = json.load(open(os.path.join(D, "hulls", layer, "catalog.json")))["provenance"]
    check("hull catalog %s states the upstream energy corrections" % layer,
          pv.get("correction_model") != "none" and bool(pv.get("upstream_energy_corrections")), pv.get("correction_model"))
pv = json.load(open(os.path.join(D, "hulls", "dft-mp", "catalog.json")))["provenance"]
check("MP GGA_GGA+U hulls: MP2020 named", pv["thermo_type"] != "GGA_GGA+U" or "MP2020" in pv["upstream_energy_corrections"])

# ------------------------------------------------------------------ quotes: at most 500 characters, printed value kept
import quote_clip as QC  # noqa: E402

n_q = n_clip = 0
bad_q, bad_clip = [], []
for line in gzip.open(os.path.join(D, "literature", "orchid-literature.jsonl.gz"), "rt", encoding="utf-8"):
    rec = json.loads(line)
    q = rec.get("quote")
    if not q:
        continue
    n_q += 1
    if len(q) > QC.MAX_QUOTE_CHARS:
        bad_q.append((rec["id"], len(q)))
    fl = rec.get("flags") or {}
    if fl.get("quote_clipped"):
        n_clip += 1
        req = QC.numbers(rec.get("value_printed")) or QC.numbers(rec.get("temperature_printed"))
        if (not all(QC._spans(q, t) for t in req) or rec.get("quote_page") is None
                or fl.get("quote_chars_original", 0) <= QC.MAX_QUOTE_CHARS or not fl.get("quote_original_sha256")):
            bad_clip.append(rec["id"])
    elif len(q) <= QC.MAX_QUOTE_CHARS and fl.get("quote_chars_original"):
        bad_clip.append(rec["id"])
check("quotes: none of %d literature quotes is longer than %d characters" % (n_q, QC.MAX_QUOTE_CHARS), not bad_q, bad_q[:3])
check("quotes: %d clipped quotes carry quote_clipped, keep the printed value and quote_page" % n_clip, not bad_clip, bad_clip[:3])


def long_quotes(o, path=""):
    if isinstance(o, dict):
        for k, v in o.items():
            if k == "quote" and isinstance(v, str) and len(v) > QC.MAX_QUOTE_CHARS:
                yield path, len(v)
            else:
                yield from long_quotes(v, path + "/" + k)
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from long_quotes(v, path)


bad_files = []
for pat in ("literature/contexts.json", "literature/systems/*.json", "literature/dft_parity.json", "espei/**/*.json",
            "assessed-phase-assemblages/**/*.json", "phase-diagrams/**/*.json", "bibliography/*.json"):
    for f in glob.glob(os.path.join(D, pat), recursive=True):
        try:
            doc = json.load(open(f, encoding="utf-8"))
        except ValueError:
            continue
        if next(long_quotes(doc), None):
            bad_files.append(os.path.relpath(f, D))
check("quotes: no JSON file of the release has a `quote` field above %d characters" % QC.MAX_QUOTE_CHARS, not bad_files,
      bad_files[:3])
if pd is not None:
    for tbl, col in (("literature_values", "quote"), ("experimental_context", "quote"), ("dft_parity", "quote")):
        t = pd.read_parquet(os.path.join(D, "parquet", tbl + ".parquet"), columns=[col])[col].dropna()
        check("quotes: Parquet %s.%s has no value above %d characters (max %d)" % (tbl, col, QC.MAX_QUOTE_CHARS,
              int(t.str.len().max())), int(t.str.len().max()) <= QC.MAX_QUOTE_CHARS)

print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
