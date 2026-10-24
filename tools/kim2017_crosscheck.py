#!/usr/bin/env python3
"""Cross-check ORCHID experimental formation enthalpies against Kim et al. 2017 (data/literature/kim2017_crosscheck.*).

Reference: G. Kim, S. V. Meschel, P. Nash, W. Chen, Experimental formation enthalpies for intermetallic phases and
other inorganic compounds, Scientific Data 4 (2017) 170162, doi:10.1038/sdata.2017.162 (1,276 calorimetric values;
Figshare dataset doi:10.6084/m9.figshare.7289723, MIT).

The reference table is read from a LOCAL copy only (no download). Looked for, in order: --kim-csv, then the
paths in CANDIDATES. When none exists the output records that and the search, and no comparison is made.

Comparison: ORCHID records with property HM_FORM, unit kJ/mol-atom, evidence class experimental, extraction
accepted, not unsigned-in-print; formula parsed from the phase label as in tools/build_dft_parity.py; reduced
composition (4-decimal mole fractions) matched to the reference formula. Per formula, the median of the ORCHID
values and the median of the reference values are compared: n overlap, MAE, median |difference|, sign agreement,
and every formula differing by more than 10 kJ/mol-atom (with the ORCHID record ids). Deterministic; standard
library only; no network.

  python3 tools/kim2017_crosscheck.py [--kim-csv PATH]
  python3 tools/kim2017_crosscheck.py --check
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import os
import re
import sys

sys.dont_write_bytecode = True
SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LIT = os.path.join(SITE, "data", "literature")
OUT_JSON = os.path.join(LIT, "kim2017_crosscheck.json")
OUT_MD = os.path.join(LIT, "kim2017_crosscheck.md")
HOME = os.path.expanduser("~")
CANDIDATES = [  # relative to $HOME; recorded in the output so the search is reproducible
    "projects/open-tdb/build/kim2017/formation_enthalpy_expt.csv",
    "projects/open-tdb/examples/thermo-extraction-benchmark/real-paper/formation_enthalpy_expt.csv",
    "data/kim2017/formation_enthalpy_expt.csv",
    "data/PD-extraction/datasets/kim2017/formation_enthalpy_expt.csv",
    "data/kim2017-expt-formation-enthalpy/formation_enthalpy_expt.csv",
]
SEARCH_NOTE = ("Searched on 2026-10-22 without a match: the CANDIDATES paths; every file named *formation_enthalpy* "
               "or *expt_formation* under $HOME (excluding ~/Library; only an unrelated HEACalculator pair-enthalpy "
               "table was found); every *.csv / *.csv.gz / *.json.gz under ~/projects and ~/data (depth 6) for the "
               "SHA-256 prefix 694ba240c2b45da9 of formation_enthalpy_expt.csv recorded in open-tdb's benchmark "
               "intake (examples/thermo-extraction-benchmark/real-paper/kim2018-figshare-formation-enthalpy-intake-v1.json); "
               "no matminer dataset cache. The intake names the source file but the bytes are not on this machine.")
DIFF = 10.0
REF = {"citation": "G. Kim, S. V. Meschel, P. Nash, W. Chen, Experimental formation enthalpies for intermetallic phases "
                   "and other inorganic compounds, Scientific Data 4 (2017) 170162.",
       "doi": "10.1038/sdata.2017.162", "dataset_doi": "10.6084/m9.figshare.7289723",
       "dataset_file": "formation_enthalpy_expt.csv", "dataset_sha256_expected": "694ba240c2b45da9a11b772cb2be30554941f9e04ee1a51f783416f4cd6fe240",
       "license": "MIT (Figshare dataset)", "n_values_reported": 1276}


def load_parity():
    spec = importlib.util.spec_from_file_location("build_dft_parity_k", os.path.join(SITE, "tools", "build_dft_parity.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


def read_kim(path, P):
    """-> {reduced key: [kJ/mol-atom values]}, column info. Tolerant to the column naming."""
    raw = open(path, "rb").read()
    text = (gzip.decompress(raw) if path.endswith(".gz") else raw).decode("utf-8-sig")
    rows = list(csv.reader(io.StringIO(text)))
    head = [h.strip() for h in rows[0]]
    low = [h.lower() for h in head]
    fcol = next(i for i, h in enumerate(low) if "formula" in h)
    ecol = next(i for i, h in enumerate(low) if re.search(r"exp", h) and re.search(r"form|enth|h_?f|e_?f", h))
    unit = "ev" if "ev" in low[ecol] else "kj"
    if unit == "kj":
        # Kim et al. 2017 publish e_form in eV/atom without saying so in the header; formation enthalpies of
        # compounds never fall within +-10 kJ/mol-atom as a whole column, so a column bounded by 10 is eV/atom.
        vals = [float(r[ecol]) for r in rows[1:] if len(r) > ecol and re.fullmatch(r"\s*-?[0-9.eE+-]+\s*", r[ecol] or "")]
        if vals and max(abs(v) for v in vals) < 10:
            unit = "ev"
    out = {}
    for row in rows[1:]:
        if len(row) <= max(fcol, ecol):
            continue
        try:
            v = float(row[ecol])
        except ValueError:
            continue
        c, _ = P.parse_phase(row[fcol].strip(), None)
        if not c:
            continue
        out.setdefault(P.reduced_key(c), []).append(v * 96.485 if unit == "ev" else v)
    return out, {"formula_column": head[fcol], "value_column": head[ecol], "unit_read": unit,
                 "rows": len(rows) - 1, "sha256": hashlib.sha256(raw).hexdigest()}


def build(kim_csv):
    P = load_parity()
    path = None
    for c in ([kim_csv] if kim_csv else []) + [os.path.join(HOME, x) for x in CANDIDATES]:
        if c and os.path.exists(c):
            path = c
            break
    lit_path = os.path.join(LIT, "orchid-literature.jsonl.gz")
    lit_sha = hashlib.sha256(open(lit_path, "rb").read()).hexdigest()
    out = {"schema": "orchid-kim2017-crosscheck-v1", "reference": REF, "unit": "kJ/mol-atom",
           "inputs": {"literature": {"path": "data/literature/orchid-literature.jsonl.gz", "sha256": lit_sha}},
           "candidates_searched": ["~/" + x for x in CANDIDATES], "search_note": SEARCH_NOTE}
    ours = {}
    with gzip.open(lit_path, "rt", encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            fl = r.get("flags") or {}
            if (r.get("property") != "HM_FORM" or r.get("unit") != "kJ/mol-atom" or r.get("evidence_class") != "experimental"
                    or not fl.get("extraction_accepted") or fl.get("unsigned_in_print")
                    or not isinstance(r.get("value"), (int, float))):
                continue
            sysel = [e.capitalize() for e in r.get("elements") or []]
            found = {}
            for lab in [p.get("label") for p in r.get("phases") or [] if p.get("label")]:
                c, _ = P.parse_phase(lab, sysel)
                if c:
                    found[P.reduced_key(c)] = c
            if len(found) != 1:
                continue
            k = next(iter(found))
            ours.setdefault(k, []).append((r["value"], r["id"], P.formula_str(found[k])))
    out["orchid_candidates"] = {"formulas": len(ours), "records": sum(len(v) for v in ours.values())}
    if path is None:
        out["status"] = "reference_not_available_locally"
        out["result"] = None
        return out
    kim, meta = read_kim(path, P)
    out["status"] = "compared"
    out["reference_file"] = dict(meta, path=path.replace(HOME, "~"),
                                 complete=meta["sha256"] == REF["dataset_sha256_expected"])
    pairs, big = [], []
    for k in sorted(set(ours) & set(kim)):
        o = median([v for v, _, _ in ours[k]])
        ref = median(kim[k])
        pairs.append((o, ref))
        if abs(o - ref) > DIFF:
            big.append({"formula": ours[k][0][2], "orchid_median": round(o, 3), "kim_median": round(ref, 3),
                        "difference": round(o - ref, 3), "orchid_ids": sorted(i for _, i, _ in ours[k]),
                        "kim_values": sorted(round(x, 3) for x in kim[k])})
    n = len(pairs)
    d = [o - k for o, k in pairs]
    out["result"] = {"n_overlap_formulas": n, "mae": round(sum(abs(x) for x in d) / n, 3) if n else None,
                     "median_abs_difference": round(median([abs(x) for x in d]), 3) if n else None,
                     "sign_agreement": round(sum(1 for o, k in pairs if (o < 0) == (k < 0)) / n, 4) if n else None,
                     "n_differ_more_than_10": len(big),
                     "differ_more_than_10": sorted(big, key=lambda x: (-abs(x["difference"]), x["formula"]))}
    return out


def md(d):
    L = ["# Cross-check with Kim et al. 2017", "",
         "Reference: %s doi:%s (Figshare dataset doi:%s, %s)." % (d["reference"]["citation"], d["reference"]["doi"],
                                                                  d["reference"]["dataset_doi"], d["reference"]["license"]), "",
         "ORCHID candidates: %d experimental, accepted, signed HM_FORM records (kJ/mol-atom) over %d formulas."
         % (d["orchid_candidates"]["records"], d["orchid_candidates"]["formulas"]), ""]
    if d["status"] != "compared":
        L += ["**Status: not compared.** No local copy of the reference table (`%s`) was found and nothing was "
              "downloaded (no-network policy). %s" % (d["reference"]["dataset_file"], d["search_note"]), "",
              "To run the comparison, place the Figshare file at one of these paths (or pass `--kim-csv`) and rerun "
              "`python3 tools/kim2017_crosscheck.py`:", ""] + ["* `%s`" % c for c in d["candidates_searched"]]
        return "\n".join(L) + "\n"
    r = d["result"]
    rf = d["reference_file"]
    L += ["Reference file: `%s` (%d rows; %s)." % (rf["path"], rf["rows"], "complete file, SHA-256 matches"
                                                     if rf["complete"] else "NOT the complete file: a slice or another version"), "",
          "| n formulas in common | MAE | median abs diff | sign agreement | differ > 10 |", "|---:|---:|---:|---:|---:|",
          "| %d | %s | %s | %s | %d |" % (r["n_overlap_formulas"], r["mae"], r["median_abs_difference"], r["sign_agreement"],
                                       r["n_differ_more_than_10"]), "",
          "| formula | ORCHID median | Kim median | difference | ORCHID ids |", "|---|---:|---:|---:|---|"]
    for x in r["differ_more_than_10"]:
        L.append("| %s | %s | %s | %s | %s |" % (x["formula"], x["orchid_median"], x["kim_median"], x["difference"],
                                                 ", ".join(x["orchid_ids"][:6])))
    return "\n".join(L) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kim-csv", help="local copy of formation_enthalpy_expt.csv (never downloaded)")
    ap.add_argument("--check", action="store_true", help="exit 1 if the files on disk are stale")
    a = ap.parse_args(argv)
    d = build(a.kim_csv)
    js = (json.dumps(d, sort_keys=True, ensure_ascii=False, indent=1) + "\n").encode()
    mdb = md(d).encode()
    if a.check:
        bad = [p for p, b in ((OUT_JSON, js), (OUT_MD, mdb)) if not os.path.exists(p) or open(p, "rb").read() != b]
        if bad:
            print("stale: " + ", ".join(os.path.relpath(p, SITE) for p in bad), file=sys.stderr)
            return 1
        print("kim2017 cross-check is current (%s)" % d["status"])
        return 0
    for p, b in ((OUT_JSON, js), (OUT_MD, mdb)):
        with open(p, "wb") as fh:
            fh.write(b)
    print(json.dumps({"status": d["status"], "orchid_candidates": d["orchid_candidates"], "result": d["result"] and
                      {k: v for k, v in d["result"].items() if k != "differ_more_than_10"}}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
