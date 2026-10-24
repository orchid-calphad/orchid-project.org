#!/usr/bin/env python3
"""Permanent invariants from the ORCHID audit round 2, A4 (DFT, hulls, parity; 2026-10-24) and A1 C5 (one ML fold per
record in literature_values and dft_parity). Unit checks always run; release checks run when the files exist.

  E1  a point whose record names a structure is paired with an entry of that space group whenever the snapshot has
      one at the composition (else structure_unmatched, never in the default subset)
  E2  ehull = Ef - hull(x) of ORCHID's released hull (|d| < 1e-3 kJ/mol-atom), for MP and OQMD
  E4  'NiSi_2' parses as NiSi2
  E5  a measured composition more than 0.02 from the compound keeps the point out of the default subset
  E6  experimental values are not stored at 0 K
  E7  catalogs name the provider's correction model and source_release; liquid/amorphous labels are not paired;
      non-integer formula strings name the integer compound
  E8  every calibration summary carries a system-cluster bootstrap CI
  C5  dft_parity chemical_system / fold_lso / fold_lpo equal the literature record's

usage: ~/orchid-py313/bin/python tools/tests/test_r2_dft.py [--data ~/projects/orchid-data] [--open-tdb ~/projects/open-tdb]
"""
import argparse
import gzip
import json
import os
import sys
from fractions import Fraction

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import build_dft_parity as P  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--data", default=os.path.expanduser("~/projects/orchid-data"))
ap.add_argument("--open-tdb", default=os.path.expanduser("~/projects/open-tdb"))
ap.add_argument("--parity", default=None, help="dft_parity.json to test (default: <data>/literature/dft_parity.json)")
ap.add_argument("--hulls", default=None, help="hulls root (default: <data>/hulls)")
A = ap.parse_args()
fails = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + (("  " + str(detail)[:300]) if not cond and detail else ""))
    if not cond:
        fails.append(name)


# ------------------------------------------------------------------ unit checks
c, _ = P.parse_phase("αNiSi_2", ["Ni", "Si"])
check("E4 'αNiSi_2' parses as NiSi2", c == {"Ni": 1, "Si": 2}, c)
c, _ = P.parse_phase("NiSi_{2}", ["Ni", "Si"])
check("E4 'NiSi_{2}' parses as NiSi2", c == {"Ni": 1, "Si": 2}, c)
c, why = P.parse_phase("NiTi_B2", ["Ni", "Ti"])
check("E4 underscore before a designator still separates it", c == {"Ni": 1, "Ti": 1}, (c, why))
c, why = P.parse_phase("CaO(liq-am)", ["Ca", "O"])
check("E7 CaO(liq-am) is not paired with a crystalline DFT entry", c is None and "liquid" in why, why)
check("E7 U0.5C0.5 is the compound UC", P.formula_str({"U": Fraction(1, 2), "C": Fraction(1, 2)}) == "UC")
check("E7 Si0.3333Ta0.6667 is SiTa2", P.formula_str({"Si": Fraction("0.3333"), "Ta": Fraction("0.6667")}) == "SiTa2",
      P.formula_str({"Si": Fraction("0.3333"), "Ta": Fraction("0.6667")}))
D = P.structure_designators
check("E1 'CuZr-B2' names B2", set(D(["CuZr-B2"], [], {"Cu": 1, "Zr": 1})) == {"B2"})
check("E1 \"TiNi B19'\" names B19'", set(D(["TiNi B19'"], [], {"Ni": 1, "Ti": 1})) == {"B19'"})
check("E1 'Al3Li(D03)' names D03", set(D(["Al3Li(D03)"], [], {"Al": 3, "Li": 1})) == {"D03"})
check("E1 'MgZn2-type' names C14", set(D(["Fe2Ti (MgZn2-type)"], [], {"Fe": 2, "Ti": 1})) == {"C14"})
check("E1 model phase C14_LAVES names C14", set(D(["Fe2Ti"], ["C14_LAVES"], {"Fe": 2, "Ti": 1})) == {"C14"})
check("E1 model phase of another compound (NISI_B31 for NiSi2) names nothing",
      not D(["αNiSi2"], ["NISI_B31"], {"Ni": 1, "Si": 2}))
check("E1 a formula token is not a designator (B2 in a B system)", not D(["TaB2", "B2"], [], {"Ta": 1, "B": 2}))
check("E1 solution designators A1/A2/A3 are ignored", not D(["Cu3Sn"], ["BCC_A2"], {"Cu": 3, "Sn": 1}))
check("E1 B2 -> Pm-3m (221), B19' -> P2_1/m (11), C14 -> 194, C15 -> 227",
      P.STRUKTUR_SG["B2"] == {221} and P.STRUKTUR_SG["B19'"] == {11} and P.STRUKTUR_SG["C14"] == {194}
      and P.STRUKTUR_SG["C15"] == {227})

# ------------------------------------------------------------------ release checks
par_path = A.parity or os.path.join(A.data, "literature", "dft_parity.json")
hull_root = A.hulls or os.path.join(A.data, "hulls")
if os.path.exists(par_path):
    dp = json.load(open(par_path, encoding="utf-8"))
    R, pts = dp["records"], dp["points"]
    check("parity points carry smatch / default / ehull_upstream",
          all(k in p for p in pts for k in ("smatch", "default", "ehull_upstream", "excluded_by")))
    # E1 against the snapshots
    snaps = {}
    for src, rel in (("mp", "hulls/snapshot/snapshot.json"), ("oqmd", "hulls/snapshot-oqmd/snapshot.json")):
        p = os.path.join(A.open_tdb, rel)
        if os.path.exists(p):
            idx = {}
            for e in json.load(open(p, encoding="utf-8"))["entries"]:
                if e.get("deprecated"):
                    continue
                comp = {k: float(v) for k, v in e["composition_reduced"].items()}
                tot = sum(comp.values())
                idx.setdefault(tuple(sorted((k, round(v / tot, 4)) for k, v in comp.items())), []).append(e)
            snaps[src] = idx
    if snaps:
        bad_sg, bad_unm = [], []
        for p in pts:
            r = R[p["rid"]]
            if not r.get("designators") or p["src"] not in snaps:
                continue
            sgs = set().union(*(P.STRUKTUR_SG[d] for d in r["designators"]))
            cands = snaps[p["src"]].get(tuple(sorted((k, round(v, 4)) for k, v in r["comp"].items()))) or []
            have = [e for e in cands if e.get("spacegroup_number") in sgs]
            match = next((e for e in cands if e["material_id"] == p["mid"]), None)
            if have and (match is None or match.get("spacegroup_number") not in sgs):
                bad_sg.append((p["rid"], p["src"], p["mid"]))
            if not have and (p["smatch"] != "structure_unmatched" or p["default"]):
                bad_unm.append((p["rid"], p["src"]))
        check("E1 sg(match) in SG(designator) whenever such an entry exists", not bad_sg, bad_sg[:5])
        check("E1 no same-SG entry -> structure_unmatched and not default", not bad_unm, bad_unm[:5])
    # E2 against the released hulls
    if os.path.isdir(os.path.join(hull_root, "dft-mp")):
        cache, bad = {}, []
        for p in pts:
            comp = R[p["rid"]]["comp"]
            sid = "-".join(sorted(e.upper() for e in comp))
            key = (p["src"], sid)
            if key not in cache:
                f = os.path.join(hull_root, "dft-" + p["src"], sid.split("-")[0], sid + ".json")
                cache[key] = json.load(open(f, encoding="utf-8")) if os.path.exists(f) else None
            d = cache[key]
            if d is None:
                bad.append((p["rid"], p["src"], "no hull"))
                continue
            x = comp.get(d["x_element"].capitalize(), 0.0)
            hs = [it["intercept"] + it["slope"] * x for it in d["isotherms"][0]["intervals"]
                  if it["x_lo"] - 1e-3 <= x <= it["x_hi"] + 1e-3]
            eh = (p["dft_eV"] - min(hs)) * P.KJ_PER_EV
            # x in `comp` is rounded to 4 decimals: allow the slope x rounding term
            slope = max(abs(it["slope"]) for it in d["isotherms"][0]["intervals"]) * P.KJ_PER_EV
            if p["ehull"] is None or abs(eh - p["ehull"]) > 1e-3 + slope * 5e-5:
                bad.append((p["rid"], p["src"], p["ehull"], round(eh, 4)))
        check("E2 ehull = Ef - hull(x) of the released hull (MP and OQMD)", not bad, bad[:5])
    # E5 / E6 / default
    off = [p["rid"] for p in pts if p.get("default") and (R[p["rid"]].get("dx_measured") or 0) > P.OFF_STOICH_MAX]
    check("E5 no default point measured more than 0.02 off the compound", not off, off[:5])
    zero = sorted({rid for rid, r in R.items() if r["ev"] == "experimental" and r["T"] == 0})
    check("E6 experimental parity records are not stored at 0 K", not zero, "%d e.g. %s" % (len(zero), zero[:5]))
    check("E6 an experimental 0 K value is never in the default subset",
          not [p for p in pts if p.get("default") and R[p["rid"]]["ev"] == "experimental" and R[p["rid"]]["T"] == 0])
    # C5 against the literature records
    lit = os.path.join(os.path.dirname(par_path), "orchid-literature.jsonl.gz")
    if not os.path.exists(lit):
        lit = os.path.join(A.data, "literature", "orchid-literature.jsonl.gz")
    if os.path.exists(lit):
        L = {}
        with gzip.open(lit, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                if r["id"] in R:
                    L[r["id"]] = (r.get("chemical_system"), r.get("fold_lso"), r.get("fold_lpo"))
        bad = [rid for rid, r in R.items() if rid in L and (r.get("chemical_system"), r.get("fold_lso"),
                                                            r.get("fold_lpo")) != L[rid]]
        check("C5 dft_parity chemical_system / fold_lso / fold_lpo = literature record's", not bad and len(L) == len(R),
              "%d differ, %d missing; e.g. %s" % (len(bad), len(R) - len(L), bad[:3]))
pq = os.path.join(A.data, "parquet")
if os.path.exists(os.path.join(pq, "dft_parity.parquet")) and os.path.exists(os.path.join(pq, "literature_values.parquet")):
    try:
        import pyarrow.parquet as ppq
        a = ppq.read_table(os.path.join(pq, "dft_parity.parquet"), columns=["record_id", "fold_lso", "chemical_system"]).to_pylist()
        b = {r["id"]: (r["fold_lso"], r["chemical_system"]) for r in ppq.read_table(
            os.path.join(pq, "literature_values.parquet"), columns=["id", "fold_lso", "chemical_system"]).to_pylist()}
        bad = [r["record_id"] for r in a if b.get(r["record_id"]) != (r["fold_lso"], r["chemical_system"])]
        check("C5 Parquet dft_parity fold_lso / chemical_system = literature_values (joined on record_id)", not bad,
              "%d e.g. %s" % (len(bad), bad[:3]))
    except ImportError:
        print("SKIP Parquet checks (pyarrow missing)")
# E7 catalogs
for layer in ("dft-mp", "dft-oqmd"):
    f = os.path.join(hull_root, layer, "catalog.json")
    if os.path.exists(f):
        pv = json.load(open(f, encoding="utf-8"))["provenance"]
        check("E7 %s catalog correction_model names the provider scheme" % layer,
              pv.get("correction_model") not in (None, "none") and (
                  pv["correction_model"].startswith("MP2020") if layer == "dft-mp" else
                  pv["correction_model"].startswith("OQMD fit-partial")), pv.get("correction_model"))
        check("E7 %s provenance uses source_release (no mp_release)" % layer,
              "source_release" in pv and "mp_release" not in pv)
# E8
cal = os.path.join(A.data, "espei", "dft", "CALIBRATION.json")
if os.path.exists(cal):
    C = json.load(open(cal, encoding="utf-8"))
    sets = [s for src in C["sources"].values() for s in [src["all"]] + list(src["strata"].values()) if s]
    check("E8 every calibration set carries ci95_cluster_system", sets and all("ci95_cluster_system" in s for s in sets))

print("%d failed" % len(fails))
sys.exit(1 if fails else 0)
