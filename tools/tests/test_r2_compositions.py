#!/usr/bin/env python3
"""Audit round 2 (A1: elements, chemical systems, compositions): unit tests of tools/literature_r2_compositions.py and
permanent invariant tests on an exported literature release.

  python3 tools/tests/test_r2_compositions.py                       # unit tests + invariants on the default export
  ORCHID_LIT_JSONL=<path/orchid-literature.jsonl.gz> python3 ...    # invariants on another export
  ORCHID_LIT_JSONL=none python3 ...                                 # unit tests only

The default export is <site>/data/literature/orchid-literature.jsonl.gz (where the exporter writes), else the
orchid-data checkout. One invariant per error class C1-C7 of the A1 audit (ORCHID_AUDIT_2026-10-24/A1.md)."""
import copy
import gzip
import json
import os
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)
import literature_r2_compositions as M  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


N = [0]


def rec(asset="doi:10.1/x/p1/tab-1", system="NI-ZR", comps=(), prop="INVARIANT_TEMPERATURE", basis=None, **kw):
    N[0] += 1
    r = {"id": "ORC-U%04d" % N[0], "system": system, "elements": system.split("-"), "property": prop,
         "record_kind": "invariant_table", "paper": "doi:10.1/x", "flags": {"extraction_accepted": True, "espei_ready": True},
         "locator": {"asset": asset, "page": 1, "row": N[0], "caption": kw.pop("caption", "")},
         "phases": [{"label": p} for p, _ in comps],
         "phase_compositions": [{"phase": p, "composition": c, "composition_source": "printed",
                                 "basis": basis or "mole_fraction (converted from atomic percent)"} for p, c in comps]}
    r.update(kw)
    return r


def vc(key, verdict, **kw):
    d = {"key": key, "verdict": verdict, "page": 1, "image_sha256": "0" * 64, "pdf_sha256_prefix": "0" * 16,
         "seen": "fixture", "checked": "2026-10-24"}
    d.update(kw)
    return d


def codes(r):
    return {g["code"] for g in r.get("gate_failures") or []}


# ------------------------------------------------------------------------------------------------ helpers
check("pure_element parses terminal labels", [M.pure_element(x) for x in ("(αFe)", "βSn", "α(Fe)", "Fcc(Ni)", "Mg-hcp",
                                                                        "L", "Liquid", "fcc_A1")]
      == ["FE", "SN", "FE", "NI", "MG", None, None, None])
check("expected_x of a formula and of a terminal phase",
      abs(M.expected_x("Ti3Pt", "PT", ["PT", "TI"]) - 0.25) < 1e-9 and M.expected_x("(Ti)", "PT", ["PT", "TI"]) == 0.0)

# ------------------------------------------------------------------------------------------------ C1
a = "doi:10.1/c1/p2/tab-1"
t = [rec(asset=a, system="PT-TI", comps=[("L", {"PT": 0.0016}), ("Ti3Pt", {"PT": 0.0025})]),
     rec(asset=a, system="PT-TI", comps=[("βTiPt", {"PT": 0.005}), ("TiPt8", {"PT": 0.00889})])]
check("C1 table detected", a in M.c1_tables(t))
u = copy.deepcopy(t)
rep = M._Rep()
M.c1_double_division(u, {}, rep)
check("C1 unchecked table: flagged, not ESPEI-ready, values unchanged",
      all("composition_double_division_suspect" in codes(r) and not r["flags"]["espei_ready"] for r in u)
      and u[0]["phase_compositions"][1]["composition"] == {"PT": 0.0025})
u = copy.deepcopy(t)
M.c1_double_division(u, {"a1_c1:" + a: vc("a1_c1:" + a, "confirmed")}, M._Rep())
check("C1 confirmed table: x100, stored value kept, correction recorded",
      abs(u[0]["phase_compositions"][1]["composition"]["PT"] - 0.25) < 1e-9
      and u[0]["phase_compositions"][1]["composition_stored"] == {"PT": 0.0025}
      and any(c["rule"] == "fraction_printed_under_atomic_percent_header" for c in u[0]["corrections"]))
u = copy.deepcopy(t)
M.c1_double_division(u, {"a1_c1:" + a: vc("a1_c1:" + a, "rejected")}, M._Rep())
check("C1 page-rejected (genuinely dilute) table: untouched", not u[0].get("corrections") and not u[0].get("gate_failures"))

# ------------------------------------------------------------------------------------------------ C2
a2 = "doi:10.1/c2/p9/tab-7"
t = [rec(asset=a2, system="BA-GE", comps=[("Liq", {"GE": 0.96}), ("Ba", {"GE": 1.0}), ("Ba2Ge", {"GE": 0.667})]),
     rec(asset=a2, system="BA-GE", comps=[("Liq", {"GE": 0.6}), ("Ba2Ge", {"GE": 0.667}), ("Ba5Ge3", {"GE": 0.625})])]
u = copy.deepcopy(t)
M.c2_axis(u, {"a1_c2:" + a2: vc("a1_c2:" + a2, "confirmed")}, M._Rep())
check("C2 confirmed table flipped (Ba at x(Ge)=0, Ba2Ge at 1/3)",
      u[0]["phase_compositions"][1]["composition"] == {"GE": 0.0}
      and abs(u[0]["phase_compositions"][2]["composition"]["GE"] - 0.333) < 1e-9
      and u[0]["flags"].get("axis_header_contradicted"))
u = copy.deepcopy(t)
M.c2_axis(u, {}, M._Rep())
check("C2 unchecked reversed table: flagged, not ESPEI-ready", all("composition_axis_suspect" in codes(r) for r in u))

# ------------------------------------------------------------------------------------------------ C3
r3 = rec(system="MG-NI", comps=[("L", {"NI": 0.291}), ("MgNi2", {"NI": 0.333}), ("Mg2Ni", {"NI": 0.665})])
u = [copy.deepcopy(r3)]
M.c3_permutation(u, {}, M._Rep())
check("C3 swapped stoichiometric values re-assigned (one transposition)",
      u[0]["phase_compositions"][1]["composition"] == {"NI": 0.665}
      and u[0]["phase_compositions"][2]["composition"] == {"NI": 0.333} and u[0]["phase_compositions"][0]["composition"] == {"NI": 0.291})
u = [copy.deepcopy(r3)]
M.c3_permutation(u, {"a1_c3:" + r3["id"]: vc("a1_c3:" + r3["id"], "rejected")}, M._Rep())
check("C3 page-rejected swap keeps the stored order", u[0]["phase_compositions"][1]["composition"] == {"NI": 0.333})
r3b = rec(system="AL-SC", comps=[("L", {"SC": 0.43}), ("ScAl2", {"SC": 0.033}), ("ScAl", {"SC": 0.5})])
u = [r3b]
M.c3_permutation(u, {}, M._Rep())
check("C3 label rejected by its composition and no swap: label_conflicts_composition, not ESPEI-ready",
      "label_conflicts_composition" in codes(u[0]) and not u[0]["flags"]["espei_ready"])

# ------------------------------------------------------------------------------------------------ C4
r4 = rec(system="AL-MG-O", prop="HM_FORM", comps=[], record_kind="general_table", overall_composition={"AL": 0.4, "O": 0.6, "MG": 0.0},
         overall_composition_basis="from_formula", phases=[{"label": "Al2O3"}])
r4b = rec(system="C-CE-LA-PR", comps=[("(βLa)", {"C": 0.05}), ("La2C3", {"C": 0.6}), ("(αLa)", {"C": 0.01})])
r4c = rec(system="AG-CR-ZR", prop="PHASE_COMPOSITION", comps=[], phases=[{"label": "AgZr"}],
          caption="Primary phases of the ternary Ag-Cr-Zr alloys.")
u = [r4, r4b, r4c]
M.c4_state_system(u, M._Rep())
check("C4 binary compound in a ternary table filed under its binary; source_system kept; zero keys dropped",
      u[0]["system"] == "AL-O" and u[0]["elements"] == ["AL", "O"] and u[0]["source_system"] == "AL-MG-O"
      and "MG" not in u[0]["overall_composition"])
check("C4 binary invariant (terminal solutions + compound) in a C-RE table filed under C-LA", u[1]["system"] == "C-LA")
check("C4 phase observed in a ternary alloy (caption does not name the subsystem) stays ternary", u[2]["system"] == "AG-CR-ZR")

# ------------------------------------------------------------------------------------------------ C6 / C7
r6 = rec(system="CE-O", prop="LATTICE_A", comps=[], overall_composition={"CE": 0.36, "F": 0.0, "O": 0.64})
r6b = rec(asset="sha256:85fd/p149/tab-i", system="C-SI", comps=[("Liquid", {"S": 0.0}), ("diamond", {"S": 0.0}), ("SiC", {"S": 0.5})])
u = [r6, r6b]
M.c6_elements(u, {"a1_c6:" + r6b["id"]: vc("a1_c6:" + r6b["id"], "corrected", set={"Liquid": {"SI": 1.0}, "diamond": {"SI": 1.0},
                                                                                     "SiC": {"SI": 0.5}})}, M._Rep())
check("C6 stale zero element (fluorite 'F') dropped", "F" not in u[0]["overall_composition"])
check("C6 element symbol read from the page (S -> SI)", u[1]["phase_compositions"][0]["composition"] == {"SI": 1.0})
r7 = rec(system="CA-MG", comps=[("L", {"MG": 28.6})], basis="unresolved", asset="doi:10.1/c7/p5/tab-ii")
r7b = rec(system="AG-TI", comps=[("AgTi2", {"TI": 667.0})], basis="unresolved")
u = [r7, r7b]
M.c7_ranges(u, {"a1_c7:" + r7["locator"]["asset"]: vc("a1_c7:x", "corrected", factor=0.01)}, M._Rep())
check("C7 table whose header the page shows as mol%: resolved, still accepted",
      abs(u[0]["phase_compositions"][0]["composition"]["MG"] - 0.286) < 1e-9 and u[0]["flags"]["extraction_accepted"])
check("C7 raw percentage not resolved: not accepted", not u[1]["flags"]["extraction_accepted"]
      and "composition_outside_unit_interval" in codes(u[1]))


# ------------------------------------------------------------------------------------------------ invariants
def default_export():
    p = os.environ.get("ORCHID_LIT_JSONL")
    if p:
        return None if p.lower() == "none" else p
    for c in (os.path.join(os.path.dirname(HERE), "data", "literature", "orchid-literature.jsonl.gz"),
              os.path.expanduser("~/projects/orchid-data/literature/orchid-literature.jsonl.gz")):
        if os.path.exists(c):
            return c
    return None


path = default_export()
if path:
    print("invariants on", path)
    R = [json.loads(l) for l in gzip.open(path, "rt", encoding="utf-8")]
    visual = {}
    vf = os.path.join(HERE, "literature_visual_checks.json")
    if os.path.exists(vf):
        visual = {c["key"]: c for c in json.load(open(vf, encoding="utf-8"))["checks"]}
    # C1: no table (not page-rejected as genuinely dilute) whose every at.%-converted composition is <= 0.0101
    bad = [a for a in M.c1_tables(R) if (visual.get("a1_c1:" + a) or {}).get("verdict") != "rejected"]
    check("INV C1 no table of fractions divided by 100 twice (%d)" % len(bad), not bad)
    # C2: no binary table whose terminal / stoichiometric phases mostly fit 1 - x
    T2 = M.c2_tables(R)
    bad = [a for a, (rows, st) in T2.items() if M._c2_majority(rows, st)
           and (visual.get("a1_c2:" + a) or {}).get("verdict") != "rejected"]
    check("INV C2 no table with the composition axis reversed (%d)" % len(bad), not bad)
    bad = [r["id"] for r in R if r["flags"].get("espei_ready")
           and any(M.c2_cell_fit(r, pc, e, x) in ("rev", "c1") for pc, e, x in M.printed_cells(r))]
    check("INV C2/C1 no ESPEI-ready record with a reversed or 1/100 cell (%d)" % len(bad), not bad)
    # C3: no ESPEI-ready binary record with a stoichiometric phase > 0.15 from its formula and no calculated range
    bad = []
    for r in R:
        if not r["flags"].get("espei_ready") or len(r["elements"]) != 2:
            continue
        for pc, e, x in M.printed_cells(r):
            fx = M.expected_x(pc.get("phase"), e, r["elements"]) if len(pc.get("composition") or {}) == 1 else None
            if fx is not None and 0 < fx < 1 and abs(x - fx) > 0.15 and "label_conflicts_composition" not in \
                    {g["code"] for g in r.get("gate_failures") or []}:
                bad.append(r["id"])
    check("INV C3 no ESPEI-ready stoichiometric phase > 0.15 off its formula without a flag (%d)" % len(bad), len(bad) <= 5)
    # C4: compound-intrinsic records: the system is the support of the state
    bad = [r["id"] for r in R if M.state_support(r)[0]]
    check("INV C4 no record filed under a system larger than its own state (%d)" % len(bad), not bad)
    bad = [r["id"] for r in R if r.get("source_system") and not set(r["elements"]) <= set(r["source_system"].split("-"))]
    check("INV C4 the state system is a subsystem of source_system (%d)" % len(bad), not bad)
    bad = [r["id"] for r in R if r.get("chemical_system") != r["system"] or sorted(r["elements"]) != r["system"].split("-")]
    check("INV C4/C5 chemical_system == system == sorted elements (%d)" % len(bad), not bad)
    import ml_folds  # noqa: E402
    bad = [r["id"] for r in R if r.get("fold_lso") != ml_folds.fold(r["system"])]
    check("INV C5 fold_lso recomputed from the corrected system (%d)" % len(bad), not bad)
    # C6: every composition element is a system element
    bad = [r["id"] for r in R if any(e not in r["elements"] for e in (r.get("overall_composition") or {}))
           or any(e not in r["elements"] for pc in r.get("phase_compositions") or [] for e in (pc.get("composition") or {}))]
    check("INV C6 every composition element is in the system (%d)" % len(bad), not bad)
    # C7: accepted mole fractions inside [0, 1]; ESPEI-ready never with an unresolved basis
    bad = [r["id"] for r in R if r["flags"].get("extraction_accepted") and any(
        isinstance(v, (int, float)) and not -1e-9 <= v <= 1 + 1e-9 for pc in r.get("phase_compositions") or []
        if str(pc.get("basis", "")).startswith(("mole_fraction", "unresolved")) for v in (pc.get("composition") or {}).values())]
    check("INV C7 accepted compositions inside [0, 1] (%d)" % len(bad), not bad)
    bad = [r["id"] for r in R if r["flags"].get("espei_ready") and any(
        str(pc.get("basis", "")).startswith("unresolved") for pc in r.get("phase_compositions") or [])]
    check("INV C7 no ESPEI-ready record with an unresolved composition basis (%d)" % len(bad), not bad)

print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
