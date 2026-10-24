#!/usr/bin/env python3
"""Permanent invariant tests for the phase / phase-mapping error classes of the ORCHID audit round 2 (A2 classes 1-8,
A1 C8; 2026-10-24). Unit probes of tools/phase_mapper.py and tools/phase_label_check.py, then invariants of a built
release (phase registry, assessed phase assemblages, literature model phases).

usage: ~/orchid-py313/bin/python tools/tests/test_r2_phases.py [--data ~/projects/orchid-data] [--units-only]
"""
import argparse
import collections
import glob
import gzip
import json
import os
import re
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.dirname(HERE)
sys.path.insert(0, TOOLS)
import phase_mapper as PM  # noqa: E402
import phase_label_check as PLC  # noqa: E402

HOME = os.path.expanduser("~")
ap = argparse.ArgumentParser()
ap.add_argument("--data", default=os.path.join(HOME, "projects", "orchid-data"))
ap.add_argument("--tdb-dir", default=os.path.join(HOME, "projects", "open-tdb", "databases", "sgte-binary"))
ap.add_argument("--units-only", action="store_true", help="skip the release invariants")
A = ap.parse_args()
fails = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + ((" -- " + str(detail)[:300]) if not cond and detail else ""))
    if not cond:
        fails.append(name)


def m(tdb_file, label, info, axis):
    p = os.path.join(A.tdb_dir, tdb_file)
    return PM.map_phase(label, info, PM.tdb_phases(p), PM.tdb_compositions(p, list(axis)), list(axis))


# ------------------------------------------------------------------ class 3: mapper rule defects (a)-(f)
check("(a) Greek table covers theta..omega", all(g in PM.GREEK for g in "θκλμσστχω"))
check("(b) eta selects by name token, not substring of BETA/ZETA",
      "ETA" not in PM.name_tokens("CU6SN5_BETA") and "ETA" in PM.name_tokens("ALCU_ETA_1"))
syn = {"GRAPHITE_A9": {"pearson": "hP4", "constituents": [["C"]]}, "DHCP": {"constituents": [["PR"], ["C", "VA"]]},
       "LIQUID": {}}
check("(c) solution (alphaPr) hP4 never GRAPHITE_A9 (Pr not a constituent)",
      PM.map_phase("(αPr)", {"pearson_symbol": "hP4"}, syn, None, ["C", "PR"])[0] != "GRAPHITE_A9")
check("(c) solution (alphaPr) hP4 -> DHCP holding Pr",
      PM.map_phase("(αPr)", {"pearson_symbol": "hP4"}, syn, None, ["C", "PR"])[0] == "DHCP")
lav = {"C14_LAVES": {"pearson": "hP12"}, "FE2ZR": {"pearson": "cF24", "constituents": [["FE"], ["ZR"]]}}
check("(d) Pearson-only standard name (hP12 -> C14) does not run before the formula",
      PM.map_phase("Fe2Zr", {"pearson_symbol": "hP12"}, lav, None, ["FE", "ZR"])[0] != "C14_LAVES")
check("(e) a '?' label is never mapped", PM.map_phase("?", {}, {"LIQUID": {}, "FCC_A1": {}}, None, ["A", "B"])[0] is None)
check("(f) ION is a liquid for L", PM.map_phase("L", {}, {"ION": {}, "FCC_A1": {}}, None, ["A", "B"])[0] == "ION")
check("liquid / glass TDB phases never match a solid formula",
      PM.map_phase("SiO2", {}, {"SIO2_L": {"comment": "SiO2 liquid and glass"}, "LIQUID": {}}, None, ["O", "SI"])[0] is None)
check("Pearson equality: hR counts x3 (hexagonal vs rhombohedral cell)",
      PM.pearson_equal("hR21", "hR7") and PM.pearson_equal("hR3", "hR9") and PM.pearson_equal("cP2", "oC8") is False)

# ------------------------------------------------------------------ class 2: polymorph collapse (real TDBs)
if os.path.isdir(A.tdb_dir):
    B2 = {"pearson_symbol": "cP2", "strukturbericht": "B2", "prototype": "CsCl"}
    check("beta-PtZr B2 never PTZR_B33", m("PtZr-10Gao.tdb", "βPtZr", B2, ("PT", "ZR"))[0] != "PTZR_B33")
    check("alpha-PtZr B33 -> PTZR_B33 with +pearson",
          m("PtZr-10Gao.tdb", "αPtZr", {"pearson_symbol": "oC8", "strukturbericht": "B33"}, ("PT", "ZR"))
          == ("PTZR_B33", "compound_by_formula+pearson", ["PTZR_B33"]))
    check("beta-MnTe B4 never MNTE_B81",
          m("MnTe-96Mam.tdb", "βMnTe", {"pearson_symbol": "hP4", "strukturbericht": "B4"}, ("MN", "TE"))[0] != "MNTE_B81")
    check("beta-MoSi2 C40 never MOSI2_C11B",
          m("MoSi-00Liu.tdb", "βMoSi2", {"pearson_symbol": "hP9", "strukturbericht": "C40"}, ("MO", "SI"))[0] != "MOSI2_C11B")
    check("beta-Fe2Pu C36 never FE2PU_C15",
          m("FePu-01Kur.tdb", "βFe2Pu", {"pearson_symbol": "hP24", "strukturbericht": "C36"}, ("FE", "PU"))[0] != "FE2PU_C15")
    check("solid SiO2 never SIO2_L", m("SiO-04Sch.tdb", "αSiO2", {"pearson_symbol": "hP9"}, ("O", "SI"))[0] != "SIO2_L")
    check("B20 CoSi still COSI_B20 (Strukturbericht token decides over a misread Pearson)",
          m("CoSi-06Zha.tdb", "CoSi", {"pearson_symbol": "cF8", "strukturbericht": "B20", "prototype": "FeSi"}, ("CO", "SI"))[0] == "COSI_B20")

# ------------------------------------------------------------------ class 1: literature label -> model phase
E = ["FE", "TA"]
check("eps -> BCC_A2 rejected (Fe-Ta Laves)", PLC.check("ε", "BCC_A2", "FE-TA", E, None, []) is not None)
check("beta-Ni3Ge -> BCC_A2 rejected", PLC.check("βNi3Ge", "BCC_A2", "GE-NI", ["GE", "NI"], None, []) is not None)
check("(alphaTi) -> HCP_A3 kept", PLC.check("(αTi)", "HCP_A3", "AL-TI", ["AL", "TI"], None, []) is None)
check("bare beta in Cu-Zr (beta-Zr) -> BCC_A2 kept", PLC.check("β", "BCC_A2", "CU-ZR", ["CU", "ZR"], None, []) is None)
check("L(phi2) -> GRAPHITE_A9 rejected", PLC.check("L(φ2)", "GRAPHITE_A9", "C-U", ["C", "U"], None, []) is not None)
check("phi1 -> GRAPHITE_A9 rejected", PLC.check("φ1", "GRAPHITE_A9", "C-U", ["C", "U"], None, []) is not None)
check("bcc-La -> DHCP rejected", PLC.check("bcc-La", "DHCP", "LA-ZN", ["LA", "ZN"], None, []) is not None)
check("beta(B) -> BCC_A2 rejected (beta-B is rhombohedral)",
      PLC.check("β(B)", "BCC_A2", "B-SC", ["B", "SC"], None, []) is not None)
tic = {"models": [{"sublattices": [{"constituents": ["TI"]}, {"constituents": ["C", "VA"]}]}], "stability": []}
check("TiC -> FCC_A1 (M)(C,VA) kept (SGTE interstitial convention)",
      PLC.check("TiC", "FCC_A1", "C-TI", ["C", "TI"], tic, []) is None)
check("composition outside the model phase range rejects the mapping",
      PLC.check("(Bi)", "RHOMBOHEDRAL_A7", "BI-LU", ["BI", "LU"],
                {"stability": [{"stable": True, "x_element": "LU", "x_min": 0.0, "x_max": 0.01, "tdb": "x.tdb"}]},
                [{"LU": 1.0}]) == "composition_outside_model_phase_range")

if A.units_only:
    print("%d failed" % len(fails) if fails else "all passed")
    sys.exit(1 if fails else 0)

D = A.data
# ------------------------------------------------------------------ registry (classes 4, 5, 6, 7, 8)
reg = {}
for f in sorted(glob.glob(os.path.join(D, "phases", "systems", "*.json"))):
    for p in json.load(open(f, encoding="utf-8"))["phases"]:
        reg[p["id"]] = p
bad = [(p["id"], mo["tdb"]) for p in reg.values() for mo in p.get("models") or []
       if ((mo.get("tdb_definition") or {}).get("parsed") or {}).get("n_sublattices") is not None
       and len(mo.get("sublattices") or []) != mo["tdb_definition"]["parsed"]["n_sublattices"]]
check("class 4: every model has as many sublattices as its PHASE statement", not bad, bad[:5])
bad = [p["id"] for p in reg.values() if any(mo.get("disordered_part_of") and mo["phase_name"] == p["preferred_name"]
                                            for mo in p.get("models") or [])
       and any(not mo.get("disordered_part_of") for mo in p.get("models") or [])]
check("class 7: preferred_name is never a DIS_PART phase when an independent model exists", not bad, bad[:5])
bad = [p["id"] for p in reg.values() if not p.get("models") and not any(n.get("namespace") == "paper" for n in p["names"])]
check("class 8: every registry phase has a TDB model or a printed paper label", not bad, bad[:5])
bad = []
for p in reg.values():
    st = p.get("structure") or {}
    if st.get("strukturbericht") != "A1":
        continue
    X = {c for mo in p["models"] for sl in (mo.get("sublattices") or [])[1:2] for c in sl.get("constituents") or []
         if c in ("C", "N", "H", "O", "B")} & set(p["elements"])
    rows = [r for r in p.get("stability") or [] if r.get("stable") and r.get("x_min") is not None]
    if len(X) == 1 and rows:
        x = X.pop()
        lo = min(r["x_min"] if r["x_element"] == x else 1 - r["x_max"] for r in rows)
        if lo > 0.25:
            bad.append(p["id"])
check("class 5: an FCC_A1 (M)(X,VA) phase stable only as MX is not labelled A1 Cu cF4", not bad, bad[:5])
for pid, sb in (("phase:ni-ti:niti-b19", "B19'"), ("phase:c-ti:fcc-a1", "B1"), ("phase:h-zr:fcc-a1", "C1")):
    if pid in reg:
        check("structure of %s is %s" % (pid, sb), (reg[pid].get("structure") or {}).get("strukturbericht") == sb,
              reg[pid].get("structure"))
st = json.load(open(os.path.join(D, "phases", "structure_types.json"), encoding="utf-8"))["structure_types"]
check("Al4C3 D7_1 is hR7", any(t.get("prototype") == "Al4C3" and t.get("pearson") == "hR7" for t in st))
liq = {("CO-GD", "KAP_LIQ"), ("LU-PB", "ASSOCI"), ("MO-O", "OXIDE_L")}
bad = [p["id"] for p in reg.values() for mo in p["models"] if (p["system"], mo["phase_name"]) in liq
       and (p.get("structure") or {}).get("structure_type") != "st:liquid"]
check("class 8: KAP_LIQ / ASSOCI / OXIDE_L are liquids", not bad, bad)
bad = [p["id"] for p in reg.values() if (p.get("structure") or {}).get("strukturbericht") in ("B2", "L1_2", "L12", "D0_3")
       and any(r.get("stable") and r.get("x_min") is not None and (r["x_min"] < 0.03 or r["x_max"] > 0.97)
               and not r.get("ordered_structure_range_not_resolved") for r in p.get("stability") or [])
       and any(mo.get("disordered_part") for mo in p["models"])]
check("class 6: an ordered order/disorder phase reaching a pure end carries the unresolved-order caveat", not bad, bad[:5])
bad = [p["id"] for p in reg.values() if len({r["tdb"] for r in p.get("stability") or []}) != len(p.get("stability") or [])]
check("class 6: one stability row per assessment", not bad, bad[:5])
bad = [(p["id"], n["label"]) for p in reg.values() for n in p["names"] if n.get("namespace") == "paper"
       and PLC.ELEMENT_MODEL.match(p["preferred_name"]) and PLC.check(n["label"], p["preferred_name"], p["system"],
                                                                         p["elements"], p, []) not in (None,)]
check("class 1/8: no paper alias a label check would reject on an element-structure phase", not bad, bad[:8])

# ------------------------------------------------------------------ Dataset 6 (classes 2, 3, 7; A1 C8)
pa = os.path.join(D, "assessed-phase-assemblages", "records.jsonl.gz")
if os.path.exists(pa):
    recs = [json.loads(l) for l in gzip.open(pa, "rt", encoding="utf-8")]
    check("C8: every point has T_K > 0", all(r["T_K"] > 0 for r in recs))
    bad = [r["id"] for r in recs if "schematic_field_unbounded_below" in r.get("flags", []) and r["T_K"] < 298.15]
    check("C8: no point below room temperature in a schematic field without a lower boundary", not bad, bad[:5])
    bad = [(r["id"], a["label"]) for r in recs for a in r["assemblage"] if "?" in a["label"] and a.get("phase_id")]
    check("class 3e: no '?' label mapped", not bad, bad[:5])
    tp = {}
    bad, groups = [], collections.defaultdict(set)
    for r in recs:
        for a in r["assemblage"]:
            mm = re.search(r"(\S+\.tdb):(\S+) \(([^)]*)\)", a.get("mapping_basis") or "")
            if mm and "+pearson" in mm.group(3) and a.get("pearson_symbol"):
                f = mm.group(1)
                if f not in tp:
                    tp[f] = PM.tdb_phases(os.path.join(A.tdb_dir, f)) if os.path.exists(os.path.join(A.tdb_dir, f)) else {}
                t = tp[f].get(mm.group(2).upper(), {})
                if t.get("pearson") and PM.pearson_equal(t["pearson"], a["pearson_symbol"]) is False:
                    bad.append((r["id"], a["label"], a["pearson_symbol"], t["pearson"]))
            if a.get("phase_id") and a.get("strukturbericht"):
                groups[(r["diagram"], a["phase_id"])].add(re.sub(r"[_\s]", "", a["strukturbericht"]))
    check("class 2: basis '+pearson' only when the printed and TDB Pearson symbols agree", not bad, bad[:5])
    bad = [k for k, v in groups.items() if len(v) >= 2]
    check("class 2: within a diagram, labels with different printed Strukturbericht never share a registry phase",
          not bad, bad[:5])
    regc = {pid: {c for mo in p["models"] for sl in mo.get("sublattices") or [] for c in sl.get("constituents") or []}
            for pid, p in reg.items() if p.get("models")}
    bad = []
    for r in recs:
        for a in r["assemblage"]:
            lab = a["label"]
            if lab.startswith("(") and lab.endswith(")") and a.get("phase_id") in regc:
                el = re.sub("[%s]" % PM.GREEK_CHARS, "", lab[1:-1]).strip(" '")
                if re.fullmatch(r"[A-Z][a-z]?", el) and el.upper() not in regc[a["phase_id"]]:
                    bad.append((r["id"], lab, a["phase_id"]))
    check("class 3c: the element of an '(El)' label is a constituent of the mapped phase", not bad, bad[:5])
    indep = {(p["system"], mo["phase_name"]) for p in reg.values() for mo in p["models"] if not mo.get("disordered_part_of")}
    bad = []
    for f in glob.glob(os.path.join(D, "assessed-phase-assemblages", "espei-zpf", "*.json")):
        d = json.load(open(f, encoding="utf-8"))
        sysn = "-".join(sorted(e for e in d["components"] if e != "VA"))
        bad += [(os.path.basename(f), v[0]) for reg_ in d["values"] for v in reg_
                if v[0] != "__HYPERPLANE__" and (sysn, v[0]) not in indep]
    check("class 7: every espei-zpf phase name is an independent PHASE in a TDB of the system", not bad, bad[:5])

# ------------------------------------------------------------------ literature model phases (class 1)
pl = os.path.join(D, "literature", "orchid-literature.jsonl.gz")
if os.path.exists(pl):
    bad, unattributed = [], 0
    for line in gzip.open(pl, "rt", encoding="utf-8"):
        r = json.loads(line)
        els = [e.upper() for e in r.get("elements") or []]
        for p in r.get("phases") or []:
            mp = str(p.get("model_phase") or "").upper().split("#")[0]
            if (p.get("pearson_symbol") or p.get("prototype")) and p.get("structure_basis") != "mapped_model_phase":
                unattributed += 1
            if not mp or not PLC.ELEMENT_MODEL.match(mp):
                continue
            k = PLC.label_kind(p.get("label"), els)
            if k in ("greek_variant", "compound"):
                bad.append((r["id"], p.get("label"), mp))
            elif k == "greek":
                fams = set().union(*[PLC.POLYMORPH[e].get(PLC.normalize_label(p.get("label")), set())
                                     for e in els if e in PLC.POLYMORPH] or [set()])
                if PLC.family(mp) not in fams:
                    bad.append((r["id"], p.get("label"), mp))
    check("class 1: no Greek / compound label on an element-structure model phase (except element polymorph names)",
          not bad, "%d e.g. %s" % (len(bad), bad[:5]))
    check("class 1: every literature Pearson symbol / prototype says it comes from the mapped model phase",
          unattributed == 0, unattributed)

print("%d failed" % len(fails) if fails else "all passed")
sys.exit(1 if fails else 0)
