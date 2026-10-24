#!/usr/bin/env python3
"""Audit round 3 (2026-10-24; POSTFIX.md re-audit + PI additions): permanent tests, one per error class.

Unit part (always): the rule helpers. Release part: invariants on the exported literature records (ORCHID_LITERATURE,
default the orchid-data checkout's literature/orchid-literature.jsonl.gz), on literature/dft_parity.json
(ORCHID_PARITY) and on literature/papers.json.
usage: python3 tools/tests/test_r3.py"""
import gzip
import json
import os
import re
import sys

sys.dont_write_bytecode = True
TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, TOOLS)
import build_dft_parity as P  # noqa: E402
import literature_quality as Q  # noqa: E402
import literature_r2_compositions as C  # noqa: E402
import literature_r3_invariants as I  # noqa: E402
import literature_r3_thermo as T  # noqa: E402
import orchid_cite  # noqa: E402
import uncertainty_budget as U  # noqa: E402

fails = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + ("" if cond or not detail else "  [%s]" % (detail,)))
    if not cond:
        fails.append(name)


# ------------------------------------------------------------------------------------------------ unit part
mgb = {"value": -0.01855, "value_printed": -18.55, "unit_printed": "J/mole of atoms", "unit": "kJ/mol-atom",
       "locator": {"caption": "Heat of formation of the compounds (k J/mole of atoms)"}}
check("A4-N1 unit prefix: 'k J' read as J is detected (value = printed / 1000)", T.unit_prefix_lost(mgb) == 1000.0)
alta = dict(mgb, value=-29.954, value_printed=-29954, unit_printed="J mol^-1 per mol-atom",
            locator={"caption": "Standard properties of the Al-Ta phases (kJ/mol)"})
check("A4-N1 unit prefix: a value printed in J under a kJ caption is not this class", T.unit_prefix_lost(alta) is None)

base = {"elements": ["DY", "GD", "LA", "ND", "RU"], "property": "HM_FORM", "overall_composition_basis": "from_formula"}
r = dict(base, phases=[{"label": "Nd3Ru"}], overall_composition={"GD": 0.0, "LA": 0.0, "ND": 1.0, "RU": 0.0})
sup, _ = C.state_support(r)
check("A1-N1 C4: a corrupt overall composition (Nd = 1) cannot narrow Nd3Ru to a unary", sup == {"ND", "RU"}, sup)
r = {"elements": ["CE", "RU", "SM"], "property": "HM_FORM", "phases": [{"label": "Ln3Ru"}],
     "overall_composition": {"RU": 0.0, "SM": 1.0}, "overall_composition_basis": "from_formula"}
sup, basis = C.state_support(r)
check("A1-N1 C4: a generic label (Ln3Ru) is read with the row's element (Sm3Ru -> RU-SM)", sup == {"RU", "SM"}, (sup, basis))
r = {"elements": ["BI", "CA", "CU", "O", "SR"], "property": "TRANSITION_TEMPERATURE",
     "phases": [{"label": "014x24"}, {"label": "liquid"}], "overall_composition": {"CA": 0, "CU": 0, "O": 1.0, "SR": 0},
     "overall_composition_basis": "from_formula"}
check("A1-N1 C4: no unary oxygen state for an oxide phase", C.state_support(r)[0] is None)
r = {"elements": ["AL", "CO", "FE"], "property": "HM_FORM", "phases": [{"label": "Fe3Al"}], "overall_composition": {"AL": 0.25}}
check("A1 C4 residual: a partial printed composition consistent with the formula moves Fe3Al to AL-FE",
      C.state_support(r)[0] == {"AL", "FE"})
recs = [{"id": "x", "system": "NI-SI", "property": "HM_FORM", "elements": ["NI", "SI"], "phases": [{"label": "Ni5Si2"}],
         "overall_composition": {"NI": 0.75, "SI": 0.25}, "overall_composition_basis": "from_formula", "flags": {}}]
C.r3_overall_composition(recs, C._Rep())
check("A1-N2 split row: a from-formula composition is the printed formula's (Ni5Si2 Si = 0.2857)",
      abs(recs[0]["overall_composition"]["SI"] - 2 / 7) < 1e-4)
c, _ = P.parse_phase("PuO1,52", None)
check("A3-N6 comma decimal: PuO1,52 has 2.52 atoms", c and abs(float(sum(c.values())) - 2.52) < 1e-9, c)
check("A3-N2 normalised formula: Ru0.4Sn0.6 is per atom", U.label_is_normalised(
    {"elements": ["RU", "SN"], "phases": [{"label": "Ru0.4Sn0.6"}]}) and not U.label_is_normalised(
    {"elements": ["RU", "SN"], "phases": [{"label": "RuSn2"}]}))
check("A3 sec.10: a = 1 of the pure component has no reading uncertainty", U.reading_component(
    {"property": "ACR_SN", "value": 1.0, "value_printed": 1, "overall_composition": {"RU": 0.0, "SN": 1.0}},
    "value")["value"] == 0.0)
check("A3-N4 EMF intercept / slope labels", bool(T.EMF_INTERCEPT.match("E°")) and bool(T.EMF_INTERCEPT.match("a"))
      and bool(T.EMF_SLOPE.match("E = a + bT: b")) and not T.EMF_INTERCEPT.match("EMF"))
check("A3 sec.6 three-term coefficients", bool(T.EXPR3_A.search("ΔG° (A of A + B·T + C·T lnT)")) and
      bool(T.EXPR3_A.search("ΔmixGAl = x²(A + Bx) + 8.314T ln(1-x): A")))
check("A3-N7 heat content base temperature", T.HBASE.search("HT-H700").group(1) == "700" and
      T.HBASE2.search("Hm(373K)-Hm(273K)").group(2) == "273")
check("A5 F2: a page range after 'p' is not a year (Kimizuka 1974, p 1801-1802)", Q.cited_year(
    "N. Kimizuka and T. Katsura, Bull. Chem. Soc. Jpn., 1974, 47(7), p 1801-1802", "1802", 2000)[0] == 1974)
check("A5 F2: volume:page is not a year (Helv Chim Acta 1953;36:1820)",
      Q.cited_year("Treadwell WD, Stecher J. Helv Chim Acta 1953;36:1820.", "1820", 2000)[0] == 1953)
L, R, a = I.parse_reaction("L -> Cu10Zr7 + CuZr")
check("invariants: reaction parse and class", (L, R) == (["L"], ["Cu10Zr7", "CuZr"]) and I.derive_class(L, R, a) == "eutectic")
L, R, a = I.parse_reaction("Liquid + bcc(W) -> σ(RuW)")
check("invariants: peritectic / metatectic are not decided by the arrow direction",
      I.derive_class(L, R, a) == "peritectic|metatectic/catatectic" and I.compatible("peritectic", I.derive_class(L, R, a)))
check("invariants: L1 / L2 are liquids (monotectic)", I.derive_class(*I.parse_reaction("L1 -> L2+(Ni)")) == "monotectic")
check("invariants: comma lists and labelled invariants are not parsed", I.parse_reaction("β3 ↔ β2, (Ni), γ") is None
      and I.parse_reaction("p1: L+(Ag) -> ξAg") is None)
ref = P.reference_state_mismatch({"reference_state": {"PT": {"phase": "LIQUID"}, "SN": {"phase": "LIQUID"}},
                                  "temperature_K": 1659}, {"Pt": 3, "Sn": 1})
check("A4-N2 parity: liquid element references are not the 0 K DFT reference", bool(ref), ref)
check("A4-N2 parity: SER solid references pass", P.reference_state_mismatch(
    {"reference_state": {"PT": {"phase": "FCC_A1"}, "SN": {"phase": "BCT_A5"}}}, {"Pt": 3, "Sn": 1}) is None)

# ------------------------------------------------------------------------------------------------ release part
data = orchid_cite.DEFAULT_DATA
path = os.environ.get("ORCHID_LITERATURE") or os.path.join(data, "literature", "orchid-literature.jsonl.gz")
if os.path.exists(path):
    recs = [json.loads(x) for x in gzip.open(path, "rt", encoding="utf-8")]
    acc = [r for r in recs if r["flags"].get("extraction_accepted")]
    by = {r["id"]: r for r in recs}
    tiny = [r["id"] for r in acc if r["property"] in ("HM_FORM", "GM_FORM") and isinstance(r.get("value"), (int, float))
            and 0 < abs(r["value"]) < 0.1 and str(r.get("unit", "")).startswith("kJ") and r["flags"].get("espei_ready")]
    check("A4-N1: no ESPEI-ready formation value below 0.1 kJ/mol-atom without a page check", not tiny, tiny[:5])
    lost = [r["id"] for r in acc if T.unit_prefix_lost(r)]
    check("A4-N1: no accepted value stored as printed/1000 under a kJ caption", not lost, lost[:5])
    mgb = by.get("ORC-ITEAB5IK2I")
    check("A4-N1: Mg-B Table 2 MgB2 is -18.55 kJ/mol-atom", mgb is None or abs(mgb["value"] + 18.55) < 1e-6,
          mgb and mgb["value"])
    bad = []
    for r in recs:
        els = set(r.get("elements") or r["system"].split("-"))
        if len(els) == 1 and any((C.formula_fractions(p.get("label"), list(r.get("source_system", r["system"]).split("-")))
                                  or {}).keys() - els for p in r.get("phases") or [] if p.get("label")):
            bad.append(r["id"])
    check("A1-N1: no record filed under a unary system carries a compound formula label", not bad, bad[:5])
    bad = [r["id"] for r in recs if isinstance(r.get("reference_state"), dict)
           and set(r["reference_state"]) - set(r.get("elements") or r["system"].split("-"))]
    check("A3-N8: no reference-state key outside the record's system", not bad, bad[:5])
    bad = []
    for r in acc:
        ph = [p.get("label") for p in r.get("phases") or [] if p.get("label")]
        if len(ph) != 1 or r.get("overall_composition_basis") != "from_formula" or not C.INTRINSIC.match(r["property"]):
            continue
        f = C.formula_fractions(ph[0], r.get("elements"))
        oc = r.get("overall_composition") or {}
        if f and len(f) >= 2 and abs(sum(oc.values()) - 1) <= 0.01 and max(abs(oc.get(e, 0) - f.get(e, 0)) for e in f) > 0.02:
            bad.append(r["id"])
    check("A1-N2: every from-formula overall composition equals its printed formula", not bad, bad[:5])
    pu = by.get("ORC-EBMTXVKTB3.dc77b1")
    check("A3-N6: PuO1,52 -838.0 per formula = -332.5 kJ/mol-atom", pu is None or abs(pu["value"] + 332.54) < 0.05,
          pu and pu["value"])
    req = ("INVARIANT_TEMPERATURE", "TRANSITION_TEMPERATURE", "PHASE_COMPOSITION", "PHASE_ASSEMBLAGE", "HM_FORM",
           "GM_FORM", "SM_FORM")
    bad = [r["id"] for r in acc if r["property"] in req and not r.get("phases")]
    check("state rule: no accepted invariant / transition / phase composition / formation value without phases",
          not bad, (len(bad), bad[:5]))
    bad = []
    for r in acc:
        if r["property"] != "INVARIANT_TEMPERATURE" or r.get("record_kind") == "running_text":
            continue
        pr = I.parse_reaction((r.get("reaction") or {}).get("printed"))
        if pr and len({x for x in pr[0] + pr[1]}) > len(r.get("phases") or []) and \
                not any(g["gate"] == "R3I" for g in r.get("gate_failures") or []):
            ph = [p.get("label") for p in r.get("phases") or []]
            if any(not any(I._same(x, h) for h in ph) for x in pr[0] + pr[1]):
                bad.append(r["id"])
    check("invariants (a): no table invariant misses a phase its printed reaction names (unless flagged)", not bad, bad[:5])
    bad = [r["id"] for r in acc if r["property"] == "INVARIANT_TEMPERATURE" and (r.get("reaction") or {}).get("printed")
           and I.parse_reaction(r["reaction"]["printed"]) and not (r.get("reaction") or {}).get("type")
           and "|" not in str(I.derive_class(*I.parse_reaction(r["reaction"]["printed"]))) and
           I.derive_class(*I.parse_reaction(r["reaction"]["printed"]))]
    check("invariants (b): a printed reaction that decides its class has a class", not bad, bad[:5])
    bad = []
    for r in acc:
        rx = r.get("reaction") or {}
        pr = I.parse_reaction(rx.get("printed"))
        if r["property"] == "INVARIANT_TEMPERATURE" and pr and rx.get("type") and \
                not I.compatible(rx["type"], I.derive_class(*pr)) and \
                not any(g["code"] == "reaction_class_conflicts_printed_reaction" for g in r.get("gate_failures") or []):
            bad.append(r["id"])
    check("invariants (d): every class agrees with its printed reaction or is flagged", not bad, bad[:5])
    bad = [r["id"] for r in acc if r["property"] == "INVARIANT_TEMPERATURE" and r["flags"].get("espei_ready")
           and isinstance(r.get("value"), (int, float)) and (r["value"] > 6000 or r["value"] < 20)]
    check("invariants (f): no ESPEI-ready invariant temperature below 20 K or above 6000 K", not bad, bad[:5])
    bad = [r["id"] for r in acc if r["evidence_class"] == "model" and T.DFT_ROW.search(
        str(r.get("method") or "") + " " + str(r.get("property_printed") or "")) and not T.NOT_DFT_ROW.search(
        str(r.get("method") or "") + " " + str(r.get("property_printed") or ""))]
    check("A6-N6: no model-class value whose own method / label names a first-principles method", not bad, bad[:5])
    bad = [r["id"] for r in acc if r["property"] == "EMF" and (T.EMF_INTERCEPT.match(str(r.get("property_printed") or ""))
                                                              or T.EMF_SLOPE.match(str(r.get("property_printed") or "")))]
    check("A3-N4: no EMF intercept / slope stored as an EMF value", not bad, bad[:5])
    bad = [r["id"] for r in acc if r["property"] == "HEAT_CONTENT" and T.HBASE.search(str(r.get("property_printed") or ""))
           and abs(float(T.HBASE.search(r["property_printed"]).group(1)) - 298.15) > 1]
    check("A3-N7: no heat content relative to a base T other than 298 K stored as HEAT_CONTENT", not bad, bad[:5])
    bad = [r["id"] for r in recs if r["property"] in ("HM_FORM", "GM_FORM") and r["flags"].get("espei_ready")
           and r.get("unit") not in ("kJ/mol-atom", "kJ/mol (basis unstated)")]
    check("H9: no ESPEI-ready formation value in a unit other than kJ/mol-atom / kJ/mol (basis unstated)", not bad, bad[:5])
    bad = [r["id"] for r in acc for p in r.get("phases") or [] if p.get("model_phase") and str(p["model_phase"]).upper()
           .split("#")[0] in ("BCC_A2", "FCC_A1", "HCP_A3") and re.search(r"[A-Z][a-z]?\d*[A-Z][a-z]?\d*O\d|SiC|Y2O3",
                                                                       str(p.get("label") or ""))]
    check("A2-N2: no oxide / carbide compound label mapped to an element structure (betaSiC(Cubic), beta-Y2O3)", not bad, bad[:5])
    bad = [r["id"] for r in acc for p in r.get("phases") or [] if p.get("model_phase") and
           re.search(r"(?<![a-z])[A-Z][a-z]?x(?![a-z])", re.sub("[α-ω]", "", str(p.get("label") or "")))
           and not re.match(r"^(FCC|BCC|HCP)", str(p["model_phase"]).upper())]
    check("A2-N1: no variable-formula label (InxLa, LixSi) mapped to a fixed compound", not bad, bad[:5])
    pp = os.path.join(os.path.dirname(path), "papers.json")
    if os.path.exists(pp):
        papers = json.load(open(pp, encoding="utf-8"))
        mgzr = [r["id"] for r in recs if r.get("paper") == "doi:10.1016/s0364-5916(98)00037-6"
                and ((r.get("locator") or {}).get("page") or 0) > 6]
        check("A6-N5: no value read from the COST 507 book is filed under the Mg-Zr article DOI", not mgzr, mgzr[:3])
        cost = papers.get("sha256:85fd53a1a176b694") or {}
        check("A6-N5: the COST 507 book is registered with its own citation", "COST 507" in str(cost.get("citation")))
ppath = os.environ.get("ORCHID_PARITY") or os.path.join(data, "literature", "dft_parity.json")
if os.path.exists(ppath):
    d = json.load(open(ppath, encoding="utf-8"))
    recs_p = d["records"]
    bad = [p["rid"] for p in d["points"] if p["default"] and (recs_p[p["rid"]].get("ref_mismatch")
                                                              or recs_p[p["rid"]].get("gate_failures"))]
    check("A4-N2 / A5-N1: no default parity point with a non-SER reference state or a gate failure", not bad, bad[:5])
    oc = d["counts"].get("outlier_classification") or {}
    check("outlier classification: the three classes sum to the examined total",
          oc and oc["per_formula"] + oc["unsigned"] + oc["unexplained"] == oc["examined"], oc)
    check("A6-N1: dft_parity.json no longer describes the lowest-energy pairing only or the withdrawn conversion",
          "converted (basis_corrected_from_formula)" not in json.dumps(d["rules"]) and
          "against the lowest-energy DFT entry at the same composition" not in d["description"])

print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
