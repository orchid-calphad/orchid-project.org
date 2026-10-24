#!/usr/bin/env python3
"""Audit round 2 (2026-10-24), audits A3 and A4: permanent invariant tests, one per error class found.

Unit part (always): rule helpers (method classes, reading resolution, basis screen, expression terms, label
temperatures). Release part: every invariant on the exported literature records (ORCHID_LITERATURE, default the
orchid-data checkout's literature/orchid-literature.jsonl.gz) and on its UNCERTAINTY_CALIBRATION.json.
usage: python3 tools/tests/test_r2_thermo.py"""
import gzip
import json
import os
import re
import sys

sys.dont_write_bytecode = True
TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, TOOLS)
import literature_gates as G  # noqa: E402
import literature_r2_thermo as T  # noqa: E402
import orchid_cite  # noqa: E402
import uncertainty_budget as U  # noqa: E402

fails = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + ("" if cond or not detail else "  [%s]" % detail))
    if not cond:
        fails.append(name)


# ------------------------------------------------------------------------------------------------ unit part
mc = lambda m: U.method_class({"method": m, "evidence_class": "experimental"}, "HM_FORM")   # noqa: E731
check("A3 §1 method class: Knudsen cell / misprints / VPTM / HTMS are KEMS, never EMF",
      all(mc(m) == "kems_vapour" for m in ("Knudsen cell", "Knusden cell", "Kundsen eff. method", "VPTM", "HTMS",
                                           "Knudsen cell mass spectrometry", "knudsen effusion and emf")))
check("A3 §1 method class: galvanic cells and 'emf 3rd law' stay EMF", mc("galvanic cell") == "emf" and
      mc("emf 3rd law") == "emf")
check("A3 §1 method class: evaluated / estimated / reviewed values are not an experimental method",
      all(mc(m) == "evaluated_compilation" for m in ("Assessment", "estimated", "critical review", "evaluated")))
check("A3 §10 reading resolution: 1.0 has one decimal, '1.00' two, 2 none",
      U.decimals(1.0) == 1 and U.decimals("1.00") == 2 and U.decimals(2) == 0 and U.decimals(0.25) == 2)
key = ("HM_FORM", (("MG", 0.667), ("SI", 0.333)), "std", ((), "crystalline"))
vals = [("a", -26.0, "emf", "1"), ("b", -25.5, "emf", "2"), ("c", -79.1, "emf", "3"), ("d", -27.0, "emf", "4")]
check("A3 §1 calibration: a per-formula value among per-atom values of one state is screened out",
      [x[0] for x in U.basis_screen(key, vals)] == ["a", "b", "d"])
r_liq = {"reference_state": {"SB": {"phase": "LIQUID"}, "ZN": {"phase": "HCP_A3"}}, "phases": [{"label": "SbZn"}]}
r_sol = {"reference_state": {"SB": {"phase": "RHOMBOHEDRAL_A7"}, "ZN": {"phase": "HCP_A3"}}, "phases": [{"label": "SbZn"}]}
check("A3 §1 calibration: the state key separates liquid and solid reference states",
      U.ref_signature(r_liq) != U.ref_signature(r_sol))
r_split = {"elements": ["CU", "O"], "system": "CU-O", "overall_composition": {"O": 0.333}, "phases": [{"label": "CuO"}]}
check("A3 §1 calibration: the state composition comes from the printed compound, not a split row's overall composition",
      U.state_comp(r_split) == (("CU", 0.5), ("O", 0.5)))
r_cs = {"cited_sources": [{"first_author": "Lukashenko", "year": "1966", "reference": "G.M. Lukashenko and V.N. Eremenko"}]}
r_cs2 = {"cited_sources": [{"first_author": "Lukashenko", "year": "1966", "reference": "G.M. Lukashenko, V.N. Eremenko, Russ. J."}]}
check("A3 §1 calibration: one cited work printed two ways is one source", U.source_identity(r_cs) == U.source_identity(r_cs2))
check("A3 §6 expression terms: A / B / Cp coefficients recognised",
      T.expression_term({"property_printed": "Δ°G (constant term A of A + B·T)"})[0] == "A" and
      T.expression_term({"property_printed": "ΔG° (coefficient of T)"})[0] == "B" and
      T.expression_term({"property_printed": "cp coefficient d"})[0] == "d" and
      T.expression_term({"property_printed": "ΔfH°298"}) is None)
check("A3 §8 temperature printed in the quantity label", T.T_IN_LABEL.findall("G^ex_Li (633 K)") == ["633"] and
      T.T_IN_LABEL.findall("ΔGf at 1000 K") == ["1000"] and not T.T_IN_LABEL.findall("S°298"))
sp, _ = U.printed_component({"uncertainty": {"basis": "printed", "value": 3570.3}, "value": 147433.0,
                              "locator": {"caption": "The uncertainty corresponds to a 95 % confidence interval."}},
                             "value", "HEAT_CONTENT")
check("A3 §10 printed sigma kind: a caption-stated 95 % interval enters as half the printed value",
      sp["kind"] == "two_sigma_or_95_percent" and abs(sp["value"] - 1785.15) < 1.0)
check("A4 E3 the DFT-driven basis conversion is gone from literature_quality",
      "basis_ok and big and len(cands) == 1" not in open(os.path.join(TOOLS, "literature_quality.py")).read())

# ------------------------------------------------------------------------------------------------ release part
data = orchid_cite.DEFAULT_DATA
path = os.environ.get("ORCHID_LITERATURE") or os.path.join(data, "literature", "orchid-literature.jsonl.gz")
if not os.path.exists(path):
    print("SKIP release part (no %s)" % path)
else:
    R = [json.loads(l) for l in gzip.open(path, "rt", encoding="utf-8")]
    visual = G.load_visual()
    num = T.num

    check("A4 E3: no value carries the circular basis_corrected_from_formula conversion",
          not any((r.get("flags") or {}).get("basis_corrected_from_formula") for r in R))

    # per-column page decisions: one basis per column, value = printed x unit factor / atoms
    bad = []
    for key, chk in visual.items():
        if not key.startswith("table_basis:") or chk.get("verdict") != "decided":
            continue
        asset = key[len("table_basis:"):]
        for r in R:
            if G.asset_of(r) != asset or r.get("property") not in (chk.get("columns") or {}):
                continue
            col = chk["columns"][r["property"]]
            if r["id"] in (chk.get("rows_excluded") or []) or any(g["code"].startswith("table_basis_") for g in r.get("gate_failures") or []):
                continue
            n = T.record_atoms(r) if col["basis"] == "per_formula_unit" else 1.0
            f = T.unit_factor(r["property"], col.get("unit"))
            want = num(r.get("value_printed")) * f / n
            if r.get("unit_basis") != "page_checked_table" or abs(r["value"] - want) > 1e-6 * max(1, abs(want)):
                bad.append(r["id"])
    check("A3 §2 / A4 E3: one page-read basis per printed column (value = printed / atoms)", not bad, bad[:5])

    bad = []
    for key, chk in visual.items():
        if key.startswith("table_sign:") and chk.get("verdict") == "decided":
            bad += [r["id"] for r in R if G.asset_of(r) == key[len("table_sign:"):] and r.get("property") == chk["property"]
                    and re.search(chk.get("column_printed") or "", str(r.get("property_printed") or ""))
                    and isinstance(r.get("value"), (int, float)) and r["value"] > 0]
    check("A3 §7: a column printed without signs is negative as a whole (Co-Y Co17Y2, Co5Y)", not bad, bad[:5])

    bad = []
    for r in R:
        d, u = r.get("per_formula_unit_divisor"), r.get("uncertainty") or {}
        if isinstance(d, (int, float)) and d > 1 and u.get("basis") in ("printed", "printed_range") \
                and r.get("unit_basis") != "page_checked_table" and not u.get("rescaled_with_value"):
            bad.append(r["id"])
    check("A3 §4: a printed sigma is divided by the same atom count as its value", not bad, bad[:5])

    bad = [r["id"] for r in R if r.get("property") == "SM_FORM" and isinstance(r.get("value"), (int, float))
           and (r.get("flags") or {}).get("extraction_accepted")
           and r["value"] > (T.SABS_LIMIT if re.search(r"atom|mol-at", r.get("unit") or "") else T.SABS_LIMIT_FORMULA)]
    check("A3 §3: no absolute entropy is labelled a formation entropy (SM_FORM <= 15 J/K/mol-atom)", not bad, bad[:5])

    bad = [r["id"] for r in R if r.get("property") == "HEAT_CONTENT" and r.get("value")
           and isinstance(r.get("temperature_K"), (int, float)) and abs(r["temperature_K"] - 298.15) <= 2
           and not (T.HREF_FIRST.search(str(r.get("property_printed") or "")) and
                    not re.search(r"H\s*°?\s*[_(]?\s*T", str(r.get("property_printed") or "")))
           and (r.get("flags") or {}).get("extraction_accepted")]
    check("A3 §5: an accepted H(T) - H(298) with a non-zero value is never stored at 298 K", not bad, bad[:5])

    bad = [r["id"] for r in R if r.get("property") in T.EXPR_PROPS and T.expression_term(r) and r.get("temperature_K") is not None]
    check("A3 §6: no coefficient of a printed expression in T is stored as a value at one temperature", not bad, bad[:5])

    bad = [r["id"] for r in R if r.get("property") not in ("HEAT_CONTENT", "INVARIANT_TEMPERATURE", "TRANSITION_TEMPERATURE")
           and not str(r.get("property")).endswith("_EXPRESSION_TERM")
           and isinstance(r.get("temperature_K"), (int, float)) and abs(r["temperature_K"] - 298.15) > 2
           and T.T298_LABEL.search(str(r.get("property_printed") or "")) and not T.HEAT_INCREMENT.search(str(r.get("property_printed") or ""))
           and r.get("temperature_basis") not in ("quantity_label", "row_temperature")]
    check("A3 §8: a quantity printed with a 298 subscript is stored at 298.15 K", not bad, bad[:5])

    bad = [r["id"] for r in R if r.get("evidence_class") == "experimental" and r.get("temperature_K") == 0
           and r.get("property") not in ("INVARIANT_TEMPERATURE", "TRANSITION_TEMPERATURE")]
    check("A4 E6: no experimental value carries 0 K", not bad, bad[:5])

    bad = [r["id"] for r in R if r.get("property") == "HM_FORM" and r.get("evidence_class") == "experimental"
           and T.KLEPPA.search(str(r.get("method") or "")) and isinstance(r.get("temperature_K"), (int, float))
           and r["temperature_K"] >= 600 and not T.AT_T_CAPTION.search(str((r.get("locator") or {}).get("caption") or ""))
           and not re.search(r"oxide|f,ox", str(r.get("property_printed") or "") + str((r.get("locator") or {}).get("caption") or ""), re.I)]
    check("A3 §8: direct-synthesis / solute-solvent drop standard enthalpies are stored at 298.15 K", not bad, bad[:5])

    bad = [r["id"] for r in R if str(r.get("property")).endswith("_FORM") and r.get("property") != "GM_FORM_EXPRESSION_TERM"
           and (T._ref_compounds(r) or T.FROM_COMPOUNDS.search(" ".join(str(x or "") for x in (r.get("property_printed"), (r.get("locator") or {}).get("caption")))))
           and not re.search(r"f\s*,\s*elements|from (?:the )?elements", " ".join(str(x or "") for x in (r.get("property_printed"), (r.get("locator") or {}).get("caption"))), re.I)]
    check("A3 §9: no _FORM quantity is referred to compounds (oxides)", not bad, bad[:5])

    bad = [r["id"] for r in R if isinstance(r.get("reference_state"), dict) and any(
        isinstance(v, dict) and isinstance(v.get("phase"), str) and T.REF_UNIT.match(v["phase"]) for v in r["reference_state"].values())]
    check("A3 §9: no reference state is a unit string", not bad, bad[:5])

    bad = [r["id"] for r in R if "basis unstated" in (r.get("unit") or "") and T.PER_ATOM_SLASH.search(str(r.get("unit_printed") or ""))]
    check("A3 §9: a per-atom unit printed 'kJ/mol/at' keeps its per-atom basis", not bad, bad[:5])

    bad = [r["id"] for r in R if r.get("temperature_K") is None and r.get("value") is not None and not r.get("temperature_range_K")
           and len(set(T.T_IN_LABEL.findall(str(r.get("property_printed") or "")))) == 1]
    check("A3 §8: a temperature printed in the quantity label is stored", not bad, bad[:5])

    bad = [r["id"] for r in R if r.get("evidence_class") == "model" and r.get("property") in ("HM_FORM", "SM_FORM", "GM_FORM")
           and T.FP_CAPTION.search(str((r.get("locator") or {}).get("caption") or "") + " " + str(r.get("property_printed") or ""))
           and T.DFT_METHOD.search(str(r.get("method") or "")) and not T.NOT_DFT.search(str(r.get("method") or ""))]
    check("A4 E7: first-principles rows of mixed tables are evidence class dft, not model", not bad, bad[:5])

    bad = [r["id"] for r in R if any(T.NONCRYST.search(str(p.get("label") or "")) for p in r.get("phases") or [])
           and not (r.get("flags") or {}).get("non_crystalline_state")]
    check("A4 E7: liquid / amorphous states printed at room temperature are flagged non_crystalline_state", not bad, bad[:5])

    RELPH = re.compile(r"(?:relative to|referred to|with respect to)\s+(?:the\s+)?(?:B2|B19|L1_?0|L1_?2|D0_?3|martensit\w*|austenit\w*)\b", re.I)
    bad = [r["id"] for r in R if str(r.get("property")).endswith("_FORM") and RELPH.search(
        " ".join(str(x or "") for x in (r.get("property_printed"), (r.get("locator") or {}).get("caption"))))]
    check("A3 follow-up: energies referred to another phase (e.g. B2) are not formation quantities", not bad, bad[:5])

    calp = os.path.join(os.path.dirname(path), "UNCERTAINTY_CALIBRATION.json")
    if os.path.exists(calp):
        cal = json.load(open(calp))
        e = cal["value_types"]["HM_FORM"]["classes"].get("emf") or {}
        check("A3 §1: the EMF formation-enthalpy floor is not inflated by basis errors (calibrated <= 4 or prior)",
              e.get("basis") != "calibrated" or e.get("floor", 0) <= 4.0, e.get("floor"))
        k = cal["value_types"]["HM_FORM"]["classes"].get("kems_vapour") or {}
        check("A3 §1: no spurious KEMS floor from per-formula values (calibrated <= 8 or prior)",
              k.get("basis") != "calibrated" or k.get("floor", 0) <= 8.0, k.get("floor"))
        check("A3 §1: evaluated / estimated values are never calibration pairs",
              not (cal["value_types"]["HM_FORM"]["classes"].get("evaluated_compilation") or {}).get("n_pairs"))

print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
