#!/usr/bin/env python3
"""Fixture tests of tools/literature_gates.py (G1-G12). No corpus, no PDFs, no network: records are built here, the
SGTE liquid bounds come from a stub, visual verdicts are passed in. Checks every gate's correction / demotion / flag,
that printed values and quotes are never altered and that no record is deleted except by the G7 duplicate merge.
usage: python3 tools/tests/test_literature_gates.py"""
import copy
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import literature_gates as G  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


class Rule:
    """Stub of phase_rule.PhaseRule: lowest SGTE liquid temperatures and invariants."""
    liq_min = {"AG-PT": (1240.0, "AgPt.tdb"), "AG-CA": (736.7, "AgCa.tdb"), "AL-NI": (913.0, "AlNi.tdb")}
    inv_T = {"AG-PT": [1459.0], "AG-CU": [1052.0], "AG-SN": [494.0, 997.0], "CU-SN": [500.0, 862.0, 1029.0]}

    def liquid_bound(self, system, factor=1.0):
        t, tdb = self.liq_min[system]
        return t - factor * max(100.0, 0.10 * t), {"tdb": tdb, "lowest_liquid_T_K": t}


N = [0]


def rec(**kw):
    N[0] += 1
    r = {"id": "ORC-T%04d" % N[0], "pid": "doi:10.1/x/p1/tab-1/r%02d#x" % N[0], "system": "AL-NI", "elements": ["AL", "NI"],
         "record_kind": "general_table", "observation_type": "equilibrium_property", "evidence_class": "experimental",
         "paper": "doi:10.1/x", "locator": {"artifact": "table", "asset": "doi:10.1/x/p1/tab-1", "page": 1, "row": N[0]},
         "pressure_Pa": 101325, "pressure_basis": "stated_default_1_atm",
         "flags": {"extraction_accepted": True, "espei_ready": True, "espei_certified": True}}
    r.update(kw)
    if "system" in kw and "elements" not in kw:
        r["elements"] = kw["system"].split("-")
    return r


def vc(**kw):
    d = {"page": 1, "pdf_sha256_prefix": "0" * 16, "image_sha256": "1" * 64, "seen": "fixture", "checked": "2026-10-23"}
    d.update(kw)
    return d


def run(records, visual=None, units=None, plaus=None, describe=None):
    printed = [(r["id"], r.get("value_printed"), r.get("quote"), r.get("temperature_printed")) for r in records]
    out, rep = G.apply(records, Rule(), visual=visual or {}, store_units=units or {}, plaus_store=plaus or {}, describe=describe)
    after = {r["id"]: (r["id"], r.get("value_printed"), r.get("quote"), r.get("temperature_printed")) for r in out}
    check("printed values / quotes / printed temperatures never altered",
          all(after[x[0]] == x for x in printed if x[0] in after))
    return {r["id"]: r for r in out}, rep


def codes(r):
    return {g["code"] for g in r.get("gate_failures") or []}


# ------------------------------------------------------------------ G1a table temperature unit (Ag-Pt Table 6 shape)
def agpt_table(asset):
    rows = []
    for t, x in [(961.9, 0.0), (1000, 0.02), (1186, 0.201), (1769, 1.0)]:
        comp = [{"phase": "L", "composition": {"PT": x}, "composition_source": "printed", "basis": "mole_fraction"}]
        if x in (0.0, 1.0):
            comp.append({"phase": "(Ag)" if x == 0 else "(Pt)", "composition": {"PT": x}, "composition_source": "printed",
                         "basis": "mole_fraction"})
        rows.append(rec(system="AG-PT", property="PHASE_COMPOSITION", observation_type="phase_region_vertex_set",
                        temperature_K=float(t), temperature_printed=str(t), phases=[{"label": "L"}, {"label": "(Ag)"}],
                        phase_compositions=comp, locator={"artifact": "table", "asset": asset, "page": 7, "row": len(rows) + 1},
                        flags={"extraction_accepted": True, "espei_ready": True, "plausibility": {"result": "fail"}}))
    return rows


rows = agpt_table("doi:x/p7/tab-6")
out, rep = run(rows)
check("G1a °C table stored as K without a visual verdict: flagged, not ESPEI-ready, value kept",
      all("temperature_unit_table_suspect" in codes(r) and not r["flags"]["espei_ready"] and r["temperature_K"] < 1800
          for r in out.values()))
rows = agpt_table("doi:x/p7/tab-6b")
out, rep = run(rows, {"table_T_unit:doi:x/p7/tab-6b": vc(key="table_T_unit:doi:x/p7/tab-6b", verdict="confirmed", unit="C")})
r0 = sorted(out.values(), key=lambda r: r["temperature_K_stored"])[0]
check("G1a confirmed °C table: whole table converted from the printed numbers",
      abs(r0["temperature_K"] - 1235.05) < 1e-6 and r0["temperature_K_stored"] == 961.9 and r0["temperature_printed"] == "961.9")
check("G1a correction carries rule, evidence and visual_check",
      r0["corrections"][0]["rule"] == "table_temperature_unit_C" and r0["corrections"][0]["visual_check"]["verdict"] == "confirmed"
      and "pure_element_rows_match_celsius" in r0["corrections"][0]["evidence"])
check("G1a stale CALPHAD check moved aside", "plausibility" not in r0["flags"] and "plausibility_before_correction" in r0["flags"])

# ------------------------------------------------------------------ G1 sibling row, non-temperatures, ranges, liquids
sib_a = rec(system="RU-U", property="INVARIANT_TEMPERATURE", observation_type="invariant_reaction", value=1588.0, unit="K",
            temperature_K=1588.0, phases=[{"label": "L"}, {"label": "Ru3U"}, {"label": "hcp"}])
sib_b = rec(system="RU-U", property="INVARIANT_TEMPERATURE", observation_type="invariant_reaction", value=1861.0, unit="K",
            temperature_K=1861.15, temperature_printed="1588 ± 20 °C / 1861 ± 20 K",
            phases=[{"label": "Liquid"}, {"label": "Ru3U"}, {"label": "hcp"}])
dash = rec(property="INVARIANT_TEMPERATURE", observation_type="invariant_reaction", value=0.0, unit="K", temperature_K=0.0,
           temperature_printed="–", phases=[{"label": "L"}, {"label": "NiAl"}])
compT = rec(property="INVARIANT_TEMPERATURE", observation_type="invariant_reaction", value=25.0, unit="K", temperature_K=25.0,
            phases=[{"label": "L"}, {"label": "SrAl4"}],
            phase_compositions=[{"phase": "L", "composition": {"NI": 0.25}, "composition_source": "printed", "basis": "mole_fraction"}])
rng = rec(property="PHASE_COMPOSITION", observation_type="phase_region_vertex_set", temperature_K=1480.15,
          temperature_printed="418.9 to 1207", phases=[{"label": "liquid"}])
supercond = rec(system="GA-ZR", property="TRANSITION_TEMPERATURE", observation_type="phase_transition", value=0.38, unit="K",
                temperature_K=0.38, temperature_printed="0.38 K", record_kind="running_text", quote="Tc of 0.38 K",
                phases=[{"label": "GaZr2"}])
lowliq = rec(system="AG-CA", property="TRANSITION_TEMPERATURE", observation_type="phase_transition", value=100.0, unit="K",
             temperature_K=100.0, temperature_printed="100", phases=[{"label": "L"}, {"label": "Ca2Ag7"}])
nearliq = rec(property="PHASE_COMPOSITION", observation_type="phase_region_vertex_set", temperature_K=650.0,
              temperature_printed="650", phases=[{"label": "liquid"}, {"label": "(Al)"}])
out, rep = run([sib_a, sib_b, dash, compT, rng, supercond, lowliq, nearliq])
check("G1 sibling row: °C number stored as K is not accepted", "temperature_unit_sibling_row" in codes(out[sib_a["id"]])
      and not out[sib_a["id"]]["flags"]["extraction_accepted"] and not codes(out[sib_b["id"]]) & {"temperature_unit_sibling_row"})
check("G1b dash temperature not accepted", "temperature_not_a_temperature" in codes(out[dash["id"]]))
check("G1b temperature equal to a printed at.% not accepted", "temperature_not_a_temperature" in codes(out[compT["id"]]))
check("G1b superconducting 0.38 K printed with its unit is kept", "temperature_not_a_temperature" not in codes(out[supercond["id"]]))
check("G1d printed range stored as a point flagged", "temperature_range_stored_as_point" in codes(out[rng["id"]]))
check("G1c liquid far below the SGTE liquid: not accepted", "liquid_below_physical_minimum" in codes(out[lowliq["id"]])
      and not out[lowliq["id"]]["flags"]["extraction_accepted"])
check("G1c liquid inside the margin: flagged, still accepted", "liquid_below_physical_minimum" in codes(out[nearliq["id"]])
      and out[nearliq["id"]]["flags"]["extraction_accepted"] and not out[nearliq["id"]]["flags"]["espei_ready"])

# ------------------------------------------------------------------ G2 magnitudes
act = [rec(system="FE-TI", property="ACR_TI", unit="1", value=v, value_printed=v, property_printed="a_Ti·10^3",
           locator={"artifact": "table", "asset": "doi:x/p2/tab-iii", "page": 2, "row": i}, temperature_K=1873.0,
           phases=[{"label": "L"}]) for i, v in enumerate([3.16, 1.84, 4.37])]
out, rep = run(copy.deepcopy(act))
check("G2 x10^3 activity header without verdict: not accepted (table scale)", all("table_unit_scale_suspect" in codes(r) for r in out.values()))
out, rep = run(copy.deepcopy(act), {"table_scale:doi:x/p2/tab-iii:ACR_TI": vc(verdict="confirmed", factor=0.001)})
r = sorted(out.values(), key=lambda r: r["value_stored"])[0]
check("G2 x10^3 activity header confirmed: value / 1000, printed kept", abs(r["value"] - 0.00184) < 1e-12 and r["value_printed"] == 1.84)
jk = rec(property="HM_FORM", unit="kJ/mol-atom", value=-8210.0, value_printed=-8210.0, unit_printed="J/mol", phases=[{"label": "SbZn"}],
         system="SB-ZN")
mass = rec(system="AG-PB", property="CPM", unit="J/mol (basis unstated)/K", value=0.153, value_printed=0.153, unit_printed="J/g K",
           phases=[{"label": "UNKNOWN"}])
gat = rec(system="AG-PB", property="CPM", unit="J/mol-atom/K", value=26.0, value_printed=26.0, unit_printed="J/g at. K",
          phases=[{"label": "UNKNOWN"}])
oor = rec(property="HM_FORM", unit="kJ/mol-atom", value=496.7, value_printed=496.7, unit_printed="kJ/mol", phases=[{"label": "Ho (fcc)"}],
          system="AU-HO", record_kind="general_table", evidence_class="model",
          locator={"artifact": "table", "asset": "doi:x/p2/tab-9", "page": 2, "row": 1})
oor2 = rec(property="HM_FORM", unit="kJ/mol-atom", value=-40.0, value_printed=-40.0, unit_printed="kJ/mol", phases=[{"label": "AuHo"}],
           system="AU-HO", locator={"artifact": "table", "asset": "doi:x/p2/tab-9", "page": 2, "row": 2})
oor3 = rec(property="HM_FORM", unit="kJ/mol-atom", value=-35.0, value_printed=-35.0, unit_printed="kJ/mol", phases=[{"label": "AuHo2"}],
           system="AU-HO", locator={"artifact": "table", "asset": "doi:x/p2/tab-9", "page": 2, "row": 3})
z = rec(system="CU-ZR", property="ACR_CU", unit="1", value=0.0, value_printed=0.0, phases=[{"label": "L"}])
out, rep = run([jk, mass, gat, oor, oor2, oor3, z])
check("G2 J value under a printed J unit stored as kJ: divided by 1000", abs(out[jk["id"]]["value"] + 8.21) < 1e-9
      and out[jk["id"]]["value_stored"] == -8210.0)
check("G2 per-gram heat capacity relabelled J/g/K and not ESPEI-ready", out[mass["id"]]["unit"] == "J/g/K"
      and "mass_basis_unit" in codes(out[mass["id"]]))
check("G2 per gram-atom is per mole of atoms: untouched", out[gat["id"]]["unit"] == "J/mol-atom/K" and not codes(out[gat["id"]]))
check("G2 unexplained out-of-range value not accepted, magnitude_suspect set", "value_out_of_physical_range" in codes(out[oor["id"]])
      and not out[oor["id"]]["flags"]["extraction_accepted"] and out[oor["id"]]["flags"]["magnitude_suspect"])
check("G2/G9 activity stored as 0 not ESPEI-ready", "activity_rounded_to_zero" in codes(out[z["id"]]))
out, rep = run([copy.deepcopy(oor)], {"record:" + oor["id"]: vc(verdict="corrected", value=-49.67, value_printed_page="-49.67")})
check("G2 value read from the page when unambiguous", out[oor["id"]]["value"] == -49.67 and out[oor["id"]]["value_stored"] == 496.7)

# ------------------------------------------------------------------ G3 signs
sg = [rec(system="AL-PT", property="HM_FORM", unit="kJ/mol-atom", value=v, value_printed=v, property_printed="(−ΔfH)",
          locator={"artifact": "table", "asset": "doi:x/p4/tab-4", "page": 4, "row": i}, phases=[{"label": "AlPt"}])
      for i, v in enumerate([100.0, 80.0])]
out, rep = run(copy.deepcopy(sg))
check("G3 positive values under a (-Δ) heading without verdict: flagged", all("sign_heading_not_applied" in codes(r) for r in out.values()))
out, rep = run(copy.deepcopy(sg), {"table_sign:doi:x/p4/tab-4": vc(verdict="confirmed")})
check("G3 confirmed (-Δ) heading: negated, sign_from_heading", all(r["value"] < 0 and r["flags"]["sign_from_heading"] for r in out.values()))
pos = rec(system="LI-SB", property="GM_FORM", unit="kJ/mol-atom", value=25.0, value_printed=25.0, phases=[{"label": "Li3Sb"}])
uns = rec(property="HM_FORM", unit="kJ/mol-atom", value=-10.0, phases=[{"label": "AlNi"}],
          flags={"extraction_accepted": True, "espei_ready": True, "unsigned_in_print": True})
out, rep = run([pos, uns])
check("G3 positive formation of a compound flagged", "formation_positive_check_sign" in codes(out[pos["id"]]))
check("G3 unsigned_in_print blocks ESPEI", not out[uns["id"]]["flags"]["espei_ready"])

# ------------------------------------------------------------------ G4 system
ox = rec(system="B-CE", property="PHASE_COMPOSITION", phases=[{"label": "γ(Ce)"}, {"label": "B4Ce"}, {"label": "CeO2"}])
unk = rec(system="N-O-U-W", property="HM_FORM", value=-1.0, unit="kJ/mol-atom", phases=[{"label": "UNKNOWN"}])
typo = rec(system="FE-HF", property="INVARIANT_TEMPERATURE", value=1533.0, unit="K", temperature_K=1533.0,
           phases=[{"label": "(βHf)"}, {"label": "FeHF2"}, {"label": "TiNi B19"}])
proto = rec(system="AL-SR", property="HM_FORM", value=-10.0, unit="kJ/mol-atom", evidence_class="dft",
            phases=[{"label": "Ca8In3-type"}, {"label": "Al2Sr CeCu2-type"}])
cem = rec(system="C-S", property="HM_FORM", value=-2932.0, unit="kJ/mol (basis unstated)",
          phases=[{"label": "Monoclinic-Ca3SiO5"}], locator={"artifact": "table", "asset": "doi:x/p8/tab-7", "page": 8, "row": 1,
                                                              "caption": "H 298 and S 298 of the compounds in the C-S system"})
out, rep = run([ox, unk, typo, proto, copy.deepcopy(cem)])
check("G4 oxide phase in a metal-boron system: system_mismatch, not accepted", "system_mismatch" in codes(out[ox["id"]])
      and not out[ox["id"]]["flags"]["extraction_accepted"])
check("G4 'UNKNOWN' is not a formula", not codes(out[unk["id"]]) & {"system_mismatch"})
check("G4 capitalisation slip FeHF2 and structure symbol B19 are not foreign elements", "system_mismatch" not in codes(out[typo["id"]]))
check("G4 prototype '-type' labels are not phase formulas", "system_mismatch" not in codes(out[proto["id"]]))
check("G4 cement C-S without verdict: not accepted", "system_mismatch" in codes(out[cem["id"]]))
out, rep = run([copy.deepcopy(cem)], {"system:doi:x/p8/tab-7": vc(verdict="confirmed", system="CA-O-SI")})
r = out[cem["id"]]
check("G4 cement C-S confirmed: reassigned to CA-O-SI, original kept", r["system"] == "CA-O-SI" and r["system_stored"] == "C-S"
      and r["elements"] == ["CA", "O", "SI"] and r["flags"]["system_reassigned"])

# ------------------------------------------------------------------ G5 semantics
gm = rec(property="GM_MIX", unit="kJ/mol-atom", value=-27.0, property_printed="Gibbs free energy of formation", phases=[{"label": "YZn3"}],
         system="Y-ZN")
fus = rec(system="AL", property="HM_FORM", unit="kJ/mol-atom", value=10.7, property_printed="ΔH fusion", phases=[{"label": "liquid"}])
uni = rec(system="MO", property="SM_FORM", unit="J/mol-atom/K", value=-1.3, property_printed="ΔfS^298", phases=[{"label": "ε"}])
inc = rec(property="SM", unit="J/mol-atom/K", value=30.0, property_printed="S(T)−S(298)", phases=[{"label": "AlNi"}])
rk = rec(property="SM", unit="J/mol-atom/K", value=-6049.0, property_printed="0L", phases=[{"label": "liquid"}])
cf = rec(system="BA-O", property="CPM", unit="J/mol-atom/K", value=-227826.0, property_printed="cp coefficient d", phases=[{"label": "BaO"}])
mg = rec(property="SM", unit="J/mol-atom/K", value=-0.89, property_printed="0S298.15 mag", phases=[{"label": "Fe17Lu2"}])
out, rep = run([gm, fus, uni, inc, rk, cf, mg])
check("G5 'Gibbs free energy of formation' stored as GM_MIX -> GM_FORM", out[gm["id"]]["property"] == "GM_FORM"
      and out[gm["id"]]["property_code_before_gates"] == "GM_MIX")
check("G5 unary fusion quantity -> HM_FUSION", out[fus["id"]]["property"] == "HM_FUSION")
check("G5 formation code in a unary system without fusion/lattice wording: not accepted",
      "formation_code_in_unary_system" in codes(out[uni["id"]]) and not out[uni["id"]]["flags"]["extraction_accepted"])
check("G5 S(T)-S(298) -> SM_INCREMENT", out[inc["id"]]["property"] == "SM_INCREMENT")
check("G5 Redlich-Kister parameter not accepted", "model_parameter_not_property" in codes(out[rk["id"]]))
check("G5 Cp polynomial coefficient not accepted", "model_parameter_not_property" in codes(out[cf["id"]]))
check("G5 magnetic entropy contribution -> SM_MAGNETIC", out[mg["id"]]["property"] == "SM_MAGNETIC")

# ------------------------------------------------------------------ G6 invariants
dupL = rec(property="INVARIANT_TEMPERATURE", observation_type="invariant_reaction", value=1000.0, unit="K", temperature_K=1000.0,
           transition_kind="eutectic", phases=[{"label": "L"}, {"label": "Liquid", "model_phase": "LIQUID"}, {"label": "NiAl"},
                                                {"label": "Ni3Al"}, {"label": "eutectic"}])
noL = rec(property="INVARIANT_TEMPERATURE", observation_type="invariant_reaction", value=1000.0, unit="K", temperature_K=1000.0,
          transition_kind="eutectic", phases=[{"label": "NiAl"}, {"label": "Ni3Al"}, {"label": "Ni5Al3"}])
four = rec(property="INVARIANT_TEMPERATURE", observation_type="invariant_reaction", value=1100.0, unit="K", temperature_K=1100.0,
           transition_kind="eutectic", phases=[{"label": "L"}, {"label": "NiAl"}, {"label": "Ni3Al"}, {"label": "Ni2Al3"}])
out, rep = run([dupL, noL, four])
r = out[dupL["id"]]
check("G6 repeated liquid and non-phase label dropped; kept the mapped liquid",
      [p["label"] for p in r["phases"]] == ["Liquid", "NiAl", "Ni3Al"] and len(r["phases_dropped"]) == 2 and not codes(r))
check("G6 eutectic without a liquid flagged", "reaction_liquid_count_inconsistent" in codes(out[noL["id"]]))
check("G6 binary row with four phases flagged", "binary_more_than_three_phases" in codes(out[four["id"]]))

# ------------------------------------------------------------------ G7 duplicates
a = rec(property="HM_FORM", unit="kJ/mol-atom", value=-41.7, phases=[{"label": "NiAl"}], locator={"artifact": "table",
        "asset": "doi:x/p3/tab-1", "page": 3, "row": 1}, readers=["a", "b"])
b = copy.deepcopy(a)
b["id"], b["pid"], b["phases"], b["readers"] = "ORC-DUP", "pid-dup", [{"label": "UNKNOWN"}], ["a"]
c = copy.deepcopy(a)
c["id"], c["phases"] = "ORC-DIFF", [{"label": "Ni3Al"}]
out, rep = run([copy.deepcopy(a), copy.deepcopy(b)])
check("G7 UNKNOWN twin merged into the named record", "ORC-DUP" not in out and "ORC-DUP" in out[a["id"]]["merged_ids"])
out, rep = run([a, b, c])
check("G7 a different named phase is not a duplicate; an UNKNOWN twin of two candidates is kept",
      "ORC-DIFF" in out and "ORC-DUP" in out)

# ------------------------------------------------------------------ G8 evidence class
t = dict(record_kind="running_text", evidence_class="model", property="HM_FORM", unit="kJ/mol-atom", value=-20.0,
         locator={"artifact": "text_passage", "asset": "doi:x/p2/text-1", "page": 2}, phases=[])
d1 = rec(quote="The enthalpy of formation obtained from first-principles (VASP, PBE) is -20.0 kJ/mol-atom.", **t)
d2 = rec(quote="The value -20.0 kJ/mol-atom from the present CALPHAD parameters agrees with DFT.", **t)
e1 = rec(quote="Smith et al. [12] measured -20.0 kJ/mol-atom by solution calorimetry.", cited_sources=[{"label": "12"}], **t)
e2 = rec(quote="The calculated -20.0 kJ/mol-atom is compared with the value measured by [12].", cited_sources=[{"label": "12"}], **t)
asr = rec(evidence_class="experimental", origin={"own_work_of_paper": True}, property="INVARIANT_TEMPERATURE", value=900.0, unit="K",
          temperature_K=900.0, phases=[], paper="doi:10.1/eval")
d1b, e1b = copy.deepcopy(d1), copy.deepcopy(e1)
out, rep = run([d1b, e1b])
check("G8 wording alone (not checked on the page): class kept, flagged for review",
      out[d1b["id"]]["evidence_class"] == "model" and "evidence_class_review" in codes(out[e1b["id"]]))
out, rep = run([d1, d2, e1, e2, asr], {"evidence:" + d1["id"]: vc(verdict="confirmed"), "evidence:" + e1["id"]: vc(verdict="confirmed"),
                                       "evidence:" + e2["id"]: vc(verdict="confirmed")},
               describe=lambda p: {"title": "Ag-Pt (Silver-Platinum)", "container": "Journal of Phase Equilibria"} if p == "doi:10.1/eval" else {})
check("G8 first-principles wording confirmed on the page -> dft", out[d1["id"]]["evidence_class"] == "dft")
check("G8 CALPHAD wording with a DFT comparison stays model", out[d2["id"]]["evidence_class"] == "model")
check("G8 cited measurement wording confirmed on the page -> experimental", out[e1["id"]]["evidence_class"] == "experimental")
check("G8 comparison sentence stays model (rule does not fire even with a verdict)", out[e2["id"]]["evidence_class"] == "model")
check("G8 evaluation (assessed) source is not experimental", out[asr["id"]]["evidence_class"] == "unspecified"
      and out[asr["id"]]["flags"]["assessed_source"])

# ------------------------------------------------------------------ G9 compositions, G11, G12
s1 = rec(property="PHASE_COMPOSITION", observation_type="phase_region_vertex_set", temperature_K=1200.0, temperature_printed="1200",
         phases=[{"label": "NiAl"}], phase_compositions=[{"phase": "NiAl", "composition": {"AL": 0.5, "NI": 0.3},
                                                          "composition_source": "printed", "basis": "mole_fraction"}])
pl = rec(property="INVARIANT_TEMPERATURE", observation_type="invariant_reaction", value=1000.0, unit="K", temperature_K=1000.0,
         phases=[{"label": "Au"}, {"label": "L"}],
         flags={"extraction_accepted": True, "plausibility": {"result": "fail", "dT_K": -1.1, "model": "x.tdb"}})
out, rep = run([s1, pl], plaus={pl["pid"]: {"dx": {"Au": -0.747, "L": -0.014}, "calculation": "pycalphad binary topology"}})
check("G9 composition not summing to 1 flagged", "phase_composition_sum_not_one" in codes(out[s1["id"]]))
p = out[pl["id"]]["flags"]["plausibility"]
check("G11 dx, tolerance and fail_reason exported", p["dx"]["Au"] == -0.747 and p["tolerance"]["dx"] == 0.05
      and p["fail_reason"] == "composition_pure_element_label")
check("G12 1 atm default relabelled, stored basis kept", out[s1["id"]]["pressure_basis"] == "assumed_default_1_atm_not_printed"
      and out[s1["id"]]["pressure_basis_stored"] == "stated_default_1_atm")
check("report has before/after counts and markdown", rep["records_in"] >= rep["records_out"] and "## 7." in G.gates_md(rep))

# ------------------------------------------------------------------ cleanup 2026-10-24: page-checked phase lists etc.
check("phase_key: glyph / LaTeX / Greek-name / parenthesis / separator spellings coincide",
      G.phase_key('α"') == G.phase_key("α″") and G.phase_key("\\beta Yb") == G.phase_key("βYb")
      and G.phase_key("gamma") == G.phase_key("γ") and G.phase_key("(εSb)") == G.phase_key("εSb")
      and G.phase_key("In3Ir-LT") == G.phase_key("In3 Ir_LT") and G.phase_key("C₁Pu") == G.phase_key("C1Pu")
      and G.phase_key("betaS") == G.phase_key("βS") and G.phase_key("Laves_C15") != G.phase_key("Laves_C36"))
check("numbered liquids are liquids in G6; Strukturbericht L12 / La are not",
      all(G.is_liq_g6({"label": l}) for l in ("L1", "L2", "L#1", "Liquid#2", "LIQUID#1", "L'"))
      and not any(G.is_liq_g6({"label": l}) for l in ("L12", "La", "Lu")))
inv = dict(property="INVARIANT_TEMPERATURE", observation_type="invariant_reaction", unit="K", record_kind="invariant_table")
loc = {"artifact": "table", "asset": "doi:10.1/t/p5/tab-2", "page": 5}
hdr = rec(value=739.0, temperature_K=739.0, transition_kind="eutectic", locator=dict(loc, row=1),
          phases=[{"label": "L"}, {"label": "Mg2Pb"}, {"label": "Phase 2"}, {"label": "Phase 3"}, {"label": "HCP_A3 (Mg)"}],
          phase_compositions=[{"phase": "Phase 2", "composition": {"PB": 0.0975}, "composition_source": "printed",
                               "basis": "mole_fraction"}], **inv)
mono = rec(value=1100.0, temperature_K=1100.0, transition_kind="monotectic", locator=dict(loc, row=2),
           phases=[{"label": "Liquid#1"}, {"label": "Liquid#2"}, {"label": "Ni3Al"}], **inv)
miss = rec(value=1000.0, temperature_K=1000.0, transition_kind="eutectic", locator=dict(loc, row=3),
           phases=[{"label": "L"}, {"label": "NiAl"}], reaction={"printed": "L -> NiAl + Ni3Al"}, **inv)
eo = rec(value=900.0, temperature_K=900.0, transition_kind="eutectic", locator=dict(loc, row=4),
         phases=[{"label": "NiAl"}, {"label": "Ni3Al"}, {"label": "Ni5Al3"}], **inv)
ter = rec(value=880.0, temperature_K=880.0, transition_kind="other", locator=dict(loc, row=5),
          phases=[{"label": "L"}, {"label": "NiAl"}, {"label": "Ni3Al"}, {"label": "τ"}], **inv)
deg = rec(value=850.0, temperature_K=850.0, transition_kind="polymorphic", locator=dict(loc, row=6),
          phases=[{"label": "B27"}, {"label": "B33"}, {"label": "NiAl"}, {"label": "Ni3Al"}], **inv)
old_mis = rec(value=740.0, temperature_K=740.0, transition_kind="eutectic", locator=dict(loc, row=7),
              phases=[{"label": "L"}, {"label": "Mg2Pb"}, {"label": "Phase 2"}, {"label": "HCP_A3 (Mg)"}], **inv)
tchk = vc(kind="invariant_table", verdict="checked", drop_labels={"Phase 2": "composition column header", "Phase 3": "composition column header"},
          records={miss["id"]: {"action": "add_labels", "labels": ["Ni3Al"]}, eo["id"]: {"action": "reaction_type", "type": "eutectoid", "printed": "EO"},
                   ter["id"]: {"action": "higher_order", "system": "AL-CR-NI"}, deg["id"]: {"action": "as_printed"}})
out, rep = run([hdr, mono, miss, eo, ter, deg, old_mis], {"invariant_table:doi:10.1/t/p5/tab-2": tchk,
                                                          "record:" + old_mis["id"]: vc(kind="record_phases", verdict="misread")})
h = out[hdr["id"]]
check("G6 page check: header labels dropped, stored list and composition kept aside, no flag",
      [p["label"] for p in h["phases"]] == ["L", "Mg2Pb", "HCP_A3 (Mg)"] and h["phases_stored"][2] == "Phase 2"
      and h["phase_compositions_dropped"][0]["phase"] == "Phase 2" and not codes(h)
      and any(c["rule"] == "non_phase_labels_dropped_on_page" and c.get("visual_check") for c in h["corrections"]))
check("G6 monotectic with Liquid#1 and Liquid#2 is consistent (both kept)", len(out[mono["id"]]["phases"]) == 3 and not codes(out[mono["id"]]))
check("G6 page check: missing printed phase added", [p["label"] for p in out[miss["id"]]["phases"]] == ["L", "NiAl", "Ni3Al"]
      and not codes(out[miss["id"]]))
check("G6 page check: reaction type set to the printed one, stored kept",
      out[eo["id"]]["transition_kind"] == "eutectoid" and out[eo["id"]]["transition_kind_stored"] == "eutectic" and not codes(out[eo["id"]]))
check("G6 page check: higher-order row reassigned to the named system", out[ter["id"]]["system"] == "AL-CR-NI"
      and out[ter["id"]]["flags"]["extraction_accepted"])
check("G6 page check: degenerate four-phase row as printed stays accepted, not ESPEI-ready",
      "binary_degenerate_reaction_as_printed" in codes(out[deg["id"]]) and out[deg["id"]]["flags"]["extraction_accepted"]
      and not out[deg["id"]]["flags"]["espei_ready"])
check("page: a record-level 'misread' of the phase list is not demoted once the table check corrects it",
      out[old_mis["id"]]["flags"]["extraction_accepted"] and not codes(out[old_mis["id"]]))
pf = rec(property="HM_FORM", unit="kJ/mol-atom", value=196.18, value_printed=1373.27, phases=[{"label": "Mn3O4"}], system="MN-O",
         temperature_K=298.0, temperature_printed="298 K", flags={"extraction_accepted": True, "espei_ready": True, "magnitude_suspect": True})
st = rec(property="SM", unit="J/mol-atom/K", value=30.0, temperature_K=298.0, temperature_printed=None, phases=[{"label": "NiAl"}])
asp = rec(property="HM_FORM", unit="kJ/mol-atom", value=-260.0, phases=[{"label": "NiAl"}])
out, rep = run([pf, st, asp], {"page_fix:" + pf["id"]: vc(kind="page_fix", verdict="corrected", fields={"value": -196.18}, printed="-1373.27"),
                               "record:" + asp["id"]: vc(kind="record_value", verdict="as_stored")})
check("page_fix: value read on the page replaces the stored one (kept), magnitude flag cleared",
      out[pf["id"]]["value"] == -196.18 and out[pf["id"]]["value_stored"] == 196.18 and not out[pf["id"]]["flags"]["magnitude_suspect"]
      and out[pf["id"]]["value_printed"] == 1373.27)
check("298 K is 298.15 K (printed string kept, correction entry)", out[pf["id"]]["temperature_K"] == 298.15
      and out[st["id"]]["temperature_K"] == 298.15 and out[pf["id"]]["temperature_printed"] == "298 K"
      and any(c["rule"] == "standard_reference_temperature" for c in out[st["id"]]["corrections"]))
check("G2 out-of-range value printed so: accepted, flagged, not ESPEI-ready",
      out[asp["id"]]["flags"]["extraction_accepted"] and "value_outside_typical_range_as_printed" in codes(out[asp["id"]])
      and not out[asp["id"]]["flags"]["espei_ready"])
g7a = rec(property="HM_FORM", unit="kJ/mol-atom", value=-12.0, phases=[{"label": "\\sigma"}], locator={"artifact": "table", "asset": "doi:x/p9/tab-1", "page": 9, "row": 1})
g7b = copy.deepcopy(g7a); g7b["id"], g7b["phases"] = "ORC-SIG2", [{"label": "σ"}]
g7c = copy.deepcopy(g7a); g7c["id"], g7c["phases"], g7c["locator"] = "ORC-FE3O4", [{"label": "Magnetite (Fe3O4)"}], dict(g7a["locator"], row=2)
g7d = copy.deepcopy(g7c); g7d["id"], g7d["phases"] = "ORC-FE3O4B", [{"label": "Fe3O4"}]
g7e = copy.deepcopy(g7a); g7e["id"], g7e["phases"], g7e["locator"] = "ORC-PD7", [{"label": "Pd7Er"}], dict(g7a["locator"], row=3)
g7f = copy.deepcopy(g7e); g7f["id"], g7f["phases"] = "ORC-PD3", [{"label": "Pd3Er"}]
out, rep = run([g7a, g7b, g7c, g7d, g7e, g7f])
check("G7 spelling variants of one label merged (\\sigma / σ; Fe3O4 / Magnetite (Fe3O4))",
      ("ORC-SIG2" not in out or g7a["id"] not in out) and ("ORC-FE3O4" not in out or "ORC-FE3O4B" not in out))
check("G7 two different readings of one printed value kept apart and flagged",
      "phase_label_readings_disagree" in codes(out["ORC-PD7"]) and "phase_label_readings_disagree" in codes(out["ORC-PD3"]))
cm = rec(property="PHASE_COMPOSITION", observation_type="phase_region_vertex_set", temperature_K=1200.0, record_kind="general_table",
         phases=[{"label": "NiAl"}], phase_compositions=[{"phase": "NiAl", "composition": {"AL": 0.30, "NI": 0.70},
                                                          "composition_source": "printed", "basis": "mole_fraction"}])
out, rep = run([cm])
check("G9 composition deviating from the formula (printed solubility) is flag-only",
      "composition_label_mismatch" in codes(out[cm["id"]]) and out[cm["id"]]["flags"]["espei_ready"])
lq = rec(property="CPM", unit="J/mol-atom/K", value=42.0, temperature_K=298.15, system="CA-O", phases=[{"label": "CaO(liq-am)"}])
out, rep = run([lq])
check("G5 liquid / amorphous property at room temperature not ESPEI-ready", "liquid_property_at_room_temperature" in codes(out[lq["id"]]))

cz = rec(property="HM_FORM", unit="kJ/mol-atom", value=-24.4, temperature_K=1123.0, system="CU-ZR", phases=[{"label": "CuZr"}],
         uncertainty={"value": 0.96, "basis": "printed"})
out, rep = run([cz], {"page_fix:" + cz["id"]: vc(kind="page_fix", verdict="corrected", fields={
    "temperature_K": 298.15, "uncertainty": {"value": 8.749, "basis": "scatter_between_sources"}}),
    "flag:" + cz["id"]: vc(kind="page_flag", verdict="confirmed", related_records=["ORC-X"], flags=[
        {"code": "measured_at_other_stoichiometry", "effect": "not_espei_ready", "detail": "x"},
        {"code": "same_measurement_listed_in_other_paper", "effect": "not_espei_ready", "detail": "y"}])})
o = out[cz["id"]]
check("page_fix: temperature and scatter-aware sigma replace the stored ones (kept); several page flags apply",
      o["temperature_K"] == 298.15 and o["temperature_K_stored"] == 1123.0 and o["uncertainty"]["value"] == 8.749
      and o["uncertainty_stored"]["value"] == 0.96 and {"measured_at_other_stoichiometry", "same_measurement_listed_in_other_paper"} <= codes(o)
      and o["related_records"] == ["ORC-X"] and not o["flags"]["espei_ready"])

# PI fix 2026-10-24: a value printed against a non-SER reference phase is converted on the page's evidence (printed value
# kept, conversion in the correction); a historical formula of the same phase is flag-only with an inflated sigma;
# primary works resolved from the reference list
z6 = rec(property="HM_FORM", unit="kJ/mol-atom", value=-10.389, value_printed=-10389, value_stored=-10389.0,
         temperature_K=1123.0, system="CU-ZR", phases=[{"label": "Cu5Zr"}], uncertainty={"value": 0.077, "basis": "printed"},
         reference_state={"CU": {"phase": "FCC_A1"}, "ZR": {"phase": "HCP_A3"}}, cited_sources=[{"label": "46", "status": "unresolved"}])
conv = {"printed_value": -10389.0, "x_Zr": 0.181818, "lattice_stability": {"value_J_per_mol": 4206.77}}
out, rep = run([z6], {"page_fix:" + z6["id"]: vc(kind="page_fix", verdict="corrected", fields={
    "value": -9.624, "reference_state": {"CU": {"phase": "FCC_A1", "fixed_state_variables": {"T": 1123.0}},
                                         "ZR": {"phase": "HCP_A3", "fixed_state_variables": {"T": 1123.0}}},
    "uncertainty": {"value": 0.788, "basis": "printed_with_composition_offset"},
    "cited_sources": [{"label": "46", "doi": "10.1023/a:1022372503921", "status": "resolved"}]},
    rules={"value": "converted_to_ser_reference_state"}, conversions={"value": conv}),
    "flag:" + z6["id"]: vc(kind="page_flag", verdict="confirmed", flags=[
        {"code": "stoichiometry_mismatch", "effect": "flag_only", "detail": "historical formula"},
        {"code": "reference_state_converted_to_ser", "effect": "flag_only", "detail": "converted"}])})
o = out[z6["id"]]
cv = [c for c in o.get("corrections") or [] if c["field"] == "value"]
check("page_fix: converted value with its conversion in the correction; printed value kept; flag-only flags keep "
      "espei_ready; cited source resolved (unresolved entry kept in the correction)",
      o["value"] == -9.624 and o["value_printed"] == -10389 and o["value_stored"] == -10389.0
      and cv and cv[0]["before"] == -10.389 and cv[0]["rule"] == "converted_to_ser_reference_state"
      and cv[0]["evidence"]["conversion"] == conv and o["flags"]["espei_ready"]
      and o["uncertainty"]["basis"] == "printed_with_composition_offset" and o["uncertainty_stored"]["value"] == 0.077
      and o["cited_sources"][0]["doi"] == "10.1023/a:1022372503921" and "cited_sources_stored" not in o
      and [c for c in o["corrections"] if c["field"] == "cited_sources"][0]["before"][0]["status"] == "unresolved"
      and {"stoichiometry_mismatch", "reference_state_converted_to_ser"} <= codes(o))

# ------------------------------------------------------------------ G4 units / cover gases as elements (James Saal review 2026-10-24)
lani = rec(system="LA-NI-PA", property="TRANSITION_TEMPERATURE", record_kind="text_passage",
           quote="LaNi5H7 decomposed into LaNi5 and H2 at 347 and 379 K under 1 and 3 MPa, respectively.")
argas = rec(system="AR-PU", property="TRANSITION_TEMPERATURE", record_kind="text_passage",
            quote="Under argon (in green), a liquid temperature was observed in an Ar atmosphere at 3050 K.")
paox = rec(system="O-PA", property="TRANSITION_TEMPERATURE", record_kind="text_passage",
           quote="PaO2 melts at 2830 K at 101 325 Pa.")
cute = rec(system="CU-HE-TE", property="INVARIANT_TEMPERATURE", record_kind="text_passage",
           quote="the invariant reaction L + HeCu4Te3 <-> HeCuTe")
fixed = rec(system="CU-HE-TE", property="INVARIANT_TEMPERATURE", record_kind="text_passage",
            quote="the invariant reaction L + HeCu4Te3 <-> HeCuTe")
out, rep = run([lani, argas, paox, cute, fixed],
               visual={"record:" + fixed["id"]: vc(key="record:" + fixed["id"], verdict="misread", system="CU-TE")})
check("G4u Pa from MPa: not accepted", "system_element_from_unit_or_gas" in codes(out[lani["id"]])
      and not out[lani["id"]]["flags"]["extraction_accepted"])
check("G4u Ar cover gas: not accepted", "system_element_from_unit_or_gas" in codes(out[argas["id"]]))
check("G4u real protactinium (PaO2) with a Pa unit elsewhere: kept", "system_element_from_unit_or_gas" not in codes(out[paox["id"]]))
check("G4u garbled He prefix without a page check: not accepted", "system_element_from_unit_or_gas" in codes(out[cute["id"]]))
check("G4u page check reassigns He record to Cu-Te and the unit gate passes",
      out[fixed["id"]]["system"] == "CU-TE" and out[fixed["id"]]["system_stored"] == "CU-HE-TE"
      and "system_element_from_unit_or_gas" not in codes(out[fixed["id"]]))

print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
