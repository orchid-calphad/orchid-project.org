#!/usr/bin/env python3
"""Uncertainty budget (tools/uncertainty_budget.py; PI decision 2026-10-24: method-based floors, components kept
separate). Unit part always (fixture records): the floor is applied, never below ESPEI's experimental default,
printed sigma and printed value preserved, the rule combines the components as documented. Release part when the
orchid-data checkout is present: every released value has a budget made with the published calibration, every
quantity carries every component field, sigma_effective >= every max-combined component, the ESPEI copy of the
calibration is byte-identical and every ESPEI thermochemical dataset weight follows sigma_effective.
usage: python3 tools/tests/test_uncertainty_budget.py"""
import copy
import glob
import gzip
import json
import math
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import uncertainty_budget as U  # noqa: E402
import orchid_cite  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


def rec(i, src, v, method="calorimetry", unit="kJ/mol-atom", sigma=None, comp=None, prop="HM_FORM", T=298.15):
    r = {"id": "R%d" % i, "property": prop, "value": v, "unit": unit, "evidence_class": "experimental",
         "flags": {"extraction_accepted": True}, "elements": ["AL", "NI"], "system": "AL-NI",
         "overall_composition": comp or {"AL": 0.5, "NI": 0.5}, "temperature_K": T, "method": method,
         "origin": {"measured_in": src, "own_work_of_paper": True}, "paper": src, "phases": [{"label": "B2"}]}
    if sigma is not None:
        r["uncertainty"] = {"basis": "printed", "value": sigma}
    return r


# fixture: 12 compounds, 3 calorimetric sources each, scatter about 2 kJ/mol-atom
recs, k = [], 0
for c in range(12):
    x = round(0.2 + 0.05 * c, 3)
    for s, off in (("A", 0.0), ("B", 2.0 + 0.1 * c), ("C", -1.5)):
        k += 1
        recs.append(rec(k, "paper%s%d" % (s, c), -40.0 + c + off, sigma=0.08,
                        comp={"AL": x, "NI": round(1 - x, 3)}))
recs.append(rec(999, "paperX", -30.0, method=None, unit="kJ/mol-atom"))           # no method, no sigma
recs.append(rec(998, "paperY", -30.0, method="calorimetry", unit="J/g/K", prop="CPM"))  # mass basis
recs.append(rec(997, "paperZ", -0.1, method="calorimetry", sigma=0.01, comp={"AL": 0.9, "NI": 0.1}))           # tiny printed sigma
orig = copy.deepcopy(recs)
cal = U.calibrate(recs, None)
U.annotate(recs, cal)
e = cal["value_types"]["HM_FORM"]["classes"]["solution_reaction_calorimetry"]
check("calorimetry class calibrated from 36 same-class pairs of 12 states", e["basis"] == "calibrated"
      and e["n_pairs"] == 36 and e["n_states"] == 12)
check("floor = robust sigma rounded up to 0.1, >= ESPEI default 0.5", e["floor"] >= 0.5 and
      abs(e["floor"] - U.ceil_step(e["sigma_robust"], 0.1)) < 1e-9 and e["ci95"][0] <= e["sigma_robust"] <= e["ci95"][1])
q = recs[0]["uncertainty_budget"]["quantities"]["value"]
check("printed sigma kept as printed (0.08) and the floor applied", q["sigma_printed"]["value"] == 0.08 and
      q["sigma_method"]["value"] == e["floor"] and q["sigma_effective"] == e["floor"])
check("printed values and printed uncertainties unchanged", all(
    r["value"] == o["value"] and r.get("uncertainty") == o.get("uncertainty") for r, o in zip(recs, orig)))
check("every record has a budget with every component field", all(
    set(("sigma_printed", "sigma_reading", "sigma_method", "sigma_composition", "sigma_scatter", "sigma_conversion",
         "sigma_effective")) <= set(qq) for r in recs for qq in r["uncertainty_budget"]["quantities"].values()))
qn = recs[-3]["uncertainty_budget"]["quantities"]["value"]
check("no method printed: class unknown, prior (no pairs), never below 0.5 kJ/mol-atom",
      qn["sigma_method"]["class"] == "unknown" and qn["sigma_method"]["basis"] == "prior" and qn["sigma_effective"] >= 0.5)
qm = recs[-2]["uncertainty_budget"]["quantities"]["value"]
check("mass-basis unit: no sigma_effective, reason given", qm["sigma_effective"] is None and
      "not convertible" in qm["sigma_effective_reason"])
qt = recs[-1]["uncertainty_budget"]["quantities"]["value"]
check("tiny printed sigma (0.01) raised to the floor", qt["sigma_effective"] >= 0.5 and qt["sigma_printed"]["value"] == 0.01)
# rule
qq = {"sigma_printed": {"value": 1.0}, "sigma_reading": {"value": 0.05}, "sigma_method": {"value": 3.0},
      "sigma_composition": {"value": 4.0}, "sigma_scatter": None, "sigma_conversion": None}
check("rule: max(printed, reading, method) in quadrature with composition", abs(U.combine(qq) - 5.0) < 1e-9)
qq["sigma_scatter"] = {"value": 8.0}
check("rule: scatter replaces when larger", U.combine(qq) == 8.0)
qq["sigma_printed"] = {"value": 9.0}
qq["sigma_composition"] = None
check("rule: printed above the floor is kept (no double counting)", U.combine(qq) == 9.0)
check("classes: KEMS second law, EMF, combustion, DTA temperature, EPMA composition",
      U.method_class({"method": "Vapor-pressure measurements (2nd law)"}, "HM_FORM") == "kems_vapour_second_law" and
      U.method_class({"method": "emf 3rd law"}, "HM_FORM") == "emf" and
      U.method_class({"method": "Fluorine bomb calorimetry"}, "HM_FORM") == "combustion_calorimetry" and
      U.method_class({"method": "DTA, XRD, optical microscopy, SEM, EPMA"}, "T_INVARIANT") == "thermal_analysis" and
      U.method_class({"method": "DTA, XRD, optical microscopy, SEM, EPMA"}, "X_PHASE") == "microanalysis")
check("basis-error pair (ratio = formula atom count) excluded",
      U.basis_error_pair(("HM_FORM", (("MG", 0.667), ("SI", 0.333)), "std"), -79.1, -26.0))

# ---------------------------------------------------------------------------------------------- release --
data = orchid_cite.DEFAULT_DATA
lit = os.path.join(data, "literature", "orchid-literature.jsonl.gz")
if os.path.exists(lit):
    cb = open(os.path.join(data, "literature", "UNCERTAINTY_CALIBRATION.json"), "rb").read()
    calr = json.loads(cb)
    cid = U.calibration_id(calr)
    n = bad = nobud = below = 0
    for line in gzip.open(lit, "rt", encoding="utf-8"):
        r = json.loads(line)
        n += 1
        b = r.get("uncertainty_budget")
        if not b or b.get("calibration") != cid or "quantities" not in b:
            nobud += 1
            continue
        if U.quantities(r) and set(b["quantities"]) != {qn for qn, _ in U.quantities(r)}:
            bad += 1
        for qn_, q in b["quantities"].items():
            comps = [q.get(k) for k in ("sigma_printed", "sigma_reading", "sigma_method")]
            vals = [c["value"] for c in comps if c and isinstance(c.get("value"), (int, float))]
            se = q.get("sigma_effective")
            if se is not None and vals and se < max(vals) * (1 - 1e-3) - 1e-12:
                below += 1
            if se is None and not q.get("sigma_effective_reason"):
                bad += 1
            m = q.get("sigma_method") or {}
            if (r.get("evidence_class") == "experimental" and q.get("value_type") in ("HM_FORM", "GM_FORM") and
                    m.get("value") is not None and U.energy_factor(r.get("unit")) and
                    m["value"] * U.energy_factor(r.get("unit")) < 0.5 - 1e-9):
                below += 1
    check("release: every one of %d records has a budget of calibration %s" % (n, cid), nobud == 0)
    check("release: budget quantities match the record's measured quantities; reasons where no sigma", bad == 0)
    check("release: sigma_effective >= max(printed, reading, method); formation floors >= 0.5 kJ/mol-atom", below == 0)
    ec = os.path.join(data, "espei", "experimental", "CALIBRATION_EXPERIMENT.json")
    if os.path.exists(ec):
        check("ESPEI calibration copy byte-identical to the literature calibration", open(ec, "rb").read() == cb)
        wbad = nds = 0
        for f in glob.glob(os.path.join(data, "espei", "experimental", "*", "*.json")) + \
                glob.glob(os.path.join(data, "espei", "model", "*", "*.json")):
            d = json.load(open(f, encoding="utf-8"))
            if d.get("output") == "ZPF":
                continue
            u = d.get("orchid_uncertainty") or {}
            se = [e.get("sigma_effective_used") for e in u.get("values") or []]
            if not se or any(x is None for x in se) or "weight" not in d:
                wbad += 1
                continue
            nds += 1
            s = sorted(se)
            med = s[len(s) // 2] if len(s) % 2 else 0.5 * (s[len(s) // 2 - 1] + s[len(s) // 2])
            if abs(d["weight"] - round(u["default_sigma"]["value"] / med, 6)) > 1e-6 or \
                    any("budget" not in e for e in u["values"]):
                wbad += 1
        check("ESPEI: all %d thermochemical datasets weighted by sigma_effective, budgets per value" % nds, wbad == 0)
print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
