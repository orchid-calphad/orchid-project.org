"""Placeholder recovery lifts the R3 placeholder-phase gate only where it was the record's sole acceptance failure."""
import sys, copy
from pathlib import Path
TOOLS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS))
import importlib.util
spec = importlib.util.spec_from_file_location("exp", TOOLS / "export_literature_data.py")
exp = importlib.util.module_from_spec(spec); spec.loader.exec_module(exp)
import literature_r3_thermo as r3t
fails = []
def check(name, ok):
    print(("PASS " if ok else "FAIL ") + name)
    if not ok: fails.append(name)
def rec(rid, phases, failures, pr=True):
    return {"id": rid, "property": "T_INVARIANT", "phases": [{"label": x} for x in phases],
            "corrections": [{"gate": "PR20261027", "field": "phases"}] if pr else [],
            "gate_failures": [{"gate": "R3", "code": c, "effect": e} for c, e in failures],
            "flags": {"extraction_accepted": False, "espei_ready": False, "espei_certified": False,
                      "gate_failures": [c for c, _ in failures]}}
recs = [rec("A", ["L", "Ni3Ti"], [("state_phases_missing", "not_accepted"), ("human_review_residual_concern", "not_espei_ready")]),
        rec("B", ["L"], [("state_phases_missing", "not_accepted"), ("overall_composition_sum_not_one", "not_accepted")]),
        rec("C", ["unknown"], [("state_phases_missing", "not_accepted")]),
        rec("D", ["L"], [("state_phases_missing", "not_accepted")], pr=False)]
before = copy.deepcopy(recs)
n = exp.placeholder_reaccept(recs, r3t.has_real_phases)
check("only the sole-failure record is accepted again", n == 1 and recs[0]["flags"]["extraction_accepted"] is True)
check("ESPEI readiness stays off", not recs[0]["flags"]["espei_ready"] and not recs[0]["flags"]["espei_certified"])
check("the lifted failure is recorded as a correction", recs[0]["corrections"][-1]["rule"] == "placeholder_recovery_lifts_state_phases_missing")
check("a record with another not_accepted failure stays unaccepted", recs[1]["flags"]["extraction_accepted"] is False)
check("a placeholder label is not a recovered phase", recs[2] == before[2])
check("records without the recovery gate are untouched", recs[3] == before[3])
check("replay changes nothing", exp.placeholder_reaccept(recs, r3t.has_real_phases) == 0)
print("%d failed" % len(fails) if fails else "all passed")
sys.exit(1 if fails else 0)
