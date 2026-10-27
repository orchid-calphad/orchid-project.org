#!/usr/bin/env python3
"""Audit round 4 (2026-10-27): a state needs real phases. R3 tested only for an empty phase list, so placeholder lists
(['UNKNOWN'], ['other'], ['unspecified']) passed as 'phases recorded'. Unit part: the helper and the gate on fixtures.
Release part: no accepted record of a state-phase property family carries only placeholder phases.
usage: python3 tools/tests/test_r4_placeholder_phases.py"""
import gzip
import json
import os
import sys

sys.dont_write_bytecode = True
TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, TOOLS)
import literature_r3_thermo as T  # noqa: E402

fails = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + ("" if cond or not detail else "  [%s]" % (detail,)))
    if not cond:
        fails.append(name)


def rec(prop, labels, kind="general_table", quote=None):
    return {"id": "ORC-TEST", "property": prop, "record_kind": kind, "quote": quote, "system": "AL-NI",
            "elements": ["AL", "NI"], "phases": [{"label": x} for x in labels],
            "flags": {"extraction_accepted": True, "espei_ready": False}, "gate_failures": []}


class Rep:
    def __init__(self):
        self.rows = []

    def add(self, *a, **k):
        self.rows.append((a, k))


for labs in (["UNKNOWN"], ["other"], ["unspecified"], ["unknown", "Other"], []):
    check("no real phases: %r" % (labs,), not T.has_real_phases(rec("HM_FORM", labs)))
for labs in (["AlNi"], ["L", "UNKNOWN"], ["FCC_A1"], ["Al3Ni2", "other"]):
    check("real phase present: %r" % (labs,), T.has_real_phases(rec("HM_FORM", labs)))

check("sodium is a real phase label, not a placeholder ('Na', 'NA', 'na')", all(T.has_real_phases(rec("HM_FORM", [x]))
      for x in ("Na", "NA", "na")))
for labs in (["N/A"], [" unknown "], ["-"], [None]):
    check("placeholder variants: %r" % (labs,), not T.has_real_phases(rec("HM_FORM", labs)))
for prop in sorted(T.STATE_PHASE_PROPS):
    r = rec(prop, ["UNKNOWN"])
    T.r3_state_phases([r], {}, Rep())
    check("state-phase family %s: placeholder-only is not accepted" % prop, r["flags"]["extraction_accepted"] is False)
r = rec("HM_FORM", ["Na"])
T.r3_state_phases([r], {}, Rep())
check("an elemental sodium record is not demoted", r["flags"]["extraction_accepted"] is True and not r["gate_failures"])
r = rec("HM_FORM", ["UNKNOWN"]); r["flags"]["extraction_accepted"] = False
T.r3_state_phases([r], {}, Rep())
check("an already rejected record keeps a single gate entry", sum(g["code"] == "state_phases_missing" for g in r["gate_failures"]) <= 1)
for prop in ("TRANSITION_TEMPERATURE", "HM_FORM", "INVARIANT_TEMPERATURE"):
    r = rec(prop, ["UNKNOWN"])
    T.r3_state_phases([r], {}, Rep())
    check("%s with placeholder phases only is not accepted" % prop,
          r["flags"]["extraction_accepted"] is False and any(g["code"] == "state_phases_missing" for g in r["gate_failures"]))
r = rec("LATTICE_A", ["UNKNOWN"])
T.r3_state_phases([r], {}, Rep())
check("a property outside the state-phase families is untouched (lattice constant)", r["flags"]["extraction_accepted"] is True)
r = rec("HM_FORM", ["AlNi"])
T.r3_state_phases([r], {}, Rep())
check("a real phase list is untouched", r["flags"]["extraction_accepted"] is True and not r["gate_failures"])

# release part
path = os.environ.get("ORCHID_LITERATURE") or os.path.join(os.path.expanduser("~"), "projects", "orchid-data", "literature",
                                                           "orchid-literature.jsonl.gz")
if os.path.exists(path):
    bad = []
    n = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            x = json.loads(line)
            n += 1
            if x.get("property") in T.STATE_PHASE_PROPS and (x.get("flags") or {}).get("extraction_accepted") \
                    and not T.has_real_phases(x):
                bad.append(x["id"])
    check("release: no accepted state-phase record has placeholder-only phases (%d records read)" % n, not bad,
          "%d, e.g. %s" % (len(bad), bad[:3]))
else:
    print("SKIP release part (no literature file at %s)" % path)
print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
