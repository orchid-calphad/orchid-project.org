"""Independent source check of ESPEI-ready records (ORCHID-IC-20261027-v1).

A second reviewer checked 386 ESPEI-ready records on the printed pages. Records with a confirmed defect, and members of
the systematic classes that the check found (screening candidates), keep their values and stay in the release but are
not ESPEI-ready until corrected or adjudicated. This module only adds a gate failure; it never changes a value."""
import json
import os

import literature_gates as G

GATE = "IC20261027"
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "independent_check_20261027.json"), encoding="utf-8") as fh:
    REVIEW = json.load(fh)


def annotate(records):
    n = 0
    for r in records:
        e = REVIEW["records"].get(r["id"])
        if not e:
            continue
        for code in e["codes"]:
            if G.fail(r, GATE, code, "not_espei_ready", "; ".join(e["notes"])[:400]):
                n += 1
        # G.fail returns early when the code is already present: enforce the flag effect on a replay too
        fl = r.setdefault("flags", {})
        if fl.get("espei_ready"):
            fl["espei_ready_before_gates"] = True
        fl["espei_ready"] = False
        fl["espei_certified"] = False
        mirror = fl.setdefault("gate_failures", [])
        for g in r.get("gate_failures") or []:
            if g.get("code") and g["code"] not in mirror:
                mirror.append(g["code"])
    return n
