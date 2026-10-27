#!/usr/bin/env python3
"""Independent source check ledger (ORCHID-IC-20261027-v1): listed records are not ESPEI-ready, values untouched.
usage: python3 tools/tests/test_independent_check.py"""
import gzip
import json
import os
import sys

sys.dont_write_bytecode = True
TOOLS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, TOOLS)
import independent_check_review as IC  # noqa: E402

fails = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + ("" if cond or not detail else "  [%s]" % (detail,)))
    if not cond:
        fails.append(name)


rec = {"id": "ORC-TESTIC", "value": 1.0, "flags": {"extraction_accepted": True, "espei_ready": True, "espei_certified": True},
       "gate_failures": []}
ids = list(IC.REVIEW["records"])
rec["id"] = ids[0]
before = dict(rec, flags=dict(rec["flags"]))
n = IC.annotate([rec])
check("a listed record is not ESPEI-ready or certified after annotation", rec["flags"]["espei_ready"] is False
      and rec["flags"]["espei_certified"] is False and n >= 1)
check("the record stays accepted and its value is untouched", rec["flags"]["extraction_accepted"] is True
      and rec["value"] == before["value"])
other = {"id": "ORC-NOT-LISTED", "flags": {"extraction_accepted": True, "espei_ready": True}, "gate_failures": []}
IC.annotate([other])
rec2 = {"id": ids[0], "flags": {"extraction_accepted": True, "espei_ready": True, "espei_certified": True},
        "gate_failures": []}
IC.annotate([rec2])
rec2["flags"]["espei_ready"] = True          # upstream recomputation resets the readiness flag
rec2["flags"]["espei_certified"] = True
IC.annotate([rec2])
check("repeated annotation re-enforces readiness and certification", rec2["flags"]["espei_ready"] is False
      and rec2["flags"]["espei_certified"] is False)
check("the flag mirror of gate codes is synchronized", set(g["code"] for g in rec2["gate_failures"]) <= set(rec2["flags"]["gate_failures"]))
for cleared in ("ORC-4OAJ3KNMAB", "ORC-INPXMODQVD", "ORC-YY4Y7EGCLN"):
    check("record cleared on the page is not in the ledger as a confirmed defect: " + cleared,
          "independent_check_confirmed_defect" not in IC.REVIEW["records"].get(cleared, {}).get("codes", []))
check("an unlisted record is untouched", other["flags"]["espei_ready"] is True and not other["gate_failures"])
check("ledger has the expected size", len(ids) >= 500, str(len(ids)))
path = os.environ.get("ORCHID_LITERATURE") or os.path.join(os.path.expanduser("~"), "projects", "orchid-data", "literature",
                                                           "orchid-literature.jsonl.gz")
if os.path.exists(path):
    bad = []
    seen = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            x = json.loads(line)
            if x["id"] in IC.REVIEW["records"]:
                seen += 1
                if (x.get("flags") or {}).get("espei_ready") or (x.get("flags") or {}).get("espei_certified"):
                    bad.append(x["id"])
    check("release: every ledger id is present (%d of %d)" % (seen, len(IC.REVIEW["records"])), seen == len(IC.REVIEW["records"]))
    check("release: no listed record is ESPEI-ready or certified", not bad, "%d e.g. %s" % (len(bad), bad[:3]))
print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
