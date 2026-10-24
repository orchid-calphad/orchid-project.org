#!/usr/bin/env python3
"""ESPEI bundle inclusion gate (tools/build_espei_bundle.py, cleanup 2026-10-24): every value the literature release
marks not ESPEI-ready is excluded; 298 K is the standard temperature 298.15 K; no solver dataset lies below ESPEI's
298.15 K temperature filter; tiny printed sigmas use the resolution floor. Unit part always; release part when the
orchid-data checkout is present (no written dataset carries a not_espei_ready value; no solver dataset below 298.15 K).
usage: python3 tools/tests/test_espei_bundle_gate.py"""
import glob, gzip, json, os, re, sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import build_espei_bundle as B  # noqa: E402
import orchid_cite  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


ok = {"evidence_class": "experimental", "flags": {"espei_ready": True}}
g6 = {"evidence_class": "experimental", "flags": {"espei_ready": False},
      "gate_failures": [{"gate": "G6", "code": "reaction_phase_count_inconsistent", "effect": "not_espei_ready"}]}
mod = {"evidence_class": "model", "flags": {"espei_ready": False}}
modg = {"evidence_class": "model", "flags": {}, "gate_failures": [{"gate": "G1", "code": "temperature_range_stored_as_point",
                                                                     "effect": "not_accepted"}]}
flag = {"evidence_class": "model", "flags": {}, "gate_failures": [{"gate": "G9", "code": "composition_label_mismatch",
                                                                     "effect": "flag_only"}]}
nr = {"evidence_class": "experimental", "flags": {"espei_ready": False}}
check("espei_ready experimental value passes", B.gate_blocks(ok) is None)
check("G6-flagged experimental value excluded", B.gate_blocks(g6) == "not_espei_ready:literature_gate")
check("experimental value not espei_ready excluded", B.gate_blocks(nr) == "phase_rule:not_espei_ready")
check("model value without blocking gate passes; with one is excluded",
      B.gate_blocks(mod) is None and B.gate_blocks(modg) == "not_espei_ready:literature_gate")
check("flag-only gate does not exclude", B.gate_blocks(flag) is None)
check("every gate drop reason is documented", all(k in B.DROP_REASONS for k in (
    "not_espei_ready:literature_gate", "phase_rule:not_espei_ready", "below_espei_temperature_floor", "duplicate_of_other_dataset")))
check("298 K is 298.15 K", B.std_T(298) == 298.15 and B.std_T([298.0, 500.0]) == [298.15, 500.0] and B.std_T(298.5) == 298.5)
check("298 K dataset matches a 298.15 K published value", B.check_temperature(298.15, {"temperature_K": 298.15}, None)
      and B.check_temperature(298.0, {"temperature_K": 298.15}, None))

D = orchid_cite.DEFAULT_DATA
E = os.path.join(D, "espei")
lit = os.path.join(D, "literature", "orchid-literature.jsonl.gz")
if os.path.isdir(E) and os.path.exists(lit):
    pub = {}
    for line in gzip.open(lit, "rt", encoding="utf-8"):
        r = json.loads(line)
        pub[r["id"]] = r
    bad, low, n = [], [], 0
    for f in glob.glob(os.path.join(E, "*", "*", "*.json")):
        rel = os.path.relpath(f, E)
        if rel.split(os.sep)[0] not in ("experimental", "model", "invariants"):
            continue
        d = json.load(open(f, encoding="utf-8"))
        n += 1
        for o in (d.get("orchid") or {}).get("observations") or []:
            r = pub.get(o.get("id"))
            if r is not None and B.gate_blocks(r):
                bad.append((rel, o["id"]))
        if d.get("solver"):
            T = d["conditions"]["T"]
            if any(t < 298.15 for t in (T if isinstance(T, list) else [T])):
                low.append(rel)
    check("release: no ESPEI dataset value comes from a record the release marks not ESPEI-ready (%d datasets; %d bad)"
          % (n, len(bad)), not bad)
    check("release: no solver dataset below ESPEI's 298.15 K filter (%d)" % len(low), not low)
print("%d failures" % len(fails))
sys.exit(1 if fails else 0)
