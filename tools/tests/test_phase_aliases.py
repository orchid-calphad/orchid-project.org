#!/usr/bin/env python3
"""Tests of the per-system phase alias table (tools/phase_aliases.json, tools/phase_aliases.py; 2026-10-24):
label parsing, alias uniqueness, structure / composition consistency of every alias with the registry and the TDBs,
the Cu3Sn eps case (eps in Cu-Sn -> the Cu3Sn phase, never BCC_A2), the integration in tools/phase_label_check.py
(alias before the Greek rule, provenance, corrections, twin rule), the ESPEI ZPF relabel, and release invariants.

usage: python3 tools/tests/test_phase_aliases.py [--data ~/projects/orchid-data] [--units-only]
"""
import argparse
import copy
import gzip
import json
import os
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
TOOLS = os.path.dirname(HERE)
sys.path.insert(0, TOOLS)
import phase_aliases as PA  # noqa: E402
import phase_label_check as PLC  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--data", default=os.path.expanduser("~/projects/orchid-data"))
ap.add_argument("--tdb-dir", default=PA.DEFAULT_TDB_DIR)
ap.add_argument("--units-only", action="store_true")
A = ap.parse_args()
fails = []


def check(name, cond, detail=""):
    print(("PASS " if cond else "FAIL ") + name + ((" -- " + str(detail)[:300]) if not cond and detail else ""))
    if not cond:
        fails.append(name)


T = PA.Aliases()

# ------------------------------------------------------------------ label parsing
cases = {"ε": ("ε", None), "epsilon": ("ε", None), "ϵ": ("ε", None), "ε (Cu3Sn)": ("ε", "Cu3Sn"),
         "εCu3Sn": ("ε", "Cu3Sn"), "ε-Cu3Sn": ("ε", "Cu3Sn"), "Cu3Sn(ε)": ("ε", "Cu3Sn"), "γ'": ("γ′", None),
         "γ’": ("γ′", None), "γ′": ("γ′", None), "β1": ("β1", None), "β phase": ("β", None), "Γ1": ("Γ1", None),
         "(Cu)": None, "L": None, "Cu3Sn": None, "bcc": None}
bad = {k: (PA.label_key(k), v) for k, v in cases.items() if PA.label_key(k) != v}
check("label_key parses Greek core / prime / index / formula forms", not bad, bad)

# ------------------------------------------------------------------ uniqueness and curated consistency
probs = PA.validate(T)
check("validate(): alias keys unique, labels canonical, ranges sane, shared TDB phases only for disjoint ranges",
      not probs, probs[:5])
keys = [(a["system"], lab) for a in T.aliases for lab in [a["label"]] + list(a.get("label_variants") or [])]
check("alias uniqueness: one alias per (system, printed label)", len(keys) == len(set(keys)))
check("every alias has a basis, a model phase, >= 1 TDB phase and >= 1 literature record",
      all(a.get("basis") and a.get("model_phase") and a.get("tdb_phases") and
          a["evidence"]["literature_records_total"] >= 1 for a in T.aliases),
      [a["id"] for a in T.aliases if not a["evidence"]["literature_records_total"]])
check("no global Greek -> structure rule: every alias names one system", all(a["system"].count("-") == 1 for a in T.aliases))

# ------------------------------------------------------------------ structure / composition consistency (registry, TDBs)
text, probs = PA.refresh(data=A.data, tdb_dir=A.tdb_dir, write=False)
check("evidence current and consistent (structure, composition range, calculated field, one phase per TDB)",
      not probs and text == open(PA.ALIAS_FILE, encoding="utf-8").read(), probs[:5] or "stale: run --refresh")
reg = PA.load_registry(A.data, sorted({a["system"] for a in T.aliases}))
bad = []
for a in T.aliases:
    phases = {p["id"]: p for p in (reg.get(a["system"]) or {}).get("phases") or []}
    lo, hi = a["x_range"]
    for tdb, name in a["tdb_phases"].items():
        ms = [(rid, m) for rid in a["registry_ids"] for m in (phases.get(rid) or {}).get("models") or []
              if m["tdb"] == tdb and m["phase_name"] == name and not m.get("disordered_part_of")]
        if len(ms) != 1:
            bad.append((a["id"], tdb, name, "models %d" % len(ms)))
            continue
        rid, m = ms[0]
        r = PA.model_range(m, a["x_element"])
        if r and (r[1] < lo - PA.X_TOL or r[0] > hi + PA.X_TOL):
            bad.append((a["id"], tdb, name, "model range", r))
        st = phases[rid].get("structure") or {}
        if PA.pearson_equal((a.get("structure") or {}).get("pearson"), st.get("pearson")) is False:
            bad.append((a["id"], rid, "pearson", st.get("pearson")))
check("structure / composition: each TDB phase is one independent model reaching the alias range, Pearson agrees",
      not bad, bad[:5])

# ------------------------------------------------------------------ the Cu3Sn eps case
eps = T.by_id["cu-sn:epsilon"]
h = T.lookup("CU-SN", "ε", [{"SN": 0.25}], 723.0)
cu3sn = {"CU3SN_L", "CU3SN", "CU3SN_LT"}
check("Cu-Sn eps -> the Cu3Sn phase (CU3SN_L by default), not BCC_A2",
      h and h.get("target") == "CU3SN_L" and h["target"] != "BCC_A2", h)
check("Cu-Sn eps: every Cu-Sn TDB names a Cu3Sn phase; none is BCC_A2",
      set(eps["tdb_phases"].values()) <= cu3sn and "BCC_A2" not in eps["tdb_phases"].values()
      and len(eps["tdb_phases"]) == len(PA.tdb_files_of(A.tdb_dir, "CU-SN")), eps["tdb_phases"])
check("Cu-Sn eps in CuSn-12Gie -> CU3SN, in CuSn-14Don -> CU3SN_LT",
      T.lookup("CU-SN", "ε", [], 600, "CuSn-12Gie.tdb")["target"] == "CU3SN"
      and T.lookup("CU-SN", "ε", [], 600, "CuSn-14Don.tdb")["target"] == "CU3SN_LT")
check("Cu-Sn 'ε (Cu3Sn)' -> Cu3Sn phase; 'ε' at x(Sn) = 0.45 not mapped",
      T.lookup("CU-SN", "ε (Cu3Sn)", [], 723)["target"] == "CU3SN_L"
      and "miss" in T.lookup("CU-SN", "ε", [{"SN": 0.45}], 723))
check("aliases are per system: eps in Cu-Zn is HCP_IM_A3, eps in Fe-Ta is C14_LAVES, eps in Al-Fe is not BCC_A2",
      T.lookup("CU-ZN", "ε", [{"ZN": 0.8}], 700)["target"] == "HCP_IM_A3"
      and T.lookup("FE-TA", "ε", [{"TA": 0.3333}], 1070)["target"] == "C14_LAVES"
      and T.lookup("AL-FE", "ε", [{"FE": 0.391}], 1400)["target"] in ("AL8FE5_D82", "AL5FE4"))
check("a system without an alias is untouched (beta in Cu-Zr: no alias)", T.lookup("CU-ZR", "β", [], 1200) is None)

# ------------------------------------------------------------------ integration in phase_label_check (before the Greek rule)
rec = {"id": "TEST-1", "system": "CU-SN", "elements": ["CU", "SN"], "temperature_K": 949.0,
       "phases": [{"label": "γ", "model_phase": "CU3SN_D03"}, {"label": "ε", "model_phase": "BCC_A2",
                                                                "pearson_symbol": "cI2", "prototype": "W"}],
       "phase_compositions": [{"phase": "ε", "composition": {"SN": 0.25}, "composition_source": "printed",
                               "basis": "mole_fraction"}]}
r1 = copy.deepcopy(rec)
rep = PLC.apply([r1], None)
p = r1["phases"][1]
check("label check: eps BCC_A2 -> CU3SN_L with model_phase_basis system_alias:cu-sn:epsilon and a PL correction",
      p.get("model_phase") == "CU3SN_L" and p.get("model_phase_basis") == "system_alias:cu-sn:epsilon"
      and not p.get("model_phase_rejected") and p.get("structure_basis") == "system_alias"
      and any(c.get("rule") == "system_alias:cu-sn:epsilon" and c.get("before") == "BCC_A2" for c in r1["corrections"]),
      (p, r1.get("corrections")))
r2 = copy.deepcopy(rec)
r2["phase_compositions"][0]["composition"] = {"SN": 0.45}
PLC.apply([r2], None)
check("label check: eps outside the alias range -> alias not applied, Greek rule rejects BCC_A2",
      not r2["phases"][1].get("model_phase") and r2["phases"][1]["model_phase_rejected"]["model_phase"] == "BCC_A2")
r3 = {"id": "TEST-3", "system": "CU-ZN", "elements": ["CU", "ZN"], "temperature_K": 1176.0,
      "phases": [{"label": "(Cu)", "model_phase": "FCC_A1"}, {"label": "β", "model_phase": "BCC_A2"},
                 {"label": "L", "model_phase": "LIQUID"}],
      "phase_compositions": [{"phase": "β", "composition": {"ZN": 0.361}, "composition_source": "printed",
                              "basis": "mole_fraction"}]}
PLC.apply([r3], None)
check("label check: beta brass (Cu-Zn beta -> BCC_A2) kept, provenance system_alias:cu-zn:beta, no correction",
      r3["phases"][1].get("model_phase") == "BCC_A2" and r3["phases"][1].get("model_phase_basis") ==
      "system_alias:cu-zn:beta" and not r3.get("corrections"))
r4 = {"id": "TEST-4", "system": "CU-FE", "elements": ["CU", "FE"], "temperature_K": 1300.0,
      "phases": [{"label": "γ", "model_phase": "FCC_A1"}, {"label": "ε", "model_phase": "BCC_A2"}],
      "phase_compositions": [{"phase": "ε", "composition": {"FE": 0.03}, "composition_source": "printed",
                              "basis": "mole_fraction"}]}
PLC.apply([r4], None)
check("label check: an alias never puts two printed phases of one row on one model phase (Cu-Fe gamma + eps)",
      r4["phases"][1].get("model_phase") != "FCC_A1")

# ------------------------------------------------------------------ ESPEI ZPF relabel
try:
    import build_espei_bundle as EB  # noqa: E402
    region = [["CU3SN_D03", [["SN"], [0.2]]], ["BCC_A2", [["SN"], [0.25]]]]
    reg2, why = EB.alias_relabel(region, r1, "CuSn-12Gie")
    reg3, why3 = EB.alias_relabel(region, r1, "CuSn-96Shi-LB")
    check("ESPEI ZPF: the corpus vertex BCC_A2 of the aliased eps becomes the target TDB's Cu3Sn phase",
          why is None and reg2[1][0] == "CU3SN" and reg2[0][0] == "CU3SN_D03" and reg3[1][0] == "CU3SN_L", (reg2, why))
except ImportError as e:
    check("ESPEI bundle importable", False, e)

# ------------------------------------------------------------------ release invariants
pl = os.path.join(A.data, "literature", "orchid-literature.jsonl.gz")
if not A.units_only and os.path.exists(pl):
    bad, cusn, n = [], [], 0
    for line in gzip.open(pl, "rt", encoding="utf-8"):
        r = json.loads(line)
        for p in r.get("phases") or []:
            b = str(p.get("model_phase_basis") or "")
            if r["system"] == "CU-SN" and PA.label_key(p.get("label")) and PA.label_key(p.get("label"))[0] == "ε" \
                    and str(p.get("model_phase") or "").upper().startswith("BCC"):
                cusn.append(r["id"])
            if not b.startswith("system_alias:"):
                continue
            n += 1
            a = T.by_id.get(b.split(":", 1)[1])
            if a is None or p.get("model_phase") not in PA.Aliases.names(a):
                bad.append((r["id"], p.get("label"), p.get("model_phase"), b))
                continue
            for c in PA.record_compositions(r, p.get("label")):
                x = PA.x_of(c, a["x_element"])
                if x is not None and not (a["x_range"][0] - PA.X_TOL <= x <= a["x_range"][1] + PA.X_TOL):
                    bad.append((r["id"], p.get("label"), x))
    check("release: every system_alias mapping names a TDB phase of its alias at a consistent composition (%d)" % n,
          n > 0 and not bad, bad[:5])
    check("release: no Cu-Sn eps on a bcc model phase", not cusn, cusn[:5])

print("%d failed" % len(fails) if fails else "all passed")
sys.exit(1 if fails else 0)
