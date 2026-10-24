#!/usr/bin/env python3
"""Physical phase rule (tools/phase_rule.py), ESPEI phase assignment (tools/build_espei_bundle.py) and the
espei_ready / espei_certified / calphad_check fields of the literature release.

Unit part (always): liquid labels and captions, no default phase, compound-property rule, temperature gates, value
gates. Release part (when the orchid-data checkout is present): the PI's examples (Fe-Y Saenko 2017 table 5,
Fe-Hf, Al-Ca, Ni-Sc, Ru-Si) are no longer LIQUID; no non-ZPF dataset names LIQUID without a printed liquid; no LIQUID
value lies below the liquid bound; compound formation values are `solver` endmember datasets; every value of every
ESPEI dataset is published with the phase rule passed; espei_certified implies espei_ready and plausibility pass.
usage: python3 tools/tests/test_phase_rule.py"""
import glob, gzip, json, os, sys
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import phase_rule as PR  # noqa: E402
import build_espei_bundle as B  # noqa: E402
import orchid_cite  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


# ---- unit
for lab in ("liquid", "L", "Liq", "(Bi,U)liq", "melt", "liquid solution"):
    check("liquid label %r" % lab, PR.is_liquid_label(lab))
for lab in ("L12", "γ'-L12", "Fe17Y2", "UNKNOWN", "liquidus", "B2", "L1_0", "solid"):
    check("not a liquid label %r" % lab, not PR.is_liquid_label(lab))
check("model phase LIQUID is liquid", PR.is_liquid_label("x", "LIQUID"))
check("caption: Sb dropped into liquid Li", PR.caption_liquid(
    "Integral and partial enthalpies of mixing of Sb dropped into liquid Li at 979 K with respect to the pure liquid components"))
check("caption: Si-Ge liquid alloys", PR.caption_liquid("Activity of Si and Ge of the Si-Ge Liquid Alloys at 1723 K"))
check("caption: reference state only is not liquid", not PR.caption_liquid(
    "The calculated boron activities (T = 1850 K; and reference states: pure liquid components)."))
check("caption: solid Li3Sb is not liquid", not PR.caption_liquid(
    "Integral and partial enthalpies of mixing of Sb dropped into solid Li 3 Sb at 979 K"))
check("caption: intermetallic compounds is not liquid", not PR.caption_liquid(
    "Enthalpies of formation of the intermetallic compounds relative to liquid Al"))
check("activity 2.1 rejected", B.value_gate("ACR_TI", 2.1, {}, 1500) == "activity_out_of_range")
check("activity 0 rejected", B.value_gate("ACR_CU", 0.0, {}, 1500) == "activity_out_of_range")
check("activity 0.5 kept", B.value_gate("ACR_CU", 0.5, {}, 1500) is None)
check("negative SM rejected", B.value_gate("SM", -3.0, {"property": "SM"}, 298) == "absolute_entropy_negative")
check("per-gram Cp rejected", B.value_gate("CPM", 0.15, {"property": "CPM"}, 700) == "heat_capacity_unit_implausible")
check("increment rejected", B.value_gate("SM", 30.0, {"property": "SM", "property_printed": "S(T) − S(298)"}, 800)
      == "increment_not_absolute")
check("HM_FORM +196 kJ rejected", B.value_gate("HMR", 196e3, {"property": "HM_FORM"}, 298) == "formation_magnitude_implausible")
check("HM_FORM -40 kJ kept", B.value_gate("HMR", -40e3, {"property": "HM_FORM"}, 298) is None)

data = orchid_cite.DEFAULT_DATA
if not os.path.isdir(os.path.join(data, "espei")):
    print("orchid-data not found: release part skipped")
    sys.exit(1 if fails else 0)

R = PR.PhaseRule(data)
rec = {"observation_type": "equilibrium_property", "system": "FE-Y", "property": "HM_FORM", "temperature_K": 973.0,
       "phases": [{"label": "UNKNOWN"}], "locator": {"caption": "Enthalpy of formation"}, "paper": "doi:x"}
check("no printed phase -> fail (no default LIQUID)", R.classify(rec)["reason"] == "phase_not_printed")
rec2 = dict(rec, phases=[{"label": "liquid"}], property="HM_MIX", temperature_K=800.0)
check("Fe-Y liquid at 800 K impossible", R.classify(rec2)["reason"] == "liquid_impossible_at_T")
rec3 = dict(rec2, temperature_K=1900.0)
check("Fe-Y liquid at 1900 K possible", R.classify(rec3)["result"] == "pass")
rec4 = dict(rec, observation_type="invariant_reaction", temperature_K=18.0, system="B-TA", phases=[{"label": "L"}])
check("invariant at 18 K implausible", R.classify(rec4)["reason"] == "temperature_implausible")
rec5 = dict(rec, temperature_K=0.0, phases=[{"label": "Fe17Y2"}], paper="doi:10.1007/s11669-017-0574-3")
check("T = 0 implausible", R.classify(rec5)["reason"] == "temperature_implausible")

pub = {}
with gzip.open(os.path.join(data, "literature", "orchid-literature.jsonl.gz"), "rt", encoding="utf-8") as fh:
    for line in fh:
        r = json.loads(line)
        pub[r["id"]] = r
E = os.path.join(data, "espei")
DS = {}
for f in glob.glob(os.path.join(E, "*", "*", "*.json")):
    DS[os.path.relpath(f, E)] = json.load(open(f, encoding="utf-8"))

for sysn, stem in (("FE-Y", "FeY-17Sae.hmr-03"), ("FE-HF", "FeHf-08Idb.hmr-01"), ("AL-CA", "AlCa-05Ozt.hmr-01"),
                   ("NI-SC", "NiSc-13Kar.hmr-00"), ("RU-SI", "RuSi-01Liu.hmr-01")):
    hits = [rel for rel in DS if rel.startswith("experimental/%s/%s" % (sysn, stem))]
    check("%s %s: no LIQUID dataset" % (sysn, stem), all("LIQUID" not in DS[r]["phases"] for r in hits))
fy = [r for r in DS if r.startswith("experimental/FE-Y/FeY-17Sae.hmr-03")]
check("Fe-Y Saenko table 5: four datasets, one compound phase each (endmember `solver` where stoichiometric in "
      "FeY-17Sae.tdb)", len(fy) == 4 and all(len(DS[r]["phases"]) == 1 and DS[r]["phases"] != ["LIQUID"] for r in fy)
      and all(("solver" in DS[r]) == B.ec.Registry(data).stoichiometric("FeY-17Sae", DS[r]["phases"][0]) for r in fy))
check("Fe17Y2 endmember FE:Y", any(DS[r]["phases"] == ["FE17Y2"] and DS[r]["solver"]["sublattice_configurations"] == [["FE", "Y"]]
                                   for r in fy))

bad_liq, bad_T, not_pub, bad_rule = [], [], [], []
for rel, d in DS.items():
    if rel.split("/")[0] not in ("experimental", "model", "invariants"):
        continue
    ids = [o["id"] for o in d["orchid"]["observations"]]
    for oid in ids:
        r = pub.get(oid)
        if r is None:
            not_pub.append((rel, oid))
            continue
        c = R.classify(r)
        if c["result"] != "pass":
            bad_rule.append((rel, oid, c.get("reason")))
        if d["output"] != "ZPF" and "LIQUID" in d["phases"]:
            if not any(p["kind"] == "liquid" for p in c.get("phases") or []):
                bad_liq.append((rel, oid))
            tb, _ = R.liquid_bound(r["system"])
            if tb is not None and r.get("temperature_K") is not None and r["temperature_K"] < tb:
                bad_T.append((rel, oid))
check("every ESPEI value is published (%d missing)" % len(not_pub), not not_pub)
check("every ESPEI value passes the phase rule (%d fail)" % len(bad_rule), not bad_rule)
check("no LIQUID without a printed liquid (%d)" % len(bad_liq), not bad_liq)
check("no LIQUID below the liquid bound (%d)" % len(bad_T), not bad_T)
check("no dataset names LIQUID beside a compound by default",
      not [r for r, d in DS.items() if d["output"] != "ZPF" and "solver" not in d and "LIQUID" in d["phases"]
           and len(d["phases"]) > 1 and not d.get("orchid_phase_assignment")])

n_ready = n_cert = 0
okc = True
for r in pub.values():
    fl = r["flags"]
    n_ready += bool(fl.get("espei_ready"))
    n_cert += bool(fl.get("espei_certified"))
    pl = (fl.get("plausibility") or {}).get("result")
    if fl.get("espei_certified") and not (fl.get("espei_ready") and pl == "pass"):
        okc = False
    if fl.get("calphad_check") not in ("agrees", "disagrees", "not_available"):
        okc = False
check("espei_certified => espei_ready and plausibility pass; calphad_check vocabulary", okc)
idx = json.load(open(os.path.join(data, "literature", "index.json"), encoding="utf-8"))
check("index totals.espei agree with the records",
      idx["totals"]["espei"]["espei_ready"] == n_ready and idx["totals"]["espei"]["espei_certified"] == n_cert)
print("%d failed" % len(fails) if fails else "all passed")
sys.exit(1 if fails else 0)
