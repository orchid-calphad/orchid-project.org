#!/usr/bin/env python3
"""Tests of the supplementary precision sample (build_supplement.py, score_supplement.py). Temp files only.

Fixture: a small synthetic release (phase registry, ESPEI datasets of every folder kind, literature records, assessed
assemblages with run geometry, an earlier blind key). Checks: allocation (minimums, caps, sums); ESPEI value
flattening (multi-composition, multi-temperature, solver, ZPF); dashed-boundary detection; the draw is deterministic,
respects the strata, takes one S-ASM field per diagram, excludes the earlier blind and expert diagrams, and keeps
S-ASM blind (no record id, machine name, basis or stratum in sample.json or the pages); draw refuses an existing
sample; packets partition the items with ~10 % overlap; `repartition` re-cuts the same sample into 5 balanced packets
(sample.json, key.json, index.html unchanged; old packets kept as *_3pk) and the scorer follows it; `check` passes on an intact kit, fails on a changed or extra
file and only reports release drift; the scorer refuses a foreign result file, prefers the owner's answer, computes
Wilson intervals and agreement, and scores S-ASM names with the blind sample's equivalences (identical to open-tdb's
blind_score.py when that file is present); the kit's .gitignore keeps restricted files out of git.
usage: python3 tools/test_supplement.py"""
import contextlib
import gzip
import io
import json
import math
import random
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import build_supplement as B  # noqa: E402
import score_supplement as S  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


def quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        try:
            return fn(*a, **k)
        except SystemExit as e:
            return ("exit", str(e))


TMP = Path(tempfile.mkdtemp(prefix="orchid-supp-test-"))


# ---------------------------------------------------------------- fixture

def jdump(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(obj, ensure_ascii=False), encoding="utf-8")


def make_release(root):
    data = root / "orchid-data"
    rnd = random.Random(1)
    systems = ["AL-NI", "CU-ZR", "B-FE", "AG-SI"]
    # phase registry
    phases, full = [], {}
    tiers = [("assigned", {"tier": "A", "basis": "standard_name"}), ("assigned", {"tier": "B", "basis": "dft_unique_match"}),
             ("assigned", {"tier": "C", "basis": "assessment_declared"}), ("assigned", {"basis": "liquid"}),
             ("assigned", {"basis": "gas"}), ("needs_review", None), ("unresolved_no_candidates", None)]
    for s in systems:
        for j in range(30):
            status, st = tiers[j % len(tiers)]
            pid = "phase:%s:p%02d" % (s.lower(), j)
            el = s.split("-")
            st2 = dict(st, prototype="Cu", pearson="cF4", strukturbericht="A1", space_group={"number": 225, "symbol": "Fm-3m"}) if st else None
            phases.append({"id": pid, "system": s, "elements": el, "preferred_name": "P%02d" % j, "status": status,
                           "structure": st2, "detail": "systems/%s.json" % s})
            full.setdefault(s, []).append({
                "id": pid, "system": s, "elements": el, "preferred_name": "P%02d" % j, "status": status, "structure": st2,
                "names": [{"name": "tdb:X.tdb:P%02d" % j, "namespace": "tdb", "label": "P%02d" % j}],
                "models": [{"tdb": "X%d.tdb" % m, "path": "databases/sgte-binary/X%d.tdb" % m, "phase_name": "P%02d" % j,
                            "tdb_definition": {"verbatim": {"phase": {"line_start": 3, "line_end": 3, "text": "PHASE P%02d %% 1 1 !" % j}},
                                               "provenance": {"website_url": "https://orchid-project.org/data/x.tdb"}}}
                           for m in range(2)],
                "references": [{"doi": "10.1000/r%d" % j, "citation": "R%d" % j, "roles": ["reports values for this phase"]}]})
    jdump(data / "phases" / "registry.json", {"phases": phases})
    for s, ph in full.items():
        jdump(data / "phases" / "systems" / ("%s.json" % s), {"system": s, "phases": ph})
    # literature + ESPEI
    lit, rows, papers = [], [], {}
    n = 0

    def rec(sysname):
        nonlocal n
        n += 1
        rid = "ORC-%06d" % n
        lit.append({"id": rid, "system": sysname, "paper": "doi:10.1000/p%d" % (n % 5), "value": n, "value_printed": n,
                    "unit": "J/mol", "evidence_class": "experimental", "record_kind": "general_table",
                    "locator": {"artifact": "table", "page": 1 + n % 4, "label": "1", "row": n},
                    "flags": {"extraction_accepted": True, "espei_ready": True}})
        return rid

    kinds = [("experimental", "HM_FORM", "experimental"), ("experimental", "ZPF", "experimental"),
             ("model", "HMR", "model"), ("model", "CPM", "model"), ("experimental", "ACR_CU", "experimental"),
             ("invariants", "ZPF", "experimental"), ("dft/literature", "HM_FORM", "dft"), ("experimental", "SM", "experimental")]
    for s in systems:
        for folder, out, ev in kinds:
            for f in range(3):
                path = "%s/%s/%s.%s-%02d.json" % (folder, s, s.replace("-", ""), out.lower(), f)
                k = 2 + f
                ids = [rec(s) for _ in range(k)]
                d = {"components": s.split("-") + ["VA"], "phases": ["LIQUID"], "output": out, "tags": ["orchid"],
                     "orchid": {"evidence_class": ev},
                     "orchid_uncertainty": {"unit": "J/mol-atom", "values": [{"sigma": None, "budget": {"sigma_effective": 1.0 * i}} for i in range(k)]}}
                if out == "ZPF":
                    d["conditions"] = {"P": 101325, "T": [1000.0 + i for i in range(k)]}
                    d["values"] = [[["LIQUID", [s.split("-")[1]], [0.1 * (i + 1)]]] for i in range(k)]
                elif folder == "dft/literature":
                    d["conditions"] = {"P": 101325, "T": 298.15}
                    d["values"] = [[[-1000.0 * (i + 1) for i in range(k)]]]
                    d["solver"] = {"mode": "manual", "sublattice_site_ratios": [1, 1], "sublattice_configurations": [[s.split("-")[0], s.split("-")[1]]] * k}
                    d["comment"] = "0 K DFT. eV/atom x 96485.33212 = J/mol-atom."
                else:
                    d["conditions"] = {"P": 101325, "T": 1273.0, "X_" + s.split("-")[1]: [0.1 * (i + 1) for i in range(k)]}
                    d["values"] = [[[100.0 * (i + 1) for i in range(k)]]]
                if folder == "dft/literature":
                    d["orchid"]["mapping"] = [{"id": i} for i in ids]
                else:
                    d["orchid"]["observations"] = [{"id": i, "paper": "doi:10.1000/p1"} for i in ids]
                jdump(data / "espei" / path, d)
                rows.append({"path": path, "folder": folder, "output": out, "evidence_class": ev})
    rows.append({"path": "assessed/X/y.json", "folder": "assessed", "output": "ZPF", "evidence_class": "assessed_phase_diagram"})
    jdump(data / "espei" / "index.json", {"datasets": rows})
    for j in range(5):
        papers["doi:10.1000/p%d" % j] = {"doi": "10.1000/p%d" % j, "citation": "P%d" % j, "doi_url": "https://doi.org/10.1000/p%d" % j}
    (data / "literature").mkdir(parents=True, exist_ok=True)
    (data / "literature" / "orchid-literature.jsonl.gz").write_bytes(
        gzip.compress("".join(json.dumps(r) + "\n" for r in lit).encode(), mtime=0))
    jdump(data / "literature" / "papers.json", papers)
    # assessed assemblages + run geometry
    import numpy as np
    annot = root / "annot"
    run = "run-test"
    recs = []
    bases = list(B.BASIS_ORDER[:5])
    diagrams = ["%s-p%d" % (s, 10 + i) for i, s in enumerate(systems * 6)]
    for di, d in enumerate(diagrams):
        sd = annot / run / "systems" / d
        sd.mkdir(parents=True, exist_ok=True)
        xs = np.linspace(0, 1, 100)
        Ts = np.linspace(0, 1000, 100)
        lab = np.zeros((100, 100), dtype=np.int32)
        lab[:, :50] = 1
        lab[:, 50:] = 2
        np.savez(sd / "grid.npz", lab=lab, x=xs, T=Ts)
        joins = [{"points": [[0.5, 500.0], [0.5, 510.0]], "dashed": True}] if di % 3 == 0 else [{"points": [[0.5, 500.0]]}]
        jdump(sd / "fields.json", {"endpoint_joins": joins, "phase_table": {"(A)": {}, "AB": {}}, "diagram_labels": ["L+AB"]})
        jdump(sd / "calibrated.json", {"page": 10 + di, "component": d.split("-")[1]})
        for f in (1, 2):
            recs.append({"id": "aspa-%s-%d" % (d, f), "system": d.rsplit("-p", 1)[0], "diagram": d, "field": f,
                         "assemblage_label": "L+AB" if f == 1 else "(A)", "assemblage": [{"canonical_name": "L"}],
                         "basis": bases[(di + f) % 5], "confidence": "high",
                         "review_state": "expert decision (tier 5)" if di == 1 else "not reviewed by the expert",
                         "T_K": 900.0, "T_C": 626.85, "x": 0.25 if f == 1 else 0.75, "x_component": d.split("-")[1],
                         "method": {"annotation_run": run}})
    (data / "assessed-phase-assemblages").mkdir(parents=True, exist_ok=True)
    (data / "assessed-phase-assemblages" / "records.jsonl.gz").write_bytes(
        gzip.compress("".join(json.dumps(r) + "\n" for r in recs).encode(), mtime=0))
    blind = root / "blind-key.json"
    jdump(blind, {"items": {"B001": {"diagram": diagrams[2]}, "B002": {"diagram": diagrams[3]}}})
    return data, annot, blind, diagrams


data, annot, blind, diagrams = make_release(TMP)

# ---------------------------------------------------------------- unit checks

a = B.allocate({"a": 1000, "b": 10, "c": 3, "d": 0}, 50, 6)
check("allocate: sums to n, minimum, cap at population, empty stratum 0",
      sum(a.values()) == 50 and a["b"] == 6 and a["c"] == 3 and "d" not in a and a["a"] == 41)
a2 = B.allocate2({("p", "x"): 100, ("p", "y"): 5, ("q", "x"): 50, ("q", "y"): 0}, 20, 6)
check("allocate2: secondary levels get at least one when there is room", a2[("p", "y")] >= 1 and sum(a2.values()) == 20
      and a2[("q", "y")] == 0)

d_multi = {"output": "HMR", "conditions": {"P": 101325, "T": [300.0, 400.0], "X_B": [0.1, 0.2, 0.3]},
           "values": [[[1, 2, 3], [4, 5, 6]]]}
v, c, _ = B.espei_value(d_multi, 4)
check("espei_value: P x T x X flattening", v == 5 and c["T"] == 400.0 and c["X_B"] == 0.2)
d_zpf = {"output": "ZPF", "conditions": {"P": 101325, "T": [900.0, 950.0]}, "values": [[["L", ["B"], [0.1]]], [["L", ["B"], [0.2]]]]}
v, c, _ = B.espei_value(d_zpf, 1)
check("espei_value: ZPF row and its temperature", v == [["L", ["B"], [0.2]]] and c["T"] == 950.0)
d_sol = {"output": "HM_FORM", "conditions": {"P": 101325, "T": 298.15}, "values": [[[-1, -2]]],
         "solver": {"mode": "manual", "sublattice_site_ratios": [2, 1], "sublattice_configurations": [["A", "B"], ["A", "C"]]}}
v, c, sol = B.espei_value(d_sol, 1)
check("espei_value: solver configuration of the value", v == -2 and sol["configuration"] == ["A", "C"])

import numpy as np  # noqa: E402
lab = np.zeros((100, 100), dtype=np.int32)
lab[:, :50] = 1
lab[:, 50:] = 2
lab[:, 90:] = 3
grid = (lab, np.linspace(0, 1, 100), np.linspace(0, 1000, 100))
df = B.dashed_fields({"endpoint_joins": [{"points": [[0.5, 500.0]], "dashed": True}, {"points": [[0.95, 500.0]]}]}, grid)
check("dashed_fields: fields beside a dashed join only", df == {1, 2})

# ---------------------------------------------------------------- draw

kit = TMP / "kit"


def draw(out, *extra):
    return quiet(B.main, ["draw", "--out", str(out), "--data", str(data), "--annot", str(annot), "--blind-key", str(blind),
                          "--no-images", "--k", "3"] + list(extra))


check("draw: succeeds on the fixture", draw(kit) == 0)
check("draw: refuses a directory that holds a sample", draw(kit) == 1)
kit2 = TMP / "kit2"
draw(kit2)
check("draw: deterministic (sample.json, key.json, index.html identical)",
      all((kit / f).read_bytes() == (kit2 / f).read_bytes() for f in ("sample.json", "key.json", "index.html")))
sample = json.loads((kit / "sample.json").read_text(encoding="utf-8"))
key = json.loads((kit / "key.json").read_text(encoding="utf-8"))
man = json.loads((kit / "manifest.json").read_text(encoding="utf-8"))
items = key["items"]
for sub in ("S-REG", "S-ESPEI", "S-ASM"):
    st = key["strata"][sub]
    check("draw %s: sampled per stratum adds up to the items drawn" % sub,
          sum(s["sampled"] for s in st) == sum(1 for it in items if it["kind"] == sub) > 0)
    check("draw %s: every sampled item's stratum is a stratum cell" % sub,
          {it["stratum"] for it in items if it["kind"] == sub} <= {s["stratum"] for s in st})
asm = [it for it in items if it["kind"] == "S-ASM"]
diag = [it["source_id"].split("-", 1)[1].rsplit("-", 1)[0] for it in asm]
check("draw S-ASM: one field per diagram", len(diag) == len(set(diag)))
check("draw S-ASM: earlier blind-sample and expert diagrams excluded",
      not ({diagrams[1], diagrams[2], diagrams[3]} & set(diag)))
check("draw S-ASM: both boundary strata present", {it["stratum"].split("|")[1] for it in asm} == {"dashed", "solid"})
pub = (kit / "sample.json").read_text(encoding="utf-8")
page = (kit / "index.html").read_text(encoding="utf-8")
leak = [it["item"] for it in asm if it["source_id"] in pub or it["source_id"] in page]
asm_pub = [x for x in sample["items"] if x["kind"] == "S-ASM"]
check("blind: no S-ASM record id in sample.json or index.html", not leak)
check("blind: S-ASM items carry no stratum, basis or machine name in the page data",
      all("stratum" not in x and "secret" not in x and "basis" not in x["view"] and "machine_name" not in json.dumps(x) for x in asm_pub))
check("draw: S-REG / S-ESPEI items carry their stratum", all("stratum" in x for x in sample["items"] if x["kind"] != "S-ASM"))
own = [i for p in key["packets"] for i in p["own"]]
check("packets: own items partition the sample", sorted(own) == sorted(it["item"] for it in items) and len(own) == len(set(own)))
ov = sum(len(p["overlap"]) for p in key["packets"])
check("packets: overlap is about 10 %% of the items (%d of %d)" % (ov, len(items)), 0.05 <= ov / len(items) <= 0.2)
check("packets: overlap items come from the next packet",
      all(set(p["overlap"]) <= set(key["packets"][p["packet"] % 3]["own"]) for p in key["packets"]))
check("manifest: hashes sample.json, key.json, index.html and every packet",
      {"sample.json", "key.json", "index.html", "packets/packets.json"} <= set(man["files_sha256"])
      and sum(1 for f in man["files_sha256"] if f.startswith("packets/ORCHID_")) == 3)
check("manifest: population snapshot records the release files", "phases/registry.json" in man["population"]["files_sha256"])

# the release's own boundary flag, when present, defines the S-ASM secondary strata
rp = data / "assessed-phase-assemblages" / "records.jsonl.gz"
orig = rp.read_bytes()
rs = [json.loads(l) for l in gzip.decompress(orig).decode().splitlines() if l.strip()]
for j, r in enumerate(rs):
    r["boundary_basis"] = [None, "printed_solid", "mixed", "dashed_estimated"][j % 4]
rp.write_bytes(gzip.compress("".join(json.dumps(r) + "\n" for r in rs).encode(), mtime=0))
cand, info = B.asm_population(data, annot, blind)
check("S-ASM strata: release boundary_basis used when present (null -> frame_only)",
      {c[1] for c in cand} == {"frame_only", "printed_solid", "mixed", "dashed_estimated"}
      and info["boundary_flag"] == "release field boundary_basis")
rp.write_bytes(orig)

# ---------------------------------------------------------------- check

check("check: intact kit passes", quiet(B.main, ["check", "--out", str(kit), "--data", str(data)]) == 0)
drift = TMP / "drift.json"
esp = next(it for it in items if it["kind"] == "S-ESPEI")
path, i = esp["source_id"].rsplit("#", 1)
dd = json.loads((data / "espei" / path).read_text(encoding="utf-8"))
if dd["output"] == "ZPF":
    dd["values"][int(i)][0][2] = [0.999]
else:
    dd["values"] = [[[x + 1 for x in dd["values"][0][0]]]] if len(dd["values"][0]) == 1 else dd["values"]
(data / "espei" / path).write_text(json.dumps(dd), encoding="utf-8")
rc = quiet(B.main, ["check", "--out", str(kit), "--data", str(data), "--drift-json", str(drift)])
check("check: a release change is reported (changed), not a failure",
      rc == 0 and json.loads(drift.read_text())[esp["item"]] == "changed")
(kit / "packets" / "extra.html").write_text("x")
check("check: a file not in the freeze record fails", quiet(B.main, ["check", "--out", str(kit), "--data", str(data)]) == 1)
(kit / "packets" / "extra.html").unlink()
bak = (kit / "sample.json").read_bytes()
(kit / "sample.json").write_bytes(bak.replace(b'"seed": 20261024', b'"seed": 1'))
check("check: a changed sample.json fails", quiet(B.main, ["check", "--out", str(kit), "--data", str(data)]) == 1)
r = quiet(S.main, ["--kit", str(kit), str(drift)])
check("score: refuses a sample.json that is not the frozen one", isinstance(r, tuple) and r[0] == "exit")
(kit / "sample.json").write_bytes(bak)

# ---------------------------------------------------------------- scorer

p, lo, hi = S.wilson(8, 10)
check("wilson: 8/10 -> [0.490, 0.943]", abs(lo - 0.4902) < 1e-3 and abs(hi - 0.9433) < 1e-3)
check("asm: order, spelling and parentheses normalised", S.asm_verdict("(Cu)+L", "L+Cu") == "agree")
check("asm: occurrence numbering counts as a notation equivalence",
      S.asm_verdict("(Nb)#2+(Nb)", "(Nb)1+(Nb)2") in ("agree", "agree_notation"))
check("asm: different phases disagree", S.asm_verdict("L+AB", "L") == "disagree")
check("asm: cannot tell is not an error", S.asm_verdict("cannot_tell", "L") == "cannot_tell")
blind_py = B.ANNOT.parent / "okamoto-pilot" / "tools" / "blind_score.py"
if blind_py.exists():
    sys.path.insert(0, str(blind_py.parent))
    try:
        import blind_score as BS  # noqa: E402
        names = ["L+(Cu)", "(αFe)+AsFe2", "(Nb)#2+(Nb)", "α/βCaCu+L", "InPr1+x+L", "O2+L", "ϕ+L1", "(βNd,βPr)", "Gas+(Ag)"]
        check("asm: normalisation identical to open-tdb blind_score.py",
              all(S.phases(n) == BS.phases(n) and S.phases_equiv(n) == BS.phases_equiv(n) for n in names))
    except Exception as e:  # pragma: no cover
        check("asm: open-tdb blind_score.py importable (%s)" % e, False)

sha = man["sample_sha256"]
byitem = {it["item"]: it for it in items}


def answers_for(packet, value_answer, asm_name):
    pk = key["packets"][packet - 1]
    out = []
    for i in pk["own"] + pk["overlap"]:
        it = byitem[i]
        a = {"item": i, "n": it["n"], "kind": it["kind"], "note": ""}
        if it["kind"] == "S-ASM":
            a["decision"] = asm_name(it)
            a["confidence"] = "high"
        else:
            for q in S.QUESTIONS[it["kind"]]:
                a[q] = "yes"
            a[S.KEY_QUESTION[it["kind"]]] = value_answer
        out.append(a)
    return {"schema": B.RESULT_SCHEMA, "sample_sha256": sha, "packet": packet, "reviewer": "R%d" % packet, "answers": out}


res = []
for pk in (1, 2, 3):
    good = answers_for(pk, "yes" if pk != 2 else "no", lambda it: it["secret"]["machine_name"] if pk != 2 else "cannot_tell")
    f = TMP / ("res%d.json" % pk)
    f.write_text(json.dumps(good), encoding="utf-8")
    res.append(f)
out_json = TMP / "scores.json"
rc = quiet(S.main, ["--kit", str(kit), "--json", str(out_json)] + [str(x) for x in res])
sc = json.loads(out_json.read_text(encoding="utf-8")) if out_json.exists() else {}
n_reg = sum(1 for it in items if it["kind"] == "S-REG")
n_reg_owned_by_2 = sum(1 for it in items if it["kind"] == "S-REG" and it["owner_packet"] == 2)
check("score: runs on three packets", rc == 0 and bool(sc))
check("score: owner's answer counts (packet 2 answered 'no' only for its own items)",
      sc["S-REG"]["names"]["overall"]["judged"] == n_reg and
      sc["S-REG"]["names"]["overall"]["correct"] == n_reg - n_reg_owned_by_2)
w = sc["S-REG"]["names"]["weighted"]
check("score: population-weighted estimate with an interval", w and 0 <= w["ci95"][0] <= w["precision"] <= w["ci95"][1] <= 1)
asm_owned_not2 = sum(1 for it in asm if it["owner_packet"] != 2)
check("score: S-ASM agrees with the machine for reviewers naming it; cannot-tell counted apart",
      sc["S-ASM"]["with_notation_equivalences"]["overall"]["correct"] == asm_owned_not2 and
      sc["S-ASM"]["verdicts"].get("cannot_tell", 0) == len(asm) - asm_owned_not2)
ag = sc["agreement"]
check("score: agreement measured on the overlap items",
      ag["n_items_answered_by_two_or_more"] == ov and ag["S-REG"]["names"]["n"] > 0)
bad = json.loads(res[0].read_text())
bad["sample_sha256"] = "0" * 64
fbad = TMP / "bad.json"
fbad.write_text(json.dumps(bad))
r = quiet(S.main, ["--kit", str(kit), str(fbad)])
check("score: refuses a result file for another sample", isinstance(r, tuple) and r[0] == "exit")
k1 = S.kappa([("yes", "yes"), ("no", "no"), ("yes", "no"), ("yes", "yes")])
check("kappa: 3/4 agreement, kappa 0.5", abs(k1[0] - 0.75) < 1e-9 and abs(k1[1] - 0.5) < 1e-9)

# ---------------------------------------------------------------- repartition (same sample, new packets)

kit5 = TMP / "kit5"
shutil.copytree(kit, kit5)
before = {f: (kit5 / f).read_bytes() for f in ("sample.json", "key.json", "index.html")}
old_packets = {f.name: f.read_bytes() for f in (kit5 / "packets").iterdir()}
rc = quiet(B.main, ["repartition", "--out", str(kit5), "--k", "5", "--overlap", "4"])
check("repartition: runs", rc == 0)
check("repartition: sample.json, key.json and index.html unchanged",
      all((kit5 / f).read_bytes() == b for f, b in before.items()))
check("repartition: old packets kept with the _3pk suffix, bytes unchanged",
      all((kit5 / "packets" / ("%s_3pk%s" % (Path(n).stem, Path(n).suffix))).read_bytes() == b for n, b in old_packets.items()))
part = json.loads((kit5 / "partition.json").read_text())
man5 = json.loads((kit5 / "manifest.json").read_text())
own5 = [i for p in part["packets"] for i in p["own"]]
check("repartition: 5 packets whose own items partition the sample",
      len(part["packets"]) == 5 and sorted(own5) == sorted(it["item"] for it in items) and len(own5) == len(set(own5)))
cnt = [{k_: sum(1 for i in p["own"] if byitem[i]["kind"] == k_) for k_ in ("S-REG", "S-ESPEI", "S-ASM")} for p in part["packets"]]
check("repartition: sub-samples balanced across packets (max - min <= 1)",
      all(max(c[k_] for c in cnt) - min(c[k_] for c in cnt) <= 1 for k_ in cnt[0]))
check("repartition: 4 overlap items per packet, owned by the next packet in the ring",
      all(len(p["overlap"]) == 4 and set(p["overlap"]) <= set(part["packets"][p["packet"] % 5]["own"]) for p in part["packets"]))
check("repartition: overlap spans the three sub-samples",
      all({byitem[i]["kind"] for i in p["overlap"]} == {"S-REG", "S-ESPEI", "S-ASM"} for p in part["packets"]))
check("repartition: manifest records the partition and the new packet hashes; check passes",
      man5["partition"]["id"] == "k5" and "packets/ORCHID_supplement_packet_5.html" in man5["files_sha256"]
      and quiet(B.main, ["check", "--out", str(kit5), "--data", str(data)]) == 0)
check("repartition: refuses to repeat the same partition", quiet(B.main, ["repartition", "--out", str(kit5), "--k", "5"]) == 1)
pg = (kit5 / "packets" / "ORCHID_supplement_packet_2.html").read_text()
check("repartition: packet pages carry the partition id", '"partition": "k5"' in pg or '"partition":"k5"' in pg)
res5 = []
for p in part["packets"]:
    ans = []
    for i in p["own"] + p["overlap"]:
        it = byitem[i]
        a_ = {"item": i, "n": it["n"], "kind": it["kind"], "note": ""}
        if it["kind"] == "S-ASM":
            a_["decision"] = it["secret"]["machine_name"]
        else:
            for q in S.QUESTIONS[it["kind"]]:
                a_[q] = "yes" if p["packet"] != 3 or i in p["overlap"] else "no"
        ans.append(a_)
    f = TMP / ("r5-%d.json" % p["packet"])
    f.write_text(json.dumps({"sample_sha256": sha, "packet": p["packet"], "partition": "k5", "reviewer": "x", "answers": ans}))
    res5.append(str(f))
o5 = TMP / "s5.json"
rc = quiet(S.main, ["--kit", str(kit5), "--json", str(o5)] + res5)
s5 = json.loads(o5.read_text()) if o5.exists() else {}
n_own3 = sum(1 for i in part["packets"][2]["own"] if byitem[i]["kind"] == "S-REG")
check("score (k5): owners from the partition", rc == 0 and s5["S-REG"]["names"]["overall"]["correct"] == n_reg - n_own3)
check("score (k5): one agreement entry per overlap pair from the partition",
      len(s5["agreement"]["pairs"]) == 5 and all(pr["answered_by_both"] == 4 for pr in s5["agreement"]["pairs"]))
r = quiet(S.main, ["--kit", str(kit5)] + [str(res[0])])
check("score (k5): refuses results from superseded packets", isinstance(r, tuple) and r[0] == "exit")

# ---------------------------------------------------------------- git hygiene of the real kit

real = B.KIT
try:
    ig = subprocess.run(["git", "-C", str(real), "check-ignore", "index.html", "packets/x.html", "img/a.jpg", "key.json",
                         "answers/r.json"], capture_output=True, text=True).stdout.split()
    check(".gitignore: index.html, packets/, img/, key.json, answers/ are ignored", len(ig) == 5)
    ig2 = subprocess.run(["git", "-C", str(real), "check-ignore", "sample.json", "manifest.json", "INSTRUCTIONS.md"],
                         capture_output=True, text=True).stdout.split()
    check(".gitignore: sample.json, manifest.json, INSTRUCTIONS.md are not ignored", not ig2)
except OSError:
    pass
tpl = (HERE / "page_template.html").read_text(encoding="utf-8")
check("page: no external scripts or stylesheets", "<script src" not in tpl and "<link rel=\"stylesheet\"" not in tpl)

shutil.rmtree(TMP, ignore_errors=True)
print("\n%d failure(s)" % len(fails) if fails else "\nall passed")
sys.exit(1 if fails else 0)
