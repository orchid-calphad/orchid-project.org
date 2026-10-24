#!/usr/bin/env python3
"""Synthetic-fixture tests of tools/build_assessed_assemblages.py: the licensing guard refuses image/crop/page-text
content, point placement (centroid, pole-of-inaccessibility fallback, too thin), and tie-line side assignment.
No handbook data are read. usage: ~/orchid-py313/bin/python tools/tests/test_assessed_assemblages.py"""
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import build_assessed_assemblages as B  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


def refused(obj):
    try:
        B.guard(obj, "t")
    except SystemExit:
        return True
    return False


# ---- licensing guard
check("guard: clean record passes", not refused({"assemblage_label": "L+(Al)", "T_K": 900.0, "citations": [{"citation": "A. B., J. Phase Equilib. 1 (1990) 1"}]}))
check("guard: image key", refused({"page_image": "x"}))
check("guard: crop key", refused({"crop_box": [1, 2]}))
check("guard: page_text key", refused({"page_text": "abc"}))
check("guard: png path value", refused({"note": "see p0046/page.png"}))
check("guard: pdf path value", refused({"src": "handbook.pdf"}))
check("guard: base64 blob", refused({"x": "iVBORw0KGgo" + "A" * 200}))
check("guard: data URI", refused({"x": "data:image/png;base64,AAAA"}))
check("guard: long free text", refused({"label": "word " * 120}))
check("guard: build path", refused({"x": "/Users/me/projects/open-tdb/build/okamoto-annotation/x"}))
try:
    B.guard_text("a,b\n1,crop.png\n", "csv")
    check("guard_text refuses an image path in a cell", False)
except SystemExit:
    check("guard_text refuses an image path in a cell", True)

# ---- point placement
nT, nx = 201, 301
xs, Ts = np.linspace(0, 1, nx), np.linspace(0, 1000, nT)
m = np.zeros((nT, nx), bool)
m[50:150, 50:250] = True
p = B.field_point(m, xs, Ts)
check("rectangle -> centroid", p["method"] == "centroid" and abs(xs[p["col"]] - 0.5) < 0.01 and abs(Ts[p["row"]] - 500) < 6)
yy, xx = np.mgrid[0:nT, 0:nx]
ring = ((yy - 100) ** 2 + (xx - 150) ** 2 < 90 ** 2) & ((yy - 100) ** 2 + (xx - 150) ** 2 > 60 ** 2) & (xx < 150)
p = B.field_point(ring, xs, Ts)
check("crescent -> pole of inaccessibility, centroid outside", p["method"] == "pole_of_inaccessibility" and ring[p["row"], p["col"]]
      and p["centroid"]["inside"] is False)
thin = np.zeros((nT, nx), bool)
thin[100, 20:280] = True
check("one-cell strip -> no point", B.field_point(thin, xs, Ts)["method"] is None)

# ---- tie line: field 2 (a+b) between single-phase fields 1 (a) and 3 (b)
row = np.array([1] * 100 + [0] + [2] * 100 + [0] + [3] * 100)
xr = np.linspace(0, 1, len(row))
names = {1: ["a"], 2: ["a", "b"], 3: ["b"]}
left, right = B.run_extent(row, 150, 2, -1), B.run_extent(row, 150, 2, +1)
tl, why = B.tie_line(["a", "b"], left, right, xr, names, {}, 0.001, float(xr[150]), 900.0, False, False)
check("tie line: endpoints a (left) and b (right)", tl is not None and [e["phase"] for e in tl["endpoints"]] == ["a", "b"]
      and tl["endpoints"][0]["x"] < xr[150] < tl["endpoints"][1]["x"] and "low credence" in tl["caveat"])
names_bad = {1: ["c"], 2: ["a", "b"], 3: ["b"]}
tl, why = B.tie_line(["a", "b"], left, right, xr, names_bad, {}, 0.001, float(xr[150]), 900.0, False, False)
check("tie line: unrelated neighbour -> omitted with reason", tl is None and why == "side_assignment_ambiguous")
names_open = {1: None, 2: ["a", "b"], 3: ["b"]}
tl, why = B.tie_line(["a", "b"], left, right, xr, names_open, {}, 0.001, float(xr[150]), 900.0, False, False)
check("tie line: unnamed neighbour -> omitted", tl is None and why == "neighbour_field_unnamed")
# line compound on the left: neighbour field 1 is (c + a) across the compound a, a printed at 30 at.%
row2 = np.array([1] * 100 + [0] + [2] * 100 + [0] + [3] * 100)
xr2 = np.linspace(0.0, 1.0, len(row2))
xb = (xr2[100]) * 100
names_lc = {1: ["c", "a"], 2: ["a", "b"], 3: ["b"]}
tl, why = B.tie_line(["a", "b"], B.run_extent(row2, 150, 2, -1), B.run_extent(row2, 150, 2, +1), xr2, names_lc,
                     {"a": {"composition_at_pct": [xb, xb]}}, 0.001, float(xr2[150]), 900.0, False, False)
check("tie line: line-compound endpoint takes the printed composition", tl is not None and tl["endpoints"][0].get("line_compound")
      and abs(tl["endpoints"][0]["x"] - xb / 100) < 1e-9)

# ---- uncertainty budget (PI decision 2026-10-24): components separate; synthetic record
import copy  # noqa: E402
rec = {"id": "aspa-test", "coordinate_uncertainty": {"T_K": 3.5, "x": 0.004, "calibration_max_residual_T_K": 1.0},
       "positional_margin": {"normalized": 0.01, "isotherm_mole_fraction": 0.05, "isopleth_K": 40.0}}
before = copy.deepcopy(rec)
b = B.assessed_budget(rec)
qT, qX = b["quantities"]["temperature_K"], b["quantities"]["composition"]
check("budget: same structure as the release", set(b) == {"calibration", "sigma_rule", "quantities", "boundary_basis"}
      and all(k in qT for k in ("sigma_printed", "sigma_reading", "sigma_method", "sigma_composition", "sigma_scatter",
                                "sigma_conversion", "sigma_effective")))
check("budget: digitisation = coordinate uncertainty", qT["sigma_reading"]["value"] == 3.5 and qX["sigma_reading"]["value"] == 0.004)
check("budget: handbook prior 10 K / 0.02", qT["sigma_method"]["value"] == 10.0 and qX["sigma_method"]["value"] == 0.02
      and qT["sigma_method"]["basis"] == "prior")
check("budget: margin reported, not combined", qT["sigma_boundary_margin"]["value"] == 40.0 and qT["sigma_effective"] == 10.0)
rec2 = dict(rec, coordinate_uncertainty={"T_K": 25.0, "x": 0.03})
b2 = B.assessed_budget(rec2)
check("budget: effective = max(reading, prior)", b2["quantities"]["temperature_K"]["sigma_effective"] == 25.0
      and b2["quantities"]["composition"]["sigma_effective"] == 0.03)
check("budget: coordinate fields unchanged", rec == before)
# dashed (estimated) boundaries: reading enlarged to 2x (dashed) / 1.5x (mixed) the handbook prior; solid unchanged
bd = B.assessed_budget(dict(rec, boundary_basis="dashed_estimated"))["quantities"]
bm = B.assessed_budget(dict(rec, boundary_basis="mixed"))["quantities"]
bs = B.assessed_budget(dict(rec, boundary_basis="printed_solid"))["quantities"]
check("budget: dashed reading 2x prior (20 K, 0.04), digitisation kept",
      bd["temperature_K"]["sigma_reading"]["value"] == 20.0 and bd["composition"]["sigma_reading"]["value"] == 0.04
      and bd["temperature_K"]["sigma_reading"]["digitisation"] == 3.5 and bd["temperature_K"]["sigma_effective"] == 20.0)
check("budget: mixed reading 1.5x prior", bm["temperature_K"]["sigma_reading"]["value"] == 15.0
      and abs(bm["composition"]["sigma_reading"]["value"] - 0.03) < 1e-12)
check("budget: solid keeps digitisation", bs["temperature_K"]["sigma_reading"]["value"] == 3.5)
# dashed-line re-derivation (tools/dashed_boundaries.py): a run of short collinear strokes is one dashed chain
import dashed_boundaries as DB  # noqa: E402
cal = {"calibration": {"x": {"a": 0.5, "b": 0}, "y": {"a": 5.0, "b": 0}},
       "boundaries": [{"id": "d%d" % k, "points": [[0.10 + 0.02 * k, 500.0], [0.10 + 0.02 * k + 0.012, 500.0]]} for k in range(5)]
       + [{"id": "solid", "points": [[0.5, 100.0], [0.5, 900.0]]}]}
L = DB.DiagramLines(cal, {}, (0.0, 1.0), (0.0, 1000.0))
check("dashed: 5 dashes give 8 joins (both directions, as the tracer counts)", L.n_dash_joins == 8)
check("dashed: crossing on the dashed line", L.basis_at([(0.15, 500.0)])[0] == "dashed_estimated")
check("dashed: crossing on the solid line", L.basis_at([(0.5, 400.0)])[0] == "printed_solid")
check("dashed: both", L.basis_at([(0.15, 500.0), (0.5, 400.0)])[0] == "mixed")
check("dashed: frame only", L.basis_at([(0.9, 950.0)])[0] is None)
# release part: every published record has a budget with sigma_effective >= reading and >= the handbook prior
import gzip, json  # noqa: E402,E401
rp = os.path.join(os.path.expanduser("~/projects/orchid-data"), "assessed-phase-assemblages", "records.jsonl.gz")
if os.path.exists(rp):
    recs = [json.loads(l) for l in gzip.open(rp, "rt", encoding="utf-8")]
    bad = []
    for r in recs:
        q = ((r.get("uncertainty_budget") or {}).get("quantities") or {})
        good = True
        for qn, key, pk in (("temperature_K", "T_K", "T"), ("composition", "x", "x")):
            e = q.get(qn)
            if not e or e["sigma_effective"] is None or e["sigma_effective"] < B.HANDBOOK_PRIOR[pk] - 1e-12 or \
                    (e["sigma_reading"] and e["sigma_effective"] < e["sigma_reading"]["value"] - 1e-12) or \
                    ((e["sigma_reading"] or {}).get("digitisation", (e["sigma_reading"] or {}).get("value"))
                     != r["coordinate_uncertainty"][key]):
                good = False
        if not good:
            bad.append(r["id"])
    check("release: every record has a budget, effective >= reading and >= prior (%d records)" % len(recs),
          bool(recs) and not bad)

# ---- structure fallback (James Saal review 2026-10-24): cF8 alone never names DIAMOND_A4
check("STD_PEARSON has no cF8 (diamond / rock salt / sphalerite share it)", "cF8" not in B.STD_PEARSON)

print("%d failure(s)" % len(fails))
sys.exit(1 if fails else 0)
