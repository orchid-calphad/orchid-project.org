#!/usr/bin/env python3
"""Dense field grid of tools/build_assessed_assemblages.py: synthetic grid rules, and the release's Cu-Zr dense
points reproduce the 177 points of the PhaseForge+ Cu-Zr final fit (fixtures/cuzr_dense177.json, same T bit for bit,
same phase sets, same compositions). usage: ~/orchid-py313/bin/python tools/tests/test_dense_field_grid.py"""
import gzip
import json
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import build_assessed_assemblages as B  # noqa: E402

fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


# ---- synthetic field: a 1000 x 1000 raster, field = rows 100..899, cols 200..799
nT = nx = 1000
xs = np.linspace(0.0, 1.0, nx)
Ts = np.linspace(0.0, 1000.0, nT)
mask = np.zeros((nT, nx), dtype=bool)
mask[100:900, 200:800] = True
cells = B.dense_grid_cells(mask, xs, Ts, False)
check("synthetic: points inside the field", all(mask[i, j] for i, j, _ in cells))
check("synthetic: clearance >= MARGIN", all(c >= B.MARGIN for _, _, c in cells))
check("synthetic: at most %d points" % B.DENSE_MAX_PER_FIELD, 0 < len(cells) <= B.DENSE_MAX_PER_FIELD)
check("synthetic: deterministic", cells == B.dense_grid_cells(mask, xs, Ts, False))
liq = B.dense_grid_cells(mask, xs, Ts, True)
check("liquid band: only within %d K of the lower boundary" % B.DENSE_LIQUID_BAND_K,
      bool(liq) and all(Ts[i] - Ts[100] <= B.DENSE_LIQUID_BAND_K + 1e-9 for i, _, _ in liq))
thin = np.zeros((nT, nx), dtype=bool)
thin[500:503, 100:900] = True
check("sliver thinner than the margin: no points", B.dense_grid_cells(thin, xs, Ts, False) == [])
check("validated set is Cu-Zr only", B.DENSE_VALIDATED == {"CU-ZR-p371"})

# ---- release: the 177 points of the final fit
OD = os.path.expanduser("~/projects/orchid-data/assessed-phase-assemblages")
fp = os.path.join(OD, "espei-zpf-dense", "CU-ZR-p371.json")
if os.path.exists(fp):
    fx = json.load(open(os.path.join(HERE, "fixtures", "cuzr_dense177.json")))["points"]
    d = json.load(open(fp))
    tdb = {"BCC_A2", "CU10ZR7", "CU51ZR14", "CU5ZR_C15B", "CU8ZR3", "CUZR_B2", "CUZR2", "FCC_A1", "HCP_A3", "LIQUID"}   # CuZr-94Zen-LB
    got = []
    for T, reg in zip(d["conditions"]["T"], d["values"]):
        ph = [v[0] for v in reg if v[0] != "__HYPERPLANE__"]
        if all(p in tdb for p in ph):
            got.append((T, sorted(ph), [v[2][0] for v in reg if v[2][0] is not None][-1]))
    want = [(T, sorted(ph), x) for T, ph, x in fx]
    check("release Cu-Zr: %d dense points map to CuZr-94Zen-LB (fit used %d)" % (len(got), len(want)), len(got) == len(want))
    check("release Cu-Zr: coordinates and assemblages equal the fit's 177 points, in order", got == want)
    recs = [json.loads(l) for l in gzip.open(os.path.join(OD, "dense-field-grid", "records.jsonl.gz"), "rt")]
    check("dense records: point_method dense_field_grid, flags, no tie line",
          all(r["point_method"] == "dense_field_grid" and "interior_grid" in r["flags"] and r["tie_line"] is None for r in recs))
    check("dense records: unique ids", len({r["id"] for r in recs}) == len(recs))
    byf = {}
    for r in recs:
        byf.setdefault((r["diagram"], r["field"]), []).append(r["dense_grid"]["field_weight_share"])
    check("dense records: field weight shares sum to one per field",
          all(abs(sum(v) - 1.0) < 1e-3 * len(v) for v in byf.values()))
    check("dense records: budget structure of the primary records",
          all(r["uncertainty_budget"]["quantities"]["temperature_K"]["sigma_method"]["class"] == "assessed_phase_diagram" for r in recs))
else:
    print("SKIP release checks (no espei-zpf-dense in the checkout)")

print("%d failure(s)" % len(fails))
sys.exit(1 if fails else 0)
