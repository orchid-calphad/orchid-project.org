#!/usr/bin/env python3
"""Add a stable-phase grid to diagram JSONs (2026-10-22): pycalphad equilibrium at every dx=0.01 x dT=10 K point; each
row (temperature) is run-length encoded as [label_index, count] pairs; labels are sorted phase sets joined by "+".
  add_phase_grid.py DIAGRAM_JSON_DIR TDB_DIR [names...] [--workers N]"""
import argparse, json, os, sys, time, multiprocessing as mp, warnings
warnings.filterwarnings("ignore")

def one(args):
    try:
        return _one(args)
    except ZeroDivisionError:
        try:
            return _one((args[0], args[1], 1e-3))
        except Exception as e:
            return os.path.basename(args[0]), "failed", f"{type(e).__name__}: {str(e)[:160]}"
    except Exception as e:
        return os.path.basename(args[0]), "failed", f"{type(e).__name__}: {str(e)[:160]}"

def _one(args):
    jf, tdbdir = args[:2]; EDGE = args[2] if len(args) > 2 else 1e-5
    import numpy as np
    from pycalphad import Database, equilibrium, variables as v
    d = json.load(open(jf))
    if d.get("status") != "ok" or "grid" in d: return os.path.basename(jf), "skip", 0
    t0 = time.time()
    db = Database(os.path.join(tdbdir, d["tdb"])); comps = d["components"] + (["VA"] if "VA" in db.elements else [])
    b = d["x_component"]; xs = np.round(np.linspace(0.0, 1.0, 101), 4); xs[0], xs[-1] = EDGE, 1 - EDGE
    Ts = np.arange(d["T_range"][0], d["T_range"][1] + 1, 10.0)
    eq = equilibrium(db, comps, d["phases"], {v.X(b): xs, v.T: Ts, v.P: 101325, v.N: 1})
    ph = eq.Phase.squeeze()                         # dims: T, X_b, vertex (P, N squeezed)
    ph = ph.transpose("T", f"X_{b}", "vertex").values
    labels, rows = [], []
    for i in range(ph.shape[0]):
        row = []
        for j in range(ph.shape[1]):
            s = "+".join(sorted({str(p) for p in ph[i, j] if str(p)}))
            if s not in labels: labels.append(s)
            k = labels.index(s)
            if row and row[-1][0] == k: row[-1][1] += 1
            else: row.append([k, 1])
        rows.append(row)
    d["grid"] = {"x_step": 0.01, "x_points": 101, "T0": float(Ts[0]), "dT": 10.0, "labels": labels, "rows": rows,
                 "note": f"stable phase set from pycalphad equilibrium; first/last columns evaluated at x={EDGE:g} and 1-{EDGE:g}"}
    json.dump(d, open(jf, "w"), separators=(",", ":"))
    return os.path.basename(jf), "ok", round(time.time() - t0, 1)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("dir"); ap.add_argument("tdbdir"); ap.add_argument("names", nargs="*")
    ap.add_argument("--workers", type=int, default=8); a = ap.parse_args()
    files = [os.path.join(a.dir, n) for n in (a.names or sorted(os.listdir(a.dir))) if n.endswith(".json")]
    with mp.Pool(a.workers, maxtasksperchild=4) as pool:
        for r in pool.imap_unordered(one, [(f, a.tdbdir) for f in files]): print(*r, flush=True)

if __name__ == "__main__":
    main()
