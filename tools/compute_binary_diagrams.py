#!/usr/bin/env python3
"""Compute binary phase diagrams (ZPF boundaries + invariants) from TDB files with pycalphad's BinaryStrategy and
write one compact JSON per TDB for the ORCHID phase-diagram viewer (2026-10-22).
  compute_binary_diagrams.py OUT_DIR TDB [TDB ...] [--workers N] [--tmax 2500]"""
import argparse, hashlib, json, os, re, sys, time, multiprocessing as mp
import warnings; warnings.filterwarnings("ignore")

def elements_from_name(fn):
    stem = os.path.basename(fn).split("-")[0]
    return re.findall(r"[A-Z][a-z]?", stem)

def r(a, n):  # round arrays for compact JSON
    return [round(float(x), n) for x in (a if getattr(a, "shape", None) != () else [a])]

def one(args):
    tdb, out, tmax = args
    from pycalphad import Database, variables as v
    from pycalphad.mapping import BinaryStrategy
    els = elements_from_name(tdb); rec = {"tdb": os.path.basename(tdb), "elements": els}
    try:
        h = hashlib.sha256(open(tdb, "rb").read()).hexdigest()
        db = Database(tdb); comps = [e.upper() for e in els] + (["VA"] if "VA" in db.elements else [])
        a, b = comps[0], comps[1]
        phases = sorted(db.phases.keys()); t0 = time.time()
        s = BinaryStrategy(db, comps, phases, {v.X(b): (0, 1, 0.02), v.T: (300, tmax, 10), v.P: 101325, v.N: 1})
        s.do_map()
        regions = []
        for sd in s.get_tieline_data(v.X(b), v.T):
            regions.append([{"phase": p.phase, "x": r(p.x, 5), "T": r(p.y, 2)} for p in sd.data])
        inv = []
        for sd in s.get_invariant_data(v.X(b), v.T):
            pts = [{"phase": p.phase, "x": round(float(p.x), 5)} for p in sd.data]
            inv.append({"T": round(float(sd.data[0].y), 2), "points": sorted(pts, key=lambda q: q["x"])})
        rec.update(sha256=h, x_component=b, components=[a, b], T_range=[300, tmax], phases=phases,
                   regions=regions, invariants=sorted(inv, key=lambda q: -q["T"]), seconds=round(time.time() - t0, 1),
                   engine="pycalphad " + __import__("pycalphad").__version__, status="ok")
    except Exception as e:
        rec.update(status="failed", error=f"{type(e).__name__}: {str(e)[:300]}")
    json.dump(rec, open(os.path.join(out, os.path.basename(tdb).rsplit(".", 1)[0] + ".json"), "w"), separators=(",", ":"))
    return rec["tdb"], rec["status"], rec.get("seconds")

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("tdbs", nargs="+")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--tmax", type=float, default=2500)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    with mp.Pool(a.workers, maxtasksperchild=4) as pool:
        for res in pool.imap_unordered(one, [(t, a.out, a.tmax) for t in a.tdbs]):
            print(*res, flush=True)

if __name__ == "__main__":
    main()
