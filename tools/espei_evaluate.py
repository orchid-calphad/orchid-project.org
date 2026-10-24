#!/usr/bin/env python3
"""Evaluate ESPEI datasets against their TDB with the ESPEI 0.9.1 likelihood functions (run in the ESPEI interpreter).

usage: <espei python> tools/espei_evaluate.py manifest.json out.json [workers]
manifest: [{"rel": ..., "file": <dataset json>, "tdb": <TDB path>}, ...]
out: {rel: {"status": ok|nonfinite|error|tdb_error|timeout|empty, "kind", "n", "n_nonfinite", "median_abs",
            "max_abs", "bad_regions" (ZPF: indices of regions with a non-finite driving force), "error"}}

ZPF: driving forces (J/mol) of every vertex (espei.error_functions.zpf_error); equilibrium properties: residuals
(calc_prop_differences); activities: residuals and the activity likelihood; non-equilibrium (`solver`): the
non-equilibrium likelihood. Parameters are the TDB's own (no fitting). Called by tools/build_espei_bundle.py.
"""
import json
import os
import signal
import sys
import warnings
from multiprocessing import Pool

warnings.filterwarnings("ignore")
TIMEOUT_S = 300


class _Timeout(Exception):
    pass


def _alarm(*_):
    raise _Timeout()


def _run(job):
    tdb, items = job
    import numpy as np
    from pycalphad import Database
    from espei.datasets import load_datasets
    from espei.error_functions.equilibrium_thermochemical_error import (get_equilibrium_thermochemical_data,
                                                                        calc_prop_differences)
    from espei.error_functions.non_equilibrium_thermochemical_error import (
        get_thermochemical_data, calculate_non_equilibrium_thermochemical_probability)
    from espei.error_functions.zpf_error import get_zpf_data, calculate_zpf_driving_forces
    from espei.error_functions import activity_error as A
    signal.signal(signal.SIGALRM, _alarm)
    out = {}
    try:
        dbf = Database(tdb)
    except Exception as e:  # noqa: BLE001
        msg = "%s: %s" % (type(e).__name__, str(e)[:200])
        return {it["rel"]: {"status": "tdb_error", "error": msg} for it in items}
    for it in items:
        signal.alarm(TIMEOUT_S)
        bad = []
        try:
            ds = load_datasets([it["file"]])
            d = ds.all()[0]
            comps = sorted(set(d["components"]))
            phases = sorted(dbf.phases)
            o = d["output"]
            if o == "ZPF":
                zd = get_zpf_data(dbf, comps, phases, ds, {})
                df, _ = calculate_zpf_driving_forces(zd)
                vals = [x for reg in df for x in np.asarray(reg, dtype=float).ravel().tolist()]
                kind = "zpf_driving_force_J"
                if not np.all(np.isfinite(vals)):
                    # identify the regions exactly: each region evaluated alone (indices of the file's `values`)
                    raw = json.load(open(it["file"]))
                    Ts = raw["conditions"]["T"] if isinstance(raw["conditions"]["T"], list) else \
                        [raw["conditions"]["T"]] * len(raw["values"])
                    one = it["file"] + ".one.json"
                    for k, reg in enumerate(raw["values"]):
                        r1 = dict(raw, values=[reg], conditions=dict(raw["conditions"], T=[Ts[k]]),
                                  phases=sorted({vx[0] for vx in reg if vx[0] != "__HYPERPLANE__"}))
                        with open(one, "w") as fh:
                            json.dump(r1, fh)
                        try:
                            d1 = load_datasets([one])
                            a, _ = calculate_zpf_driving_forces(get_zpf_data(dbf, comps, phases, d1, {}))
                            ok = all(np.all(np.isfinite(np.asarray(x, dtype=float))) for x in a)
                        except Exception:  # noqa: BLE001
                            ok = False
                        if not ok:
                            bad.append(k)
            elif o.startswith("ACR"):
                r, _ = A.calculate_activity_residuals(dbf, comps, phases, ds, parameters={})
                vals = list(np.asarray(r, dtype=float).ravel())
                if not np.isfinite(A.calculate_activity_error(dbf, comps, phases, ds, parameters={})):
                    vals.append(float("nan"))
                kind = "activity_residual"
            elif "solver" in d:
                td = get_thermochemical_data(dbf, comps, phases, ds)
                vals = [float(calculate_non_equilibrium_thermochemical_probability(td, dbf))]
                kind = "nonequilibrium_log_likelihood"
            else:
                vals = []
                for x in get_equilibrium_thermochemical_data(dbf, comps, phases, ds, {}):
                    r, _ = calc_prop_differences(x, np.array([]))
                    vals.extend(np.asarray(r, dtype=float).ravel().tolist())
                kind = "residual"
            a = np.asarray(vals, dtype=float)
            fin = a[np.isfinite(a)]
            st = "empty" if not len(a) else ("ok" if len(fin) == len(a) else "nonfinite")
            out[it["rel"]] = {"status": st, "kind": kind, "n": int(len(a)), "n_nonfinite": int(len(a) - len(fin)),
                              "median_abs": round(float(np.median(np.abs(fin))), 6) if len(fin) else None,
                              "max_abs": round(float(np.max(np.abs(fin))), 6) if len(fin) else None,
                              "bad_regions": bad}
        except _Timeout:
            out[it["rel"]] = {"status": "timeout"}
        except Exception as e:  # noqa: BLE001
            out[it["rel"]] = {"status": "error", "error": "%s: %s" % (type(e).__name__, str(e)[:200])}
        finally:
            signal.alarm(0)
    return out


def main():
    man, outp = sys.argv[1], sys.argv[2]
    nw = int(sys.argv[3]) if len(sys.argv) > 3 else 6
    items = json.load(open(man))
    by = {}
    for it in items:
        by.setdefault(it["tdb"], []).append(it)
    jobs = sorted(by.items(), key=lambda j: -len(j[1]))
    res = {}
    if jobs:
        with Pool(min(nw, len(jobs))) as pool:
            for r in pool.imap_unordered(_run, jobs):
                res.update(r)
    with open(outp + ".tmp", "w") as fh:
        json.dump(res, fh, sort_keys=True)
    os.replace(outp + ".tmp", outp)


if __name__ == "__main__":
    main()
