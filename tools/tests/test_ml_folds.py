#!/usr/bin/env python3
"""Tests of tools/ml_folds.py (machine-learning folds, 2026-10-24) and, when the release is present, of the fold
columns in the literature JSONL, literature_values / dft_parity Parquet and the assessed points: rule reproduced from
sha256, every record assigned, no system / paper split across folds, folds roughly balanced.
usage: ~/orchid-py313/bin/python tools/tests/test_ml_folds.py [--data ~/projects/orchid-data]"""
import argparse
import collections
import gzip
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import ml_folds as M  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--data", default=os.path.expanduser("~/projects/orchid-data"))
A = ap.parse_args()
fails = []


def check(name, cond):
    print(("PASS " if cond else "FAIL ") + name)
    if not cond:
        fails.append(name)


ref = lambda k: int(hashlib.sha256(k.encode()).hexdigest(), 16) % 5   # noqa: E731
check("fold = int(sha256(key), 16) mod 5", all(M.fold(k) == ref(k) for k in ("AL-NI", "CU-ZR", "doi:10.1016/x", "O")))
check("chemical_system sorts and upper-cases", M.chemical_system(["Ni", "al"]) == "AL-NI")
f = M.fields("AL-NI", "doi:10.1/x")
check("fields() keys", set(f) == {"chemical_system", "paper_id", "fold_lso", "fold_lpo"} and f["fold_lso"] == ref("AL-NI"))
check("missing keys give None", M.fields(None, None)["fold_lso"] is None)
check("balance counts", M.balance([{"fold_lso": 1, "fold_lpo": 2}])["fold_lso"][1] == 1)


def consistent(rows, name, sys_key="chemical_system", pap_key="paper_id"):
    bad = [r for r in rows if r.get("fold_lso") != ref(r[sys_key]) or r.get("fold_lpo") != ref(r[pap_key])]
    check("%s: every row's folds follow the rule (%d rows)" % (name, len(rows)), rows and not bad)
    by = collections.defaultdict(set)
    for r in rows:
        by[r[sys_key]].add(r["fold_lso"])
    check("%s: no system in two folds" % name, all(len(v) == 1 for v in by.values()))
    for fk in ("fold_lso", "fold_lpo"):
        c = collections.Counter(r[fk] for r in rows)
        n = len(rows)
        print("   %s %s: %s" % (name, fk, ", ".join("%d: %d (%.1f%%)" % (k, c[k], 100 * c[k] / n) for k in range(5))))
        check("%s %s: every fold holds 8-32%% of rows" % (name, fk), all(0.08 * n <= c[k] <= 0.32 * n for k in range(5)))


jl = os.path.join(A.data, "literature", "orchid-literature.jsonl.gz")
if os.path.exists(jl):
    rows = [json.loads(line) for line in gzip.open(jl, "rt", encoding="utf-8")]
    check("literature JSONL: chemical_system = system, paper_id = paper",
          all(r.get("chemical_system") == r["system"] and r.get("paper_id") == r["paper"] for r in rows))
    consistent(rows, "literature JSONL")
try:
    import pyarrow.parquet as pq
except ImportError:
    pq = None
if pq:
    for name, path, pk in (("literature_values", "parquet/literature_values.parquet", "paper_id"),
                           ("dft_parity", "parquet/dft_parity.parquet", "paper_id"),
                           ("assessed points", "assessed-phase-assemblages/records.parquet", "paper_id")):
        p = os.path.join(A.data, path)
        if os.path.exists(p):
            t = pq.read_table(p, columns=["chemical_system", pk, "fold_lso", "fold_lpo"]).to_pylist()
            consistent(t, name, pap_key=pk)
print("%d failed" % len(fails))
sys.exit(1 if fails else 0)
