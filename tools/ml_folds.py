#!/usr/bin/env python3
"""Machine-learning folds of the ORCHID release (2026-10-24, PI-approved): two deterministic 5-fold assignments
attached to every record, so that a model can be evaluated on chemical systems or papers it has not seen.

  chemical_system  sorted upper-case element symbols joined by "-" (the record's `system`, e.g. "AL-NI")
  paper_id         stable id of the source (`doi:<doi>` or `sha256:<pdf prefix>`; Dataset 6: the original evaluation's
                   `doi:<doi>`, else `handbook:<diagram>`; DFT pairs: the literature record's paper)
  fold_lso         leave-chemical-system-out fold: int(sha256(chemical_system).hexdigest(), 16) mod 5
  fold_lpo         leave-paper-out fold:           int(sha256(paper_id).hexdigest(), 16) mod 5

Hashes are of the UTF-8 string itself (no salt), so anyone can recompute them. The fold is keyed on the record's exact
chemical system: unary and lower-order subsystems are systems of their own and fall in their own folds, so a held-out
ternary A-B-C may still have A-B, A-C, B-C and unary data in training. For strict extrapolation to unseen chemistry,
drop from training every record whose system is a subsystem of a held-out system (documented in the dataset cards).
The release keeps one split ("full"); the folds are columns, not splits. Standard library only.
"""
import hashlib

N_FOLDS = 5


def fold(key):
    return int(hashlib.sha256(str(key).encode("utf-8")).hexdigest(), 16) % N_FOLDS


def chemical_system(elements):
    return "-".join(sorted(str(e).upper() for e in elements))


def fields(system, paper_id):
    """{chemical_system, paper_id, fold_lso, fold_lpo} for a record of `system` read from `paper_id`."""
    return {"chemical_system": system, "paper_id": paper_id,
            "fold_lso": fold(system) if system else None, "fold_lpo": fold(paper_id) if paper_id else None}


def annotate(r, system_key="system", paper_key="paper"):
    r.update(fields(r.get(system_key), r.get(paper_key)))
    return r


def balance(rows):
    """{fold_lso: {k: n}, fold_lpo: {k: n}} for a list of dicts with the fold fields."""
    out = {}
    for f in ("fold_lso", "fold_lpo"):
        c = {k: 0 for k in range(N_FOLDS)}
        for r in rows:
            if r.get(f) is not None:
                c[r[f]] += 1
        out[f] = c
    return out


DOC = ("`chemical_system` (sorted elements, = `system`), `paper_id` (stable source id), `fold_lso` = leave-chemical-"
       "system-out fold int(sha256(chemical_system), 16) mod 5, `fold_lpo` = leave-paper-out fold int(sha256(paper_id), "
       "16) mod 5 (tools/ml_folds.py). Subsystems are their own systems (a held-out ternary may have its binaries in "
       "training); drop subsystems of held-out systems for strict chemistry extrapolation")
