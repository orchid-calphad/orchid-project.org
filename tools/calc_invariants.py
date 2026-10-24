#!/usr/bin/env python3
"""Distinct calculated invariants of one SGTE phase-diagram file (audit A5 F1, 2026-10-24).

pycalphad's mapper finds the same invariant from several starting points, so `invariants` in
phase-diagrams/sgte/<file>.json lists the mapper's hits, which can repeat (9,481 of 24,828 hits in v2026.10).
The Parquet table phase_diagram_invariants has one row per distinct invariant.

Duplicate key (one calculated diagram): the same multiset of phases, with the mapper's '#n' suffix of a second
composition set ignored (FCC_A1 / FCC_A1#2 are the same phase name; which one the mapper calls #2 is arbitrary),
temperatures within DT_K of each other (chained, in increasing T), and the sorted phase compositions within DX of
each other. The representative is the hit with the lowest index in the file; `n_mapper_hits` counts the hits it
stands for. Two invariants that share phases and temperature but not compositions (e.g. a compound's polymorphic
transition meeting the liquid on either side) stay distinct.

Standard library only.
"""
from __future__ import annotations

import re

DT_K = 0.1
DX = 0.01


def base_phase(name):
    return re.sub(r"#\d+$", "", name or "")


def _key(inv):
    pts = inv.get("points") or []
    return tuple(sorted(base_phase(p.get("phase")) for p in pts))


def _xs(inv):
    return sorted(float(p["x"]) for p in inv.get("points") or [] if p.get("x") is not None)


def _same_x(a, b):
    return len(a) == len(b) and all(abs(u - v) <= DX for u, v in zip(a, b))


def distinct(invs):
    """[(index of the representative hit, number of hits)] in file order."""
    groups = {}
    for k, inv in enumerate(invs):
        if inv.get("T") is None:
            groups[("_noT", k)] = [k]
            continue
        groups.setdefault(_key(inv), []).append(k)
    reps = []
    for key, ks in groups.items():
        if key and key[0] == "_noT":
            reps.append((ks[0], 1))
            continue
        ks = sorted(ks, key=lambda k: (float(invs[k]["T"]), k))
        clusters = []   # each: [members]; chain on T, then compositions
        for k in ks:
            T, xs = float(invs[k]["T"]), _xs(invs[k])
            for c in clusters:
                last = c[-1]
                if T - float(invs[last]["T"]) <= DT_K and _same_x(_xs(invs[c[0]]), xs):
                    c.append(k)
                    break
            else:
                clusters.append([k])
        reps += [(min(c), len(c)) for c in clusters]
    return sorted(reps)
