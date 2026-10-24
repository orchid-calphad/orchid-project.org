#!/usr/bin/env python3
"""Solid vs dashed (estimated) boundaries of the assessed handbook diagrams (2026-10-24, James Saal review, comment 2).

The tracer (open-tdb build/okamoto-pilot/tools/pd_fields.py, l. 192-217) joins the short strokes of a dashed line into
chains: a stroke shorter than 12 pt whose end continues, within 25 degrees and across a gap of 0.3-6 pt, into another
short stroke pointing back along it. It stores only the number of such joins (`fields.json` `dash_joins`). This module
re-derives the same joins, deterministically, from the traced boundaries in `calibrated.json` (the count is checked
against the tracer's own count), and classifies, for each released point, the boundary lines that bound its field along
the point's isotherm and isopleth (the lines that define its positional margin): all printed solid, all dashed
(estimated), or mixed. Lines the tracer did not use as phase boundaries (`fields.json` `non_phase_lines`: pointer arrows,
Tc lines, cooling curves, the dashed second assessment of an overlay) are ignored. Standard library + numpy.
"""
import math

import numpy as np

SHORT_PT, GAP_MIN_PT, GAP_MAX_PT, ANGLE_DEG = 12.0, 0.3, 6.0, 25.0   # as in pd_fields.py
SNAP_NORM = 0.02     # a margin crossing is attributed to a traced line within this frame-normalised distance


def _plen(P, cx, cy):
    return sum(math.hypot((q[0] - p[0]) * 100 / cx["a"], (q[1] - p[1]) / cy["a"]) for p, q in zip(P, P[1:]))


def _dir(P, end, cx, cy):
    a, b = (P[-2], P[-1]) if end else (P[1], P[0])
    dx, dy = (b[0] - a[0]) * 100 / cx["a"], (b[1] - a[1]) / cy["a"]
    n = math.hypot(dx, dy) or 1.0
    return dx / n, dy / n


def dash_joins(cal):
    """The tracer's dash joins: [(id_a, id_b)], same rule and order as pd_fields.py."""
    cx, cy = cal["calibration"]["x"], cal["calibration"]["y"]
    B = cal["boundaries"]
    short = [(k, b) for k, b in enumerate(B) if len(b["points"]) >= 2 and _plen(b["points"], cx, cy) < SHORT_PT]
    ends = [(k, e, b["points"][-1 if e else 0], _dir(b["points"], e, cx, cy)) for k, b in short for e in (0, 1)]
    cosmax = math.cos(math.radians(ANGLE_DEG))
    out = []
    for i in range(len(ends)):
        ki, ei, pi_, di = ends[i]
        best = None
        for j in range(len(ends)):
            kj, ej, pj, dj = ends[j]
            if kj == ki:
                continue
            gx, gy = (pj[0] - pi_[0]) * 100 / cx["a"], (pj[1] - pi_[1]) / cy["a"]
            gap = math.hypot(gx, gy)
            if not (GAP_MIN_PT < gap <= GAP_MAX_PT):
                continue
            if (gx * di[0] + gy * di[1]) / gap < cosmax or -(dj[0] * di[0] + dj[1] * di[1]) < cosmax:
                continue
            if best is None or gap < best[0]:
                best = (gap, kj)
        if best:
            out.append((B[ki]["id"], B[best[1]]["id"]))
    return out


class DiagramLines:
    """Traced phase-boundary lines of one diagram in frame-normalised coordinates, each marked dashed or solid."""

    def __init__(self, cal, fields_json, x_range, T_range):
        self.joins = dash_joins(cal)
        dashed = {i for j in self.joins for i in j}
        nonphase = {d["id"] for d in (fields_json or {}).get("non_phase_lines") or []}
        (x0, x1), (T0, T1) = x_range, T_range
        sx, sT = (x1 - x0) or 1.0, (T1 - T0) or 1.0
        segs, flags = [], []
        for b in cal["boundaries"]:
            if b["id"] in nonphase or len(b["points"]) < 2:
                continue
            P = np.array([[(p[0] - x0) / sx, (p[1] - T0) / sT] for p in b["points"]], float)
            for p, q in zip(P, P[1:]):
                segs.append((p, q))
                flags.append(b["id"] in dashed)
        for a_, b_ in self.joins:   # the gaps between dashes are dashed line too
            pa = self._end(cal, a_, b_)
            if pa is not None:
                p, q = pa
                segs.append((np.array([(p[0] - x0) / sx, (p[1] - T0) / sT]), np.array([(q[0] - x0) / sx, (q[1] - T0) / sT])))
                flags.append(True)
        self.A = np.array([s[0] for s in segs]) if segs else np.zeros((0, 2))
        self.B = np.array([s[1] for s in segs]) if segs else np.zeros((0, 2))
        self.dashed = np.array(flags, bool)
        self.sx, self.sT, self.x0, self.T0 = sx, sT, x0, T0
        self.n_dash_joins = len(self.joins)

    @staticmethod
    def _end(cal, ida, idb):
        bymap = {b["id"]: b["points"] for b in cal["boundaries"]}
        Pa, Pb = bymap.get(ida), bymap.get(idb)
        if not Pa or not Pb:
            return None
        cands = [(math.hypot(p[0] - q[0], (p[1] - q[1]) / 1000.0), p, q) for p in (Pa[0], Pa[-1]) for q in (Pb[0], Pb[-1])]
        _, p, q = min(cands, key=lambda c: c[0])
        return p, q

    def nearest(self, x, T):
        """(distance, dashed) of the traced boundary nearest to (x, T), frame-normalised; (inf, None) if none."""
        if not len(self.A):
            return math.inf, None
        p = np.array([(x - self.x0) / self.sx, (T - self.T0) / self.sT])
        d = self.B - self.A
        L2 = (d * d).sum(1)
        t = np.clip(((p - self.A) * d).sum(1) / np.where(L2 > 0, L2, 1), 0, 1)
        dist = np.hypot(*(self.A + d * t[:, None] - p).T)
        k = int(np.argmin(dist))
        return float(dist[k]), bool(self.dashed[k])

    def basis_at(self, crossings):
        """crossings: [(x, T)] where the point's isotherm / isopleth leaves its field (frame edges excluded).
        Returns (boundary_basis, n_solid, n_dashed): printed_solid | dashed_estimated | mixed | None (no traced line)."""
        ns = nd = 0
        for x, T in crossings:
            dist, dashed = self.nearest(x, T)
            if dashed is None or dist > SNAP_NORM:
                continue
            if dashed:
                nd += 1
            else:
                ns += 1
        if not ns and not nd:
            return None, 0, 0
        return ("dashed_estimated" if not ns else "printed_solid" if not nd else "mixed"), ns, nd


BASIS_TEXT = {"printed_solid": "every boundary line that bounds the point's field along its isotherm and isopleth is printed solid",
              "dashed_estimated": "every such boundary line is dashed (estimated by the evaluator)",
              "mixed": "the field is bounded by solid and dashed (estimated) lines along the point's isotherm and isopleth"}
