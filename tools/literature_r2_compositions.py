#!/usr/bin/env python3
"""Round-2 audit rules for elements, chemical systems and compositions (ORCHID audit A1, 2026-10-24, PI "fix all").

Runs inside tools/literature_gates.apply() (hooks: `before_compositions` ahead of G6/G9, `after_system` after the G4
system gates). Same evidence model as the gates: a correction keeps the stored value (`composition_stored`,
`system_stored`, ...) and adds a `corrections` entry {gate "A1", field, before, after, rule, evidence[, visual_check]};
a suspicion that the page has not settled is a `gate_failures` entry. Page verdicts live in
tools/literature_visual_checks.json (keys below); no page image or text is stored.

  C1  fractions printed under an "at.%" header divided by 100 a second time.  Detection per table (asset): every
      printed composition converted from atomic percent is <= 0.0101 and the stoichiometric phases fit 100 x the
      stored value.  Page key `a1_c1:<asset>` (verdict confirmed -> x 100; rejected -> genuinely dilute, nothing
      changes); single cells in otherwise normal tables: `a1_c1_row:<record id>`.  Unchecked detections are flagged
      `composition_double_division_suspect` (not ESPEI-ready).
  C2  the printed composition header names the wrong element.  Detection per binary table: the terminal (pure
      element) and stoichiometric phases fit 1 - x better than x.  Page key `a1_c2:<asset>` (confirmed -> x -> 1 - x,
      flag `axis_header_contradicted`); single reversed cells of otherwise consistent tables are flagged
      `composition_axis_suspect` (not ESPEI-ready) unless `a1_c2_row:<record id>` confirms the flip.
  C3  compositions assigned to the wrong phase within one row: when a permutation of the row's printed values puts
      every stoichiometric phase within 0.03 of its formula while the stored order leaves one more than 0.08 off, the
      values are re-assigned by stoichiometry (rule `composition_phase_permutation`, page-checked sample); a page
      check `a1_c3:<record id>` with verdict rejected keeps the stored order.  A stoichiometric label whose printed
      composition lies outside max(registry calculated range, formula) + 0.03 and that no permutation explains is
      flagged `label_conflicts_composition` (not ESPEI-ready).
  C4  lower-order data filed under a higher-order system: the record is filed under the chemical system of its own
      state (support of its phase formulas / overall composition; unary for a pure element) and keeps the paper's
      system as `source_system`; folds follow the corrected system (tools/ml_folds.py).
  C6  element symbols: overall-composition keys outside `elements` (stale zeros dropped; non-zero ones flagged),
      record-level page checks for wrong systems / symbols (`a1_c6:<record id>`: system or element rename).
  C7  compositions outside [0, 1] or with basis `unresolved`: not accepted (outside) / not ESPEI-ready (unresolved)
      unless a page fix resolves them.

Standard library only; deterministic.
"""
from __future__ import annotations

import itertools
import json
import os
import re
import unicodedata
from collections import Counter, defaultdict

import build_dft_parity as parity

HERE = os.path.dirname(os.path.abspath(__file__))
GATE = "A1"
RELEASE = "a1-round2 2026-10-24"

# ------------------------------------------------------------------------------------------------- helpers

_STRUCT = r"(?i:hcp|fcc|bcc|dhcp|rhom\.?|hex\.?|tetr\.?|cub\.?|diamond|A1|A2|A3|A4|A5|A7|A8|A9)"
_PURE = re.compile(r"^\s*(?:[αβγδεζηθικλμ]|alpha|beta|gamma|delta|epsilon|HT|LT|ht|lt|" + _STRUCT + r"|[abhr])?[\s\-‐_]*"
                   r"\(?\s*(?:[αβγδεζηθικλμ]|alpha|beta|gamma|delta|epsilon)?[\s\-‐]*([A-Z][a-z]?)\s*\)?"
                   r"(?:[\s\-‐_]*(?:\((?:[hrl]|HT|LT|s|" + _STRUCT + r")\)|HT|LT|" + _STRUCT + r"|[αβγδε]|['′’]+|I{1,3}|IV))?\s*$")
ELEMENTS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr "
               "Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W "
               "Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf".split())


def asset_of(r):
    loc = r.get("locator") or {}
    return loc.get("asset") or loc.get("ref")


def pure_element(label, elements=None):
    """Upper-case element symbol of a pure-element / terminal-solution label ('(αFe)', 'βSn', 'Nd', 'γ-Fe'), or None."""
    if not label:
        return None
    s = unicodedata.normalize("NFKC", str(label)).strip()
    m = _PURE.match(s)
    if not m or m.group(1) not in ELEMENTS:
        return None
    e = m.group(1).upper()
    if elements is not None and e not in elements:
        return None
    return e


def formula_fractions(label, elements=None):
    """{EL: x} of a multi-element formula in the label (elements restricted to `elements` when given), or None."""
    c, _ = parity.parse_phase(label, [e.capitalize() for e in elements] if elements else None)
    if not c:
        return None
    tot = float(sum(c.values()))
    out = {e.upper(): float(n) / tot for e, n in c.items()}
    if elements is not None and not set(out) <= set(elements):
        return None
    return out


def expected_x(label, el, elements):
    """Expected mole fraction of `el` in the phase named `label` (formula or pure element) in a binary, or None."""
    f = formula_fractions(label, elements)
    if f and len(f) == 2:
        return f.get(el, 0.0)
    p = pure_element(label, elements)
    if p:
        return 1.0 if p == el else 0.0
    return None


def printed_cells(r):
    """[(pc, element, value)] of printed single-element mole-fraction compositions of a record."""
    out = []
    for pc in r.get("phase_compositions") or []:
        if pc.get("composition_source") != "printed" or not str(pc.get("basis", "")).startswith("mole_fraction"):
            continue
        comp = pc.get("composition") or {}
        for e, v in comp.items():
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                out.append((pc, e, float(v)))
    return out


def _vc(check):
    if not check:
        return None
    return {k: check[k] for k in ("page", "pdf_sha256_prefix", "image_sha256", "verdict", "seen", "checked")
            if check.get(k) is not None}


def correct(r, field, before, after, rule, evidence, check=None):
    c = {"gate": GATE, "field": field, "before": before, "after": after, "rule": rule, "evidence": evidence}
    if check:
        c["visual_check"] = _vc(check)
    r.setdefault("corrections", []).append(c)


def fail(r, code, effect, detail=None, check=None):
    gf = r.setdefault("gate_failures", [])
    if any(x["code"] == code for x in gf):
        return False
    e = {"gate": GATE, "code": code, "effect": effect}
    if detail:
        e["detail"] = detail
    if check:
        e["visual_check"] = _vc(check)
    gf.append(e)
    fl = r.setdefault("flags", {})
    if effect == "not_accepted":
        fl["extraction_accepted"] = False
        fl["espei_ready"] = False
        fl["espei_certified"] = False
    elif effect == "not_espei_ready":
        if fl.get("espei_ready"):
            fl["espei_ready_before_gates"] = True
        fl["espei_ready"] = False
        fl["espei_certified"] = False
    elif effect == "not_certified":
        fl["espei_certified"] = False
    return True


class _Rep:
    """Minimal report interface (the gates' Report is passed in by literature_gates.apply)."""

    def __init__(self):
        self.counts = Counter()
        self.items = defaultdict(list)

    def add(self, gate, code, r, action, **kw):
        self.counts[(gate, code, action)] += 1
        e = {"id": r["id"], "system": r["system"], "property": r.get("property"), "action": action}
        e.update({k: v for k, v in kw.items() if v is not None})
        self.items[(gate, code)].append(e)


def _set_comp(r, pc, new, rule, evidence, check, field_note=None):
    before = dict(pc.get("composition") or {})
    if before == new:
        return False
    pc.setdefault("composition_stored", before)
    pc["composition"] = new
    correct(r, "phase_compositions[%s]" % pc.get("phase") + (field_note or ""), before, new, rule, evidence, check)
    return True


def _rnd(v):
    return round(v, 6)


# ------------------------------------------------------------------------------------------------- C1

C1_MAX = 0.0101


def c1_tables(records):
    """{asset: rows} of tables whose every printed composition converted from atomic percent is <= 0.0101 (at least
    two values) -- the double-division signature (fractions printed under an at.% header)."""
    by = defaultdict(list)
    for r in records:
        a = asset_of(r)
        if a and r.get("phase_compositions"):
            by[a].append(r)
    out = {}
    for a, rows in by.items():
        vals, conv = [], False
        for r in rows:
            for pc, e, v in printed_cells(r):
                vals.append(v)
                conv = conv or "atomic percent" in str(pc.get("basis"))
        if conv and len(vals) >= 2 and max(vals) <= C1_MAX and max(vals) > 0:
            out[a] = rows
    return out


def _c1_scale(r, check, rule, rep, code):
    n = 0
    for pc in r.get("phase_compositions") or []:
        if pc.get("composition_source") != "printed" or "atomic percent" not in str(pc.get("basis", "")):
            continue
        comp = pc.get("composition") or {}
        if not comp or not all(isinstance(v, (int, float)) for v in comp.values()):
            continue
        if max(comp.values()) > C1_MAX:
            continue
        new = {e: _rnd(v * 100.0) for e, v in comp.items()}
        if _set_comp(r, pc, new, rule, {"printed_header": "at.%", "printed_values": "mole fractions (<= 1)",
                                        "factor": 100}, check):
            pc["basis_stored"] = pc.get("basis")
            pc["basis"] = "mole_fraction (printed as a fraction under an at.% header)"
            n += 1
    if n:
        r.setdefault("flags", {})["composition_double_division_corrected"] = True
        rep.add(GATE, code, r, "corrected", cells=n)
    return n


def c1_double_division(records, visual, rep):
    for a, rows in sorted(c1_tables(records).items()):
        chk = visual.get("a1_c1:" + a)
        v = (chk or {}).get("verdict")
        if v == "confirmed":
            for r in rows:
                _c1_scale(r, chk, "fraction_printed_under_atomic_percent_header", rep, "c1_double_division_table")
        elif v == "rejected":
            rep.counts[(GATE, "c1_genuinely_dilute_on_page", "kept")] += len(rows)
        else:
            for r in rows:
                fail(r, "composition_double_division_suspect", "not_espei_ready",
                     "every printed composition of the table is <= 0.0101 after the at.% conversion (fractions under an "
                     "at.% header?); not yet checked on the page")
                rep.add(GATE, "composition_double_division_suspect", r, "flagged")
    for k, chk in sorted(visual.items()):
        if not k.startswith("a1_c1_row:") or chk.get("verdict") != "confirmed":
            continue
    rows_checked = {k[len("a1_c1_row:"):]: c for k, c in visual.items() if k.startswith("a1_c1_row:")}
    if rows_checked:
        for r in records:
            chk = rows_checked.get(r["id"])
            if not chk or chk.get("verdict") != "confirmed":
                continue
            phases = set(chk.get("phases") or [])
            n = 0
            for pc in r.get("phase_compositions") or []:
                if phases and pc.get("phase") not in phases:
                    continue
                comp = pc.get("composition") or {}
                if pc.get("composition_source") == "printed" and comp and max(comp.values()) <= C1_MAX:
                    new = {e: _rnd(v * 100.0) for e, v in comp.items()}
                    if _set_comp(r, pc, new, "fraction_printed_under_atomic_percent_header", {"factor": 100}, chk):
                        pc["basis_stored"] = pc.get("basis")
                        pc["basis"] = "mole_fraction (printed as a fraction under an at.% header)"
                        n += 1
            if n:
                r.setdefault("flags", {})["composition_double_division_corrected"] = True
                rep.add(GATE, "c1_double_division_cell", r, "corrected", cells=n)


# ------------------------------------------------------------------------------------------------- C2

def c2_cell_fit(r, pc, e, v):
    """'ok' | 'rev' | None for one printed binary cell checked against its pure-element / stoichiometric label."""
    els = r.get("elements") or []
    if len(els) != 2 or len(pc.get("composition") or {}) != 1 or not 0 <= v <= 1:
        return None
    fx = expected_x(pc.get("phase"), e, els)
    if fx is None or abs(fx - 0.5) < 0.02:
        return None
    if 0 < v <= C1_MAX and fx > 0.05 and abs(100 * v - fx) <= 0.01:
        return "c1"      # a stray fraction under an at.% header (C1), not a reversed axis
    if abs(v - fx) <= 0.06:
        return "ok"
    if abs(v - (1 - fx)) <= 0.06:
        return "rev"
    return "other"


def c2_tables(records):
    by = defaultdict(list)
    for r in records:
        a = asset_of(r)
        if a and r.get("phase_compositions") and len(r.get("elements") or []) == 2:
            by[a].append(r)
    out = {}
    for a, rows in by.items():
        st = Counter()
        for r in rows:
            for pc, e, v in printed_cells(r):
                f = c2_cell_fit(r, pc, e, v)
                if f:
                    st[f] += 1
        out[a] = (rows, st)
    return out


def rev_labels(rows):
    return {pc.get("phase") for r in rows for pc, e, v in printed_cells(r) if c2_cell_fit(r, pc, e, v) == "rev"}


def _flip(r, pc, chk, rule, evidence):
    comp = pc.get("composition") or {}
    if len(comp) != 1:
        return False
    (e, v), = comp.items()
    return _set_comp(r, pc, {e: _rnd(1.0 - v)}, rule, evidence, chk)


def _c2_majority(rows, st):
    return st["rev"] >= 3 and st["rev"] > 2 * st["ok"] and len(rev_labels(rows)) >= 2


def c2_axis(records, visual, rep):
    """Table pass: a binary table whose terminal and stoichiometric phases fit 1 - x (>= 3 cells of >= 2 labels, more
    than twice the cells that fit x) is flipped when the page confirms the header names the wrong element."""
    for a, (rows, st) in sorted(c2_tables(records).items()):
        if not _c2_majority(rows, st):
            continue
        chk = visual.get("a1_c2:" + a)
        v = (chk or {}).get("verdict")
        if v == "confirmed":
            for r in rows:
                n = 0
                for pc, e, x in printed_cells(r):
                    if len(pc.get("composition") or {}) == 1 and _flip(
                            r, pc, chk, "composition_axis_from_terminal_and_stoichiometric_phases",
                            {"printed_header_element": e.capitalize(), "fits": dict(st)}):
                        n += 1
                if n:
                    r.setdefault("flags", {})["axis_header_contradicted"] = True
                    rep.add(GATE, "c2_axis_flipped_table", r, "corrected", cells=n)
        elif v != "rejected":
            for r in rows:
                fail(r, "composition_axis_suspect", "not_espei_ready",
                     "terminal / stoichiometric phases of the table fit 1 - x better than x (%s); not yet checked on the "
                     "page" % dict(st))
                rep.add(GATE, "composition_axis_suspect_table", r, "flagged")


def c2_rows(records, visual, rep):
    """Row pass (after C3): single cells of otherwise consistent tables.  A terminal (pure-element) phase at 1 - x is
    an axis suspect (the whole row is flipped when `a1_c2_row:<id>` confirms it); a formula label at 1 - x is a label
    whose element order the page or the reading reversed (`label_conflicts_composition`, C3 class); a cell at 1/100 of
    its formula fraction is a stray C1 fraction."""
    rows_checked = {k[len("a1_c2_row:"):]: c for k, c in visual.items() if k.startswith("a1_c2_row:")}
    for a, (rows, st) in sorted(c2_tables(records).items()):
        if _c2_majority(rows, st):
            continue
        for r in rows:
            cells = [(pc, e, x, c2_cell_fit(r, pc, e, x)) for pc, e, x in printed_cells(r)]
            c1 = [pc for pc, e, x, f in cells if f == "c1"]
            if c1 and not r.get("flags", {}).get("composition_double_division_corrected") and not (
                    visual.get("a1_c1_row:" + r["id"]) or {}).get("verdict") == "rejected":
                fail(r, "composition_double_division_suspect", "not_espei_ready",
                     "printed composition of %s is 1/100 of its formula fraction (a fraction under an at.%% header?); "
                     "not yet checked on the page" % ", ".join(sorted({str(p.get("phase")) for p in c1})))
                rep.add(GATE, "composition_double_division_suspect_row", r, "flagged")
            rev = [(pc, e, x) for pc, e, x, f in cells if f == "rev"]
            if not rev:
                continue
            rchk = rows_checked.get(r["id"])
            if rchk and rchk.get("verdict") == "confirmed":
                n = sum(_flip(r, pc, rchk, "composition_axis_from_page", {"row": "printed on the other element's axis"})
                        for pc, e, x in printed_cells(r) if len(pc.get("composition") or {}) == 1)
                if n:
                    r.setdefault("flags", {})["axis_header_contradicted"] = True
                    rep.add(GATE, "c2_axis_flipped_row", r, "corrected", cells=n)
                continue
            if rchk and rchk.get("verdict") == "rejected":
                if rchk.get("label_conflict"):
                    fail(r, "label_conflicts_composition", "not_espei_ready",
                         "page check: %s" % (rchk.get("seen") or "")[:200], rchk)
                    rep.add(GATE, "label_conflicts_composition", r, "flagged", kind="page")
                continue
            terminal = [pc for pc, e, x in rev if pure_element(pc.get("phase"), r.get("elements"))]
            formula = [pc for pc, e, x in rev if not pure_element(pc.get("phase"), r.get("elements"))]
            if terminal:
                fail(r, "composition_axis_suspect", "not_espei_ready",
                     "terminal phase %s printed at the other end of the composition axis" %
                     ", ".join(sorted({str(p.get('phase')) for p in terminal})))
                rep.add(GATE, "composition_axis_suspect_row", r, "flagged")
            if formula:
                fail(r, "label_conflicts_composition", "not_espei_ready",
                     "printed composition of %s is the formula with its elements reversed (label misprint or misread)" %
                     ", ".join(sorted({str(p.get('phase')) for p in formula})))
                rep.add(GATE, "label_conflicts_composition", r, "flagged", kind="formula_reversed")


# ------------------------------------------------------------------------------------------------- C3

C3_FIT, C3_OFF = 0.03, 0.08


def _registry_ranges(data):
    """{(SYSTEM, TDB phase name): (x_element, x_min, x_max)} from the phase registry's calculated stability."""
    out = {}
    d = os.path.join(data or "", "phases", "systems")
    if not data or not os.path.isdir(d):
        return out
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".json"):
            continue
        try:
            js = json.load(open(os.path.join(d, fn), encoding="utf-8"))
        except (OSError, ValueError):
            continue
        for p in js.get("phases") or []:
            st = [s for s in p.get("stability") or [] if s.get("stable") and s.get("x_min") is not None]
            if not st:
                continue
            xe = st[0].get("x_element")
            lo = min(s["x_min"] for s in st if s.get("x_element") == xe)
            hi = max(s["x_max"] for s in st if s.get("x_element") == xe)
            for n in p.get("tdb_names") or []:
                out[(p["system"], n.upper())] = (xe, lo, hi)
    return out


def _range_dev(r, pc, e, v, ranges):
    """Distance of the printed value from the phase's registry range (for element e), or None without a range."""
    for p in r.get("phases") or []:
        if p.get("label") == pc.get("phase") and p.get("model_phase"):
            rg = ranges.get((r["system"], str(p["model_phase"]).upper()))
            if rg:
                xe, lo, hi = rg
                if xe != e:
                    lo, hi = 1 - hi, 1 - lo
                return max(lo - v, v - hi, 0.0)
    return None


def c3_permutation(records, visual, rep, data=None):
    ranges = _registry_ranges(data)
    for r in records:
        els = r.get("elements") or []
        if len(els) != 2 or any(g["code"] in ("composition_axis_suspect", "composition_double_division_suspect")
                                for g in r.get("gate_failures") or []):
            continue
        cells = [(pc, e, v) for pc, e, v in printed_cells(r) if len(pc.get("composition") or {}) == 1 and 0 <= v <= 1]
        if not cells:
            continue
        el = cells[0][1]
        if any(e != el for _, e, _ in cells):
            continue
        chk = visual.get("a1_c3:" + r["id"])
        if chk and chk.get("verdict") == "corrected" and chk.get("set"):
            # values read on the page cell by cell (the stored row was shifted, not a clean swap)
            for pc, e, v in cells:
                nv = (chk["set"].get(pc.get("phase")) or {}).get(e)
                if nv is not None and abs(nv - v) > 1e-12:
                    _set_comp(r, pc, {e: nv}, "composition_read_from_page", {"printed_row": chk.get("printed")}, chk)
            r.setdefault("flags", {})["composition_phase_permutation_corrected"] = True
            rep.add(GATE, "c3_composition_read_from_page", r, "corrected")
            continue
        exp = [expected_x(pc.get("phase"), el, els) for pc, _, _ in cells]
        stoich = [i for i, x in enumerate(exp) if x is not None and 0 < x < 1]
        if not stoich:
            continue
        devs = {i: abs(cells[i][2] - exp[i]) for i in stoich}

        def off(i, v):
            """A stoichiometric cell is off when it is > C3_OFF from its formula and outside the calculated range."""
            if abs(v - exp[i]) <= C3_OFF:
                return False
            rd = _range_dev(r, cells[i][0], el, v, ranges)
            return rd is None or rd > C3_FIT

        def fits(i, v):
            """Within C3_FIT of the formula, or inside the phase's calculated range (+ C3_FIT)."""
            if abs(v - exp[i]) <= C3_FIT:
                return True
            rd = _range_dev(r, cells[i][0], el, v, ranges)
            return rd is not None and rd <= C3_FIT

        offs = [i for i in stoich if off(i, cells[i][2])]
        if not offs:
            continue
        chk = visual.get("a1_c3:" + r["id"])
        vals = [v for _, _, v in cells]
        best = None
        if 2 <= len(cells) <= 6 and not (chk and chk.get("verdict") == "rejected"):
            # one transposition of two printed values that puts every phase of known composition within C3_FIT
            for i, j in itertools.combinations(range(len(cells)), 2):
                if abs(vals[i] - vals[j]) < 1e-9 or not ({i, j} & set(offs)):
                    continue
                perm = list(range(len(cells)))
                perm[i], perm[j] = j, i
                if all(exp[k] is None or fits(k, vals[perm[k]]) or
                       (k not in (i, j) and not off(k, vals[k])) for k in range(len(cells))):
                    best = tuple(perm)
                    break
        if best:
            moved = []
            for i, (pc, e, v) in enumerate(cells):
                nv = vals[best[i]]
                if abs(nv - v) > 1e-12:
                    moved.append((pc, e, v, nv))
            for pc, e, v, nv in moved:
                _set_comp(r, pc, {e: nv}, "composition_phase_permutation",
                          {"rule": "values re-assigned to the phases by stoichiometry within the printed row",
                           "stored_order": [c[0].get("phase") for c in cells], "values": vals}, chk)
            r.setdefault("flags", {})["composition_phase_permutation_corrected"] = True
            rep.add(GATE, "c3_composition_phase_permutation", r, "corrected",
                    phases=[c[0].get("phase") for c in cells])
            continue
        # no permutation explains it: does the registry range (homogeneity range) explain the printed value?
        bad = []
        for i in stoich:
            if devs[i] <= C3_OFF:
                continue
            pc, e, v = cells[i]
            rd = _range_dev(r, pc, e, v, ranges)
            if rd is None:
                if devs[i] > 0.15:
                    bad.append((pc.get("phase"), v, exp[i], None))
            elif rd > C3_FIT:
                bad.append((pc.get("phase"), v, exp[i], rd))
        if chk and chk.get("verdict") == "rejected" and chk.get("label_conflict"):
            fail(r, "label_conflicts_composition", "not_espei_ready", "page check: %s" % (chk.get("seen") or "")[:200], chk)
            rep.add(GATE, "label_conflicts_composition", r, "flagged", kind="page")
        elif bad and not (chk and chk.get("verdict") == "rejected"):
            fail(r, "label_conflicts_composition", "not_espei_ready",
                 "; ".join("%s printed x(%s)=%.4g, formula %.4g%s" % (p, el, v, x, "" if rd is None else
                                                                       ", %.3g outside the calculated range" % rd)
                           for p, v, x, rd in bad), chk)
            rep.add(GATE, "label_conflicts_composition", r, "flagged")


# ------------------------------------------------------------------------------------------------- C4

INTRINSIC = re.compile(r"^(?:HM_FORM|SM_FORM|GM_FORM|LATTICE_.*|CPM|SM|ENTROPY|HEAT_CONTENT|MOLAR_VOLUME|VOLUME|"
                       r"BULK_MODULUS|SITE_FRACTIONAL_COORDINATES|TRANSITION_ENTHALPY|SM_FUSION|SM_INCREMENT|SM_MAGNETIC|"
                       r"HM_LATTICE_STABILITY|GIBBS_ENERGY|VAPOR_PRESSURE)$")
EQUILIBRIUM = {"INVARIANT_TEMPERATURE", "TRANSITION_TEMPERATURE", "PHASE_COMPOSITION", "PHASE_ASSEMBLAGE"}
_LIQ_GAS = re.compile(r"^\s*(?:l|liq\.?|liquid|melt|gas|g|vapou?r)\s*\d*\s*$", re.I)


def _phase_support(r, label, elements, allow_solution=False):
    """Element set of the state of one phase: its printed composition, else its formula / pure element; None when it
    cannot be determined (a solution label without a composition)."""
    f = formula_fractions(label, elements)
    p = pure_element(label, elements)
    lab = set(f) if f else ({p} if p and (allow_solution or not str(label).strip().startswith("(")) else None)
    for pc in r.get("phase_compositions") or []:
        if pc.get("phase") == label and str(pc.get("basis", "")).startswith("mole_fraction"):
            comp = {e: v for e, v in (pc.get("composition") or {}).items() if isinstance(v, (int, float))}
            if comp and abs(sum(comp.values()) - 1) <= 0.01:
                return {e for e, v in comp.items() if v > 1e-9}
            if comp and lab:
                # a partial composition: the balance is the label's own element(s)
                return lab | {e for e, v in comp.items() if v > 1e-9}
            return None   # a partial composition of a solution label: the balance is not determined
    return lab


def _state_support_r2(r):
    """Round-2 rule (kept as the first stage of state_support): (support, basis) or (None, reason)."""
    els = list(r.get("elements") or [])
    if len(els) < 2:
        return None, "unary"
    prop = r.get("property") or ""
    # round-3 relabels of compound-intrinsic quantities stay compound-intrinsic for the state system
    prop = "HM_FORM" if prop == "DEFECT_FORMATION_ENERGY" else re.sub(r"_(?:CONTRIBUTION|FROM_BASE_T)$", "", prop)
    phases = [p.get("label") for p in r.get("phases") or [] if p.get("label")]
    oc = {e: v for e, v in (r.get("overall_composition") or {}).items() if isinstance(v, (int, float))}
    if any(e not in els for e in oc):
        return None, "overall composition names elements outside the system"
    ocb = r.get("overall_composition_basis")
    if oc and ocb != "from_formula":
        bal = 1 - sum(oc.values())
        sup = {e for e, v in oc.items() if v > 1e-9}
        missing = [e for e in els if e not in oc]
        if bal > 1e-3:
            if len(missing) != 1:
                return None, "balance element not determined"
            sup.add(missing[0])
        if sup and len(sup) < len(els):
            if INTRINSIC.match(prop) or prop in ("HM_MIX", "GM_MIX", "SM_MIX"):
                return sup, "printed_overall_composition"
        return None, "printed overall composition spans the system"
    if oc and abs(sum(oc.values()) - 1) <= 0.01 and (
            INTRINSIC.match(prop) or (prop == "TRANSITION_TEMPERATURE" and len(phases) <= 2)):
        # every element listed with its fraction (zeros explicit): the state's own composition
        sup = {e for e, v in oc.items() if v > 1e-9}
        if sup and len(sup) < len(els):
            return sup, "overall_composition_all_elements"
    if INTRINSIC.match(prop):
        if not phases:
            return None, "no phase label"
        sups = [_phase_support(r, p, els) for p in phases]
        if any(s is None for s in sups):
            return None, "phase state not determined"
        sup = set().union(*sups)
        if oc:
            ocs = {e for e, v in oc.items() if v > 1e-9}
            if not ocs <= sup and not sup <= ocs:
                return None, "formula and overall composition disagree"
            sup |= ocs
        if sup and len(sup) < len(els):
            return sup, "phase_formula"
        return None, "spans the system"
    if prop in EQUILIBRIUM:
        if not phases or any(_LIQ_GAS.match(p) for p in phases if not any(
                pc.get("phase") == p for pc in r.get("phase_compositions") or [])):
            return None, "liquid / gas without a composition"
        # a terminal solution '(El)' is the element's own phase when every other phase of an invariant / transition
        # is fixed as well: three (or fewer) phases of one binary make a binary invariant by the phase rule
        sups = [_phase_support(r, p, els, allow_solution=prop in ("INVARIANT_TEMPERATURE", "TRANSITION_TEMPERATURE"))
                for p in phases]
        if any(s is None for s in sups):
            return None, "phase state not determined"
        sup = set().union(*sups)
        if sup and len(sup) < len(els):
            if prop in ("PHASE_COMPOSITION", "PHASE_ASSEMBLAGE"):
                # a phase observed in a multi-component alloy is a state of that alloy unless the table names the
                # subsystem itself
                import literature_gates as _g
                loc = r.get("locator") or {}
                caps = _g.caption_systems((loc.get("caption") or "") + " " + (r.get("quote") or ""))
                if not any(set(c) == sup for c in caps):
                    return None, "multi-component alloy (caption does not name the subsystem)"
                return sup, "every_phase_determined_and_caption"
            return sup, "every_phase_determined"
    return None, "not a single-state property"



# ------------------------------------------------------------------------------------------------- C4, round 3

# audit round 3 (POSTFIX A1-N1 / A3-N1): the round-2 rule trusted an overall_composition that was already corrupt (one
# element = 1.0, the others 0: Nd3Ru stored as {ND: 1, RU: 0}), so compound data were filed under a unary system. The
# printed phase formula now bounds the support from below, a generic rare-earth label (Ln3Ru, Cu5R, RMg2, Mn23RE6) is
# read with the element the row names, and a unary system needs a pure-element label.
GENERIC = re.compile(r"(?<![A-Z])(Ln|RE|R|M|Me|X|T)(?=[A-Z0-9₀-₉(\s]|$)")
_STRUCT_ONLY = re.compile(r"^\s*(?:(?:fcc|bcc|hcp|dhcp|bct|diamond|liquid|liq\.?|l|gas|vapou?r)(?:[_\-\s]?[A-Z]\d{1,2})?|"
                          r"(?:FCC|BCC|HCP|BCT|DHCP)_A\d{1,2}|A\d{1,2})\s*$", re.I)
FORMATION_LIKE = re.compile(r"^(?:HM_FORM|SM_FORM|GM_FORM|.*_MIX|.*_MIX_PARTIAL|.*_FROM_COMPOUNDS|TRANSITION_ENTHALPY)$")
EXPRESSION_TERM = re.compile(r"_EXPRESSION_TERM$")
GAS_ELEMENTS = {"H", "N", "O", "F", "CL", "BR", "I"}


def generic_substituted(label, element, elements):
    """Formula fractions of a generic label (Ln3Ru, Cu5R, RMg2) with the placeholder read as `element`, or None."""
    if not label or not GENERIC.search(str(label)):
        return None
    s = GENERIC.sub(element.capitalize(), str(label), count=1)
    if GENERIC.search(s):
        return None
    return formula_fractions(s, elements)


def label_supports(r, els):
    """(F, generic, other): F = union of the element sets of the phase labels that parse to a formula of the system's
    elements; generic = labels with a placeholder element; other = labels that are neither a formula, a generic label,
    a pure element nor a structure / liquid name."""
    F, generic, other = set(), [], []
    for p in r.get("phases") or []:
        lab = p.get("label")
        if not lab:
            continue
        f = formula_fractions(lab, els)
        if f and len(f) >= 2:
            F |= set(f)
        elif GENERIC.search(str(lab)) and not pure_element(lab, els):
            generic.append(lab)
        elif not pure_element(lab, els) and not _STRUCT_ONLY.match(str(lab)) and not _LIQ_GAS.match(str(lab)):
            other.append(lab)
    return F, generic, other


def state_support(r):
    """(support, basis) of the record's own state when it is a strict subset of the record's elements, else (None,
    reason). Round 2 decides first; round 3 checks the decision against the printed phase labels."""
    sup, basis = _state_support_r2(r)
    els = list(r.get("elements") or [])
    if len(els) < 2:
        return sup, basis
    prop = r.get("property") or ""
    prop = "HM_FORM" if prop == "DEFECT_FORMATION_ENERGY" else re.sub(r"_(?:CONTRIBUTION|FROM_BASE_T)$", "", prop)
    F, generic, other = label_supports(r, els)
    oc = {e: v for e, v in (r.get("overall_composition") or {}).items() if isinstance(v, (int, float))}
    ocs = {e for e, v in oc.items() if v > 1e-9}
    if generic:
        # the row names the element of a generic label through the stored composition (one element non-zero)
        if len(generic) == 1 and len(r.get("phases") or []) == 1 and (INTRINSIC.match(prop) or FORMATION_LIKE.match(prop)):
            for e in sorted(ocs):
                g = generic_substituted(generic[0], e, els)
                if not (g and len(g) >= 2):
                    continue
                if len(ocs) == 1 or all(abs(g.get(k, 0.0) - oc.get(k, 0.0)) <= 0.02 for k in set(g) | ocs):
                    G = set(g)
                    return (G, "generic_label_with_row_element") if G < set(els) else (None, "generic label spans the system")
        return None, "generic label: element not determined"
    if sup is None:
        if basis == "balance element not determined" and F and ocs <= F and len(r.get("phases") or []) == 1:
            f = formula_fractions(r["phases"][0]["label"], els)
            if f and all(abs(f.get(e, 0) - v) <= 0.02 for e, v in oc.items()) and F < set(els):
                return F, "phase_formula_with_partial_overall_composition"
        if EXPRESSION_TERM.search(prop) and len(r.get("phases") or []) == 1:
            lab = r["phases"][0].get("label")
            f = formula_fractions(lab, els)
            pe = pure_element(lab, els)
            S = set(f) if f and len(f) >= 2 else ({pe} if pe and not str(lab).strip().startswith("(") else None)
            if S and S < set(els) and (not ocs or ocs <= S):
                return S, "expression_term_of_one_phase"
        return None, basis
    if F and not F <= sup:
        if ocs <= F:
            sup2 = sup | F
            return (sup2, basis + "+phase_formula") if sup2 < set(els) else (None, "phase formula spans the system")
        return None, "formula and overall composition disagree"
    if len(sup) == 1:
        if F:
            return None, "compound label: not a unary state"
        (e,) = sup
        if e in GAS_ELEMENTS and not all(pure_element(p.get("label"), els) for p in r.get("phases") or [] if p.get("label")):
            # 'UNKNOWN' / '014x24' (Bi-2212) / oxide reactions stored with O = 1: never a state of pure oxygen
            return None, "gas-element unary without a pure-element label"
    return sup, basis


def c4_state_system(records, rep):
    """File each record under the chemical system of its own state; `source_system` keeps the paper's system."""
    for r in records:
        r.setdefault("source_system", r["system"])
        sup, basis = state_support(r)
        if not sup:
            continue
        new_els = sorted(sup)
        new = "-".join(new_els)
        before = r["system"]
        correct(r, "system", before, new, "system_of_the_state",
                {"basis": basis, "phases": [p.get("label") for p in r.get("phases") or []][:6],
                 "overall_composition": r.get("overall_composition")})
        r["system"] = new
        r["elements"] = new_els
        r["flags"]["system_is_state_subsystem"] = True
        # compositions keyed on elements outside the state are zeros by construction; drop them
        oc = r.get("overall_composition")
        if isinstance(oc, dict):
            drop = {e: v for e, v in oc.items() if e not in sup}
            if drop and all(not v for v in drop.values() if isinstance(v, (int, float))):
                r.setdefault("overall_composition_stored", oc)
                r["overall_composition"] = {e: v for e, v in oc.items() if e in sup}
        for pc in r.get("phase_compositions") or []:
            comp = pc.get("composition") or {}
            drop = {e: v for e, v in comp.items() if e not in sup}
            if drop and all(not v for v in drop.values() if isinstance(v, (int, float))):
                pc.setdefault("composition_stored", comp)
                pc["composition"] = {e: v for e, v in comp.items() if e in sup}
        rep.add(GATE, "c4_system_of_the_state", r, "corrected", source_system=before, n_elements=len(new_els))
    for r in records:
        # reference-state keys of elements outside the record's (narrowed) system: kept in reference_state_stored
        # (audit round 3, POSTFIX A3-N8)
        rs = r.get("reference_state")
        els = set(r.get("elements") or r["system"].split("-"))
        if isinstance(rs, dict) and any(k not in els for k in rs):
            r.setdefault("reference_state_stored", rs)
            r["reference_state"] = {k: v for k, v in rs.items() if k in els}
            rep.add(GATE, "r3_reference_state_outside_system_dropped", r, "corrected", n=len(rs) - len(r["reference_state"]))


def r3_overall_composition(records, rep):
    """Audit round 3 (POSTFIX A1-N2, A3 section 1(c)): an overall_composition whose basis is `from_formula` must be the
    composition of the record's printed formula. In split table rows the composition of one column was copied to the
    others (Ni5Si2 at Si = 0.25, the Ni3Si value; CuO at O = 1/3; Dy3Ru at Ru = 0); for a one-phase compound-intrinsic
    record the composition is recomputed from the printed phase formula (a generic label Ln3Ru / Cu5R / RMg2 with the
    element the row names), the stored one kept in overall_composition_stored. A printed (not from-formula) composition
    more than 0.1 from the formula is flagged `overall_composition_conflicts_label` (not ESPEI-ready)."""
    for r in records:
        prop = r.get("property") or ""
        if not (INTRINSIC.match(prop) or EXPRESSION_TERM.search(prop)):
            continue
        phases = [p.get("label") for p in r.get("phases") or [] if p.get("label")]
        if len(phases) != 1:
            continue
        els = list(r.get("elements") or [])
        oc = {e: v for e, v in (r.get("overall_composition") or {}).items() if isinstance(v, (int, float))}
        if not oc:
            continue
        f = formula_fractions(phases[0], els)
        basis = "printed_formula"
        if not (f and len(f) >= 2):
            ocs = [e for e, v in oc.items() if v > 1e-9]
            f = generic_substituted(phases[0], ocs[0], els) if len(ocs) == 1 else None
            basis = "generic_label_with_row_element"
        if not (f and len(f) >= 2):
            continue
        # a partial composition lists some elements (the balance is the rest): compare the listed ones
        keys = set(oc) | (set(f) if abs(sum(oc.values()) - 1) <= 0.01 else set())
        dev = max(abs(oc.get(e, 0.0) - f.get(e, 0.0)) for e in keys)
        if dev <= 0.02:
            continue
        if r.get("overall_composition_basis") == "from_formula":
            new = {e: round(f.get(e, 0.0), 6) for e in sorted(set(oc) | set(f))}  # zeros of other elements kept
            correct(r, "overall_composition", oc, new, "overall_composition_from_printed_formula",
                    {"label": phases[0], "basis": basis, "max_deviation": round(dev, 4)})
            r.setdefault("overall_composition_stored", r.get("overall_composition"))
            r["overall_composition"] = new
            rep.add(GATE, "r3_overall_composition_from_formula", r, "corrected", label=phases[0], dev=round(dev, 4))
        elif dev > 0.1 and basis == "printed_formula":
            fail(r, "overall_composition_conflicts_label", "not_espei_ready",
                 "printed overall composition %s is %.3f from the formula of %s" % (oc, dev, phases[0]))
            rep.add(GATE, "r3_overall_composition_conflicts_label", r, "flagged", label=phases[0], dev=round(dev, 4))


# ------------------------------------------------------------------------------------------------- C6 / C7

def c6_elements(records, visual, rep):
    for r in records:
        chk = visual.get("a1_c6:" + r["id"])
        if chk and chk.get("verdict") == "corrected":
            for pc in r.get("phase_compositions") or []:
                if pc.get("phase") in (chk.get("set") or {}):
                    _set_comp(r, pc, dict(chk["set"][pc["phase"]]), "composition_read_from_page",
                              {"printed_row": chk.get("printed")}, chk)
            ren = chk.get("rename_elements") or {}
            for pc in r.get("phase_compositions") or []:
                comp = pc.get("composition") or {}
                if any(k in ren for k in comp):
                    _set_comp(r, pc, {ren.get(k, k): v for k, v in comp.items()}, "element_symbol_read_from_page",
                              {"rename": ren}, chk)
            oc = r.get("overall_composition")
            if isinstance(oc, dict) and any(k in ren for k in oc):
                new = {ren.get(k, k): v for k, v in oc.items()}
                correct(r, "overall_composition", oc, new, "element_symbol_read_from_page", {"rename": ren}, chk)
                r["overall_composition_stored"], r["overall_composition"] = oc, new
            if chk.get("system"):
                new_els = sorted(set(chk["system"].upper().split("-")))
                new = "-".join(new_els)
                if new != r["system"]:
                    correct(r, "system", r["system"], new, "system_read_from_page", {"seen": chk.get("seen")}, chk)
                    r["system_stored"], r["system"], r["elements"] = r["system"], new, new_els
                    r["flags"]["system_reassigned"] = True
            rep.add(GATE, "c6_page_fix", r, "corrected")
        els = set(r.get("elements") or [])
        oc = r.get("overall_composition")
        if isinstance(oc, dict) and oc:
            extra = {e: v for e, v in oc.items() if e not in els}
            if extra:
                zero = {e: v for e, v in extra.items() if not v}
                if zero:
                    new = {e: v for e, v in oc.items() if e not in zero}
                    correct(r, "overall_composition", oc, new, "stale_zero_outside_system",
                            {"dropped": sorted(zero)})
                    r["overall_composition_stored"], r["overall_composition"] = oc, new
                    rep.add(GATE, "c6_stale_zero_dropped", r, "corrected")
                nz = sorted(e for e in extra if extra[e])
                if nz:
                    fail(r, "composition_element_outside_system", "not_espei_ready",
                         "overall composition names %s, not in the system" % ", ".join(nz))
                    rep.add(GATE, "composition_element_outside_system", r, "flagged")
        for pc in r.get("phase_compositions") or []:
            comp = pc.get("composition") or {}
            extra = sorted(e for e in comp if e not in els)
            if extra:
                zero = [e for e in extra if not comp[e]]
                if zero and len(zero) == len(extra):
                    _set_comp(r, pc, {e: v for e, v in comp.items() if e not in zero}, "stale_zero_outside_system",
                              {"dropped": zero}, None)
                    rep.add(GATE, "c6_stale_zero_dropped", r, "corrected")
                else:
                    fail(r, "composition_element_outside_system", "not_espei_ready",
                         "composition of %s names %s, not in the system" % (pc.get("phase"), ", ".join(extra)))
                    rep.add(GATE, "composition_element_outside_system", r, "flagged")
                    break


def c7_ranges(records, visual, rep):
    for r in records:
        chk = visual.get("a1_c7:" + str(asset_of(r)))
        if chk and chk.get("verdict") == "corrected" and chk.get("factor"):
            # a table whose unresolved compositions the page shows to be percentages (header 'mol%' / 'at.%')
            n = 0
            for pc in r.get("phase_compositions") or []:
                if str(pc.get("basis", "")).startswith("unresolved"):
                    comp = pc.get("composition") or {}
                    new = {e: _rnd(v * chk["factor"]) for e, v in comp.items() if isinstance(v, (int, float))}
                    if _set_comp(r, pc, new, "composition_basis_read_from_page",
                                 {"printed_header": chk.get("printed"), "factor": chk["factor"]}, chk):
                        pc["basis_stored"] = pc.get("basis")
                        pc["basis"] = chk.get("basis") or "mole_fraction (converted from atomic percent)"
                        n += 1
            if n:
                rep.add(GATE, "c7_basis_resolved_on_page", r, "corrected", cells=n)
        bad, unresolved = [], False
        for pc in r.get("phase_compositions") or []:
            b = str(pc.get("basis", ""))
            if b == "unresolved" or b.startswith("unresolved"):
                unresolved = True
            if not (b.startswith("mole_fraction") or b.startswith("unresolved")):
                continue
            for e, v in (pc.get("composition") or {}).items():
                if isinstance(v, (int, float)) and not -1e-9 <= v <= 1 + 1e-9:
                    bad.append("%s %s=%g" % (pc.get("phase"), e, v))
        oc = r.get("overall_composition") or {}
        for e, v in oc.items():
            if isinstance(v, (int, float)) and not -1e-9 <= v <= 1 + 1e-9:
                bad.append("overall %s=%g" % (e, v))
        if bad:
            fail(r, "composition_outside_unit_interval", "not_accepted",
                 "stored as a mole fraction outside [0, 1]: %s (raw percentage or mass-to-mole conversion error)" %
                 "; ".join(bad[:4]))
            rep.add(GATE, "composition_outside_unit_interval", r, "demoted")
        elif unresolved:
            fail(r, "composition_basis_unresolved", "not_espei_ready",
                 "composition basis (mole / mass / at.%) not resolved on the page")
            rep.add(GATE, "composition_basis_unresolved", r, "flagged")


# ------------------------------------------------------------------------------------------------- hooks

def before_compositions(records, visual, rep, data=None):
    """C1 -> C2 -> C3 -> C6 -> C7 (in this order: the axis and permutation tests read the rescaled values)."""
    c1_double_division(records, visual, rep)
    c2_axis(records, visual, rep)
    c3_permutation(records, visual, rep, data)
    c2_rows(records, visual, rep)
    c6_elements(records, visual, rep)
    c7_ranges(records, visual, rep)


def after_system(records, rep):
    r3_overall_composition(records, rep)
    c4_state_system(records, rep)
