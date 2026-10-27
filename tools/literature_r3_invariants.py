#!/usr/bin/env python3
"""Audit round 3 (2026-10-24; PI: "fix invariant reaction issues"): invariant-reaction QA class family.

Called from tools/literature_gates.apply() after the round-3 thermo rules (literature_r3_thermo.apply_late, which
recovers or un-accepts states without phases). Every rule keeps the stored field (`phases_stored`,
`reaction.type_stored`) and writes a `corrections` entry or a gate failure (gate "R3I"). The printed reaction string
(`reaction.printed`, the table cell read from the page) is the evidence; page checks of a sample of each class are in
tools/literature_visual_checks.json (`r3_class:<rule>`), Cu-Zr record by record (`record:<id>`).

  R3I-a  phases missing from a table record whose printed reaction names them (`fcc (Ag) + fcc (Co) -> hcp (Co)` stored
         with one phase): the phases are restored from the printed reaction (label_source printed_reaction).
  R3I-b  reaction class missing: derived from the printed reaction (reaction.type_basis derived_from_printed_reaction);
         neither a class nor a printed reaction -> flag reaction_not_recorded.
  R3I-c  phase count against the phase rule: an invariant of a C-component system has C + 1 phases (binary eutectic,
         peritectic, eutectoid, peritectoid, monotectic, metatectic, syntectic: 3); a congruent melting / congruent or
         polymorphic transformation: 2. A printed reaction with another count -> flag phase_rule_count_mismatch
         (not ESPEI-ready).
  R3I-d  class consistent with the printed reaction: liquid count and side (eutectic L -> S1 + S2; peritectic
         L + S1 -> S2; eutectoid / peritectoid without liquid; monotectic L1 -> L2 + S; metatectic S1 -> L + S2;
         syntectic L1 + L2 -> S; congruent: one phase each side). A stored class that contradicts the printed
         reaction is replaced by the derived class (reaction.type_stored kept) when the reaction is unambiguous (a
         one-way arrow, or a class that does not depend on the direction); otherwise flagged.
  R3I-f  temperatures: a value equal to a printed °C number (no +273.15) or outside a broad plausibility window of the
         system (below 0.15 x the lowest melting point of its elements or above the highest + 1500 K) -> flag
         invariant_temperature_needs_page_check (not ESPEI-ready).
  R3I-g  the same invariant twice on one page (paper, page, system, |dT| <= 0.5 K, same printed phase set, both
         accepted): the later records are flagged duplicate_invariant_same_page (not ESPEI-ready; nothing deleted).
"""
from __future__ import annotations

import re
from collections import defaultdict

import literature_gates as G
import phase_rule

GATE = "R3I"
RELEASE = "literature-r3-invariants 2026-10-24"
INV_PROPS = {"INVARIANT_TEMPERATURE"}
ARROW = re.compile(r"\s*(<=>|<->|<-->|-->|=>|->|⇌|⇄|↔|⟷|→|⟶|⇒|=)\s*")
ONE_WAY = {"->", "-->", "=>", "→", "⟶", "⇒"}
THREE_PHASE = {"eutectic", "peritectic", "eutectoid", "peritectoid", "monotectic", "monotectoid", "metatectic/catatectic",
               "syntectic"}
TWO_PHASE = {"congruent melting", "congruent transformation", "polymorphic"}


def norm_label(s):
    return re.sub(r"\s+", " ", str(s or "")).strip().strip(".").strip()


def parse_reaction(printed):
    """(left phases, right phases, arrow) of a printed reaction 'A + B -> C', or None when it is not of that form."""
    s = str(printed or "").replace("þ", "+").strip()
    if not s or "/" in s or "?" in s or re.search(r"^[A-Za-z]\d*\s*:", s) or "," in s:
        return None   # alternatives, garbled signs, invariant labels ('p1: L + ...') and comma lists are not parsed
    parts = ARROW.split(s)
    if len(parts) != 3:
        return None
    left, arrow, right = parts
    left = re.sub(r"(?i)^\s*(?:eutectic|peritectic|eutectoid|peritectoid|monotectic|congruent(?: melting)?)\s*:?\s+", "", left)
    L = [norm_label(x) for x in re.split(r"\s*\+\s*", left)]
    R = [norm_label(x) for x in re.split(r"\s*\+\s*", right)]
    if not all(L) or not all(R) or any(len(x) > 40 for x in L + R):
        return None
    return L, R, arrow


LIQ = re.compile(r"^(?:L|l|Liq|liq|Liquid|liquid|LIQUID|melt|Melt)\.?\s*(?:#?\d|[′'])?\*?$|\((?:L|l|liq|liquid)\)$|^Liquid\s*#\s*\d$")


def is_liquid(label):
    return bool(LIQ.search(norm_label(label))) or phase_rule.is_liquid_label(label, None)


def _same(a, b):
    a, b = re.sub(r"\s+", "", a).lower(), re.sub(r"\s+", "", b).lower()
    return a == b or (len(a) >= 1 and len(b) >= 1 and (a in b or b in a)) or (is_liquid(a) and is_liquid(b))


BAD_LABEL = re.compile(r"[<>:;]|^[^\w(\[]|\($|^\)|^x\b|^\d")


def derive_class(L, R, arrow):
    """Reaction class from the printed reaction, or None when the printed form does not decide it."""
    nl = [is_liquid(x) for x in L]
    nr = [is_liquid(x) for x in R]
    if len(L) == 1 and len(R) == 1:
        if nl[0] != nr[0]:
            return "congruent melting"
        return "polymorphic" if not nl[0] else None
    one_way = arrow in ONE_WAY
    if len(L) == 1 or len(R) == 1:
        single, single_liq = (L, nl[0]) if len(L) == 1 else (R, nr[0])
        other_liq = sum(nr) if len(L) == 1 else sum(nl)
        other_n = len(R) if len(L) == 1 else len(L)
        if other_n != 2:
            return None
        if single_liq and other_liq == 0:
            return "eutectic"
        if single_liq and other_liq == 1:
            return "monotectic"
        if not single_liq and other_liq == 2:
            return "syntectic"
        if not single_liq and other_liq == 1:
            # peritectic L + S1 = S2 or metatectic S1 = L + S2: the two differ only in which side is the high-T side,
            # and tables write the arrow in either the cooling or the heating direction; not decided here
            return "peritectic|metatectic/catatectic"
        if not single_liq and other_liq == 0:
            return "eutectoid|peritectoid"
    return None


def compatible(stored, derived):
    if not stored or not derived or stored in ("other", "invariant", "unknown"):
        return True
    if "|" in derived:
        return stored in derived.split("|") or stored in ("monotectoid", "quasi-peritectic") or \
            (derived.startswith("peritectic") and stored in ("peritectic", "metatectic/catatectic", "incongruent melting"))
    if stored == derived:
        return True
    pairs = {frozenset(("congruent melting", "congruent transformation")), frozenset(("polymorphic", "congruent transformation")),
             frozenset(("monotectic", "monotectoid")), frozenset(("eutectoid", "monotectoid")),
             frozenset(("order-disorder", "polymorphic")), frozenset(("magnetic", "polymorphic")),
             frozenset(("critical point", "polymorphic"))}
    return frozenset((stored, derived)) in pairs


def _labels(r):
    return [norm_label(p.get("label")) for p in r.get("phases") or [] if p.get("label")]


def r3i_reactions(records, visual, rep):
    cls = lambda k: visual.get("r3_class:" + k)  # noqa: E731
    for r in records:
        if r.get("property") not in INV_PROPS:
            continue
        rx = r.get("reaction") if isinstance(r.get("reaction"), dict) else {}
        printed = rx.get("printed")
        parsed = parse_reaction(printed)
        stored_type = rx.get("type")
        els = r.get("elements") or (r.get("system") or "").split("-")
        if not parsed:
            if not stored_type and not printed and not r.get("phases"):
                G.fail(r, GATE, "reaction_not_recorded", "flag_only", "neither a reaction class nor a printed reaction")
                rep.add(GATE, "reaction_not_recorded", r, "flagged")
            n = len(_labels(r))
            if r["flags"].get("extraction_accepted") and n and ((stored_type in THREE_PHASE and n < 3) or (n == 1 and stored_type not in ("critical point", "magnetic", "order-disorder"))):
                # (c) without a printed reaction: a three-phase class with fewer phases, or one phase only
                G.fail(r, GATE, "phase_rule_count_mismatch", "not_espei_ready",
                       "%s with %d phase(s) recorded" % (stored_type or "invariant", n))
                rep.add(GATE, "phase_rule_count_mismatch", r, "flagged", n=n)
            continue
        L, R, arrow = parsed
        printed_set = []
        for x in L + R:
            if x not in printed_set:
                printed_set.append(x)
        # (a) phases named by the printed reaction but missing from the record
        have = _labels(r)
        missing = [x for x in printed_set if not any(_same(x, h) for h in have)]
        if missing and len(have) < len(printed_set) and r.get("record_kind") != "running_text" \
                and not any(BAD_LABEL.search(x) or x.count("(") != x.count(")") for x in missing):
            r.setdefault("phases_stored", [dict(p) for p in r.get("phases") or []])
            r["phases"] = list(r.get("phases") or []) + [{"label": x, "label_source": "printed_reaction"} for x in missing]
            G.correct(r, GATE, "phases", have, _labels(r), "phases_from_printed_reaction",
                      {"printed": printed, "added": missing}, cls("phases_from_printed_reaction"))
            G.fail(r, GATE, "phases_added_from_printed_reaction", "not_espei_ready",
                   "phases %s read from the printed reaction; not mapped to model phases" % missing,
                   cls("phases_from_printed_reaction"))
            rep.add(GATE, "phases_from_printed_reaction", r, "corrected", n=len(missing))
        # (b)/(d) class: only a binary reaction written in a form that decides it fills a missing class; a printed /
        # stored class is never replaced (tables print saddle, maximum, monotectoid, degenerate ... that the reaction
        # string does not show; re-sample 2026-10-24)
        derived = derive_class(L, R, arrow) if len(els) == 2 else None
        decisive = derived in ("congruent melting", "polymorphic") or (
            derived in ("eutectic", "monotectic") and arrow in ONE_WAY and len(L) == 1 and is_liquid(L[0]))
        if not stored_type and decisive:
            rx["type"] = derived
            rx["type_basis"] = "derived_from_printed_reaction"
            r["reaction"] = rx
            rep.add(GATE, "reaction_class_from_printed_reaction", r, "corrected", type=derived)
        elif stored_type and decisive and not compatible(stored_type, derived):
            G.fail(r, GATE, "reaction_class_conflicts_printed_reaction", "not_espei_ready",
                   "class %r but the printed reaction %r reads %r" % (stored_type, printed, derived))
            rep.add(GATE, "reaction_class_conflicts_printed_reaction", r, "flagged")
        # (c) phase count
        t = rx.get("type") or derived
        n = len(printed_set)
        C = len(els)
        need = C + 1 if t in THREE_PHASE or t == "quasi-peritectic" else (2 if t in TWO_PHASE else None)
        if t in ("congruent melting", "polymorphic", "congruent transformation") and C >= 2 and n == 2:
            need = 2
        if need and n != need and not (t in THREE_PHASE and C > 2 and n == 3):
            G.fail(r, GATE, "phase_rule_count_mismatch", "not_espei_ready",
                   "%s with %d phases in a %d-component system (%r)" % (t, n, C, printed))
            rep.add(GATE, "phase_rule_count_mismatch", r, "flagged", n=n, need=need)


TM = None


def _tm():
    global TM
    if TM is None:
        from espei_constraints import ELEMENT_TM_TB
        TM = ELEMENT_TM_TB
    return TM


C_NUM = re.compile(r"(-?\d+(?:\.\d+)?)\s*(?:°|º|◦|o)\s*C\b")


def r3i_temperatures(records, visual, rep):
    tm = _tm()
    for r in records:
        if r.get("property") not in INV_PROPS or not r["flags"].get("extraction_accepted"):
            continue
        v = r.get("value")
        if not isinstance(v, (int, float)):
            continue
        why = None
        tp = str(r.get("temperature_printed") or "")
        m = C_NUM.search(tp)
        if m and abs(float(m.group(1)) - v) < 1e-6:
            why = "value equals the printed °C number (no conversion to K)"
        els = [e for e in r.get("elements") or [] if e in tm]
        if not why and els:
            lo = min(tm[e][0] for e in els)
            hi = max(tm[e][0] for e in els)
            if v < 0.15 * lo or v > max(hi + 2500, 4500):
                why = "outside %g-%g K (0.15 x lowest melting point of the elements; highest + 2500 K, at least 4500 K)" % (
                    round(0.15 * lo), round(max(hi + 2500, 4500)))
        if why:
            G.fail(r, GATE, "invariant_temperature_needs_page_check", "not_espei_ready", why)
            rep.add(GATE, "invariant_temperature_needs_page_check", r, "flagged")


def r3i_duplicates(records, rep):
    groups = defaultdict(list)
    for r in records:
        if r.get("property") not in INV_PROPS or not r["flags"].get("extraction_accepted"):
            continue
        v = r.get("value")
        if not isinstance(v, (int, float)):
            continue
        loc = r.get("locator") or {}
        ph = tuple(sorted(re.sub(r"\s+", "", x).lower() for x in _labels(r)))
        src = tuple(sorted(str(c.get("doi") or c.get("label") or c.get("reference") or "")
                           for c in r.get("cited_sources") or []))
        groups[(r.get("paper"), loc.get("page") or r.get("quote_page"), r.get("system"), ph, r.get("evidence_class"),
                src, r.get("method"))].append(r)
    for k, rs in groups.items():
        if len(rs) < 2 or not k[3]:
            continue
        rs.sort(key=lambda r: (r.get("record_kind") != "invariant_table", r["id"]))
        kept = []
        for r in rs:
            dup = next((q for q in kept if abs(q["value"] - r["value"]) <= 0.5), None)
            if dup is None:
                kept.append(r)
                continue
            G.fail(r, GATE, "duplicate_invariant_same_page", "not_espei_ready",
                   "same invariant as %s on the same page (|dT| <= 0.5 K, same phases)" % dup["id"])
            r["duplicate_of"] = dup["id"]
            rep.add(GATE, "duplicate_invariant_same_page", r, "flagged")


def apply(records, visual, rep):
    r3i_reactions(records, visual, rep)
    r3i_temperatures(records, visual, rep)
    r3i_duplicates(records, rep)
