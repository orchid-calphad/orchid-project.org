"""Deterministic data hardening of the ORCHID literature records (called by tools/export_literature_data.py).

Every rule here only re-labels, normalizes, merges or flags what the records already carry; no value is invented.
Wherever a value or label changes, the original is kept beside it (`*_printed`, `property_code_raw`,
`value_before_quality_fix`, `merged_ids`) and the change is listed in QUALITY_FIXES.md / quality_fixes.json.

  1. duplicates            identical values under several record ids -> one record + merged_ids / merged_pids
  2. controlled vocabulary  reaction / transition type -> transition_kind (raw: transition_kind_printed);
                            malformed property codes -> clean code (raw: property_code_raw)
  3. form                   garbled quote_page -> null + flags.page_unreadable; temperature_bound always a dict;
                            reference_state always a dict ({"printed": ...} when it was free text)
  4. DFT parity outliers    basis_corrected_from_formula / sign_from_heading / unsigned_in_print /
                            dft_outlier_unexplained (rules in classify_dft_outliers)
  5. energy units           formation, mixing and partial enthalpies and Gibbs energies of mixing in kJ/mol-atom
                            only when the per-atom basis is explicit, "kJ/mol (basis unstated)" otherwise

Standard library only (the DFT snapshot reader is shared with tools/build_dft_parity.py). No network.
"""
from __future__ import annotations

import importlib.util
import json
import os
import re
from collections import Counter, defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------------------------------ vocabularies

TRANSITION_KINDS = (
    "eutectic", "peritectic", "monotectic", "syntectic", "eutectoid", "peritectoid", "monotectoid",
    "quasi-peritectic", "congruent melting", "congruent transformation", "polymorphic", "magnetic",
    "order-disorder", "critical point", "metatectic/catatectic", "other")

# ordered (first match wins); applied to the lower-cased printed type with typos repaired
_KIND_RULES = [
    ("magnetic", r"magnet|curie|n[ée]el"),
    ("order-disorder", r"order[\s\-–‐]*disorder|\bordering|disordering"),
    ("critical point", r"critical|consolute"),
    ("metatectic/catatectic", r"metatect|catatect"),
    ("monotectoid", r"monotectoid"),
    ("eutectoid", r"eutectoid"),
    ("peritectoid", r"peritectoid"),
    ("monotectic", r"monotectic"),
    ("syntectic", r"syntectic"),
    ("quasi-peritectic", r"quasi[- ]?peritectic|transition reaction|\bclass ii\b|^\s*ii\s*$|^\s*u\d*\s*$"),
    ("eutectic", r"eutectic"),
    ("peritectic", r"peritectic"),
    ("congruent transformation", r"congruent transformation|congruent (solid|polymorph)"),
    ("congruent melting", r"congruent(ly)? melt|^c\.\s*m\.?$|^\s*c\.?m\.?\s*$|melting point|^\s*melting\s*$|"
                          r"^\s*m\.?p\.?\s*$|^\s*fusion\s*$"),
    ("polymorphic", r"polymorph|allotrop|solid[- ]state transition|structural transition|"
                    r"^\s*allotropy\s*$"),
]
# Rhines-type single-letter labels of ternary/binary invariants (E/e eutectic, P/p peritectic, U transition)
_LETTER = {"e": "eutectic", "p": "peritectic", "u": "quasi-peritectic", "i": "eutectic", "iii": "peritectic"}
_TYPO = [(r"peritet?ic", "peritectic"), (r"perite[ct]{1,3}i[co]d", "peritectoid"), (r"peri[ct]{1,3}toid", "peritectoid"),
         (r"peritictoid|perictectoid|peritectiod", "peritectoid"), (r"peritecitc|peritecic", "peritectic"),
         (r"eutectiod|eutecoid", "eutectoid"), (r"eutecic|eutetic", "eutectic")]
_LIQ = re.compile(r"^\s*(l|liq|liquid|melt|l\d|l['′’]+|l\s*\d)\b", re.I)


def _kind(raw, phases):
    s = (raw or "").strip()
    if not s:
        return None
    low = s.lower()
    for pat, rep in _TYPO:
        low = re.sub(pat, rep, low)
    letters = low.strip().rstrip(".").strip()
    if letters in _LETTER and len(s.strip().rstrip(".")) <= 3:
        return _LETTER[letters], "letter_label"
    hits = []
    for term, rx in _KIND_RULES:
        if re.search(rx, low):
            hits.append(term)
    # a printed alternative ("eutectic (or peritectic)") is ambiguous; a degenerate eutectic is still eutectic
    distinct = [h for h in dict.fromkeys(hits)]
    if distinct and re.search(r"\bor\b|\?", low):
        return "other", "ambiguous_printed_alternative"
    if distinct:
        return distinct[0], "keyword"
    if re.fullmatch(r"\s*(congruent|congruent (point|reaction|maximum|minimum)|c|max(imum)?|congruently)\s*\.?", low):
        if low.strip(" .") in ("c", "max", "maximum"):
            return "other", "unqualified_letter_or_extremum"
        liquid = any(_LIQ.match(p or "") for p in phases)
        return ("congruent melting" if liquid else "congruent transformation"), (
            "congruent_with_liquid_phase" if liquid else "congruent_without_liquid_phase")
    return "other", "no_rule"


PROPERTY_FIX = {
    "ACR_": "ACTIVITY", "LATTICE_A_": "LATTICE_A", "LATTICE_B_": "LATTICE_B", "LATTICE_C_": "LATTICE_C",
    "LATTICE_LATTICE_PARAMETER": "LATTICE_UNSPECIFIED_AXIS", "LATTICE_LATTICE_PARAMETERS": "LATTICE_UNSPECIFIED_AXIS",
    "LATTICE_LATTICE_PARAMETERS_ERROR_5_": "LATTICE_UNSPECIFIED_AXIS",
    "LATTICE_LATTICE_PARAMETER_A": "LATTICE_A", "LATTICE_LATTICE_PARAMETER_B": "LATTICE_B",
    "LATTICE_LATTICE_PARAMETER_C": "LATTICE_C", "LATTICE_LATTICE_PARAMETER_ALPHA": "LATTICE_ALPHA",
    "LATTICE_LATTICE_PARAMETER_BETA": "LATTICE_BETA", "LATTICE_LATTICE_PARAMETER_GAMMA": "LATTICE_GAMMA",
    "LATTICE_A_FE_": "LATTICE_A", "LATTICE_A_NI_": "LATTICE_A", "LATTICE_A_FE50NI50_": "LATTICE_A",
    "LATTICE_VAT": "LATTICE_V_ATOM", "LATTICE_V_AT": "LATTICE_V_ATOM",
    "LATTICE_4F_Z_1_": "LATTICE_4F_Z_1", "LATTICE_4F_Z_2_": "LATTICE_4F_Z_2",
    "LATTICE_4F_Z_FIRST_": "LATTICE_4F_Z_1", "LATTICE_4F_Z_SECOND_": "LATTICE_4F_Z_2",
    "LATTICE_COMMENTS_UNLABELLED_VALUE_": "LATTICE_UNLABELLED_VALUE",
}


def clean_property_code(code, unit):
    if not code:
        return code
    if code == "LATTICE__":
        return "LATTICE_UNSPECIFIED_ANGLE" if (unit or "").strip() in ("°", "deg", "degree", "degrees") \
            else "LATTICE_UNSPECIFIED_AXIS"
    if code in PROPERTY_FIX:
        return PROPERTY_FIX[code]
    new = re.sub(r"^LATTICE_LATTICE_PARAMETERS?_?", "LATTICE_", code)
    new = re.sub(r"_+$", "", re.sub(r"__+", "_", new))
    return new if new and new != "LATTICE" else "LATTICE_UNSPECIFIED_AXIS"


# ------------------------------------------------------------------------------------------ energy units

ENERGY_KJ = {"HM_FORM", "GM_FORM", "HM_MIX", "GM_MIX", "PARTIAL_ENTHALPY"}
NEW_IN_THIS_RELEASE = {"HM_MIX", "GM_MIX", "PARTIAL_ENTHALPY"}
_E = [(re.compile(r"^kcal", re.I), 4.184), (re.compile(r"^cal", re.I), 0.004184), (re.compile(r"^kj", re.I), 1.0),
      (re.compile(r"^j(\b|[\s/(.·]|$)", re.I), 0.001), (re.compile(r"^ev", re.I), 96.48533212)]
PER_ATOM = re.compile(
    r"mol(e)?s?[\s.\-·'’()−‐⁻^]*(of[\s.\-]*)?at(om|oms)?\b(?!\s*\d)|mol(e)?[\s.\-·]*at\b(?!\s*\d)|"
    r"mol(e)?\s*\^?\{?[−\-‐–⁻'’]?\s*1?\}?[\s.(]*(of\s*)?at(om)?s?\b(?!\s*\d)|"
    r"\bat(om)?s?\.?\s*\^?\{?[−\-‐–⁻]\s*1|\bat\.?[\s\-]*g\b|\bat\.g\b|g(ram)?[\s.\-Æ·]*at(om)?\b|"
    r"\bat[\s.\-−]*mol|molatom|/\s*atom|per\s+(mole?\s+(of\s+)?)?atoms?|mol-atom", re.I)
PER_FORMULA = re.compile(r"formula|f\.\s*u\.|\bf\.?u\b|of phase|of compound|per mole of compound|/mol\s*f\b", re.I)
UNSTATED = re.compile(r"basis unstated|^k?j\s*(·|\.|\s)*(/\s*mol(e)?\.?|mol(e)?\s*[\^]?\{?[−\-‐–⁻'’]?\s*1?\}?\.?|"
                      r"mol(e)?[−\-‐–⁻]\s*1)\s*$", re.I)
EXPLICIT_BASIS = {"printed_unit", "column_header", "caption", "reviewer", "state_completion_caption",
                  "state_completion_column_header", "state_completion_text", "state_completion_table_footnote",
                  "eV_per_atom_printed"}


def unit_basis_from_conversions(conversions):
    """How the builder decided the per-atom / per-formula basis (from its conversion log)."""
    for c in conversions or []:
        head, _, tail = c.partition(":")
        if head == "unit_basis_from_column_header":
            return "column_header"
        if head == "unit_basis_from_caption":
            return "caption"
        if head == "unit_basis_from_reviewer":
            return "reviewer"
        if head == "unit_basis_from_state_completion":
            return "state_completion_" + (tail or "unspecified")
        if head == "unit_basis_from_magnitude_heuristic":
            return "magnitude_heuristic"
        if head == "unit_basis_assumed_per_mole_of_atoms":
            return "assumed_per_mole_of_atoms"
        if head == "unit_basis_from_not_printed":
            return "not_printed"
        if head == "eV_per_atom_to_J_per_mol_atom":
            return "eV_per_atom_printed"
        if head in ("basis_unstated", "kJ_to_J_basis_unstated"):
            return "unstated"
    return None


def _sig(x):
    return float(f"{x:.10g}")


def normalize_energy_unit(r, fixes):
    """kJ/mol-atom when per-atom is explicit (printed unit or a printed header/caption/footnote/reviewer
    decision); 'kJ/mol (basis unstated)' when the printed unit has no basis. Formation records keep the
    builder's magnitude-heuristic per-atom decisions (flagged); for mixing / partial quantities, newly converted
    in this release, a heuristic or assumed basis counts as unstated."""
    prop = r.get("property")
    if prop not in ENERGY_KJ or isinstance(r.get("value"), bool) or not isinstance(r.get("value"), (int, float)):
        return
    u = (r.get("unit") or "").strip()
    u0 = u.lstrip("( ").strip()
    fac = next((f for pat, f in _E if pat.search(u0)), None)
    if fac is None:
        return
    printed = r.get("unit_printed") or ""
    basis = r.get("unit_basis")
    if basis is None and r.get("record_kind") == "running_text":
        basis = "printed_unit" if PER_ATOM.search(u) else ("unstated" if UNSTATED.search(u) else None)
        if basis:
            r["unit_basis"] = basis
    if basis is None and PER_ATOM.search(printed):
        basis = r["unit_basis"] = "printed_unit"
    if fac == 96.48533212 and not PER_ATOM.search(u):
        return  # eV per cell / defect / formula: not a per-atom molar value
    per_atom = bool(PER_ATOM.search(u))
    if per_atom and prop in NEW_IN_THIS_RELEASE and basis not in EXPLICIT_BASIS and not PER_ATOM.search(printed):
        per_atom = False
        fixes["mixing_basis_not_explicit_demoted"] += 1
    if per_atom:
        new = "kJ/mol-atom"
    elif UNSTATED.search(u) or PER_ATOM.search(u) or re.fullmatch(r"\(?k?j\)?", u0, re.I):
        new = "kJ/mol (basis unstated)"
    else:
        return
    if r.get("unit_stored") is None:
        r["value_stored"], r["unit_stored"] = r["value"], u
    r["value"], r["unit"] = _sig(r["value"] * fac), new
    unc = r.get("uncertainty")
    if isinstance(unc, dict) and isinstance(unc.get("value"), (int, float)):
        unc["value"] = _sig(unc["value"] * fac)
    rng = r.get("value_range")
    if isinstance(rng, dict):
        r["value_range"] = {k: (_sig(v * fac) if isinstance(v, (int, float)) else v) for k, v in rng.items()}
    if prop in NEW_IN_THIS_RELEASE:
        fixes["mixing_converted_" + ("per_atom" if per_atom else "basis_unstated")] += 1
    if prop in ("HM_FORM", "GM_FORM") and per_atom and basis == "magnitude_heuristic":
        r.setdefault("flags", {})["basis_inferred_from_magnitude"] = True


# ------------------------------------------------------------------------------------------ form

MAX_PAGE = 2000


def fix_form(r, log):
    fl = r.setdefault("flags", {})
    qp = r.get("quote_page")
    if qp is not None and (not isinstance(qp, int) or isinstance(qp, bool) or qp < 1 or qp > MAX_PAGE):
        log.append({"id": r["id"], "task": "page", "field": "quote_page", "before": str(qp), "after": None,
                    "evidence": "not a page number (> %d or not a positive integer)" % MAX_PAGE})
        r["quote_page"] = None
        fl["page_unreadable"] = True
    loc = r.get("locator") or {}
    lp = loc.get("page")
    if lp is not None and (not isinstance(lp, int) or isinstance(lp, bool) or lp < 1 or lp > MAX_PAGE):
        log.append({"id": r["id"], "task": "page", "field": "locator.page", "before": str(lp), "after": None,
                    "evidence": "not a page number"})
        loc.pop("page", None)
        fl["page_unreadable"] = True
    tb = r.get("temperature_bound")
    if isinstance(tb, str):
        s = tb.strip().lower()
        if s in ("", "none", "null"):
            r.pop("temperature_bound", None)
        else:
            op = {"upper": "<", "lower": ">", "below": "<", "above": ">"}.get(s, tb.strip())
            d = {"op": op}
            if isinstance(r.get("temperature_K"), (int, float)):
                d["value_K"] = r["temperature_K"]
            if r.get("temperature_printed"):
                d["printed"] = r["temperature_printed"]
            r["temperature_bound"] = d
    rs = r.get("reference_state")
    if isinstance(rs, str):
        s = rs.strip()
        q = r.get("quote") or ""
        if len(s) > 300 and s not in q:
            log.append({"id": r["id"], "task": "reference_state", "field": "reference_state",
                        "before": s[:120] + "... (%d characters)" % len(s), "after": None,
                        "evidence": "free text of %d characters not printed in the quote (reader commentary)" % len(s)})
            r.pop("reference_state", None)
            fl["reference_state_dropped_not_in_quote"] = True
        elif s:
            r["reference_state"] = {"printed": s}
        else:
            r.pop("reference_state", None)


# ------------------------------------------------------------------------------------------ duplicates

def _j(x):
    return json.dumps(x, sort_keys=True, ensure_ascii=False)


def _page(r):
    return r.get("quote_page") if r.get("quote_page") is not None else (r.get("locator") or {}).get("page")


def dup_key(r):
    """Same paper, page, property, value, unit, temperature, phases and quote, and nothing that tells two
    printed states apart (compositions, component, reference state, method, evidence class, cited sources,
    reaction, range/expression, structure, bound, uncertainty, measured-in)."""
    return (r.get("paper"), _page(r), r.get("system"), r.get("property"), _j(r.get("value")), r.get("unit"),
            _j(r.get("temperature_K")), _j(sorted(_j(p) for p in r.get("phases") or [])), r.get("quote"),
            _j(r.get("phase_compositions")), _j(r.get("overall_composition")), r.get("component"),
            _j(r.get("reference_state")), r.get("method"), r.get("evidence_class"),
            _j(sorted(c.get("label") or c.get("reference") or "" for c in r.get("cited_sources") or [])),
            _j(r.get("reaction")), _j(r.get("value_range")), r.get("value_expression"), _j(r.get("structure")),
            _j(r.get("temperature_bound")), _j((r.get("uncertainty") or {}).get("value")),
            r.get("observation_type"), _j(r.get("temperature_range_K")), _j(r.get("temperature_max_K")),
            (r.get("origin") or {}).get("measured_in"))


def dedupe(records, counts):
    groups = defaultdict(list)
    for r in records:
        groups[dup_key(r)].append(r)
    keep = []
    for r in records:
        g = groups[dup_key(r)]
        if len(g) == 1:
            keep.append(r)
            continue
        if g[0] is not r and "_dup_done" in g[0]:
            continue
        best = sorted(g, key=lambda x: (not (x.get("flags") or {}).get("extraction_accepted"),
                                        -len(x.get("readers") or []), x["id"]))[0]
        others = sorted((x for x in g if x is not best), key=lambda x: x["id"])
        best["merged_ids"] = [x["id"] for x in others]
        best["merged_pids"] = sorted({x.get("pid") for x in others if x.get("pid") and x.get("pid") != best.get("pid")})
        if not best["merged_pids"]:
            del best["merged_pids"]
        locs = {((x.get("locator") or {}).get("ref") or (x.get("locator") or {}).get("asset"),
                 (x.get("locator") or {}).get("label"), (x.get("locator") or {}).get("row")) for x in g}
        arts = {(a, b) for a, b, _ in locs}
        scope = "same_row" if len(locs) == 1 else ("rows_of_one_artifact" if len(arts) == 1 else "several_artifacts")
        counts["groups_" + scope] += 1
        counts["records_removed_" + scope] += len(others)
        counts["groups"] += 1
        counts["records_removed"] += len(others)
        g[0]["_dup_done"] = True
        keep.append(best)
    for r in records:
        r.pop("_dup_done", None)
    # keep the original record order
    order = {id(r): i for i, r in enumerate(records)}
    keep.sort(key=lambda r: order[id(r)])
    return keep


# ------------------------------------------------------------------------------------------ DFT outliers

OUTLIER_ABS = 50.0     # kJ/mol-atom
BASIS_TOL = 15.0       # kJ/mol-atom
SIGN_MIN_RES = 10.0    # a sign case needs |DFT - lit| above the parity band
NEG_HEADING = re.compile(r"(^|[\s(\[])[−\-–‐]\s*\(?\s*Δ|[−\-–‐]\s*Δ\s*[fH]|minus\s+(the\s+)?(enthalp|heat)|"
                         r"negative (of the )?(enthalp|heat)", re.I)
PER_COMPOUND = re.compile(r"per mol(e)? of (the )?(compound|formula|substance|phase)|per formula unit|"
                          r"mol\s*(of\s*)?formula|/\s*mol\s*f\.?\s*u|kJ/mol of compound", re.I)


def _load_parity():
    spec = importlib.util.spec_from_file_location("build_dft_parity_q", os.path.join(HERE, "build_dft_parity.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _formula_counts(P, r):
    """The record's single binary formula as printed (unreduced counts), mirroring build_dft_parity."""
    system = [e.capitalize() for e in (r.get("elements") or r.get("system", "").split("-")) if e]
    found = {}
    for lab in [p.get("label") for p in (r.get("phases") or []) if p.get("label")]:
        c, _ = P.parse_phase(lab, system)
        if c:
            found[P.reduced_key(c)] = (c, lab)
    if len(found) != 1:
        return None, None
    c, lab = next(iter(found.values()))
    return (c, lab) if len(c) == 2 else (None, None)


def _printed_negative(r):
    vp = r.get("value_printed")
    if isinstance(vp, (int, float)) and not isinstance(vp, bool):
        return vp < 0
    if isinstance(vp, str):
        return bool(re.match(r"^\s*[−\-–‐]", vp))
    q = r.get("quote")
    v = r.get("value_stored", r.get("value"))
    if q and isinstance(v, (int, float)):
        a = abs(v)
        forms = {f"{a:g}", f"{a:.1f}", f"{a:.2f}", f"{a:.0f}" if a == int(a) else f"{a:g}"}
        for f in forms:
            for m in re.finditer(re.escape(f), q):
                pre = q[max(0, m.start() - 2):m.start()]
                return bool(re.search(r"[−\-–‐]\s?$", pre))
    return None


def classify_dft_outliers(records, open_tdb, log, counts):
    P = _load_parity()
    dft = {k: P.DFT(k, open_tdb) for k in P.SOURCES}
    for r in records:
        if r.get("property") != "HM_FORM" or r.get("unit") != P.UNIT or r.get("evidence_class") not in P.CLASSES:
            continue
        v = r.get("value")
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            continue
        counts_f, label = _formula_counts(P, r)
        if counts_f is None:
            continue
        key = P.reduced_key(counts_f)
        m = {s: src.lookup(key) for s, src in dft.items()}
        m = {s: x["dft"] for s, x in m.items() if x}
        if not m:
            continue
        big = any(abs(d - v) > OUTLIER_ABS for d in m.values())
        sign = any((d < 0) != (v < 0) and abs(d - v) > SIGN_MIN_RES and v != 0 for d in m.values())
        if not (big or sign):
            continue
        counts["outliers_examined"] += 1
        fl = r.setdefault("flags", {})
        natoms = sum(counts_f.values())
        dtxt = ", ".join("%s %.1f" % ({"mp": "MP", "oqmd": "OQMD"}[s], d) for s, d in sorted(m.items()))
        loc = r.get("locator") or {}
        unit_txt = r.get("unit_printed") or r.get("unit_stored") or ""
        heading = " | ".join(x for x in (r.get("property_printed"), loc.get("caption"), r.get("unit_printed")) if x)
        context = " ".join(x for x in (loc.get("caption"), r.get("quote"), r.get("unit_printed"), r.get("unit_stored")) if x)
        entry = {"id": r["id"], "system": r.get("system"), "formula_label": label, "evidence_class": r["evidence_class"],
                 "value_before": v, "unit": r.get("unit"), "dft": m, "printed_value": r.get("value_printed"),
                 "printed_unit": unit_txt or None, "kind": "large_residual" if big else "sign"}
        # (a) per-formula value labelled per atom
        basis_ok = (not PER_ATOM.search(unit_txt)) or bool(PER_COMPOUND.search(context))
        cands = []
        if natoms != 1:
            for tag, val in (("divided by", v / float(natoms)), ("multiplied by", v * float(natoms))):
                if all(abs(d - val) <= BASIS_TOL for d in m.values()):
                    cands.append((tag, val))
        if basis_ok and big and len(cands) == 1:
            tag, val = cands[0]
            r["value_before_quality_fix"] = v
            r["value"] = _sig(round(val, 6))
            fl["basis_corrected_from_formula"] = True
            r["basis_correction"] = {"formula": label, "atoms_per_formula_unit": float(natoms),
                                     "operation": "%s %g" % (tag, float(natoms))}
            unc = r.get("uncertainty")
            if isinstance(unc, dict) and isinstance(unc.get("value"), (int, float)):
                unc["value"] = _sig(unc["value"] / float(natoms) if tag == "divided by" else unc["value"] * float(natoms))
            why = ("printed unit %r has no per-atom marker" % unit_txt) if not PER_ATOM.search(unit_txt) else \
                "caption / quote says per mole of compound"
            entry.update({"class": "basis_corrected_from_formula", "value_after": r["value"],
                          "evidence": "%s; %s has %g atoms per formula unit; %g %s %g = %.2f, within %g of DFT (%s)"
                                      % (why, label, float(natoms), v, tag, float(natoms), val, BASIS_TOL, dtxt)})
            counts["class_basis_corrected_from_formula"] += 1
            log.append(entry)
            continue
        # (b) unsigned magnitude
        neg_printed = _printed_negative(r)
        if v > 0 and neg_printed is not True and all(abs(d + v) <= BASIS_TOL < abs(d - v) - BASIS_TOL for d in m.values()) \
                and not fl.get("sign_pdf_page_verified"):
            hm = NEG_HEADING.search(heading)
            if hm:
                r["value_before_quality_fix"] = v
                r["value"] = -v
                fl["sign_from_heading"] = True
                r["sign_heading_printed"] = heading
                entry.update({"class": "sign_from_heading", "value_after": -v,
                              "evidence": "printed heading %r carries the minus sign; -value within %g of DFT (%s)"
                                          % (heading, BASIS_TOL, dtxt)})
                counts["class_sign_from_heading"] += 1
            else:
                fl["unsigned_in_print"] = True
                entry.update({"class": "unsigned_in_print", "value_after": v,
                              "evidence": "printed without a minus sign (%r) and no negated heading in the record; "
                                          "-value within %g of DFT (%s); value kept as printed, excluded from the "
                                          "default parity set" % (r.get("value_printed", v), BASIS_TOL, dtxt)})
                counts["class_unsigned_in_print"] += 1
            log.append(entry)
            continue
        fl["dft_outlier_unexplained"] = True
        note = []
        if natoms != 1 and any(all(abs(d - x) <= BASIS_TOL for d in m.values()) for x in (v / natoms, v * natoms)):
            note.append("per-formula reading would match DFT but the printed unit is explicitly per atom")
        if all(abs(d + v) <= BASIS_TOL for d in m.values()) and neg_printed:
            note.append("opposite sign would match DFT but a minus sign is printed")
        entry.update({"class": "dft_outlier_unexplained", "value_after": v,
                      "evidence": "; ".join(note) or "no deterministic rule applies (DFT %s)" % dtxt})
        counts["class_dft_outlier_unexplained"] += 1
        log.append(entry)


# ------------------------------------------------------------------------------------------ driver

def apply(records, open_tdb):
    """Mutates and returns (records, report)."""
    log, counts = [], Counter()
    kinds_map, prop_map = Counter(), Counter()
    for r in records:
        fix_form(r, log)
        rx = r.get("reaction") or {}
        raw = rx.get("type")
        if raw:
            k, rule = _kind(raw, [p.get("label") or "" for p in r.get("phases") or []])
            r["transition_kind"], r["transition_kind_printed"] = k, raw
            rx["type"] = k
            kinds_map[(raw, k, rule)] += 1
        code = r.get("property")
        new = clean_property_code(code, r.get("unit"))
        if new != code:
            r["property_code_raw"], r["property"] = code, new
            prop_map[(code, new)] += 1
        normalize_energy_unit(r, counts)
    n0 = len(records)
    dcounts = Counter()
    records = dedupe(records, dcounts)
    assert len(records) == n0 - dcounts["records_removed"]
    classify_dft_outliers(records, open_tdb, log, counts)
    report = {
        "duplicates": dict(sorted(dcounts.items())),
        "transition_kind_map": [{"printed": a, "transition_kind": b, "rule": c, "n": n}
                                for (a, b, c), n in sorted(kinds_map.items(), key=lambda x: (x[0][1], -x[1], x[0][0]))],
        "transition_kind_counts": dict(sorted(Counter(k for (_, k, _), n in kinds_map.items() for _ in range(n)).items())),
        "property_code_map": [{"raw": a, "clean": b, "n": n} for (a, b), n in sorted(prop_map.items())],
        "counts": dict(sorted(counts.items())),
        "changes": log,
    }
    report["counts"]["transition_kind_records"] = sum(kinds_map.values())
    report["counts"]["transition_kind_variants"] = len({a for a, _, _ in kinds_map})
    report["counts"]["transition_kind_relabelled"] = sum(n for (a, b, _), n in kinds_map.items() if a != b)
    report["counts"]["property_codes_cleaned"] = sum(prop_map.values())
    report["counts"]["pages_nulled"] = sum(1 for x in log if x.get("task") == "page")
    report["counts"]["reference_state_dropped"] = sum(1 for x in log if x.get("task") == "reference_state")
    return records, report


def quality_md(report, release):
    c, d = report["counts"], report["duplicates"]
    L = ["# ORCHID literature data: quality fixes", "",
         "Release: `%s`. Written by `tools/export_literature_data.py` via `tools/literature_quality.py` "
         "(deterministic). Machine-readable copy: `quality_fixes.json`." % release, "",
         "Every rule re-labels, normalizes, merges or flags; no value is invented. Where a value changed, the "
         "original is kept in the record (`value_before_quality_fix`, `transition_kind_printed`, "
         "`property_code_raw`, `value_stored`/`unit_stored`, `merged_ids`).", "",
         "## 1. Duplicates", "",
         "Records identical in paper, page, property, value, unit, temperature, phases and quote, and in every "
         "field that distinguishes two printed states (compositions, component, reference state, method, evidence "
         "class, cited sources, reaction, range, structure, bound, uncertainty), were merged into one record; the "
         "others are listed in its `merged_ids` (and `merged_pids` when they came from other rows/artifacts).", "",
         "| scope | groups | records removed |", "|---|---:|---:|"]
    for s in ("same_row", "rows_of_one_artifact", "several_artifacts"):
        L.append("| %s | %d | %d |" % (s.replace("_", " "), d.get("groups_" + s, 0), d.get("records_removed_" + s, 0)))
    L.append("| **total** | **%d** | **%d** |" % (d.get("groups", 0), d.get("records_removed", 0)))
    L += ["", "## 2. Controlled vocabularies", "",
          "`transition_kind` (%d records, %d printed variants, %d relabelled); raw value in `transition_kind_printed`; "
          "`reaction.type` carries the controlled term." % (c["transition_kind_records"], c["transition_kind_variants"],
                                                            c["transition_kind_relabelled"]), "",
          "| transition_kind | records |", "|---|---:|"]
    for k, n in report["transition_kind_counts"].items():
        L.append("| %s | %d |" % (k, n))
    L += ["", "Mapped printed variants (all; `rule` says why):", "", "| printed | transition_kind | rule | n |",
          "|---|---|---|---:|"]
    for x in report["transition_kind_map"]:
        L.append("| %s | %s | %s | %d |" % (str(x["printed"]).replace("|", "\\|"), x["transition_kind"], x["rule"], x["n"]))
    L += ["", "Property codes (%d records; raw code in `property_code_raw`):" % c["property_codes_cleaned"], "",
          "| raw | clean | n |", "|---|---|---:|"]
    for x in report["property_code_map"]:
        L.append("| `%s` | `%s` | %d |" % (x["raw"], x["clean"], x["n"]))
    pages = [x for x in report["changes"] if x.get("task") in ("page", "reference_state")]
    L += ["", "## 3. Garbled fields", "",
          "`quote_page` values that are not page numbers were set to null with `flags.page_unreadable` (%d). "
          "`temperature_bound` is always `{op, value_K, printed}` (string forms converted, \"none\" removed); "
          "`reference_state` is always an object (free text becomes `{\"printed\": ...}`); free text not printed "
          "in the quote and longer than 300 characters (reader commentary) was dropped with "
          "`flags.reference_state_dropped_not_in_quote` (%d)." % (c["pages_nulled"], c["reference_state_dropped"]), "",
          "| id | field | before | after | evidence |", "|---|---|---|---|---|"]
    for x in pages:
        L.append("| %s | %s | %s | %s | %s |" % (x["id"], x["field"], str(x["before"]).replace("|", "\\|")[:140],
                                                 x["after"], x["evidence"]))
    out = [x for x in report["changes"] if x.get("class")]
    L += ["", "## 4. DFT parity outliers", "",
          "Formation enthalpies (kJ/mol-atom, experimental / model / unspecified) whose single binary formula matches "
          "a Materials Project or OQMD entry and that differ from DFT by more than %g kJ/mol-atom, or have the opposite "
          "sign with |DFT - lit| > %g (%d records), classified in this order:" % (OUTLIER_ABS, SIGN_MIN_RES, c.get("outliers_examined", 0)), "",
          "* `basis_corrected_from_formula` (%d): a large-residual case where lit / n_atoms (or lit x n_atoms, undoing "
          "an earlier per-formula division) lies within %g kJ/mol-atom of every matched DFT source, and the printed unit has no per-atom marker or the caption/quote says per mole of "
          "compound. Value converted; original in `value_before_quality_fix`, operation in `basis_correction`."
          % (c.get("class_basis_corrected_from_formula", 0), BASIS_TOL),
          "* `sign_from_heading` (%d): printed without a minus sign, -lit within %g of DFT, and the record's printed "
          "heading carries the minus (e.g. \"-ΔfH\"): sign applied, heading in `sign_heading_printed`."
          % (c.get("class_sign_from_heading", 0), BASIS_TOL),
          "* `unsigned_in_print` (%d): printed without a minus sign, -lit within %g of DFT (and lit itself more than "
          "twice that far, so small values near zero are not classified), no negated heading in the "
          "record: value kept as printed, flagged, excluded from the default parity set." % (c.get("class_unsigned_in_print", 0), BASIS_TOL),
          "* `dft_outlier_unexplained` (%d): no rule applies; flagged `flags.dft_outlier_unexplained`, value unchanged."
          % c.get("class_dft_outlier_unexplained", 0), "",
          "| id | class | formula | before | after | evidence |", "|---|---|---|---:|---:|---|"]
    for x in sorted(out, key=lambda x: (x["class"], x["system"] or "", x["id"])):
        L.append("| %s | %s | %s | %s | %s | %s |" % (x["id"], x["class"], x.get("formula_label"), x["value_before"],
                                                      x["value_after"], str(x["evidence"]).replace("|", "\\|")))
    L += ["", "## 5. Energy units", "",
          "Formation enthalpies / Gibbs energies of formation, mixing enthalpies, Gibbs energies of mixing and partial "
          "molar enthalpies are given in kJ (`value`, `unit`) with the value before conversion in "
          "`value_stored`/`unit_stored` and the printed value/unit beside them. `unit_basis` records how the "
          "per-atom basis was established (printed_unit, column_header, caption, reviewer, state_completion_*, "
          "magnitude_heuristic, assumed_per_mole_of_atoms, unstated).", "",
          "* mixing / partial values converted to kJ/mol-atom (explicit per-atom basis): %d" % c.get("mixing_converted_per_atom", 0),
          "* mixing / partial values converted to kJ/mol (basis unstated): %d" % c.get("mixing_converted_basis_unstated", 0),
          "* of these, per-atom labels from a heuristic or assumption demoted to basis unstated: %d" % c.get("mixing_basis_not_explicit_demoted", 0),
          "* formation values whose per-atom basis comes only from the record builder's magnitude heuristic "
          "(kept, flagged `flags.basis_inferred_from_magnitude`): see `index.json` totals", ""]
    return "\n".join(L) + "\n"
