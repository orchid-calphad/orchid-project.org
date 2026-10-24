#!/usr/bin/env python3
"""Audit round 2 (2026-10-24): thermochemical values, units, temperatures and reference states.

Called from tools/literature_gates.apply() right after the record-level page fixes (g0_page_fixes), before the other
gates. Implements the fixes of audit A3 (units, magnitudes, temperatures, reference states, uncertainty budgets) and
the literature-side items of audit A4 (E3 per-formula / per-atom decisions, E6 experimental values carrying 0 K,
E7 first-principles rows coded as model). Every rule keeps the stored field beside the new one and writes a
`corrections` entry (gate "R2"); every rule that changes a value or a state has its evidence in the printed table
(header, caption, column, label) and was checked on the rendered page for a sample of each class, or for the whole
table (tools/literature_visual_checks.json, keys `table_basis:<asset>`, `table_sign:<asset>`, `r2_class:<rule>`).

Rules (codes in QUALITY_FIXES.md, section "Audit round 2"):

  R2-basis     table_basis:<asset> page decision: one per-formula / per-atom divisor for a whole column, applied to
               the printed value (value = printed x unit factor / atoms per formula unit); the printed uncertainty is
               divided by the same number. Replaces the DFT-driven `basis_corrected_from_formula` rule, which
               accepted a conversion only when it landed near DFT (circular for the DFT calibration).
  R2-sign      table_sign:<asset> page decision: a column printed without minus signs is negated as a whole.
  R2-sigma     a printed sigma divided by the same atom count as its value (per_formula_unit_divisor).
  R2-sabs      SM_FORM values that are absolute entropies (positive, > 15 J/K per mole of atoms; a condensed
               compound's formation entropy from condensed elements cannot be that large) -> SM.
  R2-hcontent  H(T) - H(298) heat contents stored at the reference temperature -> the row's temperature
               (conditions.T_measurement_K of the record store); without one, not accepted.
  R2-expr      coefficients of printed expressions (A + B T, Cp polynomials) stored as values at one temperature
               -> property <P>_EXPRESSION_TERM with `property_expression`, no single temperature.
  R2-t298      quantities printed with a 298 subscript (ΔfH°298, S°298, Cp,298) carrying another column's
               temperature -> 298.15 K.
  R2-tstd      standard formation enthalpies from direct-synthesis / solute-solvent-drop calorimetry (or labelled
               standard / at 298 K) stored at the calorimeter temperature -> 298.15 K (calorimeter_temperature_K).
  R2-t0exp     experimental values carrying 0 K from mixed first-principles / experiment tables -> 298.15 K for
               formation enthalpies (assumed standard), unstated otherwise.
  R2-dftclass  model-class rows of first-principles tables whose method is a DFT code or that are 0 K values not
               labelled CALPHAD / assessed / optimized -> evidence class dft.
  R2-refcomp   formation quantities referred to compounds (e.g. from the oxides) -> <P>_FROM_COMPOUNDS (not _FORM).
  R2-amorph    liquid / amorphous states printed at room temperature (liq-am, glass) -> flags.non_crystalline_state.
  R2-basislost per-atom units written "mol/at" kept per atom.
"""
from __future__ import annotations

import re

import build_dft_parity as parity
import literature_gates as G

GATE = "R2"
REVIEW = "ORCHID audit round 2 (A3 / A4) 2026-10-24"
RELEASE = "literature-r2-thermo 2026-10-24"

ENERGY_PROPS = {"HM_FORM", "GM_FORM", "HM_MIX", "GM_MIX", "PARTIAL_ENTHALPY"}
ENTROPY_PROPS = {"SM", "SM_FORM", "ENTROPY", "CPM", "SM_MIX"}
E_FACTOR = [(re.compile(r"^\s*kcal", re.I), 4.184), (re.compile(r"^\s*cal", re.I), 0.004184),
            (re.compile(r"^\s*kj", re.I), 1.0), (re.compile(r"^\s*j(\b|[\s/(.·]|$)", re.I), 0.001)]
S_FACTOR = [(re.compile(r"^\s*cal", re.I), 4.184), (re.compile(r"^\s*kj", re.I), 1000.0),
            (re.compile(r"^\s*j(\b|[\s/(.·]|$)", re.I), 1.0)]
STATE_SUFFIX = re.compile(r"\s*\((?:cr|s|c|cryst|solid|liq-am|am|l|liq|g)\)\s*$", re.I)


def _sig(x):
    return float("%.10g" % x)


def num(x):
    if isinstance(x, bool):
        return None
    if isinstance(x, (int, float)):
        return float(x)
    try:
        return float(str(x).replace("−", "-").replace("–", "-").replace(",", "").strip())
    except (TypeError, ValueError):
        return None


def atoms_of_label(label, elements):
    """Atoms per formula unit of a printed compound label (None when the label is not one formula)."""
    if not label:
        return None
    lab = STATE_SUFFIX.sub("", str(label))
    c, _ = parity.parse_phase(lab, [e.capitalize() for e in elements or []])
    if not c:
        return None
    n = float(sum(c.values()))
    return n if n > 0 else None


def record_atoms(r):
    found = {atoms_of_label(p.get("label"), r.get("elements")) for p in r.get("phases") or [] if p.get("label")}
    found.discard(None)
    return found.pop() if len(found) == 1 else None


def unit_factor(prop, unit):
    pats = E_FACTOR if prop in ENERGY_PROPS else S_FACTOR
    u = (unit or "").strip().lstrip("( ")
    return next((f for p, f in pats if p.search(u)), None)


def canonical_unit(prop):
    return "kJ/mol-atom" if prop in ENERGY_PROPS else "J/mol-atom/K"


def printed_sigma(r):
    """The printed sigma (record units, before any division by an atom count) or None."""
    u = r.get("uncertainty") or {}
    if u.get("basis") not in ("printed", "printed_range") or not isinstance(u.get("value"), (int, float)):
        return None
    return float(u["value"])


def set_sigma(r, new, rule, evidence, check=None):
    u = r.get("uncertainty")
    if not isinstance(u, dict) or not isinstance(u.get("value"), (int, float)):
        return
    if abs(u["value"] - new) <= 1e-12 * max(1.0, abs(new)):
        return
    r.setdefault("uncertainty_stored", dict(u))
    G.correct(r, GATE, "uncertainty.value", u["value"], _sig(new), rule, evidence, check)
    u["value"] = _sig(new)
    u["rescaled_with_value"] = True


def _checks(visual, prefix):
    return {k[len(prefix):]: v for k, v in visual.items() if k.startswith(prefix)}


# ------------------------------------------------------------------------------------------------- R2-basis

def r2_table_basis(records, visual, rep):
    """One basis decision per printed column, read on the page (table_basis:<asset>)."""
    decided = {a: c for a, c in _checks(visual, "table_basis:").items() if c.get("verdict") == "decided"}
    for r in records:
        chk = decided.get(G.asset_of(r) or "")
        if not chk:
            continue
        col = (chk.get("columns") or {}).get(r.get("property"))
        if not col or r["id"] in (chk.get("rows_excluded") or []):
            continue
        basis = col["basis"]
        vp = num(r.get("value_printed"))
        f = unit_factor(r["property"], col.get("unit") or r.get("unit_printed"))
        if vp is None or f is None:
            G.fail(r, GATE, "table_basis_value_not_convertible", "not_espei_ready",
                   "page decision %s, but the printed value / unit could not be converted" % basis, chk)
            rep.add(GATE, "table_basis_value_not_convertible", r, "flagged")
            continue
        n = record_atoms(r) if basis == "per_formula_unit" else 1.0
        if n is None:
            G.fail(r, GATE, "table_basis_formula_unparsed", "not_espei_ready",
                   "per-formula column (page) but the phase label is not one formula", chk)
            rep.add(GATE, "table_basis_formula_unparsed", r, "flagged")
            continue
        new = _sig(round(vp * f / n, 8))
        before = r.get("value")
        unit_new = canonical_unit(r["property"])
        ev = {"printed": "%s %s" % (r.get("value_printed"), col.get("unit") or r.get("unit_printed") or ""),
              "column_basis": basis, "atoms_per_formula_unit": n, "page": chk.get("page")}
        if not isinstance(before, (int, float)) or abs(before - new) > 1e-6 * max(1.0, abs(new)) or r.get("unit") != unit_new:
            r.setdefault("value_stored", before)
            r.setdefault("unit_stored", r.get("unit"))
            G.correct(r, GATE, "value", before, new, "table_basis_page", ev, chk)
            r["value"], r["unit"] = new, unit_new
            rep.add(GATE, "table_basis_value", r, "corrected", before=before, after=new)
        r["unit_basis"] = "page_checked_table"
        if basis == "per_formula_unit":
            r["basis_correction"] = {"formula_atoms": n, "operation": "divided by %g" % n, "basis": "page_checked_table"}
        else:
            r.pop("basis_correction", None)
        u0 = r.get("uncertainty") or {}
        if u0.get("basis") in ("stated_default", "reading_resolution") and isinstance(u0.get("value"), (int, float)) \
                and isinstance(before, (int, float)) and before:
            r.setdefault("uncertainty_stored", dict(u0))
            u0["value"] = _sig(abs(u0["value"] * new / before))
        s = printed_sigma(r)
        if s is not None:
            # the record's sigma is the printed one in the canonical energy unit (never divided before)
            s0 = (r.get("uncertainty_stored") or {}).get("value", s)
            set_sigma(r, s0 / n if basis == "per_formula_unit" else s0, "table_basis_page", ev, chk)
        fl = r["flags"]
        for k in ("basis_corrected_from_formula", "basis_inferred_from_magnitude", "dft_outlier_unexplained"):
            if fl.pop(k, None):
                fl[k + "_before_page_basis"] = True
        if fl.get("magnitude_suspect"):
            fl["magnitude_suspect_before_page_fix"], fl["magnitude_suspect"] = True, False
        if col.get("flag"):
            G.fail(r, GATE, col["flag"], col.get("effect") or "not_espei_ready", col.get("detail"), chk)


# ------------------------------------------------------------------------------------------------- R2-sign

def r2_table_sign(records, visual, rep):
    decided = {a: c for a, c in _checks(visual, "table_sign:").items() if c.get("verdict") == "decided"}
    for r in records:
        chk = decided.get(G.asset_of(r) or "")
        if not chk or r.get("property") != chk.get("property"):
            continue
        if chk.get("column_printed") and not re.search(chk["column_printed"], str(r.get("property_printed") or "")):
            continue
        v = r.get("value")
        if not isinstance(v, (int, float)) or v <= 0:
            continue
        r.setdefault("value_stored", v)
        G.correct(r, GATE, "value", v, -v, "column_printed_unsigned_page",
                  {"printed": r.get("value_printed"), "column": r.get("property_printed")}, chk)
        r["value"] = -v
        fl = r["flags"]
        fl.pop("unsigned_in_print", None)
        fl["sign_from_page_column"] = True
        if fl.pop("dft_outlier_unexplained", None):
            fl["dft_outlier_unexplained_before_page_sign"] = True
        rep.add(GATE, "table_sign_value", r, "corrected", before=v, after=-v)


# ------------------------------------------------------------------------------------------------- R2-sigma

def r2_sigma_with_divisor(records, rep, check=None):
    """The record builder divided a per-formula value by the atom count but kept the printed sigma."""
    for r in records:
        d = r.get("per_formula_unit_divisor")
        if not isinstance(d, (int, float)) or d <= 1 or r.get("unit_basis") == "page_checked_table":
            continue
        s = printed_sigma(r)
        if s is None or (r.get("uncertainty") or {}).get("rescaled_with_value"):
            continue
        set_sigma(r, s / d, "printed_sigma_rescaled_with_value",
                  {"divisor": d, "sigma_printed": s, "why": "value divided by the atoms per formula unit"}, check)
        rep.add(GATE, "sigma_rescaled_with_value", r, "corrected", before=s, after=_sig(s / d))


# ------------------------------------------------------------------------------------------------- R2-sabs

SABS_LIMIT = 15.0   # J/K per mole of atoms
SABS_LIMIT_FORMULA = 20.0   # J/K per mole of formula units / basis unstated
SABS_MAX = 150.0   # J/K per mole of atoms: above every absolute entropy S298 of a condensed phase


def entropy_too_large_for_formation(r):
    """A positive entropy above 15 J/K per mole of atoms (20 J/K per formula unit / basis unstated) cannot be a formation
    entropy of a condensed compound; it is an absolute entropy."""
    v, u = r.get("value"), r.get("unit") or ""
    if not isinstance(v, (int, float)) or isinstance(v, bool) or v <= 0:
        return False
    return v > (SABS_LIMIT if re.search(r"atom|mol-at", u, re.I) else SABS_LIMIT_FORMULA)


def r2_absolute_entropy(records, rep, check=None):
    for r in records:
        if r.get("property") != "SM_FORM" or not isinstance(r.get("value"), (int, float)) or r["value"] <= 0:
            continue
        u = r.get("unit") or ""
        per_atom = bool(re.search(r"atom|mol-at", u, re.I))
        v = r["value"]
        # per mole of atoms: > 15 J/K; per formula unit or basis unstated: > 20 J/K (a positive formation entropy of
        # that size is not physical for a condensed compound under either reading)
        if v <= (SABS_LIMIT if per_atom else SABS_LIMIT_FORMULA):
            continue
        n = record_atoms(r) or 1.0
        if (v if per_atom else v / n) > SABS_MAX:
            # far above any absolute entropy too: a scale error, not an S298 (left to the magnitude gate, not accepted)
            G.fail(r, GATE, "formation_entropy_out_of_range", "not_accepted",
                   "SM_FORM %g %s: above %g J/K per mole of atoms" % (v, u, SABS_MAX), check)
            rep.add(GATE, "formation_entropy_out_of_range", r, "demoted")
            continue
        r["property_code_before_gates"] = r.get("property_code_before_gates") or "SM_FORM"
        G.correct(r, GATE, "property", "SM_FORM", "SM", "absolute_entropy_labelled_formation",
                  {"printed": r.get("property_printed"), "value": v, "unit": u,
                   "why": "a positive formation entropy above %g J/K per mole of atoms is not physical for a condensed "
                          "compound; the value is the absolute (third-law) entropy S298" % SABS_LIMIT}, check)
        r["property"] = "SM"
        r["property_group"] = "heat capacity / entropy / enthalpy"
        rep.add(GATE, "absolute_entropy_labelled_formation", r, "corrected")


# ------------------------------------------------------------------------------------------------- R2-hcontent

HREF_FIRST = re.compile(r"H\s*°?\s*\(?\s*(?:298|0)\s*(?:\.15)?\s*(?:K)?\s*\)?\s*[-−–]\s*H", re.I)   # H(298)-H(0)


def r2_heat_content(records, rep, check=None):
    for r in records:
        if r.get("property") != "HEAT_CONTENT":
            continue
        T, v = r.get("temperature_K"), r.get("value")
        if not isinstance(T, (int, float)) or abs(T - 298.15) > 2 or not v:
            continue
        lab = str(r.get("property_printed") or "")
        if HREF_FIRST.search(lab) and not re.search(r"H\s*°?\s*[_(]?\s*T", lab):
            continue            # H(298) - H(0) or H(0) - H(298): a value at the reference temperature
        Tm = r.get("_T_measurement_K")
        if isinstance(Tm, (int, float)) and Tm > 0 and abs(Tm - 298.15) > 2:
            r["temperature_K_stored"] = T
            r["temperature_K"] = Tm
            r["reference_temperature_K"] = 298.15
            r["temperature_basis"] = "row_temperature"
            G.correct(r, GATE, "temperature_K", T, Tm, "heat_content_row_temperature",
                      {"printed": lab, "why": "H(T) - H(298) is a value at the row temperature T; 298 K is its reference"},
                      check)
            rep.add(GATE, "heat_content_row_temperature", r, "corrected", before=T, after=Tm)
        else:
            G.fail(r, GATE, "heat_content_temperature_missing", "not_accepted",
                   "H(T) - H(298) stored at the reference temperature and the row temperature was not read", check)
            rep.add(GATE, "heat_content_temperature_missing", r, "demoted")


# ------------------------------------------------------------------------------------------------- R2-expr

EXPR_PROPS = {"GIBBS_ENERGY", "GM_FORM", "GM_MIX", "HM_FORM", "CPM", "EMF", "SM_FORM", "HM_MIX", "CHEMICAL_POTENTIAL"}
EXPR_A = re.compile(r"constant term|coefficient A\b|\(A\)|\bterm A\b|\(constant\)", re.I)
EXPR_B = re.compile(r"coefficient of T|T[- ]coefficient|coefficient B\b|\(B\)|B\s*coefficient|B[\s\d]*coefficient of T|"
                    r"\bterm B\b", re.I)
EXPR_CP = re.compile(r"\bc_?p\b\s*coefficient\s*([a-f])\b|coefficient\s*([a-f])\b\s*\(", re.I)
EXPR_PRINTED = re.compile(r"\[printed:\s*([^\]]*?T[^\]]*)\]", re.I)


def expression_term(r):
    lab = str(r.get("property_printed") or "")
    m = EXPR_CP.search(lab)
    if m:
        return (m.group(1) or m.group(2)).lower(), "Cp = a + b T + c T^2 + d T^-2 (+ ...)", lab
    if EXPR_B.search(lab):
        return "B", "A + B T", lab
    if EXPR_A.search(lab):
        return "A", "A + B T", lab
    m = EXPR_PRINTED.search(lab)
    if m:
        return "A", "A + B T", m.group(1).strip()
    return None


def r2_expression_terms(records, rep, check=None):
    for r in records:
        p = r.get("property")
        if p not in EXPR_PROPS or p.endswith("_EXPRESSION_TERM"):
            continue
        t = expression_term(r)
        if not t:
            continue
        term, form, printed = t
        r["property_code_before_gates"] = r.get("property_code_before_gates") or p
        new = p + "_EXPRESSION_TERM"
        G.correct(r, GATE, "property", p, new, "expression_coefficient",
                  {"printed": r.get("property_printed"), "term": term, "form": form}, check)
        r["property"] = new
        r["property_group"] = "model expression"
        r["property_expression"] = {"form": form, "term": term, "printed": printed,
                                    "temperature_range_K": r.get("temperature_range_K")}
        T = r.get("temperature_K")
        if T is not None:
            r["temperature_K_stored"] = T
            r["temperature_K"] = None
            G.correct(r, GATE, "temperature_K", T, None, "expression_coefficient_no_single_temperature",
                      {"printed": r.get("property_printed"),
                       "why": "a coefficient of an expression in T has no single temperature; validity range kept"},
                      check)
        G.fail(r, GATE, "expression_coefficient_not_a_value", "not_espei_ready",
               "coefficient %s of %s" % (term, form), check)
        rep.add(GATE, "expression_coefficient", r, "corrected")


# ------------------------------------------------------------------------------------------------- R2-t298

T298_LABEL = re.compile(r"(?:H|S|C\s*_?p|G)\s*°?\s*(?:_|,)?\s*\{?\s*(?:m\s*)?\(?\s*298(?:\.15)?\s*(?:K)?\s*\)?\}?|"
                        r"298(?:\.15)?\s*K\b|\(298(?:\.15)?\s*K?\)", re.I)
HEAT_INCREMENT = re.compile(r"H\s*°?\s*\(?\s*T|H_?T\s*[-−–]|[-−–]\s*H\s*°?\s*_?\(?\s*298|S\s*°?\s*\(?\s*T\s*\)?\s*[-−–]", re.I)


def r2_t298(records, rep, check=None):
    for r in records:
        p = r.get("property")
        if p in ("HEAT_CONTENT", "INVARIANT_TEMPERATURE", "TRANSITION_TEMPERATURE") or p.endswith("_EXPRESSION_TERM"):
            continue
        lab = str(r.get("property_printed") or "")
        T = r.get("temperature_K")
        if not isinstance(T, (int, float)) or abs(T - 298.15) <= 2 or not T298_LABEL.search(lab) or HEAT_INCREMENT.search(lab):
            continue
        if r.get("temperature_basis") in ("quantity_label", "row_temperature"):
            continue
        r["temperature_K_stored"] = T
        r["temperature_K"] = 298.15
        r["temperature_basis"] = "quantity_label"
        G.correct(r, GATE, "temperature_K", T, 298.15, "quantity_printed_at_298",
                  {"printed": lab, "why": "the printed quantity carries the 298 K subscript; the stored temperature came "
                                          "from another column"}, check)
        rep.add(GATE, "quantity_printed_at_298", r, "corrected", before=T, after=298.15)


# ------------------------------------------------------------------------------------------------- R2-tstd

KLEPPA = re.compile(r"direct[- ]synthesis|solute[\s\-–]+solvent[\s\-–]+drop|high[- ]temperature (?:direct )?synthesis calorim", re.I)
STD_LABEL = re.compile(r"standard|298|T_?ref|\bTref\b", re.I)
AT_T_CAPTION = re.compile(r"\bat\s+1?\d{3,4}(?:\.\d+)?\s*K\b(?!.*298)|referred to (?:the )?(?:liquid|bcc|fcc)[^.]*at|"
                          r"at\s+\d{3,4}\s*°\s*C", re.I)
AT_298_CAPTION = re.compile(r"at\s+298(?:\.15)?\s*K|298(?:\.15)?\s*K", re.I)


def r2_standard_temperature(records, rep, check=None):
    for r in records:
        if r.get("property") != "HM_FORM" or r.get("evidence_class") != "experimental":
            continue
        T = r.get("temperature_K")
        if not isinstance(T, (int, float)) or T < 600:
            continue
        meth = str(r.get("method") or "")
        lab = str(r.get("property_printed") or "")
        cap = str((r.get("locator") or {}).get("caption") or "")
        if re.search(r"oxide|f,ox", lab + cap, re.I):
            continue
        why = None
        if KLEPPA.search(meth) and not AT_T_CAPTION.search(cap):
            why = "direct-synthesis / solute-solvent drop calorimetry gives the standard enthalpy of formation at " \
                  "298 K (the compound's heat content is measured and subtracted); %g K is the calorimeter temperature" % T
        elif re.search(r"calori", meth, re.I) and AT_298_CAPTION.search(cap) and not AT_T_CAPTION.search(cap):
            why = "caption gives the enthalpies of formation at 298 K; %g K is the calorimeter temperature" % T
        elif re.search(r"calori", meth, re.I) and STD_LABEL.search(lab) and not AT_T_CAPTION.search(cap) \
                and not re.search(r"\d{3,4}\s*K", lab):
            why = "standard enthalpy of formation (printed label %r); %g K is the calorimeter temperature" % (lab, T)
        if not why:
            continue
        r["calorimeter_temperature_K"] = T
        r["temperature_K_stored"] = T
        r["temperature_K"] = 298.15
        r["temperature_basis"] = "standard_state_of_method"
        G.correct(r, GATE, "temperature_K", T, 298.15, "standard_enthalpy_at_calorimeter_temperature",
                  {"printed": lab, "method": meth, "why": why}, check)
        rep.add(GATE, "standard_enthalpy_at_calorimeter_temperature", r, "corrected", before=T, after=298.15)


# ------------------------------------------------------------------------------------------------- R2-t0exp

def r2_experimental_zero_kelvin(records, rep, check=None):
    for r in records:
        if r.get("evidence_class") != "experimental" or r.get("temperature_K") != 0:
            continue
        p = r.get("property")
        if p in ("INVARIANT_TEMPERATURE", "TRANSITION_TEMPERATURE"):
            continue
        cap = str((r.get("locator") or {}).get("caption") or "")
        ambient = bool(re.search(r"ambient|room temperature", cap, re.I))
        new = 298.15 if p in ("HM_FORM", "SM_FORM", "GM_FORM") or ambient else None
        r["temperature_K_stored"] = 0.0
        r["temperature_K"] = new
        r["temperature_basis"] = ("caption_ambient" if ambient else "assumed_standard_298K") if new else "not_printed"
        G.correct(r, GATE, "temperature_K", 0.0, new, "experimental_value_from_0K_caption",
                  {"caption": str((r.get("locator") or {}).get("caption") or "")[:160],
                   "why": "0 K is the temperature of the first-principles columns of the table, not of this cited "
                          "experimental value; formation quantities are standard (298.15 K), others unstated"}, check)
        rep.add(GATE, "experimental_value_from_0K_caption", r, "corrected")


# ------------------------------------------------------------------------------------------------- R2-dftclass

FP_CAPTION = re.compile(r"first[- ]principles|ab[- ]initio|\bDFT\b|density functional|VASP|\bGGA\b|\bPBE\b|\bLDA\b|"
                        r"\bPAW\b|FLAPW|LMTO|\bSQS\b|0\s*K", re.I)
DFT_METHOD = re.compile(r"first[- ]principles|ab[- ]initio|\bDFT\b|density functional|VASP|\bGGA\b|\bPBE\b|\bLDA\b|"
                        r"\bPAW\b|FLAPW|LMTO|WIEN2k|CASTEP|Quantum ESPRESSO|\bSQS\b|\bΔ\s*E\s*tot", re.I)
NOT_DFT = re.compile(r"CALPHAD|assess|optimi|\bopt\b|thermo-?calc|pandat|model|fit|evaluat|estimat|miedema|exp", re.I)


def r2_dft_class(records, rep, visual):
    for r in records:
        if r.get("evidence_class") != "model" or r.get("property") not in ("HM_FORM", "SM_FORM", "GM_FORM"):
            continue
        cap = str((r.get("locator") or {}).get("caption") or "")
        meth = str(r.get("method") or "")
        if not FP_CAPTION.search(cap + " " + str(r.get("property_printed") or "")):
            continue
        if NOT_DFT.search(meth):
            continue
        dft_meth = bool(DFT_METHOD.search(meth))
        zero_k = r.get("temperature_K") == 0
        if not (dft_meth or (zero_k and not meth.strip())):
            continue
        chk = visual.get("r2_class:dft_class_from_first_principles_table")
        r["evidence_class_before_gates"] = r.get("evidence_class_before_gates") or "model"
        G.correct(r, GATE, "evidence_class", "model", "dft", "first_principles_row_of_mixed_table",
                  {"caption": cap[:160], "method": meth or None, "T_K": r.get("temperature_K")}, chk)
        r["evidence_class"] = "dft"
        r.setdefault("origin", {})["class_gate"] = "dft"
        rep.add(GATE, "first_principles_row_of_mixed_table", r, "corrected")


# ------------------------------------------------------------------------------------------------- R2-refcomp

COMPOUND_TOKEN = re.compile(r"\b(?:[A-Z][a-z]?\d*(?:\.\d+)?){2,}\b")
FROM_COMPOUNDS = re.compile(r"f\s*,\s*oxides?|from (?:the )?(?:binary )?oxides|relative to (?:the )?oxides|"
                            r"from (?:the )?(?:constituent )?compounds", re.I)


def _ref_compounds(r):
    rs = r.get("reference_state") or {}
    txt = " ".join(str(v.get("phase") or "") for v in rs.values() if isinstance(v, dict)) if isinstance(rs, dict) else str(rs)
    txt = txt.replace("printed:", " ")
    toks = [t for t in COMPOUND_TOKEN.findall(txt)
            if len(set(re.findall(r"[A-Z][a-z]?", t))) >= 2 and not re.fullmatch(r"(?:[A-Z]{2,}_?)+\w*", t)]
    return toks


def r2_reference_compounds(records, rep, check=None):
    for r in records:
        p = r.get("property")
        if p not in ("HM_FORM", "GM_FORM", "SM_FORM"):
            continue
        lab = " ".join(str(x or "") for x in (r.get("property_printed"), (r.get("locator") or {}).get("caption")))
        toks = _ref_compounds(r)
        if not toks and not FROM_COMPOUNDS.search(lab):
            continue
        if re.search(r"f\s*,\s*elements|from (?:the )?elements", lab, re.I):
            continue
        new = p.replace("_FORM", "_FROM_COMPOUNDS")
        r["property_code_before_gates"] = r.get("property_code_before_gates") or p
        G.correct(r, GATE, "property", p, new, "formation_from_compounds",
                  {"printed": r.get("property_printed"), "reference_compounds": toks or None,
                   "why": "the printed reference state is compounds (e.g. oxides), not the elements"}, check)
        r["property"] = new
        G.fail(r, GATE, "formation_from_compounds_not_from_elements", "not_espei_ready",
               "reference: %s" % (", ".join(toks) if toks else lab[:120]), check)
        rep.add(GATE, "formation_from_compounds", r, "corrected")


# ------------------------------------------------------------------------------------------------- R2-amorph

NONCRYST = re.compile(r"\((?:liq-?am|am|amorph\w*|glass\w*)\)|\bamorphous\b|\bglass(?:y)?\b|liq-?am", re.I)


def r2_noncrystalline(records, rep):
    for r in records:
        if any(NONCRYST.search(str(p.get("label") or "")) for p in r.get("phases") or []):
            if not r["flags"].get("non_crystalline_state"):
                r["flags"]["non_crystalline_state"] = True
                rep.add(GATE, "non_crystalline_state", r, "flagged")


# ------------------------------------------------------------------------------------------------- R2-property

def r2_table_property(records, visual, rep):
    """table_property:<asset> page decision: a property code read wrongly for a whole column (e.g. Gibbs energies of
    formation coded as Gibbs energies of mixing)."""
    decided = {a: c for a, c in _checks(visual, "table_property:").items() if c.get("verdict") == "decided"}
    for r in records:
        chk = decided.get(G.asset_of(r) or "")
        if not chk:
            continue
        new = (chk.get("property_map") or {}).get(r.get("property"))
        if not new:
            continue
        old = r["property"]
        r["property_code_before_gates"] = r.get("property_code_before_gates") or old
        G.correct(r, GATE, "property", old, new, "column_property_read_from_page",
                  {"printed": r.get("property_printed"), "caption": str((r.get("locator") or {}).get("caption") or "")[:160]},
                  chk)
        r["property"] = new
        r["property_group"] = "formation" if new.endswith("_FORM") else r.get("property_group")
        rep.add(GATE, "column_property_read_from_page", r, "corrected", before=old, after=new)


# ------------------------------------------------------------------------------------------------- R2-tlabel

T_IN_LABEL = re.compile(r"(?:\(\s*|\bat\s+|_\{?\s*|\(at\s+)(\d{3,4}(?:\.\d+)?)\s*K\b")


def r2_temperature_in_label(records, rep, check=None):
    """The printed quantity names its temperature (G^ex_Li (633 K), ΔGf at 1000 K) but the record has none."""
    for r in records:
        if r.get("temperature_K") is not None or r.get("value") is None or r.get("temperature_range_K"):
            continue
        lab = str(r.get("property_printed") or "")
        ts = {float(x) for x in T_IN_LABEL.findall(lab)}
        if len(ts) != 1:
            continue
        T = ts.pop()
        T = 298.15 if T in (298.0, 298.1, 298.2) else T
        r["temperature_K"] = T
        r["temperature_basis"] = "quantity_label"
        G.correct(r, GATE, "temperature_K", None, T, "temperature_printed_in_quantity_label", {"printed": lab}, check)
        rep.add(GATE, "temperature_printed_in_quantity_label", r, "corrected", after=T)


# ------------------------------------------------------------------------------------------------- R2-basislost

PER_ATOM_SLASH = re.compile(r"mol(?:e)?\s*/\s*at(?:om)?s?\.?\s*$|mol\s*/\s*at\b", re.I)


def r2_basis_lost(records, rep):
    """A per-atom unit printed as 'kJ/mol/at' was stored as 'kJ/mol (basis unstated)'."""
    for r in records:
        u = r.get("unit") or ""
        if "basis unstated" not in u or not PER_ATOM_SLASH.search(str(r.get("unit_printed") or "")):
            continue
        new = u.replace("mol (basis unstated)", "mol-atom")
        G.correct(r, GATE, "unit", u, new, "per_atom_unit_printed_mol_per_at", {"printed": r.get("unit_printed")})
        r["unit"], r["unit_basis"] = new, "printed_unit"
        rep.add(GATE, "per_atom_unit_printed_mol_per_at", r, "corrected", before=u, after=new)


# ------------------------------------------------------------------------------------------------- R2-refgarbled

REF_UNIT = re.compile(r"^\s*printed:\s*(?:k?J|cal|kcal|eV)\b", re.I)
REF_ONE_ELEMENT = re.compile(r"^\s*printed:\s*[αβγδεζηθ]?-?\(?([A-Z][a-z]?)\)?\s*$")


def r2_reference_state_text(records, rep):
    """Reference-state entries that are a unit ('printed:kJ/mol atom') or another element's phase ('printed:α-Zr' under
    Al) are not reference states: moved to reference_state_stored."""
    for r in records:
        rs = r.get("reference_state")
        if not isinstance(rs, dict):
            continue
        drop = []
        for el, v in rs.items():
            ph = v.get("phase") if isinstance(v, dict) else None
            if not isinstance(ph, str):
                continue
            m = REF_ONE_ELEMENT.match(ph)
            if REF_UNIT.match(ph) or (m and m.group(1).upper() != el.upper() and m.group(1).upper() in (r.get("elements") or [])):
                drop.append(el)
        if not drop:
            continue
        r.setdefault("reference_state_stored", {k: v for k, v in rs.items()})
        for el in drop:
            G.correct(r, GATE, "reference_state." + el, rs[el], None, "reference_state_text_not_a_phase",
                      {"printed": rs[el].get("phase")})
            del rs[el]
        if not rs:
            r.pop("reference_state")
        rep.add(GATE, "reference_state_text_not_a_phase", r, "corrected")


# ------------------------------------------------------------------------------------------------- driver

def apply(records, visual, rep):
    cls = lambda k: visual.get("r2_class:" + k)   # noqa: E731
    r2_table_basis(records, visual, rep)
    r2_table_sign(records, visual, rep)
    r2_table_property(records, visual, rep)
    r2_sigma_with_divisor(records, rep, cls("printed_sigma_rescaled_with_value"))
    r2_expression_terms(records, rep, cls("expression_coefficient"))
    r2_absolute_entropy(records, rep, cls("absolute_entropy_labelled_formation"))
    r2_heat_content(records, rep, cls("heat_content_row_temperature"))
    r2_t298(records, rep, cls("quantity_printed_at_298"))
    r2_temperature_in_label(records, rep, cls("temperature_printed_in_quantity_label"))
    r2_standard_temperature(records, rep, cls("standard_enthalpy_at_calorimeter_temperature"))
    r2_experimental_zero_kelvin(records, rep, cls("experimental_value_from_0K_caption"))
    r2_dft_class(records, rep, visual)
    r2_reference_compounds(records, rep, cls("formation_from_compounds"))
    r2_noncrystalline(records, rep)
    r2_basis_lost(records, rep)
    r2_reference_state_text(records, rep)
    for r in records:
        for k in [k for k in r if k.startswith("_T_measurement")]:
            r.pop(k)
