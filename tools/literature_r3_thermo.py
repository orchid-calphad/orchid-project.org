#!/usr/bin/env python3
"""Audit round 3 (2026-10-24): residual thermochemical classes of the post-fix re-audit (POSTFIX.md, A3 / A4).

Called from tools/literature_gates.apply() right after tools/literature_r2_thermo.apply(). Same conventions as round
2: every rule keeps the stored field beside the new one and writes a `corrections` entry (gate "R3"); a rule that
changes a value or a state has its evidence in the printed table (caption, header, label) and was checked on the
rendered page for the table or for a sample of the class (tools/literature_visual_checks.json, keys
`table_unit:<asset>`, `record:<id>`, `r3_class:<rule>`).

Rules (codes in QUALITY_FIXES.md, section "Audit round 3"):

  R3-unitprefix  a value divided by 1000 (or 4184) because the stored unit lost its prefix ("k J/mole of atoms" read
                 as J; POSTFIX A4-N1, the Mg-B table): the caption / printed unit, with spaces removed, names kJ or
                 kcal. With a page check (table_unit:<asset>, verdict corrected) the value and the printed sigma are
                 multiplied back; without one the record is not accepted.
  R3-small       accepted formation enthalpies / Gibbs energies with 0 < |value| < 0.1 kJ/mol-atom need a page check
                 (record:<id>); unchecked -> not ESPEI-ready; checked and printed so -> kept, flagged.
  R3-emf0        EMF intercepts / slopes (E°, a, A; b of E = a + bT) stored as EMF values -> EMF_EXPRESSION_TERM.
  R3-expr3       three-term expression coefficients ("A of A + B·T + C·T lnT", "...: A", "T ln(T) coefficient") ->
                 <P>_EXPRESSION_TERM (the round-2 rule knew only the two-term A + B T forms).
  R3-contrib     contribution terms (electronic, phonon / vibrational, magnetic part of S, Cp, H, G) stored as total
                 properties -> <P>_CONTRIBUTION with `contribution` (not ESPEI-ready).
  R3-hbase       heat contents relative to a base temperature other than 298.15 K (HT - H700, Hm(373K) - Hm(273K))
                 -> HEAT_CONTENT_FROM_BASE_T with `heat_content_base_temperature_K` (not ESPEI-ready, not paired).
  R3-tprinted    a temperature printed for the row (temperature_printed) but not stored: one number -> T; two numbers
                 -> temperature_range_K (no single T).
  R3-dftclass    model-class rows whose own method / quantity label names a first-principles method (GGA, Wien2k,
                 FLAPW, "ab initio", "(DFT)"), or that sit in a table whose caption names only a first-principles
                 calculation (no CALPHAD / assessed / optimised / Miedema / estimate), -> evidence class dft
                 (POSTFIX A6-N6; round 2 handled formation rows of mixed tables only).
"""
from __future__ import annotations

import re

import build_dft_parity as parity
import literature_gates as G

GATE = "R3"
RELEASE = "literature-r3-thermo 2026-10-24"


def _ns(s):
    return re.sub(r"\s+", "", str(s or "")).lower()


def _asset(r):
    loc = r.get("locator") or {}
    return loc.get("asset") or loc.get("ref")


def _num(x):
    try:
        return float(str(x).replace("−", "-").replace(",", "."))
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------------------------------------- R3-unitprefix

PREFIXED = re.compile(r"k(?:j|cal)")
BARE_UNIT = re.compile(r"^\s*(?:j|cal)\b|^\s*(?:j|cal)\s*[/·.(]", re.I)


def unit_prefix_lost(r):
    """Factor the stored value lost (1000) when the printed unit lost its k prefix, or None."""
    v, vp = r.get("value"), _num(r.get("value_printed"))
    if not isinstance(v, (int, float)) or isinstance(v, bool) or not vp:
        return None
    if abs(vp) >= 1000 or abs(v) >= 0.2:
        return None   # a value printed in J (-29954) or one whose J reading is plausible is not this class
    up = str(r.get("unit_printed") or "")
    if not BARE_UNIT.search(up):
        return None
    txt = _ns((r.get("locator") or {}).get("caption"))
    if not PREFIXED.search(txt):
        return None
    # the caption must not also print a bare J / cal unit (a table mixing J and kJ columns)
    if re.search(r"(?<![a-z])(?:j|cal)(?:/|mol|\.)", re.sub(r"k(?:j|cal)", "", txt)):
        return None
    ratio = abs(v / vp)
    for f in (1e-3, 4.184e-3):
        for n in range(1, 13):
            if abs(ratio / (f / n) - 1.0) <= 0.01:
                return 1000.0
    return None


def r3_unit_prefix(records, visual, rep):
    for r in records:
        f = unit_prefix_lost(r)
        if not f:
            continue
        chk = visual.get("table_unit:%s" % _asset(r))
        if chk and chk.get("verdict") == "corrected":
            v = r["value"]
            new = round(v * f, 10)
            r.setdefault("value_stored", v)
            G.correct(r, GATE, "value", v, new, "unit_prefix_lost_k",
                      {"unit_printed": r.get("unit_printed"), "caption": ((r.get("locator") or {}).get("caption") or "")[:120],
                       "unit_on_page": chk.get("unit_printed")}, chk)
            r["value"] = new
            if chk.get("unit_printed"):
                r["unit_printed_stored"] = r.get("unit_printed")
                r["unit_printed"] = chk["unit_printed"]
            u = r.get("uncertainty")
            if isinstance(u, dict) and isinstance(u.get("value"), (int, float)):
                r.setdefault("uncertainty_stored", dict(u))
                G.correct(r, GATE, "uncertainty.value", u["value"], round(u["value"] * f, 10), "unit_prefix_lost_k", {}, chk)
                u["value"] = round(u["value"] * f, 10)
                u["rescaled_with_value"] = True
            if r["flags"].get("magnitude_suspect"):
                r["flags"]["magnitude_suspect_before_page_fix"] = True
                r["flags"]["magnitude_suspect"] = False
            rep.add(GATE, "unit_prefix_lost_k", r, "corrected", before=v, after=new)
        else:
            G.fail(r, GATE, "unit_prefix_lost_suspect", "not_accepted",
                   "value stored as printed/1000 while the caption prints a kJ / kcal unit (%r)" % r.get("unit_printed"), chk)
            rep.add(GATE, "unit_prefix_lost_suspect", r, "demoted")


# ------------------------------------------------------------------------------------------------- R3-small

def r3_small_formation(records, visual, rep):
    for r in records:
        if r.get("property") not in ("HM_FORM", "GM_FORM") or not r["flags"].get("extraction_accepted"):
            continue
        v = r.get("value")
        if not isinstance(v, (int, float)) or isinstance(v, bool) or not (0 < abs(v) < 0.1):
            continue
        if not str(r.get("unit") or "").startswith("kJ"):
            continue
        chk = visual.get("record:" + r["id"])
        if chk and chk.get("verdict") == "as_stored":
            if chk.get("misprint"):
                G.fail(r, GATE, "formation_value_printed_misprint", "not_espei_ready",
                       "printed so, but not physical (source misprint): %s" % chk.get("seen", "")[:160], chk)
                rep.add(GATE, "formation_value_printed_misprint", r, "flagged")
            continue
        G.fail(r, GATE, "formation_below_0p1_kJ_not_page_checked", "not_espei_ready",
               "|%s| = %g %s < 0.1: a lost unit prefix is possible; needs a page check" % (r["property"], v, r.get("unit")))
        rep.add(GATE, "formation_below_0p1_kJ_not_page_checked", r, "flagged")


# ------------------------------------------------------------------------------------------------- R3-emf0 / expr3

EMF_INTERCEPT = re.compile(r"^\s*(?:E\s*[°º⁰o0]|E\^?\{?\\?circ\}?|E_?\{?0\}?|E\s*\(0\)|a|A)\s*(?:\([^)]*\))?\s*$")
EMF_SLOPE = re.compile(r"^\s*(?:b|B|dE\s*/\s*dT|\(dE/dT\).*|E\s*=\s*a\s*\+\s*bT\s*:\s*b)\s*$")
EMF_INTERCEPT2 = re.compile(r"^\s*E\s*=\s*a\s*\+\s*bT\s*:\s*a\s*$")
EXPR3_A = re.compile(r"\bA of A\s*\+\s*B|:\s*A\s*$", re.I)
EXPR3_B = re.compile(r"\bB of A\s*\+\s*B|:\s*B\s*$", re.I)
EXPR3_C = re.compile(r"\bC of A\s*\+\s*B|T\s*ln\s*\(?T\)?\s*coefficient|:\s*C\s*$", re.I)


def _relabel_term(r, term, form, rule, rep, check=None):
    p = r["property"]
    r["property_code_before_gates"] = r.get("property_code_before_gates") or p
    new = p + "_EXPRESSION_TERM"
    G.correct(r, GATE, "property", p, new, rule, {"printed": r.get("property_printed"), "term": term, "form": form}, check)
    r["property"] = new
    r["property_group"] = "model expression"
    r["property_expression"] = {"form": form, "term": term, "printed": r.get("property_printed"),
                                "temperature_range_K": r.get("temperature_range_K")}
    T = r.get("temperature_K")
    if T is not None:
        r["temperature_K_stored"] = T
        r["temperature_K"] = None
        G.correct(r, GATE, "temperature_K", T, None, rule + "_no_single_temperature",
                  {"why": "a coefficient of an expression in T has no single temperature"}, check)
    G.fail(r, GATE, "expression_coefficient_not_a_value", "not_espei_ready", "coefficient %s of %s" % (term, form), check)
    rep.add(GATE, rule, r, "corrected")


def r3_expression_terms(records, visual, rep):
    cls = visual.get("r3_class:emf_intercept_stored_as_emf")
    cls3 = visual.get("r3_class:three_term_expression_coefficient")
    for r in records:
        p = r.get("property") or ""
        if p.endswith("_EXPRESSION_TERM"):
            continue
        lab = str(r.get("property_printed") or "")
        if p == "EMF":
            if EMF_INTERCEPT.match(lab) or EMF_INTERCEPT2.match(lab):
                _relabel_term(r, "A", "E = A + B T", "emf_intercept_stored_as_emf", rep, cls)
            elif EMF_SLOPE.match(lab):
                _relabel_term(r, "B", "E = A + B T", "emf_slope_stored_as_emf", rep, cls)
            continue
        if p in ("GIBBS_ENERGY", "GM_FORM", "GM_MIX", "G_MIX_PARTIAL", "CHEMICAL_POTENTIAL", "GM_REACTION", "HM_FORM",
                 "SM_FORM"):
            for rx, term in ((EXPR3_C, "C"), (EXPR3_B, "B"), (EXPR3_A, "A")):
                if rx.search(lab):
                    form = "A + B T + C T ln T" if re.search(r"ln\s*\(?T", lab) else (
                        "printed expression: " + lab.split(":")[0].strip() if ":" in lab else "A + B T")
                    _relabel_term(r, term, form, "three_term_expression_coefficient", rep, cls3)
                    break


# ------------------------------------------------------------------------------------------------- R3-contrib

CONTRIB = [(re.compile(r"\(\s*el\.?\s*\)|electronic", re.I), "electronic"),
           (re.compile(r"\(\s*ph\.?\s*\)|phonon", re.I), "phonon"),
           (re.compile(r"\(\s*vib\.?\s*\)|vibrational", re.I), "vibrational"),
           (re.compile(r"\(\s*mag\.?\s*\)|magnetic contribution|^\s*mag|\bp\s*mag\b|_\{?mag\}?|\bmag\.?\s*$", re.I), "magnetic"),
           (re.compile(r"\(\s*conf\.?\s*\)|configurational contribution", re.I), "configurational")]
CONTRIB_PROPS = {"SM", "S_ABS", "ENTROPY", "CPM", "HEAT_CONTENT", "GIBBS_ENERGY", "HM_FORM", "SM_FORM", "GM_FORM"}


def r3_contributions(records, visual, rep):
    chk = visual.get("r3_class:contribution_stored_as_total")
    for r in records:
        p = r.get("property") or ""
        if p not in CONTRIB_PROPS:
            continue
        lab = str(r.get("property_printed") or "")
        if re.search(r"ferromagnetic|non[- ]?magnetic|paramagnetic", lab, re.I):
            continue   # a total value computed in a magnetic state, not a contribution
        kind = next((k for rx, k in CONTRIB if rx.search(lab)), None)
        if not kind or (kind == "magnetic" and p in ("SM", "ENTROPY", "S_ABS", "SM_FORM")):
            continue   # magnetic entropies are SM_MAGNETIC (gate G5)
        new = p + "_CONTRIBUTION"
        r["property_code_before_gates"] = r.get("property_code_before_gates") or p
        G.correct(r, GATE, "property", p, new, "contribution_stored_as_total", {"printed": lab, "contribution": kind}, chk)
        r["property"] = new
        r["contribution"] = kind
        G.fail(r, GATE, "contribution_not_a_total", "not_espei_ready", "%s contribution to %s" % (kind, p), chk)
        rep.add(GATE, "contribution_stored_as_total", r, "corrected", contribution=kind)


# ------------------------------------------------------------------------------------------------- R3-hbase

HBASE = re.compile(r"H\s*[°o]?\s*_?\{?\s*m?\s*_?\(?\s*T\s*\)?\}?\s*[-−–]\s*H\s*[°o]?\s*_?\{?\s*m?\s*\(?\s*(\d{3,4}(?:\.\d+)?)\s*K?\s*\)?\}?",
                   re.I)
HBASE2 = re.compile(r"H\s*m?\s*\(\s*(\d{3,4}(?:\.\d+)?)\s*K\s*\)\s*[-−–]\s*H\s*m?\s*\(\s*(\d{3,4}(?:\.\d+)?)\s*K\s*\)", re.I)


def r3_heat_content_base(records, visual, rep):
    chk = visual.get("r3_class:heat_content_base_temperature")
    for r in records:
        if r.get("property") != "HEAT_CONTENT":
            continue
        lab = str(r.get("property_printed") or "")
        T1 = None
        m2 = HBASE2.search(lab)
        if m2:
            T1, base = float(m2.group(1)), float(m2.group(2))
        else:
            m = HBASE.search(lab)
            if not m:
                continue
            base = float(m.group(1))
        if abs(base - 298.15) <= 1.0 or base == 0:
            continue
        r["property_code_before_gates"] = r.get("property_code_before_gates") or "HEAT_CONTENT"
        G.correct(r, GATE, "property", "HEAT_CONTENT", "HEAT_CONTENT_FROM_BASE_T", "heat_content_base_temperature",
                  {"printed": lab, "base_temperature_K": base}, chk)
        r["property"] = "HEAT_CONTENT_FROM_BASE_T"
        r["heat_content_base_temperature_K"] = base
        if T1 is not None and r.get("temperature_K") is None:
            r["temperature_K"] = T1
            r["temperature_basis"] = "quantity_label"
            G.correct(r, GATE, "temperature_K", None, T1, "temperature_printed_in_quantity_label", {"printed": lab}, chk)
        G.fail(r, GATE, "heat_content_base_not_298", "not_espei_ready",
               "H(T) - H(%g K): the base temperature is not 298.15 K" % base, chk)
        rep.add(GATE, "heat_content_base_temperature", r, "corrected", base=base)


# ------------------------------------------------------------------------------------------------- R3-tprinted

TNUM = re.compile(r"(\d{3,4}(?:\.\d+)?)\s*(K|°\s*C|ºC|C)?(?![\d.])")


def r3_temperature_printed(records, visual, rep):
    chk = visual.get("r3_class:temperature_printed_not_stored")
    for r in records:
        if G.temperature_withdrawn_on_page(r):
            continue
        if r.get("temperature_K") is not None or r.get("value") is None or r.get("temperature_range_K"):
            continue
        p = r.get("property") or ""
        if p.endswith("_EXPRESSION_TERM") or p in ("INVARIANT_TEMPERATURE", "TRANSITION_TEMPERATURE"):
            continue
        tp = r.get("temperature_printed")
        if not isinstance(tp, str) or not re.search(r"\d{3}", tp) or re.search(r"[<>≤≥~±]|to\b|[-–]\s*\d", tp):
            continue
        nums = TNUM.findall(tp)
        ts = []
        for x, u in nums:
            t = float(x)
            ts.append(t + 273.15 if u and "C" in u.upper() and "K" not in u.upper() else t)
        if not ts:
            continue
        if len(ts) == 1:
            T = ts[0]
            r["temperature_K"] = T
            r["temperature_basis"] = "row_temperature_printed"
            G.correct(r, GATE, "temperature_K", None, T, "temperature_printed_not_stored", {"printed": tp}, chk)
            rep.add(GATE, "temperature_printed_not_stored", r, "corrected", after=T)
        elif len(ts) == 2:
            rng = {"min": min(ts), "max": max(ts)}
            r["temperature_range_K"] = rng
            r["temperature_basis"] = "row_temperature_range_printed"
            G.correct(r, GATE, "temperature_range_K", None, rng, "temperature_range_printed_not_stored", {"printed": tp}, chk)
            G.fail(r, GATE, "no_single_temperature", "not_espei_ready",
                   "printed for a temperature range %g-%g K" % (rng["min"], rng["max"]), chk)
            rep.add(GATE, "temperature_range_printed_not_stored", r, "corrected", after=rng)


# ------------------------------------------------------------------------------------------------- R3-dftclass

DFT_ROW = re.compile(r"first[- ]principles?|1st[- ]princ|ab[- ]?initio|\bDFT\b|density functional|VASP|\bGGA\b|\bPBE\b|"
                     r"\bLDA\b|\bPAW\b|FLAPW|LMTO|WIEN2k|CASTEP|Quantum ESPRESSO", re.I)
NOT_DFT_ROW = re.compile(r"CALPHAD|assess|optimi[sz]|\bopt\b|thermo-?calc|pandat|miedema|estimat|evaluat|\bexp|"
                         r"measur|calorim|fit\b|model(?!ling with)", re.I)
NOT_DFT_CAPTION = re.compile(r"CALPHAD|assess|optimi[sz]|thermo-?calc|pandat|miedema|estimat|evaluat|present (?:work|model|"
                             r"calculation)|this work|thermodynamic (?:model|description|calculation)", re.I)


def r3_dft_class(records, visual, rep):
    chk = visual.get("r3_class:first_principles_value_coded_model")
    for r in records:
        if r.get("evidence_class") != "model":
            continue
        cap = str((r.get("locator") or {}).get("caption") or "")
        row = str(r.get("method") or "") + " " + str(r.get("property_printed") or "")
        if DFT_ROW.search(row) and not NOT_DFT_ROW.search(row):
            basis = "row_names_first_principles_method"
        elif DFT_ROW.search(cap) and not NOT_DFT_CAPTION.search(cap) and not NOT_DFT_ROW.search(row) \
                and not (r.get("cited_sources") and not DFT_ROW.search(row)):
            basis = "caption_names_only_a_first_principles_calculation"
        else:
            continue
        r["evidence_class_before_gates"] = r.get("evidence_class_before_gates") or "model"
        G.correct(r, GATE, "evidence_class", "model", "dft", "first_principles_value_coded_model",
                  {"basis": basis, "caption": cap[:160], "method": r.get("method"), "printed": r.get("property_printed")}, chk)
        r["evidence_class"] = "dft"
        r.setdefault("origin", {})["class_gate"] = "dft"
        rep.add(GATE, "first_principles_value_coded_model", r, "corrected", basis=basis)


# ------------------------------------------------------------------------------------------------- R3-defect

DEFECT = re.compile(r"/\s*defects?\b|per defect|defect formation|formation (?:energy|enthalpy) of (?:the )?(?:point )?defects?", re.I)


def r3_defect_energies(records, visual, rep):
    """Point-defect formation energies (eV/defect, kJ/mol/defects: antisites in PbTe, HfNiSn) stored as compound
    formation enthalpies -> DEFECT_FORMATION_ENERGY (not ESPEI-ready; page-checked tables, round 3)."""
    for r in records:
        p = r.get("property") or ""
        if p not in ("HM_FORM", "GM_FORM", "HM_MIX"):
            continue
        txt = str(r.get("unit_printed") or "") + " " + str(r.get("unit") or "") + " | " + \
            str((r.get("locator") or {}).get("caption") or "")
        if not DEFECT.search(txt):
            continue
        chk = visual.get("r3_table:defect_energy:%s" % _asset(r))
        r["property_code_before_gates"] = r.get("property_code_before_gates") or p
        G.correct(r, GATE, "property", p, "DEFECT_FORMATION_ENERGY", "defect_energy_stored_as_formation",
                  {"unit_printed": r.get("unit_printed"), "caption": ((r.get("locator") or {}).get("caption") or "")[:120]}, chk)
        r["property"] = "DEFECT_FORMATION_ENERGY"
        G.fail(r, GATE, "defect_energy_not_a_formation_property", "not_espei_ready",
               "point-defect formation energy, not the compound's formation enthalpy", chk)
        rep.add(GATE, "defect_energy_stored_as_formation", r, "corrected")


# ------------------------------------------------------------------------------------------------- R3-comma

COMMA_FORMULA = re.compile(r"[A-Z][a-z]?\d+,\d+")
PER_ATOM_PROPS = {"HM_FORM", "GM_FORM", "SM_FORM", "SM", "S_ABS", "CPM", "ENTROPY", "HEAT_CONTENT"}


def r3_comma_decimal(records, visual, rep):
    """A per-formula value of a comma-decimal formula (PuO1,52 = PuO1.52, 2.52 atoms) divided by the atom count of
    the formula read without its decimals (PuO: 2) -> printed value / 2.52 (POSTFIX A3-N6; page-checked)."""
    chk = visual.get("r3_class:comma_decimal_formula")
    for r in records:
        if r.get("property") not in PER_ATOM_PROPS:
            continue
        labs = [p.get("label") for p in r.get("phases") or [] if p.get("label") and COMMA_FORMULA.search(p.get("label"))]
        if len(labs) != 1:
            continue
        v, vp = r.get("value"), _num(r.get("value_printed"))
        if not isinstance(v, (int, float)) or not vp or not v:
            continue
        c, _ = parity.parse_phase(labs[0], None)
        if not c:
            continue
        n = float(sum(c.values()))
        ratio = vp / v
        if abs(ratio - n) <= 0.01 * n:
            continue   # already per atom of the decimal formula
        k = round(ratio)
        if k < 2 or abs(ratio - k) > 0.01 * k:
            continue
        new = round(vp / n, 6)
        r.setdefault("value_stored", v)
        G.correct(r, GATE, "value", v, new, "comma_decimal_formula_atoms",
                  {"label": labs[0], "atoms_per_formula": round(n, 4), "divided_before_by": k, "printed": vp}, chk)
        r["value"] = new
        u = r.get("uncertainty")
        if isinstance(u, dict) and isinstance(u.get("value"), (int, float)):
            r.setdefault("uncertainty_stored", dict(u))
            u["value"] = round(u["value"] * k / n, 6)
            u["rescaled_with_value"] = True
        r["basis_correction"] = {"atoms_per_formula_unit": round(n, 4), "operation": "divide", "gate": GATE}
        rep.add(GATE, "comma_decimal_formula_atoms", r, "corrected", before=v, after=new)


# ------------------------------------------------------------------------------------------------- R3-statephases
# Every released value is an observation of an equilibrium state (phases, T, P, composition). Running-text invariant /
# transition temperatures and formation values whose phases were not extracted (coordinator, audit round 3): the
# phases are recovered from the record's own verbatim quote (A6 checked quotes against their pages) only when the
# quote prints exactly one reaction and one temperature equal to the value (or, for a formation value, exactly one
# compound formula of the system; for a congruent melting / magnetic transition, exactly one compound); otherwise the
# record is not accepted (`state_phases_missing`).
ARROW = r"(?:<=>|<->|<-->|->|-->|⇌|⇄|↔|⟷|→|⟶|⇒|=)"
TOK = r"(?:\(?[^\s+=↔⇌⇄⟷→⟶⇒<>,;:]+\)?(?:\s\((?:[^()\s]{1,12})\))?)"
SIDE = TOK + r"(?:\s*\+\s*" + TOK + r")*"
RX = re.compile(r"(?:^|(?<=[\s(:\[]))(" + SIDE + r")\s*" + ARROW + r"\s*(" + SIDE + r")(?=\s*(?:$|[,.;:)\]]|\s(?:at|occur\w*|is|was|takes|in|with|and|has|were|are|by|of|from|reaction|temperature)\b))")
ANY_ARROW = re.compile(ARROW)
LIQ = re.compile(r"^(?:L|Liq\.?|liq\.?|Liquid|liquid|LIQUID|L\d|Liquid#?\d|liquid#?\d|L'|L′|melt)$")
STOP = {"at", "is", "the", "and", "of", "in", "a", "to", "was", "were", "which", "with", "as", "by", "for", "from"}
ELEMS = set("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf".split())
GREEK = "αβγδεζηθικλμνξπρστυφχψωΓΔΘΛΞΠΣΦΨΩ"
TEMP = re.compile(r"(?<![\d.])(\d{2,4}(?:\.\d+)?)\s*(?:±\s*\d+(?:\.\d+)?\s*)?(K|°\s*C|◦\s*C|ºC|8C|\s?C\b)")


import literature_r2_compositions as r2c
GREEK_LABEL = re.compile(r"^[%s][′'’]?\d?[′'’]?$" % GREEK)


def normalise(q):
    q = q.replace("␣", "α").replace("␤", "β").replace("␥", "γ").replace("␦", "δ").replace("þ", "+")
    for _ in range(3):
        q = re.sub(r"([A-Z][a-z]?\d+(?:\.\d+)?)\s+(?=[A-Z][a-z]?\d)", r"\1", q)   # Al3 Ni2 -> Al3Ni2
        q = re.sub(r"([A-Z][a-z]?\d+)\s+(?=[A-Z][a-z]?(?:\b|[A-Z]))", r"\1", q)       # Co17 Nd2 handled above; Al7 Cr
    q = re.sub(r"([%s])\s+(?=[A-Z][a-z]?)" % GREEK, r"\1", q)                        # α ErSb -> αErSb
    return q


def phase_ok(t, sys_els):
    t = t.strip().strip(".")
    if not t or t.lower() in STOP or len(t) > 40 or t.count("(") != t.count(")") or "[" in t or "]" in t:
        return False
    if LIQ.match(t) or GREEK_LABEL.match(t):
        return True
    els = sorted(sys_els)
    f = r2c.formula_fractions(re.sub("^[%s]-?" % GREEK, "", t), els)
    if f and len(f) >= 2:
        return True
    return bool(r2c.pure_element(t, els))


def recover_reaction(r):
    """(phases, why) recovered from the quote: one printed reaction and one printed temperature equal to the value."""
    q = normalise(str(r.get("quote") or ""))
    els = set((r.get("source_system") or r.get("system") or "").split("-")) | set(r.get("elements") or [])
    reacts = [m for m in RX.finditer(q)]
    good = []
    for m in reacts:
        sides = [[x.strip() for x in re.split(r"\s*\+\s*", m.group(i))] for i in (1, 2)]
        if all(sides[0]) and all(sides[1]) and all(phase_ok(t, els) for t in sides[0] + sides[1]) \
                and len(sides[0]) + len(sides[1]) >= 2:
            good.append(sides[0] + sides[1])
    if re.search(r"\s\?\s|[A-Za-z0-9)]\?\s", q):
        return None, "garbled_reaction_sign"
    if len(good) != 1 or len(ANY_ARROW.findall(q)) != 1:
        return None, "reactions_in_quote:%d" % len(ANY_ARROW.findall(q))
    temps = TEMP.findall(q)
    nums = re.findall(r"(?<![\d.])\d{3,4}(?:\.\d+)?(?![\d.]|\s*(?:at|wt|mol)\.?\s*%)", q)
    if len(temps) != 1 or len(nums) != 1:
        return None, "temperatures_in_quote:%d" % len(temps)
    x, u = temps[0]
    T = float(x) + (273.15 if "C" in u else 0.0)
    v = r.get("value")
    if not isinstance(v, (int, float)) or abs(T - v) > 1.0:
        return None, "temperature_not_the_value"
    out = []
    for t in good[0]:
        t = t.strip().strip(".")
        if t not in out:
            out.append(t)
    return out, "one_reaction_one_temperature"


def recover_compound(r):
    """One compound formula of the record's system printed in the quote of a formation value."""
    q = normalise(str(r.get("quote") or ""))
    q = re.sub(r"([A-Z][a-z]?\d*\.\d+)\s+(?=[A-Z][a-z]?\d*\.\d+)", r"\1", q)   # Fe0.52 Sb0.48
    els = [e.capitalize() for e in (r.get("source_system") or r.get("system") or "").split("-")]
    found = {}
    for tok in re.findall(r"[A-Z][A-Za-z0-9.()₀-₉]*", q):
        c, _ = parity.parse_phase(tok, els)
        if c and len(c) >= 2:
            found[parity.reduced_key(c)] = tok
    if len(found) != 1:
        return None, "compounds_in_quote:%d" % len(found)
    return [next(iter(found.values()))], "one_compound_in_quote"


STATE_PHASE_PROPS = {"INVARIANT_TEMPERATURE": "reaction", "TRANSITION_TEMPERATURE": "reaction",
                     "PHASE_COMPOSITION": None, "PHASE_ASSEMBLAGE": None,
                     "HM_FORM": "compound", "GM_FORM": "compound", "SM_FORM": "compound"}


def recover_single_compound(r):
    tk = str(r.get("transition_kind") or "")
    q = str(r.get("quote") or "")
    if not (tk == "magnetic" or (tk in ("congruent transformation", "congruent melting") and re.search(r"melt", q, re.I))):
        return None, "not_a_single_phase_transition"
    ph, why = recover_compound(r)
    if not ph:
        return None, why
    return (ph if tk == "magnetic" else ph + ["liquid"]), ("one_compound_magnetic" if tk == "magnetic" else "one_compound_congruent_melting")


def r3_state_phases(records, visual, rep):
    chk = visual.get("r3_class:state_phases_from_quote")
    for r in records:
        p = r.get("property") or ""
        if p not in STATE_PHASE_PROPS or has_real_phases(r):
            continue
        old_phases = [x for x in G.labels(r) if x is not None]
        ph, why = None, "no phases extracted" if not r.get("phases") else "placeholder phases only (%s)" % ", ".join(
            str(x) for x in G.labels(r))
        pchk = visual.get("state_phases:" + r["id"])
        if pchk and pchk.get("verdict") == "corrected" and pchk.get("phases"):
            r["phases"] = [{"label": x, "label_source": "page"} for x in pchk["phases"]]
            if pchk.get("reaction") and isinstance(r.get("reaction"), dict) and not r["reaction"].get("printed"):
                r["reaction"]["printed_from_page"] = pchk["reaction"]
            G.correct(r, GATE, "phases", old_phases, list(pchk["phases"]), "state_phases_read_from_page",
                      {"reaction": pchk.get("reaction")}, pchk)
            G.fail(r, GATE, "phases_read_from_page", "not_espei_ready",
                   "phases read on the rendered page; not mapped to model phases", pchk)
            rep.add(GATE, "state_phases_read_from_page", r, "corrected")
            continue
        if r.get("record_kind") == "running_text" and r.get("quote"):
            if STATE_PHASE_PROPS[p] == "reaction":
                ph, why = recover_reaction(r)
                if not ph and p == "TRANSITION_TEMPERATURE":
                    ph, why = recover_single_compound(r)
            elif STATE_PHASE_PROPS[p] == "compound":
                ph, why = recover_compound(r)
        if ph:
            r["phases"] = [{"label": x, "label_source": "quote"} for x in ph]
            G.correct(r, GATE, "phases", old_phases, [x for x in ph], "state_phases_from_quote",
                      {"basis": why, "quote_page": r.get("quote_page")}, chk)
            G.fail(r, GATE, "phases_recovered_from_quote", "not_espei_ready",
                   "phases read from the quote (%s); not mapped to model phases" % why, chk)
            rep.add(GATE, "state_phases_from_quote", r, "corrected", basis=why)
        else:
            G.fail(r, GATE, "state_phases_missing", "not_accepted",
                   "%s without the phases of its state (%s)" % (p, why))
            rep.add(GATE, "state_phases_missing", r, "demoted", reason=why.split(":")[0])


# "na" is deliberately absent: it is the element sodium (Na) in a lower-cased comparison
PLACEHOLDER_PHASE_LABELS = {"", "unknown", "other", "unspecified", "n/a", "none", "?", "-"}


def has_real_phases(r):
    """True when the record names at least one phase that is not a placeholder (audit round 4: R3 tested only for an
    empty list, so ['UNKNOWN'], ['other'] and ['unspecified'] passed as 'phases recorded')."""
    return any(str(x).strip().lower() not in PLACEHOLDER_PHASE_LABELS for x in G.labels(r) if x is not None)


def r3_assessed_database_class(records, visual, rep):
    """Values read from an assessed thermodynamic database publication (tools/pdf_identity.json `values_class`, e.g.
    the COST 507 book's calculated invariant tables) whose evidence class was left unspecified -> model."""
    import json as _json
    import os as _os
    ident = {e["key"]: e for e in _json.load(open(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)),
                                                                  "pdf_identity.json"), encoding="utf-8"))["entries"]}
    for r in records:
        asset = str(_asset(r) or "")
        e = ident.get(asset.split("/")[0])
        if not e or not e.get("values_class") or r.get("evidence_class") != "unspecified":
            continue
        chk = visual.get("r3_class:assessed_database_values_unspecified")
        r["evidence_class_before_gates"] = r.get("evidence_class_before_gates") or "unspecified"
        G.correct(r, GATE, "evidence_class", "unspecified", e["values_class"], "assessed_database_values",
                  {"source": e["citation"][:120], "why": e.get("values_class_basis")}, chk)
        r["evidence_class"] = e["values_class"]
        r.setdefault("origin", {})["class_gate"] = e["values_class"]
        rep.add(GATE, "assessed_database_values", r, "corrected")


def r3_printed_ranges(records, visual, rep):
    """A temperature printed as a range (1223-1287 K) stored as its lower bound (page-checked, ledger key range:<id>):
    the range is kept in temperature_range_K, the value stays the printed lower bound with value_basis
    printed_range_lower_bound, and the record is not ESPEI-ready (no single temperature)."""
    for r in records:
        chk = visual.get("range:" + r["id"])
        if not chk or chk.get("verdict") != "corrected" or not chk.get("range_K"):
            continue
        rng = {"min": float(min(chk["range_K"])), "max": float(max(chk["range_K"]))}
        G.correct(r, GATE, "temperature_range_K", r.get("temperature_range_K"), rng, "printed_range_read_from_page", {}, chk)
        r["temperature_range_K"] = rng
        r["value_basis"] = "printed_range_lower_bound"
        G.fail(r, GATE, "printed_as_temperature_range", "not_espei_ready",
               "the page prints %g-%g K; the value is the lower bound" % (rng["min"], rng["max"]), chk)
        rep.add(GATE, "printed_range_read_from_page", r, "corrected")


def apply(records, visual, rep):
    r3_printed_ranges(records, visual, rep)
    r3_assessed_database_class(records, visual, rep)
    r3_comma_decimal(records, visual, rep)
    r3_defect_energies(records, visual, rep)
    r3_unit_prefix(records, visual, rep)
    r3_expression_terms(records, visual, rep)
    r3_contributions(records, visual, rep)
    r3_heat_content_base(records, visual, rep)
    r3_temperature_printed(records, visual, rep)
    r3_dft_class(records, visual, rep)


def apply_late(records, visual, rep):
    """After the magnitude gates (G2): page-check requirement for tiny formation values; phases of every state."""
    r3_small_formation(records, visual, rep)
    r3_state_phases(records, visual, rep)
