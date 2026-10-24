#!/usr/bin/env python3
"""Uncertainty budget of every ORCHID literature value (PI decision 2026-10-24: method-based uncertainty floors,
components kept separate).

Printed uncertainties often reflect statistical precision only (e.g. Knudsen-effusion mass-spectrometry formation
enthalpies printed with +-0.08 kJ/mol-atom that disagree with calorimetry by 6-8 kJ/mol-atom). Every value of the
literature release therefore carries `uncertainty_budget`, one entry per measured quantity of the record
(`value`; `composition` = the printed phase compositions; `temperature_K` = the equilibration temperature of a
phase-composition / assemblage observation), with the components kept separate, in the quantity's unit:

  sigma_printed      as printed by the source (+-x, range half-width; coverage as far as the record states it)
  sigma_reading      reading resolution (half of the last printed digit) where ORCHID recorded it
  sigma_method       method-class systematic component: the robust between-study scatter of independent
                     measurements of the same state by the same method class in ORCHID (calibrated), or a
                     documented prior when fewer than MIN_PAIRS pairs exist; DFT: the calibrated sigma_DFT;
                     model-calculated values: the robust deviation of the model class from experiment
  sigma_composition  attribution / composition offset (a value measured at a historical formula, |dx| x slope)
  sigma_scatter      between-study scatter of a page-checked mutually inconsistent group (group, n)
  sigma_conversion   reference-state / unit conversion contribution (e.g. a lattice-stability conversion)
  sigma_effective    max( sqrt( max(printed, reading, method)^2 + composition^2 + conversion^2 ), scatter )

Rule (`SIGMA_RULE`): the method calibration is a between-study scatter, so it already contains each study's own
statistical error; the printed sigma and the method floor are therefore not added in quadrature (that would count
the statistical part twice) but the larger is taken. Composition offset and conversion are independent of the
measurement and are added in quadrature. A between-study scatter of a mutually inconsistent group replaces the
result when it is larger (it is itself a measure of the same disagreement the floor describes).

Calibration (`calibrate`): for each value type and method class, independent measurements of the same state from
different primary sources (origin.measured_in, else the single cited source, else the reporting paper for its own
work; one value per source and state, the median of its values; identical values from different sources are
counted once, as one measurement re-listed), all pairs within the same method class: robust sigma =
1.4826 x median|d| / sqrt(2) (pairs of zero-mean differences), 95 % percentile bootstrap over states (N_BOOT,
fixed seed). Floor = the robust sigma rounded up to the value type's step, never below the minimum (ESPEI's
experimental default sigma where ESPEI has one: 0.5 kJ/mol-atom, 0.2 J/K/mol-atom). Fewer than MIN_PAIRS pairs or
MIN_STATES states but at least MIN_PAIRS_CLEANED independent pairs (PI rule 2026-10-24): floor = max(documented prior,
cleaned robust sigma), basis cleaned_below_threshold (n, sigma and CI recorded). Fewer pairs: the documented prior
(PRIORS, basis prior). Printed values and printed uncertainties are never altered.

Standard library only; deterministic.

  python3 tools/uncertainty_budget.py [--data DIR]          # print the calibration summary for the published release
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import random
import re
import sys
from collections import Counter, defaultdict

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_dft_parity as parity  # noqa: E402
import literature_gates as gates  # noqa: E402

VERSION = "orchid-uncertainty-calibration 2026-10-24 (audit round 2)"
BUDGET_VERSION = "orchid-uncertainty-budget 2026-10-24"
CALIBRATION_FILE = "literature/UNCERTAINTY_CALIBRATION.json"
ESPEI_COPY = "espei/experimental/CALIBRATION_EXPERIMENT.json"
DFT_CALIBRATION_FILE = "espei/dft/CALIBRATION.json"
MIN_PAIRS = 20
MIN_STATES = 10            # a calibrated floor needs pairs from at least this many distinct states
MIN_PAIRS_CLEANED = 5      # PI rule (2026-10-24): below the calibration threshold but with >= 5 independent pairs, the
#                            floor is max(documented prior, cleaned robust sigma), basis cleaned_below_threshold
BASIS_RATIO_TOL = 0.03     # a pair whose ratio is within 3 % of the formula's atom count is a per-formula / per-atom
#                            basis error of one reading, not a method scatter: excluded from the calibration
N_BOOT = 1000
SEED = 20261024
R_GAS = 8.314462618
BLOCKING = ("unsigned_in_print", "magnitude_suspect", "system_mismatch", "source_disputed")
SIGMA_RULE = "max(sqrt(max(printed, reading, method)^2 + composition^2 + conversion^2), scatter)"
SIGMA_RULE_TEXT = (
    "sigma_effective = max( sqrt( max(sigma_printed, sigma_reading, sigma_method)^2 + sigma_composition^2 + "
    "sigma_conversion^2 ), sigma_scatter ); absent components drop out. The method floor is a between-study scatter "
    "and already contains each study's statistical error, so it is not added to the printed sigma in quadrature (no "
    "double counting): the larger of the two is used. The composition offset and the conversion are independent of "
    "the measurement and are added in quadrature. The scatter of a page-checked mutually inconsistent group replaces "
    "the result when larger. The reading resolution alone is not an uncertainty estimate: without a printed sigma, a "
    "method floor or a group scatter, sigma_effective is null (with the reason).")

# ------------------------------------------------------------------------------------------- value types --
# value type -> (canonical unit, kind, minimum floor, rounding step, label)
VT = {
    "HM_FORM": ("kJ/mol-atom", "energy", 0.5, 0.1, "formation enthalpy"),
    "GM_FORM": ("kJ/mol-atom", "energy", 0.5, 0.1, "formation Gibbs energy"),
    "SM_FORM": ("J/K/mol-atom", "entropy", 0.2, 0.1, "formation entropy"),
    "H_MIX_PARTIAL": ("kJ/mol", "energy", 0.5, 0.1, "mixing / partial / transformation enthalpy"),
    "G_MIX_PARTIAL": ("kJ/mol", "energy", 0.5, 0.1, "mixing / partial Gibbs energy, chemical potential"),
    "S_MIX_PARTIAL": ("J/K/mol", "entropy", 0.2, 0.1, "mixing / partial / excess entropy"),
    "S_ABS": ("J/K/mol", "entropy", 0.2, 0.1, "absolute entropy, entropy of fusion / transformation"),
    "CPM": ("J/K/mol", "entropy", 0.2, 0.1, "heat capacity"),
    "HEAT_CONTENT": ("kJ/mol", "energy", 0.5, 0.1, "heat content H(T)-H(298)"),
    "ACTIVITY": ("kJ/mol (chemical potential)", "activity", 0.5, 0.1, "activity (floor on mu = RT ln a)"),
    "EMF": ("mV", "emf", 0.1, 0.1, "electromotive force"),
    "T_INVARIANT": ("K", "temperature", 1.0, 1.0, "invariant-reaction temperature"),
    "T_TRANSITION": ("K", "temperature", 1.0, 1.0, "transition / liquidus / solidus temperature"),
    "T_CONDITION": ("K", "temperature", 1.0, 1.0, "equilibration temperature of a phase-equilibrium observation"),
    "X_PHASE": ("mole fraction", "composition", 0.001, 0.001, "phase composition (tie-line / boundary / invariant)"),
    "LATTICE_LENGTH": ("relative (fraction of the value)", "relative", 0.0002, 0.0001, "lattice parameter"),
    "LATTICE_ANGLE": ("degree", "angle", 0.01, 0.01, "lattice angle"),
    "VOLUME": ("relative (fraction of the value)", "relative", 0.0005, 0.0001, "cell / molar volume"),
}
PROP_VT = {
    "HM_FORM": "HM_FORM", "GM_FORM": "GM_FORM", "SM_FORM": "SM_FORM",
    "HM_MIX": "H_MIX_PARTIAL", "PARTIAL_ENTHALPY": "H_MIX_PARTIAL", "HM_LATTICE_STABILITY": "H_MIX_PARTIAL",
    "TRANSITION_ENTHALPY": "H_MIX_PARTIAL",
    "GM_MIX": "G_MIX_PARTIAL", "CHEMICAL_POTENTIAL": "G_MIX_PARTIAL", "GIBBS_ENERGY": "G_MIX_PARTIAL",
    "SM_MIX": "S_MIX_PARTIAL", "SM_EXCESS": "S_MIX_PARTIAL",
    "HM_FROM_COMPOUNDS": "H_MIX_PARTIAL", "GM_FROM_COMPOUNDS": "G_MIX_PARTIAL", "SM_FROM_COMPOUNDS": "S_MIX_PARTIAL",
    "SM": "S_ABS", "ENTROPY": "S_ABS", "SM_FUSION": "S_ABS", "SM_MAGNETIC": "S_ABS", "SM_INCREMENT": "S_ABS",
    "CPM": "CPM", "HEAT_CONTENT": "HEAT_CONTENT", "ACTIVITY": "ACTIVITY", "EMF": "EMF",
    "INVARIANT_TEMPERATURE": "T_INVARIANT", "TRANSITION_TEMPERATURE": "T_TRANSITION",
    "MOLAR_VOLUME": "VOLUME", "VOLUME": "VOLUME", "LATTICE_V": "VOLUME", "LATTICE_V0": "VOLUME",
    "LATTICE_V_ATOM": "VOLUME", "LATTICE_2V_ATOM": "VOLUME", "LATTICE_VOLUME": "VOLUME",
    "LATTICE_MOLAR_VOLUME": "VOLUME",
    "LATTICE_A": "LATTICE_LENGTH", "LATTICE_B": "LATTICE_LENGTH", "LATTICE_C": "LATTICE_LENGTH",
    "LATTICE_A0": "LATTICE_LENGTH", "LATTICE_C0": "LATTICE_LENGTH", "LATTICE_UNSPECIFIED_AXIS": "LATTICE_LENGTH",
    "LATTICE_A_R": "LATTICE_LENGTH", "LATTICE_A_B": "LATTICE_LENGTH", "LATTICE_C_A": "LATTICE_LENGTH",
    "LATTICE_ALPHA": "LATTICE_ANGLE", "LATTICE_BETA": "LATTICE_ANGLE", "LATTICE_GAMMA": "LATTICE_ANGLE",
    "LATTICE_UNSPECIFIED_ANGLE": "LATTICE_ANGLE",
}


def value_type(rec):
    p = rec.get("property") or ""
    if p.endswith("_EXPRESSION_TERM") or rec.get("property_expression"):
        return None                       # a coefficient of a printed expression in T is not a value of a state
    if p.startswith("ACR_"):
        return "ACTIVITY"
    if p == "GIBBS_ENERGY" and "/K" in (rec.get("unit") or "").replace(" ", ""):
        return "S_MIX_PARTIAL"            # the B (entropy-like) coefficient of a printed A + B T expression
    return PROP_VT.get(p)


def quantities(rec):
    """[(quantity name, value type or None)] measured by the record."""
    q = []
    if isinstance(rec.get("value"), (int, float)):
        q.append(("value", value_type(rec)))
    if any(isinstance(x, (int, float)) for pc in rec.get("phase_compositions") or []
           for x in (pc.get("composition") or {}).values()):
        q.append(("composition", "X_PHASE"))
    if rec.get("property") in ("PHASE_COMPOSITION", "PHASE_ASSEMBLAGE") and isinstance(rec.get("temperature_K"),
                                                                                       (int, float)):
        q.append(("temperature_K", "T_CONDITION"))
    return q


# ------------------------------------------------------------------------------------------- units --
MASS_BASIS = re.compile(r"/\s*\(?\s*g\b(?!\s*[-.·]?\s*at)|\bg\s*[-−⁻]\s*1|per\s+g(?:ram)?\b(?!\s*[-.]?\s*at)", re.I)


def energy_factor(unit):
    """Factor to kJ/mol, or None (mass basis, no energy unit, not printed)."""
    u = unit or ""
    if MASS_BASIS.search(u):
        return None
    if re.search(r"\beV\b", u):
        return 96.48533212
    if re.search(r"kcal", u, re.I):
        return 4.184
    if re.search(r"\bcal\b", u, re.I):
        return 0.004184
    if re.search(r"\bkJ\b|kJ", u):
        return 1.0
    if re.search(r"\bJ\b|J\s*/|J\s*mol|^J", u):
        return 0.001
    return None


def entropy_factor(unit):
    """Factor to J/K/mol, or None."""
    u = unit or ""
    if MASS_BASIS.search(u):
        return None
    if re.search(r"kcal", u, re.I):
        return 4184.0
    if re.search(r"\bcal\b", u, re.I):
        return 4.184
    if "kJ" in u:
        return 1000.0
    if re.search(r"J", u):
        return 1.0
    return None


def factor(rec, vt, quantity="value"):
    """Canonical-unit factor of the quantity (value x factor = canonical), 1 for relative / unitless kinds, or None."""
    kind = VT[vt][1]
    if quantity != "value" or kind in ("temperature", "composition", "relative"):
        return 1.0
    unit = rec.get("unit")
    if kind == "energy":
        return energy_factor(unit)
    if kind == "entropy":
        return entropy_factor(unit)
    if kind == "emf":
        u = (unit or "").strip()
        return 1.0 if u == "mV" else 1000.0 if u == "V" else None
    if kind == "angle":
        return 1.0 if (unit or "").strip().lower() in ("°", "deg", "degree", "degrees", "not_printed", "") else None
    if kind == "activity":
        return 1.0 if (unit or "1").strip() in ("1", "") else None
    return None


def per_atom(unit):
    return bool(re.search(r"atom|at\.|g\.?\s*at|mol-at|mol\.at", unit or "", re.I))


# ------------------------------------------------------------------------------------------- method classes --
THERMO_RULES = [   # (class, regex) on the lower-cased printed method; first match wins
    # audit round 2 (2026-10-24): evaluated / estimated / reviewed values are not independent measurements;
    # Knudsen-cell / effusion methods are tested before EMF (a bare "cell" no longer means EMF)
    ("evaluated_compilation", r"evaluat|estimat|review|interpolat|assess|selected|recommend|compil|critical|"
                              r"\bjanaf\b|tabulat"),
    ("combustion_calorimetry", r"combust|combusion|\bbomb\b|fluorine"),
    ("kems_vapour", r"knudsen|knusden|kundsen|knudson|effusion|\bkem|keml|kewl|torsion|mass[- ]spectr|transpiration|"
                    r"langmuir|isopiestic|\bvptm\b|\bhtms\b|dissociation pressure"),
    ("emf", r"\be\.?m\.?f\b|emf|galvanic|potentiomet|electrochem|coulometr|electromotive|solid[- ]electrolyte|"
            r"electrolyte cell|concentration cell|cathodic"),
    ("kems_vapour", r"vapou?r|\bvap\b|vap\.|boiling point|\bke\b|second[- ]law|third[- ]law|2nd[- ]law|3rd[- ]law|"
                    r"pressure[- ]temperature"),
    ("heat_content_calorimetry", r"adiabatic|levitation|relaxation calor|\bppms\b|low[- ]temperature (?:calor|dsc)|"
                                 r"cp measurement|cp data|from cp|heat[- ]capacity measure"),
    ("dsc_dta", r"\bdsc\b|htdsc|\bdta\b|differential scanning|differential thermal|thermal analysis"),
    ("solution_reaction_calorimetry", r"calori|\bcal\b|\bcal\.|calo\b|calo\.|solution|dissolution|soln|sol\. ?cal|"
                                      r"direct (?:reaction|synthesis)|react\.? ?cal|reaction cal|aneroid|isoperibol|"
                                      r"\bdrop\b|tian|calvet|enthalp"),
    ("equilibration_other", r"equilibri|p[-–]c[-–]t|\bpct\b|sieverts|h2/h2o|distribution|solub|derived|"
                            r"x-ray|xrd|diffraction|chemical analysis|weight loss"),
    ("experimental_unspecified", r"^\s*exp|experiment|measured|^\s*meas"),
]
NOT_INDEPENDENT = {"evaluated_compilation"}   # never paired in the experimental calibration
HEAT_CONTENT_VT = {"CPM", "HEAT_CONTENT", "S_ABS"}
PHASE_RULES = {
    "thermal_analysis": r"\bdsc\b|htdsc|\bdta\b|differential scanning|differential thermal|thermal analys|\bta\b|"
                        r"heating|cooling curve|thermal arrest|calori|thermoanaly",
    "dilatometry_resistivity": r"dilatom|dilatation|resistiv|resistance|magnetic|susceptib|thermal expansion",
    "metallography_microscopy": r"metallogr|microscop|optical|\bom\b|\bsem\b|\btem\b|microstruct|bse|"
                                r"backscatter|incipient melting|quench",
    "xrd": r"xrd|x-ray|x ray|diffraction|rietveld|xrpd|\bxrp\b|neutron|single crystal|\bscm\b|guinier|debye|lattice",
    "microanalysis": r"epma|\beds\b|\bedx\b|\bwds\b|\bwdx\b|microprobe|micro-probe|\bicp\b|chemical analysis|"
                     r"\bxrf\b|\baes\b|atom probe|\bapt\b|spectrom",
    "diffusion_couple": r"diffusion couple|diffusion-couple|diffusion multiple",
    "emf": r"\be\.?m\.?f\b|emf|electromotive|galvanic",
    "experimental_unspecified": r"^\s*exp|experiment|measured|^\s*meas",
}
PHASE_ORDER = {
    "T": ["thermal_analysis", "dilatometry_resistivity", "emf", "metallography_microscopy", "diffusion_couple",
          "xrd", "microanalysis", "experimental_unspecified"],
    "X": ["microanalysis", "diffusion_couple", "xrd", "metallography_microscopy", "thermal_analysis",
          "dilatometry_resistivity", "emf", "experimental_unspecified"],
    "L": ["xrd", "experimental_unspecified"],
}
CLASS_TEXT = {
    "solution_reaction_calorimetry": "solution / reaction calorimetry (direct synthesis, solute-solvent drop, "
                                     "high-temperature solution, isoperibolic, Calvet)",
    "combustion_calorimetry": "combustion (oxygen / fluorine bomb) calorimetry",
    "emf": "EMF (galvanic cells, solid / molten-salt electrolytes, coulometric titration)",
    "kems_vapour": "Knudsen-effusion mass spectrometry / vapour pressure (effusion, torsion, transpiration, "
                   "isopiestic, dissociation pressure)",
    "kems_vapour_second_law": "Knudsen-effusion / vapour pressure, second-law evaluation",
    "kems_vapour_third_law": "Knudsen-effusion / vapour pressure, third-law evaluation",
    "heat_content_calorimetry": "adiabatic / drop / levitation / relaxation calorimetry for Cp and H increments",
    "dsc_dta": "DSC / DTA-derived",
    "equilibration_other": "equilibration and other methods (gas equilibria, p-c-T, solubility, estimates "
                           "derived from other data)",
    "experimental_unspecified": "experimental, method not named ('Exp.', 'measured')",
    "evaluated_compilation": "evaluated / estimated / reviewed / selected values of a compilation or assessment (not an "
                             "independent measurement; never paired in the calibration)",
    "unknown": "method not printed",
    "thermal_analysis": "thermal analysis (DTA / DSC / heating-cooling curves)",
    "dilatometry_resistivity": "dilatometry, electrical resistivity, magnetic measurements",
    "metallography_microscopy": "metallography / optical, SEM, TEM microscopy of quenched samples",
    "xrd": "X-ray / neutron diffraction",
    "microanalysis": "microanalysis (EPMA, EDS / EDX, WDS, ICP, chemical analysis)",
    "diffusion_couple": "diffusion couples / multiples",
    "model_calculated": "model-calculated value (CALPHAD assessment or other model)",
    "dft": "DFT calculation (literature)",
}


def law_of(s):
    if re.search(r"second[- ]law|2nd[- ]law|2nd law|\bii[- ]law", s):
        return "kems_vapour_second_law"
    if re.search(r"third[- ]law|3rd[- ]law|\biii[- ]law", s):
        return "kems_vapour_third_law"
    return "kems_vapour"


def method_class(rec, vt):
    """Method class of one quantity of a record (value type vt)."""
    ev = rec.get("evidence_class")
    if ev == "model":
        return "model_calculated"
    if ev == "dft":
        return "dft"
    m = rec.get("method")
    s = str(m or "").strip().lower()
    if not s:
        return "unknown"
    kind = VT[vt][1] if vt in VT else None
    if vt in ("T_INVARIANT", "T_TRANSITION", "T_CONDITION", "X_PHASE", "LATTICE_LENGTH", "LATTICE_ANGLE", "VOLUME"):
        order = PHASE_ORDER["X" if vt == "X_PHASE" else "L" if vt in ("LATTICE_LENGTH", "LATTICE_ANGLE", "VOLUME")
                            else "T"]
        for c in order:
            if re.search(PHASE_RULES[c], s):
                return c
        return "unknown" if s in ("n/a", "-", "?") else "experimental_unspecified"
    for c, rx in THERMO_RULES:
        if re.search(rx, s):
            if c == "kems_vapour":
                return law_of(s)
            if c == "solution_reaction_calorimetry" and vt in HEAT_CONTENT_VT and kind in ("entropy", "energy"):
                return "heat_content_calorimetry"
            return c
    return "experimental_unspecified"


# ------------------------------------------------------------------------------------------- priors --
# Used only where ORCHID has fewer than MIN_PAIRS independent pairs of a (value type, method class). Values in the
# value type's canonical unit. Each prior names its basis; none is a fitted number.
PRIOR_SOURCES = {
    "espei": "ESPEI 0.9.1 experimental reference sigmas (500 J/mol for enthalpies and the chemical potential of "
             "activities, 0.2 J/mol/K for entropies and heat capacities; espei/error_functions, 'estimated from NIST "
             "TRC uncertainties')",
    "kubaschewski": "O. Kubaschewski, C.B. Alcock, P.J. Spencer, Materials Thermochemistry, 6th ed., Pergamon, 1993 "
                    "(chapter 2, accuracy of experimental methods: calorimetric formation enthalpies typically "
                    "+-1-2 kJ/mol-atom, EMF and vapour-pressure-derived values less accurate, entropies from the "
                    "temperature coefficient of EMF or second-law vapour-pressure evaluation the least)",
    "drowart2005": "J. Drowart, C. Chatillon, J. Hastie, D. Bonnell, High-temperature mass spectrometry: instrumental "
                   "techniques, ionization cross-sections, pressure measurements, and thermodynamic data (IUPAC "
                   "Technical Report), Pure Appl. Chem. 77 (2005) 683-737, doi:10.1351/pac200577040683 (second-law "
                   "enthalpies less accurate than third-law; systematic errors dominate the stated precision)",
    "kirklin2015": "S. Kirklin et al., The Open Quantum Materials Database (OQMD): assessing the accuracy of DFT "
                   "formation energies, npj Comput. Mater. 1 (2015) 15010, doi:10.1038/npjcompumats.2015.10 (mean "
                   "absolute difference between independent experimental formation enthalpies 0.082 eV/atom, "
                   "about 8 kJ/mol-atom, over compounds with more than one experimental value)",
    "orchid_assumed": "ORCHID assumption (PI decision 2026-10-24), stated so that it can be replaced; not fitted",
}
PRIORS = {
    # value type: {method class: (value, source key)}; "*" = any class not listed
    "HM_FORM": {"solution_reaction_calorimetry": (1.0, "kubaschewski"), "combustion_calorimetry": (1.0, "kubaschewski"),
                "emf": (1.0, "kubaschewski"), "kems_vapour_third_law": (2.0, "drowart2005"),
                "kems_vapour": (3.0, "drowart2005"), "kems_vapour_second_law": (3.0, "drowart2005"),
                "dsc_dta": (2.0, "kubaschewski"), "*": (2.0, "kubaschewski")},
    "GM_FORM": {"emf": (1.0, "kubaschewski"), "kems_vapour_third_law": (2.0, "drowart2005"),
                "kems_vapour": (2.0, "drowart2005"), "*": (2.0, "kubaschewski")},
    "SM_FORM": {"emf": (1.0, "kubaschewski"), "heat_content_calorimetry": (0.5, "kubaschewski"),
                "kems_vapour_second_law": (3.0, "drowart2005"), "kems_vapour": (2.0, "drowart2005"),
                "*": (1.0, "kubaschewski")},
    "H_MIX_PARTIAL": {"solution_reaction_calorimetry": (0.5, "espei"), "*": (1.0, "kubaschewski")},
    "G_MIX_PARTIAL": {"emf": (0.5, "espei"), "*": (1.0, "kubaschewski")},
    "S_MIX_PARTIAL": {"emf": (0.5, "kubaschewski"), "*": (1.0, "kubaschewski")},
    "S_ABS": {"heat_content_calorimetry": (0.2, "espei"), "*": (0.5, "kubaschewski")},
    "CPM": {"heat_content_calorimetry": (0.2, "espei"), "dsc_dta": (0.5, "kubaschewski"), "*": (0.5, "kubaschewski")},
    "HEAT_CONTENT": {"heat_content_calorimetry": (0.5, "espei"), "*": (1.0, "kubaschewski")},
    "ACTIVITY": {"emf": (0.5, "espei"), "*": (1.0, "kubaschewski")},
    "EMF": {"*": (1.0, "orchid_assumed")},
    "T_INVARIANT": {"thermal_analysis": (5.0, "orchid_assumed"), "*": (10.0, "orchid_assumed")},
    "T_TRANSITION": {"thermal_analysis": (5.0, "orchid_assumed"), "*": (10.0, "orchid_assumed")},
    "T_CONDITION": {"*": (5.0, "orchid_assumed")},
    "X_PHASE": {"microanalysis": (0.005, "orchid_assumed"), "*": (0.01, "orchid_assumed")},
    "LATTICE_LENGTH": {"*": (0.001, "orchid_assumed")},
    "LATTICE_ANGLE": {"*": (0.05, "orchid_assumed")},
    "VOLUME": {"*": (0.003, "orchid_assumed")},
}
MODEL_PRIORS = {   # model-calculated values, when fewer than MIN_PAIRS model-vs-experiment pairs
    "HM_FORM": 2.0, "GM_FORM": 2.0, "SM_FORM": 1.0, "H_MIX_PARTIAL": 1.0, "G_MIX_PARTIAL": 1.0, "S_MIX_PARTIAL": 1.0,
    "S_ABS": 1.0, "CPM": 1.0, "HEAT_CONTENT": 1.0, "ACTIVITY": 1.0, "EMF": 2.0, "T_INVARIANT": 10.0,
    "T_TRANSITION": 10.0, "T_CONDITION": 5.0, "X_PHASE": 0.02, "LATTICE_LENGTH": 0.002, "LATTICE_ANGLE": 0.1,
    "VOLUME": 0.005}
# model-vs-experiment calibration only where a state key identifies one state reliably (a compound, an invariant
# reaction, an invariant phase composition, a lattice parameter); mixing / partial / T-dependent quantities: prior
MODEL_CALIBRATED = {"HM_FORM", "GM_FORM", "SM_FORM", "T_INVARIANT", "X_PHASE", "LATTICE_LENGTH"}
DFT_ANALOGY = ("prior by analogy: the calibrated DFT formation-enthalpy sigma of literature DFT values (all strata), "
               "used for other DFT energies because no DFT-vs-experiment calibration exists for them")
DFT_PRIORS = {   # value type -> (floor in canonical units or "dft_all", basis text)
    "GM_FORM": ("dft_all", DFT_ANALOGY), "H_MIX_PARTIAL": ("dft_all", DFT_ANALOGY),
    "G_MIX_PARTIAL": ("dft_all", DFT_ANALOGY), "HEAT_CONTENT": ("dft_all", DFT_ANALOGY),
    "SM_FORM": (1.0, "DFT (phonon / quasi-harmonic) entropies: ORCHID assumption, about 1 J/K/mol-atom"),
    "S_ABS": (1.0, "DFT (phonon / quasi-harmonic) entropies: ORCHID assumption, about 1 J/K/mol-atom"),
    "S_MIX_PARTIAL": (1.0, "DFT entropies: ORCHID assumption, about 1 J/K/mol"),
    "CPM": (1.0, "DFT (quasi-harmonic) heat capacities: ORCHID assumption, about 1 J/K/mol-atom"),
    "T_INVARIANT": (50.0, "temperatures calculated from DFT energies: ORCHID assumption, about 50 K"),
    "T_TRANSITION": (50.0, "temperatures calculated from DFT energies: ORCHID assumption, about 50 K"),
    "X_PHASE": (0.02, "compositions calculated from DFT energies: ORCHID assumption, about 0.02"),
    "LATTICE_LENGTH": (0.01, "DFT lattice parameters: typically within about 1 % of experiment (semi-local "
                             "functionals; ORCHID assumption)"),
    "VOLUME": (0.03, "DFT volumes: typically within about 3 % of experiment (ORCHID assumption)"),
    "LATTICE_ANGLE": (0.1, "ORCHID assumption, 0.1 degree"),
    "ACTIVITY": ("dft_all", DFT_ANALOGY), "EMF": (20.0, "ORCHID assumption, 20 mV")}
CONVERSION_PRIOR_FRACTION = 0.10   # lattice-stability conversion: 10 % of the lattice stability (assumed)


def ceil_step(x, step):
    return round(math.ceil(x / step - 1e-9) * step, 6)


# ------------------------------------------------------------------------------------------- state keys --
def comp_of(rec, nd=3):
    oc = {k.upper(): v for k, v in (rec.get("overall_composition") or {}).items() if isinstance(v, (int, float))}
    els = [e.upper() for e in rec.get("elements") or (rec.get("system") or "").split("-") if e]
    if oc:
        miss = [e for e in els if e not in oc]
        if len(miss) == 1:
            oc[miss[0]] = 1.0 - sum(oc.values())
        elif miss:
            return None
    else:
        lab = next((p.get("label") for p in rec.get("phases") or [] if p.get("label")), None)
        if not lab:
            return None
        c, _ = parity.parse_phase(lab, [e.capitalize() for e in els])
        if not c:
            return None
        tot = float(sum(c.values()))
        oc = {k.upper(): float(v) / tot for k, v in c.items()}
    out = tuple(sorted((e, round(x, nd)) for e, x in oc.items() if x > 1e-6))
    return out or None


def label_comp(rec, nd=3):
    """Composition of the record's single compound label (all labels parse to the same formula), or None."""
    els = [e.capitalize() for e in rec.get("elements") or (rec.get("system") or "").split("-") if e]
    found = set()
    for p in rec.get("phases") or []:
        lab = p.get("label")
        if not lab:
            continue
        c, _ = parity.parse_phase(re.sub(r"\s*\((?:cr|s|c|l|liq|liq-am|am|g)\)\s*$", "", str(lab), flags=re.I), els)
        if c:
            tot = float(sum(c.values()))
            found.add(tuple(sorted((k.upper(), round(float(v) / tot, nd)) for k, v in c.items())))
    return found.pop() if len(found) == 1 else None


def state_comp(rec, nd=3):
    """Composition of a compound state: the printed label's formula when it disagrees with overall_composition by more
    than 0.02 (split-row records carry the overall composition of another phase), else overall_composition."""
    c = comp_of(rec, nd)
    lc = label_comp(rec, nd)
    if lc and (not c or {e for e, _ in lc} != {e for e, _ in c} or
               max(abs(dict(lc).get(e, 0) - dict(c).get(e, 0)) for e in {e for e, _ in lc} | {e for e, _ in c}) > 0.02):
        return lc
    return c


LIQ_REF = re.compile(r"liq|melt|\bl\b|liquid", re.I)
NONCRYST = re.compile(r"liq|amorph|glass|\(am\)|\(l\)|liq-am", re.I)


def ref_signature(rec):
    """Reference-state part of a state key: elements referred to their liquid (a liquid reference shifts a formation
    enthalpy by x * the enthalpy of fusion), and a non-crystalline phase state."""
    rs = rec.get("reference_state")
    liq = []
    if isinstance(rs, dict):
        for el, v in sorted(rs.items()):
            ph = v.get("phase") if isinstance(v, dict) else v
            if isinstance(ph, str) and LIQ_REF.search(ph):
                liq.append(el.upper())
    state = "noncrystalline" if any(NONCRYST.search(str(p.get("label") or "")) for p in rec.get("phases") or []) \
        else "crystalline"
    return (tuple(liq), state)


def tbin(T, step):
    if not isinstance(T, (int, float)):
        return "na"
    if abs(T - 298.15) <= 2.0:
        return "std"
    return int(round(T / step) * step)


def pkeys(rec):
    ks = {gates.phase_key(p.get("label")) for p in rec.get("phases") or [] if p.get("label")}
    return frozenset(k for k in ks if k and k.upper() not in ("UNKNOWN", "NONE", "?"))


def first_pkey(rec):
    for p in rec.get("phases") or []:
        k = gates.phase_key(p.get("label"))
        if k and k.upper() not in ("UNKNOWN", "NONE"):
            return k
    return ""


def state_keys(rec, vt):
    """[(state key, value)] of a quantity for pairing (value in canonical units; relative types: raw value)."""
    T = rec.get("temperature_K")
    if vt == "X_PHASE":
        out = []
        if rec.get("property") in ("INVARIANT_TEMPERATURE",):
            base = ("inv", rec.get("system"), tuple(sorted(pkeys(rec))))
            if len(base[2]) < 2:
                return []
        else:
            base = ("pc", rec.get("system"), tbin(T, 5), tuple(sorted(pkeys(rec))))
        for pc in rec.get("phase_compositions") or []:
            comp = {k.upper(): v for k, v in (pc.get("composition") or {}).items() if isinstance(v, (int, float))}
            if not comp or pc.get("composition_source") not in (None, "printed"):
                continue
            ph = gates.phase_key(pc.get("phase"))
            for el in sorted(comp)[:max(1, len(comp) - 1)] if len(comp) > 1 else sorted(comp):
                out.append((base + (ph, el), comp[el]))
        return out
    if vt == "T_CONDITION":
        return []
    v = rec.get("value")
    if not isinstance(v, (int, float)):
        return []
    f = factor(rec, vt)
    if f is None:
        return []
    v = v * f
    p = rec.get("property")
    if vt == "T_INVARIANT":
        ks = tuple(sorted(pkeys(rec)))
        return [(("inv", rec.get("system"), ks), v)] if len(ks) >= 2 else []
    if vt == "T_TRANSITION":
        c, tk, ks = comp_of(rec, 2), rec.get("transition_kind"), tuple(sorted(pkeys(rec)))
        if not c or not ks or tk in (None, "", "other"):
            return []           # a transition is one state only with its kind and phases printed
        return [(("tr", rec.get("system"), c, tk, ks), v)]
    if vt in ("HM_FORM", "GM_FORM", "SM_FORM"):
        c = state_comp(rec, 3)
    else:
        c = comp_of(rec, 3 if vt in ("LATTICE_LENGTH", "LATTICE_ANGLE", "VOLUME") else 2)
    if not c:
        return []
    if vt in ("HM_FORM", "GM_FORM", "SM_FORM") and len(c) < 2:
        return []                       # a formation quantity of a compound needs two or more elements
    if vt in ("HM_FORM", "SM_FORM"):
        if not per_atom(rec.get("unit")):
            return []
        return [((vt, c, tbin(T, 25), ref_signature(rec)), v)]
    if vt == "GM_FORM":
        if not per_atom(rec.get("unit")):
            return []
        return [((vt, c, tbin(T, 10), ref_signature(rec)), v)]
    if vt in ("LATTICE_LENGTH", "LATTICE_ANGLE", "VOLUME"):
        return [((p, c, first_pkey(rec), tbin(T, 25)), v)]
    if vt in ("S_ABS", "CPM", "HEAT_CONTENT"):
        return [((p, c, first_pkey(rec), tbin(T, 5)), v)]
    if vt == "ACTIVITY":
        if not isinstance(T, (int, float)) or not (0 < v <= 1):
            return []
        return [((p, c, tbin(T, 10)), R_GAS * T * math.log(v) / 1000.0)]   # mu - mu_ref, kJ/mol
    if vt == "EMF":
        return []
    return [((p, c, tbin(T, 10), str(rec.get("component") or ""), first_pkey(rec)), v)]


def _ascii(s):
    import unicodedata
    return unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode().lower().strip()


def source_identity(rec):
    """One independent source = (first-author surname, year) when the record names exactly one cited work with both
    (free-text citations of one work printed differently in two papers are one source; audit round 2), else the
    primary source."""
    cs = [c for c in rec.get("cited_sources") or [] if c.get("first_author") and c.get("year")]
    if len(cs) == 1 and len(rec.get("cited_sources") or []) == 1:
        sur = re.sub(r"[^a-z]", "", _ascii(cs[0]["first_author"]).split(",")[0].split()[-1] if
                     _ascii(cs[0]["first_author"]).split() else _ascii(cs[0]["first_author"]))
        if sur:
            return "author:%s:%s" % (sur, str(cs[0]["year"])[:4])
    return primary_source(rec)


def primary_source(rec):
    o = rec.get("origin") or {}
    if o.get("measured_in"):
        return o["measured_in"]
    cs = sorted({c.get("paper") or c.get("doi") or c.get("reference") or "" for c in rec.get("cited_sources") or []}
                - {""})
    if len(cs) == 1:
        return cs[0]
    if cs:
        return "several:" + "|".join(cs)
    if o.get("own_work_of_paper"):
        return rec.get("paper")
    return "unattributed:" + str(rec.get("paper"))


# values whose basis, unit or reliability the page leaves open (audit round 2): kept as printed, never calibration pairs
CALIBRATION_EXCLUDED_CODES = {"row_basis_contradicts_column_header", "unit_inconsistent_with_figure_on_page",
                              "energy_basis_not_printed_rows_inconsistent", "source_marks_value_unreliable",
                              "table_basis_value_not_convertible", "table_basis_formula_unparsed",
                              "formation_from_compounds_not_from_elements", "expression_coefficient_not_a_value"}


def usable(rec):
    fl = rec.get("flags") or {}
    if any(g.get("code") in CALIBRATION_EXCLUDED_CODES for g in rec.get("gate_failures") or []):
        return False
    return fl.get("extraction_accepted") and not any(fl.get(b) for b in BLOCKING)


# ------------------------------------------------------------------------------------------- statistics --
def median(xs):
    s = sorted(xs)
    n = len(s)
    if not n:
        return None
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def wmedian(pairs):
    """Weighted median of [(x, w)]."""
    s = sorted(pairs)
    tot = sum(w for _, w in s)
    if not s or tot <= 0:
        return None
    acc = 0.0
    for i, (x, w) in enumerate(s):
        acc += w
        if acc >= 0.5 * tot - 1e-12:
            if abs(acc - 0.5 * tot) <= 1e-12 and i + 1 < len(s):
                return 0.5 * (x + s[i + 1][0])
            return x
    return s[-1][0]


def state_weighted(clusters):
    """[(|d|, w)]: every state (cluster) carries the same total weight, whatever its number of pairs."""
    return [(abs(x), 1.0 / len(d)) for d in clusters.values() for x in d]


def robust_pair_sigma(d):
    """1.4826 x median|d| / sqrt(2): sigma of one measurement from zero-mean pair differences."""
    return 1.4826 * median([abs(x) for x in d]) / math.sqrt(2.0) if d else None


def robust_state_sigma(clusters, one_sided=False):
    """1.4826 x state-weighted median|d| (/ sqrt(2) for pairs of two measurements)."""
    m = wmedian(state_weighted(clusters))
    if m is None:
        return None
    return 1.4826 * m / (1.0 if one_sided else math.sqrt(2.0))


def boot_ci(clusters, est, seed):
    """95 % percentile bootstrap over clusters (states); est takes a {state: [d]} dict."""
    rng = random.Random(seed)
    keys = sorted(clusters, key=repr)
    vals = []
    for _ in range(N_BOOT):
        cl = {}
        for i in range(len(keys)):
            k = keys[rng.randrange(len(keys))]
            cl[(i, k)] = clusters[k]
        v = est(cl)
        if v is not None:
            vals.append(v)
    vals.sort()
    if not vals:
        return None
    lo = vals[int(0.025 * (len(vals) - 1))]
    hi = vals[int(math.ceil(0.975 * (len(vals) - 1)))]
    return [round(lo, 6), round(hi, 6)]


def relative_diff(a, b):
    m = 0.5 * (abs(a) + abs(b))
    return (a - b) / m if m else None


def same_value(a, b):
    return abs(a - b) <= 1e-9 + 1e-6 * max(abs(a), abs(b))


def atoms_in_formula(comp):
    """Atom count of the smallest integer formula of a composition ((el, x), ...), or None."""
    xs = [x for _, x in comp]
    for n in range(2, 13):
        if all(abs(x * n - round(x * n)) <= 0.02 * n for x in xs) and all(round(x * n) >= 1 for x in xs):
            return n
    return None


def basis_error_pair(key, a, b):
    """True when the two values differ by the formula's atom count (one reading per formula, one per atom)."""
    comp = key[1] if len(key) > 1 and isinstance(key[1], tuple) else None
    if not comp or not a or not b or (a > 0) != (b > 0):
        return False
    n = atoms_in_formula(comp)
    q = max(abs(a), abs(b)) / min(abs(a), abs(b))
    return bool(n) and abs(q / n - 1.0) <= BASIS_RATIO_TOL


BASIS_SCREEN_TOL = 0.15   # a value within 15 % of k x (or 1/k x) the state median, k an atom count >= 2


def basis_factors(comp):
    n = atoms_in_formula(comp)
    if not n:
        return []
    ks = {n} | {round(x * n) for _, x in comp}
    return sorted(k for k in ks if k >= 2)


def basis_screen(key, vals):
    """Values of one state that differ from the state median by an atom-count factor (one reading per formula, the
    other per atom; also a per-formula and a per-atom value of one cited work) are left out of the pairs."""
    comp = key[1] if len(key) > 1 and isinstance(key[1], tuple) else None
    ks = basis_factors(comp) if comp else []
    if not ks or len(vals) < 2:
        return vals
    m = median([v for _, v, _, _ in vals])
    if not m:
        return vals
    out = []
    for x in vals:
        v = x[1]
        q = v / m if m else None
        bad = q is not None and q > 0 and any(abs(q - k) < BASIS_SCREEN_TOL * k or abs(q - 1.0 / k) < BASIS_SCREEN_TOL / k
                                               for k in ks)
        if not bad:
            out.append(x)
    return out


# ------------------------------------------------------------------------------------------- calibration --
def calibrate(records, dft_calibration=None):
    """Calibration of the method floors from the records (experimental between-study scatter; model vs
    experiment). Returns a JSON-able dict (deterministic)."""
    exp = defaultdict(lambda: defaultdict(list))     # vt -> key -> [(source, value, class)]
    mod = defaultdict(lambda: defaultdict(list))
    for r in records:
        if not usable(r):
            continue
        ev = r.get("evidence_class")
        if ev not in ("experimental", "model"):
            continue
        for qn, vt in quantities(r):
            if vt is None:
                continue
            mc = method_class(r, vt)
            if ev == "experimental" and mc in NOT_INDEPENDENT:
                continue
            for key, v in state_keys(r, vt):
                (exp if ev == "experimental" else mod)[vt][key].append((source_identity(r), v, mc, r["id"]))
    out = {"version": VERSION, "min_pairs": MIN_PAIRS, "min_states": MIN_STATES, "min_pairs_cleaned": MIN_PAIRS_CLEANED, "n_boot": N_BOOT, "seed": SEED,
           "estimator": "robust sigma = 1.4826 x median|d| / sqrt(2) over pairs of independent same-class "
                        "measurements of one state (zero-mean differences); 95 % percentile bootstrap over states",
           "independence": "one value per primary source and state (origin.measured_in, else the single cited source, "
                           "else the reporting paper for its own work, else the reporting paper marked unattributed; "
                           "median of the source's values); identical values from different sources are counted once "
                           "(one measurement re-listed)",
           "sample": "experimental, extraction_accepted, no blocking flag (%s); values convertible to the value "
                     "type's canonical unit; formation quantities on a per-atom basis only" % ", ".join(BLOCKING),
           "rule": SIGMA_RULE_TEXT, "value_types": {}, "priors": PRIOR_SOURCES, "class_text": CLASS_TEXT,
           "model": {}, "dft": None}
    for vt in VT:
        unit, kind, mn, step, label = VT[vt]
        groups = exp.get(vt, {})
        within = defaultdict(lambda: defaultdict(list))   # class -> key -> [d]
        pooled = defaultdict(list)
        n_states = Counter()
        relisted = 0
        basis_pairs = 0
        basis_screened = 0
        src_values = {}
        for key in sorted(groups, key=repr):
            bysrc = defaultdict(list)
            vals = groups[key]
            if vt in ("HM_FORM", "GM_FORM", "SM_FORM"):
                kept = basis_screen(key, vals)
                basis_screened += len(vals) - len(kept)
                vals = kept
            for s, v, mc, _ in vals:
                bysrc[s].append((v, mc))
            items = []
            for s in sorted(bysrc):
                v = median([x for x, _ in bysrc[s]])
                mc = Counter(c for _, c in bysrc[s]).most_common(1)[0][0]
                if any(same_value(v, w) for _, w, _ in items):
                    relisted += 1
                    continue
                items.append((s, v, mc))
            src_values[key] = items
            if len(items) < 2:
                continue
            for i in range(len(items)):
                for j in range(i + 1, len(items)):
                    a, b = items[i], items[j]
                    d = relative_diff(a[1], b[1]) if kind == "relative" else a[1] - b[1]
                    if d is None:
                        continue
                    if vt in ("HM_FORM", "GM_FORM", "SM_FORM") and basis_error_pair(key, a[1], b[1]):
                        basis_pairs += 1
                        continue
                    pooled[key].append(d)
                    if a[2] == b[2]:
                        within[a[2]][key].append(d)
        classes = {}
        allc = sorted(set(within) | set(PRIORS.get(vt, {})) - {"*"} | {"unknown", "experimental_unspecified"})
        for c in allc:
            cl = within.get(c, {})
            d = [x for k in cl for x in cl[k]]
            e = {"class": c, "text": CLASS_TEXT.get(c, c), "n_pairs": len(d), "n_states": len(cl)}  # noqa
            if len(d) >= MIN_PAIRS and len(cl) >= MIN_STATES:
                s = robust_state_sigma(cl)
                ci = boot_ci(cl, robust_state_sigma, SEED + len(d))
                e.update({"basis": "calibrated", "sigma_robust": round(s, 6), "ci95": ci,
                          "sigma_robust_unweighted": round(robust_pair_sigma(d), 6),
                          "median_abs_diff": round(median([abs(x) for x in d]), 6),
                          "floor": max(ceil_step(s, step), mn)})
            elif len(d) >= MIN_PAIRS_CLEANED:
                pv, src = PRIORS[vt].get(c) or PRIORS[vt]["*"]
                s = robust_state_sigma(cl)
                ci = boot_ci(cl, robust_state_sigma, SEED + len(d))
                e.update({"basis": "cleaned_below_threshold", "prior_source": src, "prior": pv,
                          "why": "%d pairs from %d states (below the >= %d pairs / >= %d states calibration rule, "
                                 ">= %d pairs): floor = max(documented prior, cleaned robust sigma)"
                                 % (len(d), len(cl), MIN_PAIRS, MIN_STATES, MIN_PAIRS_CLEANED),
                          "sigma_robust": round(s, 6), "ci95": ci,
                          "floor": max(pv, ceil_step(s, step), mn)})
            else:
                pv, src = PRIORS[vt].get(c) or PRIORS[vt]["*"]
                e.update({"basis": "prior", "prior_source": src,
                          "why_prior": "%d pairs from %d states (calibration needs >= %d pairs and >= %d states)"
                                       % (len(d), len(cl), MIN_PAIRS, MIN_STATES),
                          "sigma_robust_insufficient_data": round(robust_state_sigma(cl), 6) if d else None,
                          "floor": max(pv, mn)})
            e["floor_at_minimum"] = e["floor"] == mn
            classes[c] = e
        pd = [x for k in pooled for x in pooled[k]]
        out["value_types"][vt] = {
            "label": label, "canonical_unit": unit, "minimum_floor": mn, "rounding_step": step,
            "minimum_basis": PRIOR_SOURCES["espei"] if kind in ("energy", "entropy", "activity") else
            "smallest meaningful floor for the quantity (ORCHID choice)",
            "states_with_two_or_more_sources": sum(1 for k in src_values if len(src_values[k]) >= 2),
            "relisted_values_counted_once": relisted,
            "basis_error_pairs_excluded": basis_pairs,
            "basis_screened_values_excluded": basis_screened,
            "pooled_all_classes": {"n_pairs": len(pd), "sigma_robust": round(robust_pair_sigma(pd), 6) if pd else None},
            "classes": classes}
        # model-calculated values vs experiment of the same state
        md = defaultdict(list)
        for key, items in mod.get(vt, {}).items():
            ex = src_values.get(key) or []
            if not ex:
                continue
            for s, v, _, _ in items:
                for _, w, _ in ex:
                    d = relative_diff(v, w) if kind == "relative" else v - w
                    if d is not None:
                        md[key].append(d)
        d = [x for k in md for x in md[k]]
        est = lambda cl: robust_state_sigma(cl, one_sided=True)   # noqa: E731  (model vs one measurement)
        if vt in MODEL_CALIBRATED and len(d) >= MIN_PAIRS and len(md) >= MIN_STATES:
            s = est(md)
            out["model"][vt] = {"basis": "calibrated", "n_pairs": len(d), "n_states": len(md),
                                "sigma_robust": round(s, 6), "ci95": boot_ci(md, est, SEED + 7 + len(d)),
                                "floor": max(ceil_step(s, step), mn),
                                "estimator": "1.4826 x median|model - experiment| over states with both (includes the "
                                             "experimental error; conservative); the models were often fitted to "
                                             "these experiments, so it is a lower bound of the model's error"}
        elif vt in MODEL_CALIBRATED and len(d) >= MIN_PAIRS_CLEANED:
            s = est(md)
            out["model"][vt] = {"basis": "cleaned_below_threshold", "n_pairs": len(d), "n_states": len(md),
                                "sigma_robust": round(s, 6), "ci95": boot_ci(md, est, SEED + 7 + len(d)),
                                "prior": MODEL_PRIORS[vt], "prior_source": "orchid_assumed",
                                "floor": max(MODEL_PRIORS[vt], ceil_step(s, step), mn),
                                "why": "%d pairs from %d states (below the calibration rule, >= %d pairs): floor = "
                                       "max(documented prior, cleaned robust sigma)" % (len(d), len(md), MIN_PAIRS_CLEANED)}
        else:
            out["model"][vt] = {"basis": "prior", "n_pairs": len(d), "n_states": len(md),
                                "floor": max(MODEL_PRIORS[vt], mn), "prior_source": "orchid_assumed",
                                "why_prior": "value type not calibrated against experiment (state keys of this type "
                                             "do not identify one state reliably enough)" if vt not in
                                MODEL_CALIBRATED else "%d pairs from %d states" % (len(d), len(md))}
    if dft_calibration:
        out["dft"] = {"file": DFT_CALIBRATION_FILE,
                      "note": "sources.literature of the ESPEI DFT calibration (values only; its input hash is "
                              "left out so that this file depends on values, not on the release bytes)",
                      "literature": (dft_calibration.get("sources") or {}).get("literature"),
                      "priors": {}}
        L = out["dft"]["literature"] or {}
        allsig = (L.get("all") or {}).get("sigma_robust")
        for k, (v, t) in DFT_PRIORS.items():
            if v == "dft_all":
                v = ceil_step(allsig, VT[k][3]) if allsig else None
            out["dft"]["priors"][k] = {"floor": v, "text": t}
    return out


def calibration_bytes(cal):
    return (json.dumps(cal, indent=1, sort_keys=True, ensure_ascii=False) + "\n").encode()


def calibration_id(cal):
    return "%s sha256 %s" % (VERSION, hashlib.sha256(calibration_bytes(cal)).hexdigest()[:12])


# ------------------------------------------------------------------------------------------- budgets --
def _r(x):
    if x is None:
        return None
    if x == 0:
        return 0.0
    return float("%.4g" % x)


TWO_SIGMA = re.compile(r"2\s*σ|2\s*sigma|two standard deviations|twice the standard|95\s*%|95 per ?cent|"
                       r"expanded uncertaint|\bk\s*=\s*2\b|coverage factor (?:of )?2", re.I)
NUMTXT = re.compile(r"[-+−]?\d+(?:[.,](\d+))?")


def decimals(x):
    """Decimal places of a printed number (string or float), or None."""
    if isinstance(x, bool) or x is None:
        return None
    if isinstance(x, (int, float)):
        t = repr(float(x)) if isinstance(x, float) else str(x)
        if "e" in t or "E" in t:
            return None
        if isinstance(x, int):
            return 0
        d = t.split(".")[1] if "." in t else ""
        # a JSON float such as 1.0 was a token with a decimal point (an integer is stored as int): one decimal at
        # least (audit round 2: "1.0" has one decimal, not zero)
        return len(d) if d != "0" else 1
    m = NUMTXT.search(str(x))
    if not m:
        return None
    return len(m.group(1) or "")


def reading_component(rec, quantity):
    """Half of the last printed digit, in the quantity's unit (sigma_reading)."""
    if quantity == "value":
        v, vp = rec.get("value"), rec.get("value_printed")
        src = vp if vp is not None else v
        d = decimals(src)
        if d is None or not isinstance(v, (int, float)):
            return None
        res = 0.5 * 10.0 ** (-d)
        try:
            vpn = float(str(vp).replace("−", "-").replace(",", ".")) if vp is not None else None
        except ValueError:
            vpn = None
        scale = abs(v / vpn) if vpn else 1.0
        return {"value": _r(res * scale), "basis": "last_printed_digit"}
    if quantity == "temperature_K":
        tp = rec.get("temperature_printed")
        d = decimals(tp) if tp is not None else decimals(rec.get("temperature_K"))
        return {"value": _r(0.5 * 10.0 ** (-d)), "basis": "last_printed_digit"} if d is not None else None
    if quantity == "composition":
        ds = [decimals(x) for pc in rec.get("phase_compositions") or [] for x in (pc.get("composition") or {}).values()
              if isinstance(x, (int, float))]
        ds = [d for d in ds if d is not None]
        if not ds:
            return None
        return {"value": _r(0.5 * 10.0 ** (-max(ds))), "basis": "last_stored_digit"}
    return None


def printed_component(rec, quantity, vt):
    """sigma_printed and sigma_reading of a quantity (record units)."""
    if quantity == "composition":
        return None, reading_component(rec, quantity)
    if quantity == "temperature_K" or vt in ("T_INVARIANT", "T_TRANSITION") and not rec.get("uncertainty"):
        u = rec.get("temperature_uncertainty") or {}
    else:
        u = rec.get("uncertainty") or {}
    b = u.get("basis")
    st = rec.get("uncertainty_stored") or {}
    if b in ("printed", "printed_range") and isinstance(u.get("value"), (int, float)):
        note = " ".join(str(x or "") for x in (u.get("note"), (rec.get("locator") or {}).get("caption")))
        kind = ("range_half_width" if b == "printed_range" else
                "two_sigma_or_95_percent" if TWO_SIGMA.search(note) else
                "standard_deviation" if re.search(r"standard deviation|std\.? dev|standard uncertaint", note, re.I) else
                "standard_error" if re.search(r"standard error", note, re.I) else "plus_minus_coverage_not_stated")
        if kind == "two_sigma_or_95_percent":   # audit round 2: the printed coverage is resolved from the caption
            return {"value": _r(abs(u["value"]) / 2.0), "kind": kind, "printed": _r(abs(u["value"])),
                    "taken_as": "1 sigma = half of the printed 2-sigma / 95 % interval"}, reading_component(rec, quantity)
        return {"value": _r(abs(u["value"])), "kind": kind,
                "taken_as": "half-width" if kind == "range_half_width" else "1 sigma"}, reading_component(rec, quantity)
    if b in ("scatter_between_sources", "printed_with_composition_offset"):
        pv = u.get("printed")
        if not isinstance(pv, (int, float)) and st.get("basis") in ("printed", "printed_range"):
            pv = st.get("value")
        if isinstance(pv, (int, float)):
            return {"value": _r(abs(pv)), "kind": "plus_minus_coverage_not_stated", "taken_as": "1 sigma",
                    "where": "uncertainty_stored" if not isinstance(u.get("printed"), (int, float)) else
                    "uncertainty.printed"}, reading_component(rec, quantity)
    return None, reading_component(rec, quantity)


def method_component(rec, quantity, vt, cal):
    if vt is None:
        return {"value": None, "reason": "no method floor defined for property %s" % rec.get("property")}
    unit, kind, mn, step, label = VT[vt]
    ev = rec.get("evidence_class")
    mc = method_class(rec, vt)
    f = factor(rec, vt, quantity)
    ref = None
    if ev == "dft":
        if vt == "HM_FORM" and (cal.get("dft") or {}).get("literature"):
            L = cal["dft"]["literature"]
            els = {e.capitalize() for e in rec.get("elements") or []}
            P_BLOCK = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Se", "Te", "As", "Sb", "Cl", "Br", "I", "At"}
            chem = "pblock_anion" if els & P_BLOCK else "magnetic_3d" if els & {"Cr", "Mn", "Fe", "Co", "Ni"} \
                else "metal"
            st = (L.get("strata") or {}).get(chem) or L.get("all")
            if (L.get("strata") or {}).get(chem) is None:
                chem = "all"
            e = {"class": "dft", "stratum": "literature/" + chem, "floor": _r(st["sigma_robust"]),
                 "basis": "calibrated", "n_pairs": st["n"], "ci95": (st.get("ci95") or {}).get("sigma_robust"),
                 "entry": "dft/literature/" + ("all" if chem == "all" else "strata/" + chem)}
        elif vt in DFT_PRIORS and (cal.get("dft") or {}).get("priors", {}).get(vt, {}).get("floor") is not None:
            pr = cal["dft"]["priors"][vt]
            e = {"class": "dft", "floor": pr["floor"], "basis": "prior", "entry": "dft/priors/" + vt}
        else:
            return {"value": None, "class": "dft",
                    "reason": "no DFT calibration or prior for %s (the DFT calibration covers formation "
                              "enthalpies)" % label}
    elif ev == "model":
        m = cal["model"][vt]
        e = {"class": "model_calculated", "floor": m["floor"], "basis": m["basis"], "n_pairs": m["n_pairs"],
             "ci95": m.get("ci95"), "entry": "model/" + vt}
    else:
        c = cal["value_types"][vt]["classes"].get(mc) or cal["value_types"][vt]["classes"]["unknown"]
        e = {"class": mc, "floor": c["floor"], "basis": c["basis"], "n_pairs": c["n_pairs"],
             "ci95": c.get("ci95"), "entry": "value_types/%s/classes/%s" % (vt, c["class"])}
        if c["basis"] in ("prior", "cleaned_below_threshold"):
            e["prior_source"] = c["prior_source"]
    e["floor_unit"] = unit
    # in the quantity's unit
    v = rec.get("value")
    if kind == "relative":
        if not isinstance(v, (int, float)):
            return dict(e, value=None, reason="relative floor needs a numeric value")
        e["value"] = _r(e["floor"] * abs(v))
    elif kind == "activity":
        T = rec.get("temperature_K")
        if not isinstance(v, (int, float)) or not isinstance(T, (int, float)) or T <= 0 or not 0 < v <= 1:
            return dict(e, value=None, reason="activity floor (on mu) needs T > 0 and 0 < a <= 1")
        e["value"] = _r(abs(v) * e["floor"] * 1000.0 / (R_GAS * T))   # sigma_a = a sigma_mu / (R T)
    elif f is None:
        return dict(e, value=None, reason="unit %r not convertible to %s" % (rec.get("unit"), unit))
    else:
        e["value"] = _r(e["floor"] / f)
        if vt in ("HM_FORM", "GM_FORM", "SM_FORM") and not per_atom(rec.get("unit")):
            # basis unstated: the value may be per formula unit, so the per-atom floor is scaled by the atoms per
            # formula unit (audit round 2; conservative when the value is in fact per atom)
            lc = label_comp(rec)
            n = atoms_in_formula(lc) if lc else None
            if n and n > 1:
                e["value"] = _r(e["value"] * n)
                e["scaled_for_unstated_basis"] = n
    return e


def composition_component(rec, quantity):
    if quantity != "value":
        return None
    u = rec.get("uncertainty") or {}
    note = str(u.get("note") or "")
    if u.get("basis") == "printed_with_composition_offset" and isinstance(u.get("sigma_composition_offset"),
                                                                          (int, float)):
        m = re.search(r"\|dx\|\s*([\d.]+)\s*x\s*slope\s*([\d.]+)", note)
        return {"value": _r(u["sigma_composition_offset"]), "dx": float(m.group(1)) if m else None,
                "slope_kJ_per_mol_per_x": float(m.group(2)) if m else None, "basis": "page_checked",
                "note": "uncertainty.note"}
    m = re.search(r"own sigma ([\d.]+) \(printed_with_composition_offset\)", note)
    st = rec.get("uncertainty_stored") or {}
    if u.get("basis") == "scatter_between_sources" and m and isinstance(st.get("value"), (int, float)):
        own = float(m.group(1))
        return {"value": _r(math.sqrt(max(own ** 2 - st["value"] ** 2, 0.0))), "dx": None,
                "slope_kJ_per_mol_per_x": None, "basis": "page_checked",
                "note": "sqrt(own %g^2 - printed %g^2) of uncertainty.note" % (own, st["value"])}
    return None


def scatter_component(rec, quantity):
    if quantity != "value":
        return None
    u = rec.get("uncertainty") or {}
    if u.get("basis") != "scatter_between_sources" or not isinstance(u.get("value"), (int, float)):
        return None
    ids = sorted(set(re.findall(r"ORC-[A-Z0-9]+(?:\.[a-f0-9]{6})?", str(u.get("note") or ""))) | {rec["id"]})
    return {"value": _r(u["value"]), "n_sources": u.get("n_sources"), "group": ids,
            "basis": "page_checked mutually_inconsistent group", "note": "uncertainty.note"}


def conversion_component(rec, quantity):
    if quantity != "value":
        return None
    for c in rec.get("corrections") or []:
        conv = (c.get("evidence") or {}).get("conversion")
        if c.get("field") == "value" and conv and isinstance((conv.get("lattice_stability") or {}).get(
                "value_J_per_mol"), (int, float)):
            ls = conv["lattice_stability"]
            x = conv.get("x_Zr") or 0.0
            sJ = CONVERSION_PRIOR_FRACTION * abs(ls["value_J_per_mol"]) * x
            f = energy_factor(rec.get("unit")) or 1.0
            return {"value": _r(sJ / 1000.0 / f), "basis": "prior",
                    "note": "%g x |lattice stability %g J/mol| x x %g (corrections[].evidence.conversion)"
                            % (CONVERSION_PRIOR_FRACTION, ls["value_J_per_mol"], x)}
    return None


def combine(q):
    base = [c["value"] for c in (q.get("sigma_printed"), q.get("sigma_reading"), q.get("sigma_method"))
            if c and isinstance(c.get("value"), (int, float))]
    quad = [c["value"] for c in (q.get("sigma_composition"), q.get("sigma_conversion"))
            if c and isinstance(c.get("value"), (int, float))]
    sc = q.get("sigma_scatter")
    anchor = [c for c in (q.get("sigma_printed"), q.get("sigma_method"), sc)
              if c and isinstance(c.get("value"), (int, float))]
    if not anchor:
        return None          # the reading resolution alone is not an uncertainty estimate
    s = math.sqrt((max(base) if base else 0.0) ** 2 + sum(x * x for x in quad))
    if sc and isinstance(sc.get("value"), (int, float)):
        s = max(s, sc["value"])
    return _r(s)


def quantity_unit(rec, quantity, vt):
    if quantity == "composition":
        return "mole fraction"
    if quantity == "temperature_K":
        return "K"
    return rec.get("unit")


def budget(rec, cal, cal_id):
    qs = {}
    for qn, vt in quantities(rec):
        sp, sr = printed_component(rec, qn, vt)
        q = {"unit": quantity_unit(rec, qn, vt), "value_type": vt, "sigma_printed": sp, "sigma_reading": sr,
             "sigma_method": method_component(rec, qn, vt, cal), "sigma_composition": composition_component(rec, qn),
             "sigma_scatter": scatter_component(rec, qn), "sigma_conversion": conversion_component(rec, qn)}
        q["sigma_effective"] = combine(q)
        if q["sigma_effective"] is None:
            q["sigma_effective_reason"] = (q["sigma_method"] or {}).get("reason") or "no component applies"
        qs[qn] = q
    b = {"calibration": cal_id, "sigma_rule": SIGMA_RULE, "quantities": qs}
    if not qs:
        b["reason"] = "no numeric quantity (no value, no printed phase composition)"
    return b


def annotate(records, cal):
    cid = calibration_id(cal)
    for r in records:
        r["uncertainty_budget"] = budget(r, cal, cid)
    return records


def summary(records, cal):
    """Counts per value type x method class: values, with printed sigma, sigma_effective > printed."""
    rows = defaultdict(Counter)
    for r in records:
        for qn, q in (r.get("uncertainty_budget") or {}).get("quantities", {}).items():
            m = q.get("sigma_method") or {}
            k = (r.get("evidence_class"), q.get("value_type") or "none", m.get("class") or "none")
            rows[k]["n"] += 1
            sp = (q.get("sigma_printed") or {}).get("value")
            se = q.get("sigma_effective")
            if sp is not None:
                rows[k]["printed"] += 1
                if se is not None and se > sp * (1 + 1e-9):
                    rows[k]["raised"] += 1
            if se is None:
                rows[k]["no_effective"] += 1
    return rows


def totals(records, cal):
    """Release-level counts of the budgets (index.json totals)."""
    c = Counter()
    for r in records:
        b = r.get("uncertainty_budget") or {}
        c["records"] += 1
        c["records_with_budget"] += 1 if b else 0
        for qn, q in (b.get("quantities") or {}).items():
            c["quantities"] += 1
            c["quantity:" + qn] += 1
            c["method_basis:" + str((q.get("sigma_method") or {}).get("basis"))] += 1
            if q.get("sigma_effective") is not None:
                c["with_sigma_effective"] += 1
            sp = (q.get("sigma_printed") or {}).get("value")
            if sp is not None:
                c["with_sigma_printed"] += 1
                if q.get("sigma_effective") is not None and q["sigma_effective"] > sp * (1 + 1e-9):
                    c["sigma_effective_above_printed"] += 1
            for k in ("sigma_composition", "sigma_scatter", "sigma_conversion"):
                if q.get(k):
                    c["with_" + k] += 1
    return dict(sorted(c.items()), calibration=calibration_id(cal), rule=SIGMA_RULE)


def fmt(x):
    return "" if x is None else ("%.4g" % x)


def floor_table(cal, only_calibrated=False):
    """Markdown rows: value type x method class (n pairs, n states, basis, robust sigma, 95 % CI, floor, unit)."""
    L = ["| value type | method class | pairs | states | basis | robust sigma | 95 % CI | floor | unit |",
         "|---|---|---:|---:|---|---:|---|---:|---|"]
    for vt, v in cal["value_types"].items():
        for c, e in v["classes"].items():
            if only_calibrated and e["basis"] != "calibrated":
                continue
            ci = e.get("ci95")
            L.append("| %s | %s | %d | %d | %s | %s | %s | %s | %s |" % (
                vt, c, e["n_pairs"], e["n_states"], e["basis"] + ("" if e["basis"] == "calibrated" else
                                                                 " (%s)" % e.get("prior_source")),
                fmt(e.get("sigma_robust")), "%s-%s" % (fmt(ci[0]), fmt(ci[1])) if ci else "", fmt(e["floor"]),
                v["canonical_unit"]))
        m = cal["model"][vt]
        if not only_calibrated or m["basis"] == "calibrated":
            ci = m.get("ci95")
            L.append("| %s | model_calculated (vs experiment) | %d | %d | %s | %s | %s | %s | %s |" % (
                vt, m["n_pairs"], m.get("n_states", 0), m["basis"], fmt(m.get("sigma_robust")),
                "%s-%s" % (fmt(ci[0]), fmt(ci[1])) if ci else "", fmt(m["floor"]), v["canonical_unit"]))
    return L


def count_table(records):
    rows = summary(records, None)
    L = ["| evidence class | value type | method class | quantities | with printed sigma | sigma_effective > printed "
         "| no sigma_effective |", "|---|---|---|---:|---:|---:|---:|"]
    for k in sorted(rows, key=lambda x: tuple(str(y) for y in x)):
        c = rows[k]
        L.append("| %s | %s | %s | %d | %d | %d | %d |" % (k + (c["n"], c["printed"], c["raised"], c["no_effective"])))
    return L


def calibration_md(cal, records):
    L = ["# ORCHID uncertainty budget: method-floor calibration", "",
         "Version: %s. Machine-readable: `UNCERTAINTY_CALIBRATION.json` (the same file is copied to "
         "`espei/experimental/CALIBRATION_EXPERIMENT.json`). Every record carries `uncertainty_budget` (SCHEMA.md)."
         % VERSION, "",
         "## Rule", "", SIGMA_RULE_TEXT, "",
         "## Calibration", "",
         "- estimator: " + cal["estimator"] + "; every state carries the same weight (state-weighted median), so a "
         "compound measured by many studies does not dominate;",
         "- independence: " + cal["independence"] + ";",
         "- sample: " + cal["sample"] + "; formation-quantity pairs whose ratio is within %d %% of the formula's atom "
         "count are a per-formula / per-atom basis error of one reading and are left out;" % int(BASIS_RATIO_TOL * 100),
         "- a class is calibrated with at least %d pairs from at least %d states; otherwise the documented prior "
         "(sources below) is used (basis `prior`);" % (MIN_PAIRS, MIN_STATES),
         "- floor = robust sigma rounded up to the value type's step, never below the minimum (ESPEI's experimental "
         "default for energies, entropies and activities);",
         "- model-calculated values: robust deviation of the model class from experiment for value types whose state "
         "key identifies one state (formation quantities, invariant temperatures and compositions, lattice "
         "parameters), else a prior; DFT values: the calibrated sigma_DFT of `espei/dft/CALIBRATION.json` (literature "
         "DFT, by chemistry stratum) for formation enthalpies, documented priors otherwise.", "",
         "| value type | canonical unit | minimum | label |", "|---|---|---:|---|"]
    for vt, (u, kind, mn, step, label) in VT.items():
        L.append("| %s | %s | %s | %s |" % (vt, u, fmt(mn), label))
    L += ["", "## Floors per value type and method class", ""] + floor_table(cal)
    if cal.get("dft"):
        Ld = cal["dft"]["literature"] or {}
        L += ["", "## DFT", "", "Literature DFT formation enthalpies (from `%s`): all n %s, robust sigma %s kJ/mol-atom; "
              "strata: %s." % (DFT_CALIBRATION_FILE, (Ld.get("all") or {}).get("n"), fmt((Ld.get("all") or {}).get(
                  "sigma_robust")), ", ".join("%s n %s sigma %s" % (k, v.get("n"), fmt(v.get("sigma_robust")))
                                             for k, v in sorted((Ld.get("strata") or {}).items()))), "",
              "| value type | prior floor | basis |", "|---|---:|---|"]
        for k, v in cal["dft"]["priors"].items():
            L.append("| %s | %s | %s |" % (k, fmt(v["floor"]), v["text"]))
    L += ["", "## Prior sources", ""] + ["- `%s`: %s" % (k, v) for k, v in PRIOR_SOURCES.items()]
    L += ["", "## Values per evidence class, value type and method class", ""] + count_table(records) + [""]
    return "\n".join(L)


def load_dft_calibration(data):
    p = os.path.join(data, DFT_CALIBRATION_FILE)
    return json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=os.environ.get("ORCHID_DATA", os.path.expanduser("~/projects/orchid-data")))
    a = ap.parse_args(argv)
    recs = [json.loads(l) for l in gzip.open(os.path.join(a.data, "literature", "orchid-literature.jsonl.gz"), "rt",
                                              encoding="utf-8")]
    cal = calibrate(recs, load_dft_calibration(a.data))
    for vt, v in cal["value_types"].items():
        print("%-15s states>=2src %5d relisted %5d pooled n %6d sigma %s" % (
            vt, v["states_with_two_or_more_sources"], v["relisted_values_counted_once"],
            v["pooled_all_classes"]["n_pairs"], v["pooled_all_classes"]["sigma_robust"]))
        for c, e in v["classes"].items():
            print("    %-32s n %6d %-10s sigma %-10s ci %-24s floor %s" % (
                c, e["n_pairs"], e["basis"], e.get("sigma_robust"), e.get("ci95"), e["floor"]))
        m = cal["model"][vt]
        print("    %-32s n %6d %-10s sigma %-10s ci %-24s floor %s" % ("MODEL", m["n_pairs"], m["basis"],
                                                                      m.get("sigma_robust"), m.get("ci95"), m["floor"]))
    annotate(recs, cal)
    rows = summary(recs, cal)
    for k in sorted(rows, key=lambda x: tuple(str(y) for y in x)):
        print(k, dict(rows[k]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
