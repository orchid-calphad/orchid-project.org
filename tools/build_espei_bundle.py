#!/usr/bin/env python3
"""Build the ORCHID ESPEI bundle: literature ESPEI datasets gated by the published literature values, with
uncertainties, plus DFT endmember formation enthalpies with calibrated uncertainties.

Reads (never writes)

  <corpus>/db-v2.1/general/espei/<TDB>/*.json        experimental datasets written by the corpus record builder
  <corpus>/db-v2.1/general/espei-model/<TDB>/*.json  model-calculated datasets (evidence class model, TDB-linked)
  <corpus>/db-v2.1/invariants/espei/<TDB>/*.json     ZPF datasets from invariant-reaction tables
  <corpus>/db-v2.1/invariants/espei_evaluation.json  ZPF driving forces against the assessment (self-check)
  <data>/literature/orchid-literature.jsonl.gz       the PUBLISHED literature values (inclusion gate, uncertainties)
  <data>/literature/papers.json                      citations / DOIs (nothing beyond what the release publishes)
  <data>/literature/dft_parity.json                  DFT-vs-experiment pairs (calibration of the DFT sigma)
  <data>/hulls/dft-mp, dft-oqmd/<EL>/<SYS>.json      on-hull DFT vertices (formation energies, eV/atom)
  <data>/phases/systems/<SYS>.json                   phase registry: DFT entry -> phase, paper label -> phase, TDB models

and writes <data>/espei/:

  experimental/<SYS>/*.json  model/<SYS>/*.json  invariants/<SYS>/*.json
  dft/literature/<SYS>/*.json  dft/mp/<SYS>/*.json  dft/oqmd/<SYS>/*.json   (HM_FORM endmember datasets, `solver`)
  dft/CALIBRATION.json, dft/CALIBRATION.md       DFT sigma calibration (robust sigma, bootstrap CIs, strata, LSO)
  index.json, index.csv                          one row per dataset
  references.bib                                 BibTeX of every paper a dataset cites
  validation.json, VALIDATION.md                 ESPEI check_dataset / load_datasets result, drops, exclusions
  UNCERTAINTY_SUMMARY.md                         sigma bases and distributions per folder and output
  README.md, SCHEMA.md, CITATION.cff

Inclusion (same rules as the published literature values, tools/export_literature_data.py): every value of a
corpus dataset is traced to its observation and kept only if that observation is published with
`extraction_accepted` true, in the folder's evidence class, and the published value (number, temperature, phase
compositions) agrees with the dataset value. Failing values are cut from the ESPEI arrays; a dataset left without
values is dropped. Datasets ESPEI 0.9.1 rejects for known data reasons (printed mole fractions summing above one;
a null inside a composition that ESPEI cannot convert) are excluded whole and listed; printed values are never
altered. Any other ESPEI failure fails the build.

The ESPEI checks run in the ESPEI interpreter (--espei-python, default ~/.venvs/espei/bin/python) on a temporary
copy; this script is standard library only. No network, no model calls. Deterministic (fixed bootstrap seeds).

  python3 tools/build_espei_bundle.py            # write <data>/espei/
  python3 tools/build_espei_bundle.py --check    # verify the files on disk are current (exit 1 if not)
"""
from __future__ import annotations

import argparse
import csv
import glob
import gzip
import hashlib
import io
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
import unicodedata
from collections import Counter, defaultdict
from fractions import Fraction

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import build_dft_parity as parity  # noqa: E402  (formula parser shared with the DFT benchmark)
import export_literature_data as lit  # noqa: E402
import orchid_cite  # noqa: E402
import espei_constraints as ec  # noqa: E402
import phase_rule as PR  # noqa: E402

SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_ESPEI_PY = os.path.expanduser("~/.venvs/espei/bin/python")
EVAL_WORKERS = int(os.environ.get("ORCHID_EVAL_WORKERS", "6"))
NAME = "espei"
STATUS_LINE = "Version 1, machine-extracted"

SOURCES = [("general/espei", "experimental"), ("general/espei-model", "model"), ("invariants/espei", "invariants")]
FOLDER_CLASS = {"experimental": "experimental", "model": "model", "invariants": "experimental"}
GROUPS = ("experimental", "model", "invariants", "dft/literature", "dft/mp", "dft/oqmd")
PGROUPS = ("assessed", "derived/invariant-brackets", "derived/liquid-stability")   # optional phase constraints

CAVEAT_FLAGS = ("magnitude_suspect", "system_mismatch", "source_disputed", "dft_outlier_unexplained",
                "basis_inferred_from_magnitude", "basis_corrected_from_formula", "sign_from_heading",
                "unsigned_in_print", "judge_adjudicated", "page_unreadable")
COMPOUND_X_TOL = 0.01      # printed composition vs the compound's formula composition (mole fraction)
DROP_REASONS = {
    "not_in_release": "observation not published in the literature release (added or re-read after the release)",
    "not_accepted": "published with extraction_accepted = false",
    "evidence_class_mismatch": "published evidence class differs from the folder's class (classes are never mixed)",
    "value_differs_from_release": "dataset value differs from the published value (e.g. a quality fix applied in "
                                  "the release: per-formula -> per-atom basis); the dataset value is not used",
    "temperature_differs_from_release": "dataset temperature differs from the published temperature",
    "composition_differs_from_release": "a ZPF vertex composition is not among the published phase compositions",
    "magnitude_implausible": "physically implausible magnitude for fitting (|S| > 150 J/K/mol-atom, |Cp| > 200 "
                             "J/K/mol-atom, |H| or |G| > 1000 kJ/mol-atom): a unit error in the source or the "
                             "extraction; the value stays in the literature release with its flags",
}
PHASE_DROP_TEXT = {
    "phase_rule:phase_not_printed": "physical phase rule: no phase printed for the value and the caption does not "
                                    "name the liquid (no default phase is assigned)",
    "phase_rule:phase_label_not_mapped": "physical phase rule: the printed phase label is not mapped to a model "
                                         "phase by the phase registry",
    "phase_rule:phase_label_ambiguous": "physical phase rule: the printed label maps to several registry phases",
    "phase_rule:phase_without_tdb_model": "physical phase rule: the registry phase has no TDB model",
    "phase_rule:compound_property_not_fittable": "physical phase rule: a property of a stoichiometric compound that "
                                                 "ESPEI 0.9.1 cannot fit (only HM_FORM / SM_FORM of endmembers)",
    "phase_rule:liquid_impossible_at_T": "physical phase rule: liquid printed at a temperature below the lowest "
                                         "possible liquid temperature of the system",
    "phase_rule:not_espei_ready": "the published value is not espei_ready",
    "phase_not_in_target_tdb": "the value's phase has no model in the target TDB (the assessment folder)",
    "phase_several_models_in_target_tdb": "the value's registry phase has several models in the target TDB",
    "composition_not_at_compound": "compound formation value whose printed composition differs from the compound's "
                                   "formula composition in the target TDB by more than %g" % COMPOUND_X_TOL,
    "no_unique_endmember_in_target_tdb": "no single endmember of the compound's TDB model at its composition",
    "duplicate_identical_property_dataset": "E3 (property datasets): identical to another dataset of the same "
                                            "folder and system (output, phases, conditions, values); written once, "
                                            "ids in `orchid.also_reported_as`",
    "espei_evaluation_nonfinite_region": "ZPF region with a non-finite driving force against the target TDB "
                                         "(ESPEI 0.9.1 likelihood evaluation); region not written",
    "espei_evaluation_error": "dataset raises an error in the ESPEI 0.9.1 likelihood evaluation against its TDB",
    "espei_evaluation_nonfinite": "dataset gives a non-finite value in the ESPEI 0.9.1 likelihood evaluation",
    "espei_evaluation_tdb_error": "target TDB cannot be loaded by pycalphad",
    "espei_evaluation_tdb_not_found": "target TDB not found",
    "espei_evaluation_timeout": "likelihood evaluation timed out",
    "reference_phase_not_in_target_tdb": "a reference-state phase is not a phase of the target TDB (e.g. the SER "
                                         "name 1/2_MOLE_O2(G), truncated to '1' by the record builder); ESPEI 0.9.1 "
                                         "stops with a KeyError",
    "blocked_by_flag:unsigned_in_print": "the published value is flagged unsigned_in_print (minus sign absent in "
                                         "print; sign uncertain)",
    "blocked_by_flag:magnitude_suspect": "the published value is flagged magnitude_suspect",
    "blocked_by_flag:system_mismatch": "the published value is flagged system_mismatch",
    "blocked_by_flag:source_disputed": "the published value is flagged source_disputed",
    "phase_rule:temperature_implausible": "physical phase rule: the temperature is zero, a composition, or (phase "
                                          "equilibria) below 0.25 x the lowest element melting point",
    "gas_without_printed_pressure": "ZPF region with a gas phase and no printed pressure",
    "three_phase_region_off_invariant": "binary three-phase ZPF region more than %g K from every calculated "
                                        "invariant of the system (a binary three-phase equilibrium is invariant)"
                                        % PR.INVARIANT_T_TOL,
    "zpf_region_inert": "ZPF region without any composition and without a hyperplane after E2: ESPEI 0.9.1 gives "
                        "it zero driving force, so it constrains nothing (not written, not counted)",
    "duplicate_value_same_dataset": "identical conditions and value under another id in the same dataset (another "
                                    "row or paper quoting the same value): counted once, the other id in `also_read_as`",
    "activity_out_of_range": "activity outside (0, 1] (a printed scale factor such as x10^3 not applied, or a zero "
                             "endpoint that gives -inf in ESPEI); excluded, not rescaled",
    "increment_not_absolute": "printed quantity is an increment (H(T)-H(298), S(T)-S(298)), not the absolute "
                              "property of the dataset",
    "absolute_entropy_negative": "negative absolute entropy (a formation entropy or a model parameter printed as S)",
    "heat_capacity_unit_implausible": "heat capacity below %g J/K/mol-atom above 100 K (a per-gram value)" % 5.0,
    "formation_magnitude_implausible": "formation enthalpy / Gibbs energy outside [-450, +100] kJ/mol-atom "
                                       "(per-formula value or unit error)",
    "mixing_magnitude_implausible": "mixing enthalpy / Gibbs energy outside [-200, +100] kJ/mol-atom",
    "reference_state_not_ser": "compound formation value relative to a reference phase other than the TDB's SER "
                               "phase (ESPEI _FORM is relative to SER)",
}
DROP_REASONS.update(PHASE_DROP_TEXT)
PLAUSIBLE_MAX = {"SM": 150.0, "CPM": 200.0, "HM": 1e6, "GM": 1e6}
ESPEI_EXCLUSIONS = {  # substring of the ESPEI error -> reason code
    "sum to greater than one": "mole_fractions_sum_above_one",
    "float() argument must be a string or a real number, not 'NoneType'": "null_in_composition",
}
ESPEI_EXCLUSION_TEXT = {
    "mole_fractions_sum_above_one": "printed tie-line mole fractions sum to more than one (rounding, e.g. "
                                    "0.677 + 0.333); check_dataset rejects it; printed values are not altered",
    "null_in_composition": "a composition with a null (not printed) entry next to printed ones (ternary vertex or "
                           "condition); ESPEI 0.9.1 clean_dataset cannot convert it; values are not invented",
}
DROP_REASONS["duplicate_read_same_row"] = ("E4: the same printed table row read twice (same ORCHID row id `pid`, "
                                           "conditions and value); counted once, the other id kept as `also_read_as`")
DROP_REASONS["duplicate_identical_zpf"] = ("E3: identical ZPF region (same phases, temperature and compositions) "
                                           "already written from another dataset (another assessment folder or another "
                                           "paper quoting the same invariant); kept once with every source in "
                                           "`also_reported_as`")

# E5: activity reference states, audited against the printed text of every activity table in the bundle (build fails
# for an activity table not listed here). Keyed by the table (ORCHID locator asset). `phases`: the reference phase of
# each element whose activity the table gives, as the paper states it (TDB phase names); `basis`: printed (with the
# quote and where it is printed) or not_printed_implied (the paper gives activities in the liquid at temperatures where
# the element is liquid; the reference is then the pure liquid by convention, recorded as an assumption).
ACR_REFERENCE_AUDIT = {
    "doi:10.1016/j.calphad.2006.06.002/p7/tab-3": {
        "phases": {"CU": "LIQUID", "ZR": "LIQUID"}, "basis": "printed",
        "quote": "Comparison of the measured [12] and present calculated activities of the components in the Cu-Zr melt "
                 "(reference states: liquid pure components)", "where": "table 3 caption, p. 7"},
    "doi:10.1016/j.jallcom.2009.03.044/p3/tab-3": {
        "phases": {"B": "LIQUID"}, "basis": "printed",
        "quote": "The calculated boron activities with experimental data [14] (T = 1850 K; and reference states: pure "
                 "liquid components).", "where": "table 3 caption"},
    "doi:10.1016/j.jallcom.2007.04.017/p5/tab-6": {
        "phases": {"MO": "LIQUID"}, "basis": "printed",
        "quote": "Note: The reference states are the liquid phases for the Mo and Th.", "where": "table 6 note"},
    "doi:10.1016/j.calphad.2005.06.004/p4/tab-3": {
        "phases": {"YB": "LIQUID"}, "basis": "printed", "quote": "Ytterbium reference state is Yb(liquid).",
        "where": "text next to table 3"},
    "doi:10.1007/bf02667471/p3/tab-3": {
        "phases": {"GE": "LIQUID", "SI": "LIQUID"}, "basis": "printed",
        "quote": "Fig. 4 Activities of Si and Ge at 1723 K in Si-Ge liquid alloys referred to both liquid components.",
        "where": "caption of Fig. 4 (the same activities as table 3)"},
    "doi:10.1007/s11663-998-0113-z/p2/tab-iii": {
        "phases": {"TI": "BCC_A2"}, "basis": "printed", "quote": "*In both cases, the reference state is bcc.",
        "where": "table III footnote"},
    "doi:10.1361/105497198770341923/p4/tab-2": {
        "phases": {"TI": "BCC_A2"}, "basis": "printed", "quote": "Note: The reference state is bcc for both phases.",
        "where": "table 2 note"},
    "doi:10.1007/s11669-014-0313-y/p6/tab-5": {
        "phases": {"CR": "BCC_A2"}, "basis": "printed", "quote": "Fig. 2 Activities of Cr (reference state bcc_A2) ...",
        "where": "caption of Fig. 2 (the same Cr activities as table 5)"},
    "doi:10.1007/bf02892376/p5/tab-5": {
        "phases": {"CD": "LIQUID"}, "basis": "not_printed_implied",
        "quote": "Table 5 Liquid Activities ... (a) Cd activity (Cd boiling point under external pressure).",
        "where": "table 5 heading and footnote; reference state not printed; Cd is liquid at 1373 K"},
    "doi:10.1016/0040-6031(89)87205-3/p4/tab-3": {
        "phases": {"TL": "LIQUID"}, "basis": "not_printed_implied",
        "quote": "Activity of thallium in liquid Cu-Tl alloys in the temperature range 1273-1473 K",
        "where": "table 3 caption; reference state not printed; Tl is liquid at these temperatures"},
    "doi:10.1016/0040-6031(89)87205-3/p4/tab-4": {
        "phases": {"CU": "LIQUID", "TL": "LIQUID"}, "basis": "not_printed_implied",
        "quote": "Partial Gibbs free energies and activities of copper and thallium in liquid alloys at T = 1573 K",
        "where": "table 4 caption; reference state not printed; Cu and Tl are liquid at 1573 K"},
    "doi:10.1016/0925-8388(95)02178-7/p3/tab-1": {
        "phases": {"PB": "LIQUID"}, "basis": "not_printed_implied",
        "quote": "Experimental activity of lead in [Ga + Pb] liquid alloys",
        "where": "table 1 caption; reference state not printed; Pb is liquid at 1000-1143 K"},
    "doi:10.1361/105497100770340165/p4/tab-4": {
        "phases": {"CO": "LIQUID", "CU": "LIQUID"}, "basis": "not_printed_implied",
        "quote": "... activities aCo, aCu of components Co and Cu in the liquid Co-Cu system at 1537 C",
        "where": "table 4 caption; reference state not printed; Co and Cu are liquid at 1810 K"},
}
WEIGHT_FLOOR_RULE = ("optional minimum-sigma policy (E6): sigma_eff = max(sigma, sigma_floor) with sigma_floor = "
                     "ESPEI's experimental default for the property (500 J/mol for enthalpies and for the chemical "
                     "potential of activities, 0.2 J/mol/K for entropies and heat capacities); suggested_weight = "
                     "min(weight, 1). The printed sigma is unchanged (orchid_uncertainty); use suggested_weight when "
                     "printed sigmas are much smaller than the scatter between sources (e.g. formation enthalpies "
                     "printed with +-80 J/mol-atom while other measurements differ by 1-2 kJ/mol-atom)")


def acr_reference_state(out, recs, src_rel, res):
    """E5: reference state of an activity dataset = the one the paper states (ACR_REFERENCE_AUDIT)."""
    assets = sorted({((r.get("locator") or {}).get("asset")) for _, r in recs})
    if len(assets) != 1 or assets[0] not in ACR_REFERENCE_AUDIT:
        raise SystemExit("E5: activity table not audited for its reference state: %s (%s); add it to "
                         "ACR_REFERENCE_AUDIT after reading the paper" % (assets, src_rel))
    a = ACR_REFERENCE_AUDIT[assets[0]]
    el = out["output"].split("_", 1)[1].upper()
    if el not in a["phases"]:
        raise SystemExit("E5: %s: no audited reference phase for %s" % (src_rel, el))
    ph = a["phases"][el]
    rs = out.get("reference_state") or {}
    before = list(rs.get("phases") or [])
    rs["phases"] = [ph]
    rs.setdefault("conditions", {})
    out["reference_state"] = rs
    agree = []
    for _, r in recs:
        rp = ((r.get("reference_state") or {}).get(el) or {}).get("phase")
        ok = rp == ph
        agree.append(ok)
        if not ok:
            res["upstream"]["reference_state"].append({"id": r["id"], "table": assets[0], "element": el,
                                                       "record_reference_phase": rp, "stated": ph,
                                                       "basis": a["basis"], "quote": a["quote"]})
    out["orchid_reference_state"] = {"element": el, "phase": ph, "basis": a["basis"], "quote": a["quote"],
                                     "where": a["where"], "record_reference_state_agrees": all(agree),
                                     "written_before_audit": before}
    res["e5"].append({"dataset": src_rel, "before": before, "after": [ph], "changed": before != [ph],
                      "basis": a["basis"], "records_disagree": sum(1 for x in agree if not x)})


VALUE_RTOL = 1e-4
T_TOL = 0.051          # K
X_TOL = 1e-3           # mole fraction; covers the record builder's formula-rounding rule (6e-4)

# ESPEI 0.9.1 reference standard deviations (non_equilibrium_thermochemical_error.py, equilibrium_thermochemical_
# error.py, activity_error.py: "estimated from NIST TRC uncertainties", i.e. for EXPERIMENTAL data). The effective
# sigma ESPEI uses is default / weight.
ESPEI_DEFAULT_SIGMA = {"HM": 500.0, "SM": 0.2, "CPM": 0.2, "ACR": 500.0}
ESPEI_FALLBACK_SIGMA = 1.0   # ESPEI's code uses 1.0 for any other property (e.g. GM); not an uncertainty estimate
R_GAS = 8.314462618
EV_TO_J_PER_MOL = 96485.33212  # e * N_A (CODATA 2018), J/mol per eV/atom
KJ = 1000.0

# DFT calibration
P_BLOCK_ANION = {"H", "B", "C", "N", "O", "F", "Si", "P", "S", "Se", "Te", "As", "Sb", "Cl", "Br", "I", "At"}
MAGNETIC_3D = {"Cr", "Mn", "Fe", "Co", "Ni"}
CHEM_TEXT = {"metal": "all-metal intermetallic (none of the elements below)",
             "pblock_anion": "contains a p-block / anion element (H, B, C, N, O, F, Si, P, S, Se, Te, As, Sb, halogens)",
             "magnetic_3d": "contains Cr, Mn, Fe, Co or Ni (and no p-block / anion element)"}
MATCH_TEXT = {"on_hull_single": "DFT entry on the hull and the only entry at its composition",
              "polymorph_ambiguous": "on the hull, several entries (polymorphs) at the composition",
              "above_hull": "lowest DFT entry at the composition lies above the hull"}
MIN_STRATUM = 30
N_BOOT = 1000
DFT_COMMENT = "0 K DFT used as 298.15 K enthalpy"

CITE_TITLE = "ORCHID ESPEI datasets from the literature"
CITE_KEY = "orchid_espei_2026"
CITE_ABSTRACT = ("ESPEI-format datasets (thermochemical, phase-equilibrium ZPF and invariant-reaction ZPF) built from "
                 "the ORCHID literature values, one folder per evidence class (experimental, model, DFT: literature, "
                 "Materials Project, OQMD), every value traced to its published ORCHID observation id and source DOI, "
                 "with its uncertainty in physical units (printed, ESPEI's experimental default, or calibrated for "
                 "DFT against experiment) and the ESPEI weight derived from it; plus optional, tagged phase-constraint "
                 "datasets (derived invariant-temperature brackets and liquid-stability points, assessed handbook "
                 "assemblages) in the phase names of each assessment.")
CITE_PAPERS = ("the original source of every value you use (each dataset's `reference` / `bibtex`, references.bib; "
               "index.csv lists the DOIs)")


def cite_module():
    oc = orchid_cite
    if NAME not in oc.DATASETS:
        oc.THIRD_PARTY.setdefault("papers_espei", CITE_PAPERS)
        oc.DATASETS[NAME] = (CITE_TITLE, CITE_KEY, ["papers_espei", "mp", "oqmd"])
        oc.EXTRACTED.add(NAME)
        oc.ABSTRACT[NAME] = CITE_ABSTRACT
        oc.NAMES = tuple(oc.DATASETS)
    return oc


# ------------------------------------------------------------------------------------------- helpers --
def flat_values(d):
    return [x for a in d["values"] for b in a for x in b]


def unit_factor(unit):
    return KJ if (unit or "").startswith("kJ") else 1.0


def median(xs):
    s = sorted(xs)
    n = len(s)
    if not n:
        return None
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def quantile(xs, q):
    s = sorted(xs)
    if not s:
        return None
    pos = q * (len(s) - 1)
    lo = int(math.floor(pos))
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (pos - lo)


def rnd(x, n=4):
    return None if x is None else float(f"{x:.{n}g}") if abs(x) < 1e-3 and x != 0 else round(x, n)


def prop_of(output):
    base = output.split("_")[0]
    if base.startswith("ACR"):
        return "ACR"
    if base.endswith("R") and base[:-1] in ("HM", "SM", "CPM", "GM"):
        base = base[:-1]
    return base


def value_unit(output):
    p = prop_of(output)
    return {"HM": "J/mol-atom", "GM": "J/mol-atom", "SM": "J/mol-atom/K", "CPM": "J/mol-atom/K",
            "ACR": "1 (activity)"}.get(p, "as in ESPEI")


def bib_key(pid, paper):
    fam = "anon"
    auth = (paper or {}).get("authors") or []
    if auth:
        fam = auth[0].split()[-1]
    elif (paper or {}).get("citation"):
        parts = re.findall(r"[A-Za-zÀ-ɏ'-]{2,}", paper["citation"])
        fam = parts[0] if parts else "anon"
    fam = unicodedata.normalize("NFKD", fam).encode("ascii", "ignore").decode().lower()
    fam = re.sub(r"[^a-z]", "", fam) or "anon"
    year = (paper or {}).get("year") or "nd"
    return "%s%s_%s" % (fam, year, hashlib.sha256(pid.encode()).hexdigest()[:6])


def bib_escape(s):
    return str(s).replace("\\", "\\\\").replace("{", "\\{").replace("}", "\\}")


def bib_entry(key, pid, p):
    f = []
    if p.get("authors"):
        f.append(("author", " and ".join(p["authors"])))
    if p.get("title"):
        f.append(("title", "{%s}" % bib_escape(p["title"])))
    if p.get("container"):
        f.append(("journal", p["container"]))
    if p.get("year"):
        f.append(("year", str(p["year"])))
    if p.get("doi"):
        f.append(("doi", p["doi"]))
    f.append(("note", "%s (ORCHID paper id %s)" % (p.get("citation") or "", pid)))
    kind = "article" if p.get("title") and p.get("container") else "misc"
    body = ",\n".join("  %-8s = {%s}" % (k, v if k == "title" else bib_escape(v)) for k, v in f)
    return "@%s{%s,\n%s\n}" % (kind, key, body)


def short_ref(pid, p):
    p = p or {}
    fam = (p.get("authors") or [None])[0]
    lead = "%s%s" % (fam.split()[-1] if fam else "", (" " + str(p["year"])) if p.get("year") else "")
    tail = ("doi:" + p["doi"]) if p.get("doi") else (p.get("citation") or pid)
    return ("%s, %s" % (lead.strip(), tail)) if lead.strip() else tail


def chem_class(elements):
    els = {e.capitalize() for e in elements}
    if els & P_BLOCK_ANION:
        return "pblock_anion"
    if els & MAGNETIC_3D:
        return "magnetic_3d"
    return "metal"


# ------------------------------------------------------------------------------------------- inputs --
def published(lit_dir):
    pub, alias = {}, {}
    with gzip.open(os.path.join(lit_dir, "orchid-literature.jsonl.gz"), "rt", encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            pub[r["id"]] = r
            for x in r.get("merged_ids") or []:
                alias[x] = r["id"]
    papers = json.load(open(os.path.join(lit_dir, "papers.json"), encoding="utf-8"))
    index = json.load(open(os.path.join(lit_dir, "index.json"), encoding="utf-8"))
    return pub, alias, papers, index


def store_records(corpus, pdfs):
    """Internal observation id -> {id, base, pid, reader_agreement}, with the exporter's own loaders and id rule."""
    store = lit.Store(corpus, pdfs)
    for sub in ("invariants/master", "general/records", "text"):
        for f in sorted(glob.glob(os.path.join(store.db, sub, "*.jsonl"))):
            for r in lit.jl(f):
                if r.get("record_type") == "source":
                    store.add_source(r)
    contexts = lit.load_contexts(store)
    recs = []
    lit.load_tables(store, "invariants/master", "invariant_table", contexts, recs)
    lit.load_tables(store, "general/records", "general_table", contexts, recs)
    lit.load_text(store, contexts, recs)
    seen = Counter(r["id"] for r in recs)
    out = {}
    for r in recs:
        oid, rid = r["_oid"], r["id"]
        base = rid or "ORC"
        if not rid or seen[rid] > 1:
            rid = base + "." + hashlib.sha256(oid.encode()).hexdigest()[:6]
        out[oid] = {"id": rid, "base": base, "pid": r.get("pid"), "reader_agreement": r.get("reader_agreement")}
    return out, store.db


# ------------------------------------------------------------------------------------------- checks --
def check_scalar(v, rec):
    pv = rec.get("value")
    if pv is None or v is None:
        return False
    exp = pv * unit_factor(rec.get("unit"))
    return abs(v - exp) <= VALUE_RTOL * max(1.0, abs(exp))


def check_temperature(t, rec, basis):
    pt = rec.get("temperature_K")
    if pt is None or t is None:
        return pt is None
    if abs(t - pt) <= T_TOL:
        return True
    if basis and basis.get("export_basis") == "range_midpoint" and basis.get("half_width") is not None:
        return abs(abs(t - pt) - basis["half_width"]) <= T_TOL
    return False


def check_zpf_point(region, rec):
    pubx = defaultdict(list)
    for pc in rec.get("phase_compositions") or []:
        for el, x in (pc.get("composition") or {}).items():
            if isinstance(x, (int, float)):
                pubx[el.upper()].append(x)
    for vertex in region:
        if vertex[0] == "__HYPERPLANE__":
            continue
        for el, x in zip(vertex[1], vertex[2]):
            if x is None:
                continue
            if not any(abs(x - p) <= X_TOL for p in pubx.get(el, [])):
                return False
    return True


# ------------------------------------------------------------------------------------------- uncertainty --
def printed_sigma(rec, v_espei):
    """(sigma in the ESPEI value's units, printed-info dict) from a published record, or (None, None)."""
    u = rec.get("uncertainty") or {}
    if u.get("basis") not in ("printed", "printed_range") or not isinstance(u.get("value"), (int, float)):
        return None, None
    fac = unit_factor(rec.get("unit"))
    s = abs(u["value"]) * fac
    info = {"value": u["value"], "unit": rec.get("unit")}
    if u.get("basis") == "printed":
        info["interval"] = "±x printed; coverage not stated, taken as 1 sigma"
        info["coverage"] = "1sigma_assumed"
    else:
        rng = rec.get("value_range") or {}
        info["interval"] = "printed range; sigma = half-width"
        info["coverage"] = "range_half_width"
        if rng:
            info["lower"], info["upper"] = rng.get("min"), rng.get("max")
    if u.get("note"):
        info["note"] = u["note"]
    return s, info


def uncertainty_block(out, recs, cls, zpf, T):
    """Per-value sigma (physical units), basis, weight; returns (dataset weight or None, block, tags, summary)."""
    output = out["output"]
    unit = value_unit(output)
    tags = []
    if zpf:
        vals = []
        for oid, r in recs:
            tu = r.get("temperature_uncertainty") or {}
            e = {"sigma": None, "sigma_basis": "zpf (ESPEI weighs ZPF data by a driving-force tolerance, not a value "
                                               "sigma; no weight set)"}
            if tu.get("basis") in ("printed", "printed_range") and isinstance(tu.get("value"), (int, float)):
                e["temperature_sigma_K"] = tu["value"]
                e["temperature_sigma_basis"] = tu["basis"] + (" (coverage not stated, taken as 1 sigma)"
                                                              if tu["basis"] == "printed" else " (half-width)")
            vals.append(e)
        return None, {"unit": "K (temperature only)", "values": vals,
                      "weight_rule": "none: ZPF scale is ESPEI's driving-force tolerance; weight left at 1"}, tags, \
            {"basis_counts": {"zpf": len(vals)}, "sigmas": []}
    prop = prop_of(output)
    default = ESPEI_DEFAULT_SIGMA.get(prop)
    vflat = flat_values(out)
    vals, sig_eff, basis_counts, sigmas = [], [], Counter(), []
    for (oid, r), v in zip(recs, vflat):
        if cls == "model":
            s, info = printed_sigma(r, v)
            e = {"sigma": None, "sigma_basis": "none (model-calculated; no measurement uncertainty)"}
            if info:
                e["printed_uncertainty_not_used"] = info
            basis_counts["none_model"] += 1
            vals.append(e)
            continue
        s, info = printed_sigma(r, v)
        if prop == "ACR" and s is not None:
            Tv = out["conditions"]["T"]
            s_mu = R_GAS * Tv * s / abs(v) if v else None
            e = {"sigma": rnd(s, 6), "sigma_basis": "printed", "printed": info,
                 "sigma_mu_J_per_mol": rnd(s_mu, 6) if s_mu else None}
            if s_mu:
                sig_eff.append(s_mu)
            basis_counts["printed"] += 1
            sigmas.append(s)
        elif s is not None:
            e = {"sigma": rnd(s, 6), "sigma_basis": "printed", "printed": info}
            sig_eff.append(s)
            basis_counts["printed"] += 1
            sigmas.append(s)
        else:
            e = {"sigma": None, "sigma_printed": None}
            if default is not None:
                e["sigma_assumed"] = default
                e["sigma_basis"] = "espei_default_experimental (assumed, not printed)"
                sig_eff.append(default)
            else:
                e["sigma_assumed"] = None
                e["sigma_basis"] = ("none (not printed; ESPEI has no experimental default for %s, its code uses "
                                    "sigma = %g)" % (prop, ESPEI_FALLBACK_SIGMA))
            basis_counts["none_assumed_default" if default is not None else "none"] += 1
        vals.append(e)
    weight, rule = None, None
    dsig = default if default is not None else ESPEI_FALLBACK_SIGMA
    dname = ("ESPEI experimental default sigma_%s = %g %s" % (prop, dsig, "J/mol" if prop in ("HM", "ACR") else
                                                             "J/mol/K") if default is not None else
             "ESPEI fallback sigma = %g (no default for %s)" % (dsig, prop))
    if cls == "model":
        weight = 1.0
        rule = ("weight 1 only because ESPEI needs a value; model-calculated values have no measurement "
                "uncertainty: choose your own weight or exclude them with the `model` tag")
        tags.append("model")
    elif basis_counts["printed"]:
        m = median(sig_eff)
        weight = round(dsig / m, 6) if m else None
        rule = ("dataset weight = %s / median sigma of the dataset's values (printed sigma; assumed default where "
                "none was printed%s); ESPEI 0.9.1 equilibrium%s data take one weight per dataset"
                % (dname, "; activities: sigma_mu = R T sigma_a / a" if prop == "ACR" else "",
                   "/activity" if prop == "ACR" else ""))
        for e in vals:
            ref = e.get("sigma_mu_J_per_mol") if prop == "ACR" else (e.get("sigma") or e.get("sigma_assumed"))
            e["weight_per_value_equivalent"] = round(dsig / ref, 6) if ref else None
    else:
        weight = None
        rule = "no printed uncertainty: weight left at 1 (ESPEI default)"
        tags.append("no_printed_uncertainty")
    for e in vals:
        e["weight_used"] = weight if weight is not None else 1.0
    block = {"unit": unit, "default_sigma": {"value": dsig, "name": dname}, "weight_rule": rule, "values": vals}
    return weight, block, tags, {"basis_counts": dict(basis_counts), "sigmas": sigmas}


# ------------------------------------------------------------------------------------------- corpus datasets --
# ------------------------------------------------------------------------------------------- value gates --
INCREMENT_RE = re.compile(r"(?i)[-−–]\s*[HS]\s*[°o0]?\s*\(?\s*298|[HS]\s*[°o0]?\s*\(\s*T\s*\)\s*[-−–]|increment")
FORMATION_RANGE_J = (-450e3, 100e3)   # J/mol-atom; the most stable oxides are about -410 kJ/mol-atom
MIXING_RANGE_J = (-200e3, 100e3)
CPM_MIN = 5.0                         # J/K/mol-atom above 100 K: a smaller value is a per-gram value


def value_gate(output, v, rec, T):
    """Reason a scalar value cannot be used in ESPEI (physically or by unit), or None. Values are never altered."""
    prop = rec.get("property") or ""
    if output.startswith("ACR"):
        if not (0.0 < v <= 1.0):
            return "activity_out_of_range"
        return None
    if INCREMENT_RE.search(str(rec.get("property_printed") or "")):
        return "increment_not_absolute"
    if output == "SM" and v < 0:
        return "absolute_entropy_negative"
    if output == "CPM" and v < CPM_MIN and isinstance(T, (int, float)) and T > 100:
        return "heat_capacity_unit_implausible"
    if prop.endswith("_FORM") and prop[:2] in ("HM", "GM") and not (FORMATION_RANGE_J[0] <= v <= FORMATION_RANGE_J[1]):
        return "formation_magnitude_implausible"
    if prop.endswith("_MIX") and prop[:2] in ("HM", "GM") and not (MIXING_RANGE_J[0] <= v <= MIXING_RANGE_J[1]):
        return "mixing_magnitude_implausible"
    return None


# ------------------------------------------------------------------------------------------- phases --
def assign_phases(d, src_rel, stem, keep, recs, zpf, rule, reg):
    """Phase of every kept value from the published record (tools/phase_rule.py), never from the corpus dataset's
    default. Returns (groups, drops): one group per phase set (equilibrium dataset) or per compound phase
    (non-equilibrium endmember dataset); drops [(index, record, reason)]."""
    drops = []
    models = reg.models.get(stem, {})
    if zpf:
        k2, r2 = [], []
        sysn = rec_system = None
        for i, (oid, rec) in zip(keep, recs):
            why = rule_reason(rec, rule)
            region = d["values"][i]
            T = d["conditions"]["T"][i] if isinstance(d["conditions"]["T"], list) else d["conditions"]["T"]
            names = [vx[0] for vx in region if vx[0] != "__HYPERPLANE__"]
            if not why and models and any(n not in models for n in names):
                why = "phase_not_in_target_tdb"
            if not why and "GAS" in names and rec.get("pressure_basis") != "printed":
                why = "gas_without_printed_pressure"
            if not why and len(set(names)) >= 3 and len(d["components"]) - ("VA" in d["components"]) == 2 and \
                    src_rel.startswith("general/") and not rule.near_invariant(rec.get("system"), T):
                why = "three_phase_region_off_invariant"
            if not why and ec.zpf_inert(region, stem, reg):
                why = "zpf_region_inert"
            if why:
                drops.append((i, rec, why))
            else:
                k2.append(i)
                r2.append((oid, rec))
        return ([{"kind": "zpf", "keep": k2, "recs": r2, "phases": None}] if k2 else []), drops
    G = rule.grids.get(stem)
    groups = {}
    for i, (oid, rec) in zip(keep, recs):
        why = rule_reason(rec, rule)
        if why:
            drops.append((i, rec, why))
            continue
        c = rule.classify(rec)
        names, why = [], None
        for ph in c["phases"]:
            if ph["kind"] == "liquid":
                n = "LIQUID" if "LIQUID" in models else None
            else:
                cand = sorted(set(reg.phase_tdb_names.get(ph["registry_id"], {}).get(stem, [])))
                n = cand[0] if len(cand) == 1 else None
                why = "phase_several_models_in_target_tdb" if len(cand) > 1 else None
            if n is None:
                why = why or "phase_not_in_target_tdb"
                break
            names.append(n)
        if why or not names:
            drops.append((i, rec, why or "phase_not_in_target_tdb"))
            continue
        if len(names) == 1 and (models[names[0]] or {}).get("stoichiometric"):
            if rec.get("property") not in PR.FORMATION_FITTABLE:
                drops.append((i, rec, "phase_rule:compound_property_not_fittable"))
                continue
            m = models[names[0]]
            comp = {k.upper(): v for k, v in (m.get("composition") or {}).items()}
            xr = {k.upper(): v for k, v in (rec.get("overall_composition") or {}).items()}
            if len(xr) == 1 and len(comp) == 2:
                (e1, v1), = xr.items()
                xr[[e for e in comp if e != e1][0]] = 1.0 - v1
            if xr and any(abs(xr.get(e, 0.0) - comp.get(e, 0.0)) > COMPOUND_X_TOL for e in set(comp) | set(xr)):
                drops.append((i, rec, "composition_not_at_compound"))
                continue
            ems = endmembers(m, comp)
            if len(ems) != 1:
                drops.append((i, rec, "no_unique_endmember_in_target_tdb"))
                continue
            ser = rule.ser_phases(m.get("path") or "")
            refs = rec.get("reference_state") or {}
            if any(isinstance(v, dict) and v.get("phase") and v["phase"].upper() != ser.get(el.upper())
                   for el, v in refs.items()):
                drops.append((i, rec, "reference_state_not_ser"))
                continue
            g = groups.setdefault(("compound_endmember", names[0]), {
                "kind": "compound_endmember", "keep": [], "recs": [], "phases": [names[0]], "model": m,
                "configs": [], "basis": Counter(), "flags": []})
            g["configs"].append(ems[0])
        else:
            ph = sorted(set(names))
            refp = [v.get("phase") for v in (d.get("reference_states") or {}).values()] + \
                list((d.get("reference_state") or {}).get("phases") or [])
            if models and any(p_ not in models for p_ in refp):
                drops.append((i, rec, "reference_phase_not_in_target_tdb"))
                continue
            g = groups.setdefault(("equilibrium", tuple(ph)), {"kind": "equilibrium", "keep": [], "recs": [],
                                                                "phases": ph, "basis": Counter(), "flags": []})
            if G and len(d["components"]) - ("VA" in d["components"]) == 2:
                x = PR.x_of(rec, G["x_component"])
                for n in ph:
                    if rule.in_field(stem, n, x, rec.get("temperature_K")) is False:
                        g["flags"].append({"id": rec["id"], "phase": n, "x": x, "T_K": rec.get("temperature_K")})
        g["keep"].append(i)
        g["recs"].append((oid, rec))
        g["basis"][c.get("basis") or "printed_labels"] += 1
    out = []
    for key in sorted(groups):
        g = groups[key]
        g["assignment"] = {"phases": g["phases"], "basis": dict(sorted(g.pop("basis").items())),
                           "rule": "tools/phase_rule.py; README.md, 'Phase assignment'"}
        if g["flags"]:
            g["assignment"]["outside_calculated_field"] = g["flags"]
        out.append(g)
    return out, drops


def rule_reason(rec, rule):
    """Drop reason from the published record: a blocking flag, or the physical phase rule (the record's own result;
    computed here for values that were never espei_ready, e.g. model-calculated values)."""
    fl = rec.get("flags") or {}
    bl = [f for f in PR.BLOCKING_FLAGS if fl.get(f)]
    if bl:
        return "blocked_by_flag:" + bl[0]
    pr = fl.get("espei_phase_rule") or rule.classify(rec)
    return None if pr.get("result") == "pass" else "phase_rule:" + pr["reason"]


def to_solver(out, g, recs):
    """A compound's formation enthalpies / entropies as an ESPEI non-equilibrium endmember dataset."""
    m = g["model"]
    out["output"] = recs[0][1]["property"]
    out["components"] = sorted({c for cfg in g["configs"] for c in cfg})
    out["phases"] = list(g["phases"])
    out["solver"] = {"mode": "manual", "sublattice_site_ratios": [s["site_ratio"] for s in m["sublattices"]],
                     "sublattice_configurations": g["configs"]}
    out["conditions"] = {"P": 101325, "T": out["conditions"]["T"]}
    out.pop("reference_states", None)
    out.pop("reference_state", None)
    out["orchid_phase_assignment"]["endmember"] = {
        "tdb": m["tdb"], "path": m.get("path"), "sha256": m.get("sha256"), "phase": g["phases"][0],
        "reference": "SER phases of the TDB's ELEMENT lines (ESPEI _FORM convention); the value's printed "
                     "reference state names the same phases (or none)"}
    if "compound_formation_endmember" not in out["tags"]:
        out["tags"].append("compound_formation_endmember")


def certification_block(out, recs):
    vals = [{"id": r["id"], "espei_certified": bool((r.get("flags") or {}).get("espei_certified")),
             "calphad_check": (r.get("flags") or {}).get("calphad_check") or "not_available"} for _, r in recs]
    allc = bool(vals) and all(v["espei_certified"] for v in vals)
    out["orchid_certification"] = {"basis": "pycalphad_invariant_agreement", "all_values_certified": allc,
                                   "n_certified": sum(v["espei_certified"] for v in vals), "values": vals}
    if allc and "calphad_certified" not in out["tags"]:
        out["tags"].append("calphad_certified")


def build_corpus(corpus, pdfs, data, pub, alias, papers):
    RULE, REG = PR.PhaseRule(data), ec.Registry(data)
    smap, db = store_records(corpus, pdfs)
    release_has_ra = any("reader_agreement" in r for r in pub.values())
    by_base = defaultdict(list)
    for r in pub.values():
        by_base[r["id"].split(".")[0]].append(r)

    def lookup(oid):
        s = smap.get(oid)
        if s is None:
            return None, None
        for cand in (s["id"], alias.get(s["id"])):
            if cand and cand in pub:
                return pub[cand], "id"
        same_row = [r for r in by_base.get(s["base"], []) if r.get("pid") == s["pid"]]
        if len(same_row) == 1:
            return same_row[0], "row"
        return None, None

    evals = {}
    ep = os.path.join(db, "invariants", "espei_evaluation.json")
    if os.path.exists(ep):
        for tdb, e in json.load(open(ep, encoding="utf-8")).items():
            pdl, obl = e.get("per_dataset") or [], e.get("observations") or []
            for i, fn in enumerate(e.get("files") or []):
                if i < len(pdl) and i < len(obl):
                    evals[(tdb, fn)] = (obl[i], pdl[i])

    res = {"datasets": {}, "rows": {}, "drops": [], "dropped_datasets": [], "papers": {}, "src_counts": Counter(),
           "match_kind": Counter(), "sigma_stats": [], "_ustat": {}, "release_has_ra": release_has_ra, "_recs": {},
           "upstream": {"duplicate_reads": [], "reference_state": []}, "e5": [], "e6": Counter(),
           "phase_drops": Counter(), "phase_groups": Counter()}
    for sub, default_folder in SOURCES:
        for f in sorted(glob.glob(os.path.join(db, sub, "*", "*.json"))):
            if os.path.basename(f) == "provenance.json":
                continue
            src_rel = os.path.relpath(f, db).replace(os.sep, "/")
            tdb_dir, fname = f.split(os.sep)[-2], os.path.basename(f)
            d = json.load(open(f, encoding="utf-8"))
            res["src_counts"][sub] += 1
            if sub == "general/espei" and "dft" in (d.get("tags") or []):
                res["dropped_datasets"].append({"source": src_rel, "reason": ["dft_superseded"], "values": None})
                continue
            folder = default_folder
            cls = FOLDER_CLASS[folder]
            obs = re.findall(r"obs:\S+", d.get("comment") or "")
            zpf = d["output"] == "ZPF"
            if zpf:
                T = d["conditions"]["T"] if isinstance(d["conditions"]["T"], list) else [d["conditions"]["T"]]
                npts = len(T)
            else:
                T = None
                npts = len(flat_values(d))
            if len(obs) != npts:
                raise SystemExit("unexpected: %s has %d observations for %d values" % (src_rel, len(obs), npts))
            bases = {b["observation"]: b for b in d.get("orchid_value_bases") or []}
            keep, recs, reasons = [], [], []
            fv = None if zpf else flat_values(d)
            for i, oid in enumerate(obs):
                rec, how = lookup(oid)
                why = None
                if rec is None:
                    why = "not_in_release"
                elif not rec["flags"].get("extraction_accepted"):
                    why = "not_accepted"
                elif rec.get("evidence_class") != cls:
                    why = "evidence_class_mismatch"
                elif zpf:
                    b = ((bases.get(oid) or {}).get("values") or [{}])[0]
                    if not check_temperature(T[i], rec, b):
                        why = "temperature_differs_from_release"
                    elif not check_zpf_point(d["values"][i], rec):
                        why = "composition_differs_from_release"
                else:
                    if not check_scalar(fv[i], rec):
                        why = "value_differs_from_release"
                    elif abs(fv[i]) > PLAUSIBLE_MAX.get(prop_of(d["output"]), float("inf")):
                        why = "magnitude_implausible"
                    elif value_gate(d["output"], fv[i], rec, d["conditions"]["T"]):
                        why = value_gate(d["output"], fv[i], rec, d["conditions"]["T"])
                    elif not check_temperature(d["conditions"]["T"], rec, None):
                        why = "temperature_differs_from_release"
                if why:
                    res["drops"].append({"source": src_rel, "index": i, "reason": why, "orchid_id": (rec or {}).get("id")})
                    reasons.append(why)
                else:
                    keep.append(i)
                    recs.append((oid, rec))
                    res["match_kind"][how] += 1
            if not keep:
                res["dropped_datasets"].append({"source": src_rel, "reason": sorted(set(reasons)), "values": npts})
                continue
            # E4: the same printed row read twice (same pid, conditions and value) counts once
            seen_k, seen_v, keep2, recs2, also_read = {}, {}, [], [], defaultdict(list)
            for i, (oid, rec) in zip(keep, recs):
                if zpf:
                    key = (rec.get("pid"), T[i], json.dumps(d["values"][i], sort_keys=True))
                else:
                    key = (rec.get("pid"), json.dumps({k: (v[i] if isinstance(v, list) else v)
                                                       for k, v in d["conditions"].items()}, sort_keys=True), fv[i])
                vkey = key[1:]
                if not rec.get("pid") or key not in seen_k:
                    if vkey in seen_v:          # identical conditions and value under another id (another row)
                        res["drops"].append({"source": src_rel, "index": i, "reason": "duplicate_value_same_dataset",
                                             "orchid_id": rec["id"], "same_as": seen_v[vkey]})
                        also_read[seen_v[vkey]].append(rec["id"])
                        continue
                    seen_v[vkey] = rec["id"]
                if rec.get("pid") and key in seen_k:
                    first = seen_k[key]
                    also_read[first].append(rec["id"])
                    res["drops"].append({"source": src_rel, "index": i, "reason": "duplicate_read_same_row",
                                         "orchid_id": rec["id"], "same_as": first})
                    res["upstream"]["duplicate_reads"].append({"kept_id": first, "duplicate_id": rec["id"],
                                                               "pid": rec.get("pid"), "dataset": src_rel})
                    continue
                seen_k[key] = rec["id"]
                keep2.append(i)
                recs2.append((oid, rec))
            keep, recs = keep2, recs2
            # physical phase rule (tools/phase_rule.py): phases from the paper, never a default; formation
            # quantities of compounds -> the compound endmember (non-equilibrium `solver` dataset)
            groups, pdrops = assign_phases(d, src_rel, tdb_dir, keep, recs, zpf, RULE, REG)
            for i, rec, why in pdrops:
                res["drops"].append({"source": src_rel, "index": i, "reason": why, "orchid_id": rec["id"]})
                reasons.append(why)
                res["phase_drops"][why] += 1
            if not groups:
                res["dropped_datasets"].append({"source": src_rel, "reason": sorted(set(reasons)), "values": npts})
                continue
            nkept = sum(len(g["keep"]) for g in groups)
            for g in groups:
                keep, recs = g["keep"], g["recs"]
                res["phase_groups"][g["kind"]] += len(keep)
                out = json.loads(json.dumps(d))
                if len(keep) < npts:
                    if zpf:
                        out["conditions"]["T"] = [T[i] for i in keep]
                        out["values"] = [d["values"][i] for i in keep]
                        out["phases"] = sorted({v[0] for reg in out["values"] for v in reg if v[0] != "__HYPERPLANE__"})
                    else:
                        for k, v in d["conditions"].items():
                            if k.startswith("X_"):
                                out["conditions"][k] = [v[i] for i in keep]
                        out["values"] = [[[fv[i] for i in keep]]]
                if not zpf:
                    out["phases"] = list(g["phases"])
                ids = [r["id"] for _, r in recs]
                pids, measured = [], []
                for _, r in recs:
                    for pk in [(r.get("origin") or {}).get("measured_in"), r.get("paper")] + \
                              [c.get("paper") for c in r.get("cited_sources") or []]:
                        if pk and pk not in pids:
                            pids.append(pk)
                    pk = (r.get("origin") or {}).get("measured_in") or r.get("paper")
                    if pk and pk not in measured:
                        measured.append(pk)
                for pk in pids:
                    res["papers"][pk] = papers.get(pk) or {"id": pk, "citation": pk}
                dois = sorted({(papers.get(pk) or {}).get("doi") for pk in pids} - {None})
                keys = [bib_key(pk, papers.get(pk)) for pk in measured]
                ra = [r.get("reader_agreement") if release_has_ra else (smap.get(oid) or {}).get("reader_agreement")
                      for oid, r in recs]
                ra_counts = dict(sorted(Counter(x or "not_classified" for x in ra).items()))
                caveats, per_value = [], []
                for (oid, r), a in zip(recs, ra):
                    fl = r.get("flags") or {}
                    cv = [c for c in CAVEAT_FLAGS if fl.get(c)]
                    if (fl.get("plausibility") or {}).get("result") == "fail":
                        cv.append("plausibility_fail")
                    if "basis unstated" in (r.get("unit") or ""):
                        cv.append("energy_basis_unstated")
                    caveats.extend(cv)
                    per_value.append({k: v for k, v in (("id", r["id"]), ("reader_agreement", a),
                                                        ("paper", r.get("paper")),
                                                        ("measured_in", (r.get("origin") or {}).get("measured_in")),
                                                        ("caveats", cv or None),
                                                        ("also_read_as", also_read.get(r["id"]) or None)) if v is not None})
                caveats = sorted(set(caveats))
                weight, ublock, utags, ustat = uncertainty_block(out, recs, cls, zpf, T)
                tags = list(out.get("tags") or [])
                extra = [cls] + utags
                for code in ("single_reader", "single_reader_text_layer"):
                    if ra_counts.get(code):
                        extra.append(code)
                if caveats:
                    extra.append("orchid_caveat")
                extra.append("orchid")
                for t in extra:
                    if t not in tags:
                        tags.append(t)
                out["tags"] = tags
                if weight is not None:
                    out["weight"] = weight
                out["reference"] = "; ".join(short_ref(pk, papers.get(pk)) for pk in measured)
                out["bibtex"] = keys[0] if len(keys) == 1 else " ".join(keys)
                out["comment"] = ("ORCHID literature values, %s. Evidence class: %s. ORCHID observation ids "
                                  "(orchid-literature `id`, one per value, in value order): %s"
                                  % (STATUS_LINE, cls, " ".join(ids)))
                if "orchid_value_bases" in out:
                    idof = {oid: r["id"] for oid, r in recs}
                    out["orchid_value_bases"] = [dict(b, observation=idof[b["observation"]])
                                                 for b in out["orchid_value_bases"] if b["observation"] in idof]
                me = out.get("orchid_model_evidence") or {}
                ev = None
                if (tdb_dir, fname) in evals:
                    eobs, edf = evals[(tdb_dir, fname)]
                    if eobs == obs and len(keep) == npts and edf:
                        a = sorted(abs(x) for x in edf if isinstance(x, (int, float)))
                        if a:
                            ev = {"tdb": tdb_dir, "vertices": len(a), "median_abs_driving_force_J": round(a[len(a) // 2], 2),
                                  "max_abs_driving_force_J": round(a[-1], 2)}
                out["orchid"] = {k: v for k, v in (
                    ("status", STATUS_LINE), ("evidence_class", cls), ("observations", per_value),
                    ("papers", pids), ("dois", dois), ("source_dataset", "db-v2.1/" + src_rel),
                    ("values_dropped", (npts - nkept) or None), ("espei_evaluation", ev)) if v is not None}
                out["orchid_uncertainty"] = ublock
                if d["output"].startswith("ACR"):                                                         # E5
                    acr_reference_state(out, recs, src_rel, res)
                if not zpf and cls != "model" and ESPEI_DEFAULT_SIGMA.get(prop_of(d["output"])):         # E6
                    w = weight if weight is not None else 1.0
                    out["orchid_weight_floor"] = {
                        "suggested_weight": round(min(w, 1.0), 6), "weight_written": w, "changes_weight": w > 1.0,
                        "sigma_floor": ESPEI_DEFAULT_SIGMA[prop_of(d["output"])],
                        "sigma_floor_unit": "J/mol" if prop_of(d["output"]) in ("HM", "ACR") else "J/mol/K",
                        "rule": WEIGHT_FLOOR_RULE}
                    res["e6"][w > 1.0] += 1
                if not zpf:
                    out["orchid_phase_assignment"] = g["assignment"]
                if g["kind"] == "compound_endmember":
                    to_solver(out, g, recs)
                syst = lit.system_name([c for c in d["components"] if c.upper() != "VA"])
                rel = "%s/%s/%s.%s" % (folder, syst, tdb_dir, fname)
                if len(groups) > 1 or g["kind"] == "compound_endmember":
                    rel = rel[:-5] + "-" + re.sub(r"[^A-Za-z0-9_]+", "_", "_".join(g["phases"])) + ".json"
                if rel in res["datasets"]:
                    raise SystemExit("path collision: " + rel)
                res["datasets"][rel] = out
                Ts = out["conditions"]["T"] if isinstance(out["conditions"]["T"], list) else [out["conditions"]["T"]]
                nvals = len(out["values"]) if zpf else len(flat_values(out))
                res["rows"][rel] = row(rel, folder, syst, out, nvals, npts - nkept, min(Ts), max(Ts), cls, dois, pids,
                                       ids, ra_counts, "release" if release_has_ra else "record_store", caveats,
                                       tdb=(me.get("tdb") or {}).get("name"), sha=(me.get("tdb") or {}).get("sha256"),
                                       link=me.get("link_basis"), ev=ev, src="db-v2.1/" + src_rel, weight=weight,
                                       ustat=ustat, unit=ublock.get("unit"))
                res["_ustat"][rel] = ustat
                res["_recs"][rel] = [r for _, r in recs]
    return res


def row(rel, folder, syst, out, nvals, ndrop, tmin, tmax, cls, dois, pids, ids, ra_counts, ra_src, caveats,
        tdb=None, sha=None, link=None, ev=None, src=None, weight=None, ustat=None, unit=None, calib=None):
    s = (ustat or {}).get("sigmas") or []
    return {
        "path": rel, "folder": folder, "system": syst, "components": out["components"], "output": out["output"],
        "phases": out["phases"], "n_values": nvals, "n_values_dropped": ndrop, "T_min_K": tmin, "T_max_K": tmax,
        "evidence_class": cls, "reference": out["reference"], "bibtex": out.get("bibtex"), "dois": dois,
        "papers": pids, "observation_ids": ids, "reader_agreement": ra_counts, "reader_agreement_source": ra_src,
        "caveats": caveats, "tags": out["tags"], "weight": out.get("weight") if not isinstance(out.get("weight"), list)
        else "per_value", "sigma_unit": unit, "sigma_min": rnd(min(s), 6) if s else None,
        "sigma_median": rnd(median(s), 6) if s else None, "sigma_max": rnd(max(s), 6) if s else None,
        "sigma_basis_counts": (ustat or {}).get("basis_counts") or {}, "dft_calibration_stratum": calib,
        "tdb": tdb, "tdb_sha256": sha, "model_link_basis": link,
        "espei_median_abs_driving_force_J": (ev or {}).get("median_abs_driving_force_J"), "source_dataset": src,
    }


# ------------------------------------------------------------------------------------------- DFT calibration --
def stats_of(res):
    n = len(res)
    if not n:
        return None
    med = median(res)
    mad = median([abs(x - med) for x in res])
    a = [abs(x) for x in res]
    return {"n": n, "bias_mean": sum(res) / n, "bias_median": med, "sigma_robust": 1.4826 * mad,
            "q68_abs": quantile(a, 0.68), "q95_abs": quantile(a, 0.95), "rmse": math.sqrt(sum(x * x for x in res) / n)}


def boot_ci(res, label):
    rng = random.Random(int(hashlib.sha256(label.encode()).hexdigest()[:8], 16))
    n = len(res)
    acc = defaultdict(list)
    for _ in range(N_BOOT):
        s = stats_of([res[rng.randrange(n)] for _ in range(n)])
        for k in ("bias_median", "sigma_robust", "q68_abs", "q95_abs", "rmse"):
            acc[k].append(s[k])
    return {k: [round(quantile(v, 0.025), 3), round(quantile(v, 0.975), 3)] for k, v in sorted(acc.items())}


def calibration_pairs(data, pub):
    """{source: [pair]}; pair = {res kJ/mol-atom, elements, system, match, chem, sigma_exp}."""
    dp = json.load(open(os.path.join(data, "literature", "dft_parity.json"), encoding="utf-8"))
    R = dp["records"]

    def default(r):
        return (r["ev"] == "experimental" and r["acc"] and not r["unsigned"] and not r["sign_lost"]
                and "magnitude_suspect" not in r["suspect"])

    def sig_exp(rid):
        u = (pub.get(rid) or {}).get("uncertainty") or {}
        return abs(u["value"]) if u.get("basis") == "printed" and isinstance(u.get("value"), (int, float)) else None

    pairs = defaultdict(list)
    for p in sorted(dp["points"], key=lambda p: (p["src"], p["rid"], p["mid"])):
        r = R.get(p["rid"])
        if not r or not default(r):
            continue
        els = sorted(e.upper() for e in r["comp"])
        match = "above_hull" if p["ehull"] > 0 else ("polymorph_ambiguous" if p["npoly"] > 1 else "on_hull_single")
        pairs[p["src"]].append({"res": p["res"], "system": "-".join(els), "match": match, "chem": chem_class(els),
                                "sigma_exp": sig_exp(p["rid"]), "rid": p["rid"]})
    # literature DFT vs experiment: same binary compound (reduced composition) in the default experimental subset
    exp = defaultdict(list)
    for rid, r in sorted(R.items()):
        if default(r) and len(r["comp"]) == 2:
            exp[tuple(sorted((k, round(v, 4)) for k, v in r["comp"].items()))].append((r["lit"], sig_exp(rid)))
    for r in sorted(pub.values(), key=lambda x: x["id"]):
        if not lit_dft_candidate(r):
            continue
        counts, _ = compound_counts(r)
        if counts is None:
            continue
        key = parity.reduced_key(counts)
        if key not in exp:
            continue
        ev = [x for x, _ in exp[key]]
        se = [s for _, s in exp[key] if s is not None]
        els = sorted(e.upper() for e in counts)
        pairs["literature"].append({"res": r["value"] - median(ev), "system": "-".join(els), "match": None,
                                    "chem": chem_class(els), "sigma_exp": median(se) if se else None, "rid": r["id"]})
    return pairs, dp


def lit_dft_candidate(r):
    return (r.get("evidence_class") == "dft" and r.get("property") == "HM_FORM" and r.get("unit") == "kJ/mol-atom"
            and r["flags"].get("extraction_accepted") and len(r.get("elements") or []) == 2
            and isinstance(r.get("value"), (int, float)) and not isinstance(r.get("value"), bool))


def compound_counts(r):
    """Composition of a literature record's compound (the DFT benchmark's parser); (counts, label) or (None, why)."""
    system = [e.capitalize() for e in r.get("elements") or []]
    found, reasons = {}, []
    for p in r.get("phases") or []:
        lab = p.get("label")
        if not lab:
            continue
        c, n = parity.parse_phase(lab, system)
        if c:
            found.setdefault(parity.reduced_key(c), (c, lab))
        else:
            reasons.append(n)
    if len(found) == 1:
        c, lab = next(iter(found.values()))
        return c, lab
    if len(found) > 1:
        return None, "several distinct formulas among the record's phases"
    return None, (sorted(set(reasons)) or ["no phase label"])[0]


def calibrate(pairs, needed):
    """Calibration tables; `needed` = {(source, stratum_key, system)} for leave-system-out sigmas."""
    cal = {"sources": {}, "lso": {}}
    for src in sorted(pairs):
        P = pairs[src]
        entry = {"all": None, "strata": {}}

        def summarize_set(ps, label):
            res = [p["res"] for p in ps]
            s = stats_of(res)
            if s is None:
                return None
            s = {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()}
            s["ci95"] = boot_ci(res, label)
            se = [p["sigma_exp"] for p in ps if p["sigma_exp"] is not None]
            s["sigma_exp_median"] = round(median(se), 3) if se else None
            s["n_with_printed_exp_sigma"] = len(se)
            s["sigma_deconvolved"] = (round(math.sqrt(max(s["sigma_robust"] ** 2 - s["sigma_exp_median"] ** 2, 0.0)), 3)
                                      if se else None)
            s["rmse_over_sigma_robust"] = round(s["rmse"] / s["sigma_robust"], 2) if s["sigma_robust"] else None
            return s
        entry["all"] = summarize_set(P, src + "|all")
        keys = sorted({stratum_key(p) for p in P})
        for k in keys:
            entry["strata"][k] = summarize_set([p for p in P if stratum_key(p) == k], src + "|" + k)
        for k in (sorted(CHEM_TEXT) if src != "literature" else []):
            entry["strata"].setdefault("chem:" + k, summarize_set([p for p in P if p["chem"] == k], src + "|chem:" + k))
        if src != "literature":
            for k in sorted(MATCH_TEXT):
                entry["strata"].setdefault("match:" + k, summarize_set([p for p in P if p["match"] == k],
                                                                       src + "|match:" + k))
        cal["sources"][src] = entry
    for src, sk, system in sorted(needed):
        P = pairs.get(src) or []
        own = [p for p in P if stratum_key(p) == sk and p["system"] != system]
        fallback = None
        if len(own) < MIN_STRATUM:
            own = [p for p in P if p["system"] != system]
            fallback = "source_level (stratum has < %d pairs after leaving the system out)" % MIN_STRATUM
        s = stats_of([p["res"] for p in own])
        cal["lso"]["%s|%s|%s" % (src, sk, system)] = {
            "source": src, "stratum": sk, "system": system, "n": len(own),
            "sigma_robust_kJ": round(s["sigma_robust"], 3) if s else None,
            "bias_median_kJ": round(s["bias_median"], 3) if s else None, "fallback": fallback}
    return cal


def stratum_key(p):
    return ("%s|%s" % (p["match"], p["chem"])) if p["match"] else p["chem"]


# ------------------------------------------------------------------------------------------- DFT datasets --
def registry(data):
    phases_by_sys, by_dft, by_paper = {}, {}, defaultdict(list)
    for f in sorted(glob.glob(os.path.join(data, "phases", "systems", "*.json"))):
        d = json.load(open(f, encoding="utf-8"))
        phases_by_sys[d["system"]] = d["phases"]
        for p in d["phases"]:
            for n in p.get("names") or []:
                if n.get("namespace") in ("mp", "oqmd"):
                    by_dft.setdefault(n["name"].split(":", 1)[1], []).append(p)
                elif n.get("namespace") == "paper":
                    by_paper[n["name"]].append(p)
    return phases_by_sys, by_dft, by_paper


def endmembers(model, comp):
    """Endmember configurations of a TDB phase model whose composition equals `comp` ({EL: x})."""
    subl = model.get("sublattices") or []
    if not subl:
        return []
    els = set(comp)
    choices = [[c.upper() for c in s["constituents"]] for s in subl]
    if any(not ch for ch in choices):
        return []
    out = []

    def rec(i, cfg):
        if i == len(choices):
            amt = defaultdict(float)
            for c, s in zip(cfg, subl):
                if c != "VA":
                    amt[c] += float(s["site_ratio"])
            tot = sum(amt.values())
            if tot <= 0 or set(amt) != els:
                return
            if all(abs(amt[e] / tot - comp[e]) <= 1e-4 for e in els):
                out.append(list(cfg))
            return
        for c in choices[i]:
            if c != "VA" and c not in els:
                continue
            rec(i + 1, cfg + [c])
    rec(0, [])
    return out


def map_compound(phases, comp):
    """phases (registry entries already chosen) x compound -> ([(phase, model, config)], [(phase, model, reason)])."""
    ok, bad = [], []
    for p in phases:
        for m in sorted(p.get("models") or [], key=lambda m: (m["tdb"], m["phase_name"])):
            ems = endmembers(m, comp)
            if len(ems) == 1:
                ok.append((p, m, ems[0]))
            elif not ems:
                bad.append((p, m, "no_endmember_at_composition"))
            else:
                bad.append((p, m, "several_endmembers_at_composition"))
    return ok, bad


def build_dft(data, pub, papers, pairs):
    phases_by_sys, by_dft, by_paper = registry(data)
    items = []       # (source, system, comp, value_J, provenance dict, match quality, flags)
    excluded = []    # {source, system, compound, id, reason}
    for src in ("mp", "oqmd"):
        seen = set()
        for f in sorted(glob.glob(os.path.join(data, "hulls", "dft-" + src, "*", "*.json"))):
            d = json.load(open(f, encoding="utf-8"))
            sid = d["system_id"]
            if sid in seen or len(d["elements"]) != 2:
                continue
            seen.add(sid)
            xel = d["x_element"]
            other = [e for e in d["elements"] if e != xel][0]
            for iso in d["isotherms"]:
                for v in iso["vertices"]:
                    if v["x"] in (0.0, 1.0):
                        continue
                    comp = {xel: v["x"], other: 1.0 - v["x"]}
                    cands = by_dft.get(v["phase"]) or []
                    base = {"source": src, "system": sid, "compound": v["formula"], "id": v["phase"]}
                    if not cands:
                        excluded.append(dict(base, reason="not_in_phase_registry (no registry phase carries this "
                                                          "DFT entry)"))
                        continue
                    if len(cands) > 1:
                        excluded.append(dict(base, reason="several_registry_phases"))
                        continue
                    items.append({"source": src, "system": sid, "comp": comp, "phases": cands,
                                  "value": v["energy"] * EV_TO_J_PER_MOL, "printed": {"value": v["energy"],
                                                                                      "unit": "eV/atom"},
                                  "id": v["phase"], "compound": v["formula"], "sg": v.get("spacegroup_symbol"),
                                  "match": "on_hull_single" if v.get("entries_at_x", 1) <= 1 else "polymorph_ambiguous",
                                  "snapshot": d["provenance"].get("snapshot_sha256"),
                                  "release": d["provenance"].get("mp_release"), "theoretical": v.get("theoretical")})
    for r in sorted(pub.values(), key=lambda x: x["id"]):
        if not lit_dft_candidate(r):
            continue
        sid = r["system"]
        base = {"source": "literature", "system": sid, "id": r["id"]}
        counts, lab = compound_counts(r)
        if counts is None:
            excluded.append(dict(base, compound=None, reason="composition_not_parsed (%s)" % lab))
            continue
        tot = sum(counts.values())
        comp = {e.upper(): float(Fraction(c) / Fraction(tot)) for e, c in counts.items()}
        pk = r.get("paper") or ""
        ns = ("paper:" + pk[4:]) if pk.startswith("doi:") else ("paper:sha256-" + pk[7:] if pk.startswith("sha256:")
                                                                 else None)
        labels = [p.get("label") for p in r.get("phases") or [] if p.get("label")]
        cands = {}
        for lb in labels:
            for p in by_paper.get("%s:%s" % (ns, lb), []) if ns else []:
                cands[p["id"]] = p
        basis = "paper_label"
        if not cands:
            mph = {p.get("model_phase") for p in r.get("phases") or [] if p.get("model_phase")}
            for p in phases_by_sys.get(sid, []):
                if mph & set(p.get("tdb_names") or []):
                    cands[p["id"]] = p
            basis = "model_phase_name"
        compound = parity.formula_str(counts)
        if not cands:
            excluded.append(dict(base, compound=compound, reason="label_not_mapped_to_registry_phase"))
            continue
        if len(cands) > 1:
            excluded.append(dict(base, compound=compound, reason="several_registry_phases"))
            continue
        fl = r.get("flags") or {}
        items.append({"source": "literature", "system": sid, "comp": comp, "phases": list(cands.values()),
                      "value": r["value"] * KJ, "printed": {"value": r["value"], "unit": r["unit"]}, "id": r["id"],
                      "compound": compound, "label_basis": basis, "paper": r.get("paper"), "match": None,
                      "measured_in": (r.get("origin") or {}).get("measured_in"),
                      "caveats": [c for c in CAVEAT_FLAGS if fl.get(c)]})
    # map to TDB phase endmembers; group by (source, system, TDB, phase)
    groups = defaultdict(list)
    for it in items:
        ok, bad = map_compound(it["phases"], it["comp"])
        for p, m, why in bad:
            excluded.append({"source": it["source"], "system": it["system"], "compound": it["compound"], "id": it["id"],
                             "tdb": m["tdb"], "phase": m["phase_name"], "reason": why})
        for p, m, cfg in ok:
            groups[(it["source"], it["system"], m["tdb"], m["phase_name"])].append((it, p, m, cfg))
    needed = set()
    for (src, sid, tdb, ph), lst in groups.items():
        for it, p, m, cfg in lst:
            needed.add((src, dft_stratum(it), sid))
    cal = calibrate(pairs, needed)
    datasets, rows, used_papers = {}, {}, {}
    for (src, sid, tdb, ph) in sorted(groups):
        lst = sorted(groups[(src, sid, tdb, ph)], key=lambda t: (t[3], t[0]["id"]))
        # one value per configuration: several values for one endmember stay separate rows of `values`? ESPEI
        # takes configurations x values; duplicates (same config) are allowed as separate configurations.
        m = lst[0][2]
        comps = sorted({c for _, _, _, cfg in lst for c in cfg})
        cfgs = [cfg for _, _, _, cfg in lst]
        vals = [round(it["value"], 2) for it, _, _, _ in lst]
        unc_vals, weights, strata = [], [], []
        for it, p, mm, cfg in lst:
            sk = dft_stratum(it)
            L = cal["lso"]["%s|%s|%s" % (src, sk, sid)]
            full = (cal["sources"][src]["strata"].get(sk) if not L["fallback"] else cal["sources"][src]["all"]) or {}
            sJ = L["sigma_robust_kJ"] * KJ if L["sigma_robust_kJ"] else None
            w = round(500.0 / sJ, 6) if sJ else 1.0
            weights.append(w)
            strata.append(sk if not L["fallback"] else "source_level")
            e = {"sigma": round(sJ, 1) if sJ else None, "sigma_basis": "dft_calibrated",
                 "calibration": {"source": src, "stratum": sk, "used": "source_level" if L["fallback"] else "stratum",
                                 "leave_system_out": True, "system_left_out": sid, "n": L["n"],
                                 "sigma_robust_kJ": L["sigma_robust_kJ"], "bias_median_kJ_not_applied": L["bias_median_kJ"],
                                 "ci95_sigma_robust_kJ_full_stratum": (full.get("ci95") or {}).get("sigma_robust"),
                                 "n_full_stratum": full.get("n")},
                 "weight_used": w,
                 "weight_meaning": "ESPEI's experimental reference sigma_HM (500 J/mol) / calibrated DFT sigma: the "
                                   "conversion that makes ESPEI use the calibrated sigma; 500 J/mol is not a DFT "
                                   "uncertainty",
                 "dft_value_printed": it["printed"], "id": it["id"], "compound": it["compound"]}
            unc_vals.append(e)
        tags = ["dft_" + src, "orchid"]
        if any(it.get("match") == "polymorph_ambiguous" for it, _, _, _ in lst):
            tags.append("polymorph_ambiguous")
        if any(it.get("caveats") for it, _, _, _ in lst):
            tags.append("orchid_caveat")
        if src == "literature":
            pids = []
            for it, _, _, _ in lst:
                for pk in (it.get("measured_in"), it.get("paper")):
                    if pk and pk not in pids:
                        pids.append(pk)
            for pk in pids:
                used_papers[pk] = papers.get(pk) or {"id": pk, "citation": pk}
            measured = []
            for it, _, _, _ in lst:
                pk = it.get("measured_in") or it.get("paper")
                if pk and pk not in measured:
                    measured.append(pk)
            reference = "; ".join(short_ref(pk, papers.get(pk)) for pk in measured)
            bib = " ".join(bib_key(pk, papers.get(pk)) for pk in measured)
            dois = sorted({(papers.get(pk) or {}).get("doi") for pk in pids} - {None})
        else:
            pids = []
            reference = {"mp": "Materials Project (Jain et al. 2013), doi:10.1063/1.4812323, release %s",
                         "oqmd": "OQMD (Saal et al. 2013; Kirklin et al. 2015), doi:10.1007/s11837-013-0755-4, "
                                 "release %s"}[src] % lst[0][0].get("release")
            bib = {"mp": "jain2013_mp", "oqmd": "saal2013_oqmd kirklin2015_oqmd"}[src]
            dois = {"mp": ["10.1063/1.4812323"], "oqmd": ["10.1007/s11837-013-0755-4", "10.1038/npjcompumats.2015.10"]}[src]
        ids = [it["id"] for it, _, _, _ in lst]
        ds = {
            "components": comps, "phases": [ph], "solver": {
                "mode": "manual", "sublattice_site_ratios": [s["site_ratio"] for s in m["sublattices"]],
                "sublattice_configurations": cfgs},
            "conditions": {"P": 101325, "T": 298.15}, "output": "HM_FORM", "values": [[vals]],
            "reference": reference, "bibtex": bib, "weight": weights, "tags": tags,
            "comment": ("%s. ORCHID DFT formation enthalpies (%s), %s; TDB %s (sha256 %s), phase %s. Ids, in value "
                        "order: %s. eV/atom x %.5f = J/mol-atom." % (DFT_COMMENT, src, STATUS_LINE if src == "literature"
                                                                   else "third-party DFT database", tdb, m["sha256"],
                                                                   ph, " ".join(ids), EV_TO_J_PER_MOL)),
            "orchid": {"evidence_class": "dft", "dft_source": src, "tdb": {"name": tdb, "path": m.get("path"),
                                                                             "sha256": m["sha256"], "phase": ph},
                       "registry_phase": lst[0][1]["id"],
                       "mapping": [{"id": it["id"], "compound": it["compound"],
                                    "basis": ("registry DFT name %s:%s" % (src, it["id"])) if src != "literature"
                                    else "registry %s" % it.get("label_basis"),
                                    "space_group": it.get("sg"), "match_quality": it.get("match"),
                                    **({"caveats": it["caveats"]} if it.get("caveats") else {})}
                                   for it, _, _, _ in lst],
                       "papers": pids or None, "dois": dois,
                       "status": STATUS_LINE if src == "literature" else "third-party DFT (published ORCHID hull layer)",
                       "temperature_note": DFT_COMMENT},
            "orchid_uncertainty": {"unit": "J/mol-atom", "default_sigma": {"value": 500.0, "name":
                                   "ESPEI experimental reference sigma_HM = 500 J/mol (used only to convert the "
                                   "calibrated DFT sigma into an ESPEI weight)"},
                                   "weight_rule": "weight per value = 500 J/mol / calibrated sigma_DFT (leave-system-"
                                                  "out robust sigma of the value's stratum; bias recorded, not "
                                                  "applied)", "values": unc_vals},
        }
        ds["orchid"] = {k: v for k, v in ds["orchid"].items() if v is not None}
        folder = "dft/" + src
        rel = "%s/%s/%s.%s.json" % (folder, sid, tdb[:-4] if tdb.lower().endswith(".tdb") else tdb, ph)
        datasets[rel] = ds
        sig = [e["sigma"] for e in unc_vals if e["sigma"]]
        rows[rel] = row(rel, folder, sid, ds, len(vals), 0, 298.15, 298.15, "dft", dois, pids, ids, {},
                        "not_applicable", sorted({c for it, _, _, _ in lst for c in it.get("caveats") or []}),
                        tdb=tdb, sha=m["sha256"], link="registry", src=src, ustat={"sigmas": sig,
                        "basis_counts": {"dft_calibrated": len(sig)}}, unit="J/mol-atom",
                        calib=";".join(sorted(set(strata))))
    return {"datasets": datasets, "rows": rows, "excluded": excluded, "calibration": cal, "papers": used_papers,
            "items": Counter(it["source"] for it in items)}


def dft_stratum(it):
    els = sorted(it["comp"])
    ch = chem_class(els)
    return ("%s|%s" % (it["match"], ch)) if it.get("match") else ch


# ------------------------------------------------------------------------------------------- ESPEI --
VALIDATOR = r'''
import json, os, sys, copy
import espei
from espei.datasets import check_dataset, clean_dataset, load_datasets, recursive_glob
root, mode = sys.argv[1], sys.argv[2]
groups = json.loads(sys.argv[3])
out = {"espei_version": espei.__version__, "invalid": {}, "folders": {}}
for f in recursive_glob(root, "*.json"):
    rel = os.path.relpath(f, root)
    try:
        d = json.load(open(f))
        check_dataset(copy.deepcopy(d))
        clean_dataset(copy.deepcopy(d))
    except Exception as e:
        out["invalid"][rel] = "%s: %s" % (type(e).__name__, e)
if mode == "load":
    for g in groups:
        p = os.path.join(root, g)
        if os.path.isdir(p):
            files = recursive_glob(p, "*.json")
            try:
                db = load_datasets(files)
                out["folders"][g] = {"files": len(files), "loaded": len(db.all())}
            except Exception as e:
                out["folders"][g] = {"files": len(files), "error": "%s: %s" % (type(e).__name__, e)}
print(json.dumps(out, sort_keys=True))
'''


def run_espei(espei_py, payload, mode):
    tmp = tempfile.mkdtemp(prefix="orchid-espei-")
    try:
        for rel, data in payload.items():
            p = os.path.join(tmp, rel)
            os.makedirs(os.path.dirname(p), exist_ok=True)
            with open(p, "wb") as fh:
                fh.write(data)
        r = subprocess.run([espei_py, "-I", "-c", VALIDATOR, tmp, mode, json.dumps(GROUPS + PGROUPS)], capture_output=True,
                           text=True, cwd=tmp)
        if r.returncode != 0:
            raise SystemExit("ESPEI validator failed:\n" + r.stderr[-3000:])
        return json.loads(r.stdout.strip().splitlines()[-1])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def dumps(obj):
    return (json.dumps(obj, indent=1, ensure_ascii=False) + "\n").encode("utf-8")


# ------------------------------------------------------------------------------------------- docs --
def table(rows, cols):
    o = ["| " + " | ".join(c for c, _ in cols) + " |", "|" + "|".join("---:" if r else "---" for _, r in cols) + "|"]
    for r in rows:
        o.append("| " + " | ".join("" if x is None else str(x) for x in r) + " |")
    return "\n".join(o)


def summarize(C, D, val):
    rows = list(C["rows"].values()) + list(D["rows"].values())
    by_fo = Counter((r["folder"], r["output"]) for r in rows)
    by_f, vals_f, systems = Counter(), Counter(), defaultdict(set)
    for r in rows:
        by_f[r["folder"]] += 1
        vals_f[r["folder"]] += r["n_values"]
        systems[r["folder"]].add(r["system"])
    ra = Counter()
    for r in rows:
        for k, v in r["reader_agreement"].items():
            ra[(r["folder"], k)] += v
    dft_ex = Counter((e["source"], e["reason"].split(" ")[0]) for e in D["excluded"])
    return {
        "datasets": len(rows), "values": sum(vals_f.values()),
        "by_folder": {f: {"datasets": by_f[f], "values": vals_f[f], "systems": len(systems[f])} for f in GROUPS},
        "by_folder_output": {"%s/%s" % k: v for k, v in sorted(by_fo.items())},
        "systems_total": len({r["system"] for r in rows}),
        "source_datasets": dict(sorted(C["src_counts"].items())),
        "datasets_dropped_by_gate": dict(sorted(Counter(";".join(d["reason"]) for d in C["dropped_datasets"]).items())),
        "datasets_excluded_by_espei": dict(sorted(Counter(e["reason"] for e in val["excluded"]).items())),
        "values_dropped": dict(sorted(Counter(d["reason"] for d in C["drops"]).items())),
        "reader_agreement": {"%s/%s" % k: v for k, v in sorted(ra.items())},
        "reader_agreement_source": "release" if C["release_has_ra"] else "record_store",
        "observation_match": dict(sorted(C["match_kind"].items())),
        "datasets_with_espei_evaluation": sum(1 for r in rows if r["espei_median_abs_driving_force_J"] is not None),
        "dft_candidates": dict(sorted(D["items"].items())),
        "dft_values_mapped": dict(sorted(Counter(r["folder"].split("/")[1] for r in D["rows"].values()
                                                  for _ in range(r["n_values"])).items())),
        "dft_exclusions": {"%s/%s" % k: v for k, v in sorted(dft_ex.items())},
    }


def validation_md(S, val, C, D, release):
    L = ["# ORCHID ESPEI bundle: validation", "",
         "Generated by `tools/build_espei_bundle.py` (deterministic; `--check` verifies the files). Literature release "
         "used as the inclusion gate: `%s`." % release, "",
         "## ESPEI", "",
         "Every file was checked with `espei.datasets.check_dataset` and `clean_dataset` (ESPEI %s), and each folder "
         "was loaded with `espei.datasets.load_datasets`. The build fails on any invalid file other than the "
         "exclusion classes below. Extra keys (`orchid`, `orchid_uncertainty`, `bibtex`, ...) are accepted by ESPEI."
         % val["espei_version"], "",
         table([[g, val["folders"][g].get("files"), val["folders"][g].get("loaded"),
                 "ok" if val["folders"][g].get("loaded") == val["folders"][g].get("files") else "FAIL"]
                for g in GROUPS + PGROUPS if g in val["folders"]],
               [("folder", 0), ("files", 1), ("loaded by load_datasets", 1), ("result", 0)]), "",
         "Invalid files in the final bundle: **%d**." % len(val["invalid_final"]), "",
         "### Excluded because ESPEI 0.9.1 rejects them (printed values are never altered)", ""]
    for code, txt in ESPEI_EXCLUSION_TEXT.items():
        L.append("- `%s`: %s." % (code, txt))
    L.append("")
    if val["excluded"]:
        L.append(table([[e["reason"], e["path"], e["source"]] for e in val["excluded"]],
                       [("reason", 0), ("path (not written)", 0), ("source", 0)]))
    else:
        L.append("None.")
    le = val.get("likelihood_evaluation") or {}
    L += ["", "### Likelihood evaluation against each target TDB (ESPEI 0.9.1)", "", EVAL_RULE, "",
          "Datasets evaluated in the final bundle: **%d**; all finite: **%s**." % (
              le.get("datasets_evaluated", 0), "yes" if le.get("all_finite") else "NO"), ""]
    if le.get("by_folder_kind"):
        L.append(table([[k.split("|")[0], k.split("|")[1], v["datasets"], v["values"], v["median_of_dataset_median_abs"],
                         v["p90_of_dataset_median_abs"]] for k, v in le["by_folder_kind"].items()],
                       [("folder", 0), ("quantity", 0), ("datasets", 1), ("values", 1),
                        ("median of dataset median abs", 1), ("p90 of dataset median abs", 1)]))
        L += ["", "(Residuals: J/mol-atom, J/K/mol-atom or activity; driving forces J/mol; non-equilibrium: log "
                  "likelihood. Against the TDB's own parameters, i.e. how far the published values lie from the "
                  "assessment the folder is named after.)", ""]
    ex2 = val.get("evaluation_excluded") or []
    L.append("Not written after the evaluation: %d datasets (%s); ZPF regions cut: %d in %d datasets." % (
        len(ex2), ", ".join("%s %d" % kv for kv in sorted(Counter(e["reason"] for e in ex2).items())) or "none",
        sum(c["regions"] for c in val.get("evaluation_regions_cut") or []), len(val.get("evaluation_regions_cut") or [])))
    if ex2:
        L += ["", table([[e["reason"], e["path"], (e.get("error") or "")[:120]] for e in ex2],
                        [("reason", 0), ("path (not written)", 0), ("error", 0)])]
    fx = S["fixes"]
    L += ["", "## Encoding fixes (E1-E6)", "",
          table([["E1", "VA added to `components` (target TDB has a vacancy sublattice; else system TDBs; else "
                  "metallic system)", fx["E1_va_added_datasets"]],
                 ["E2", "line-compound ZPF vertices written with composition null (datasets)", "%d vertices (%d datasets; "
                  "%d hyperplanes added)" % (fx["E2_vertices_nulled"], fx["E2_datasets"], fx["E2_hyperplanes_added"])],
                 ["E3", "identical ZPF regions written once (regions removed; datasets left empty and dropped)",
                  "%d (%d)" % (fx["E3_regions_removed"], fx["E3_datasets_emptied_and_dropped"])],
                 ["E4", "values read twice from the same printed row, counted once", fx["E4_values_removed"]],
                 ["E5", "activity datasets audited / reference state changed to the stated one / literature records "
                  "listed for correction", "%d / %d / %d" % (fx["E5_activity_datasets_audited"],
                                                              fx["E5_reference_state_changed"],
                                                              fx["E5_records_disagreeing_with_paper"])],
                 ["E6", "thermochemical datasets with a minimum-sigma weight suggestion / where it changes the weight",
                  "%d / %d" % (fx["E6_datasets_with_floor_suggestion"], fx["E6_datasets_where_floor_changes_weight"])]],
                [("fix", 0), ("what", 0), ("count", 1)]), "",
          "Records to correct in the literature release: UPSTREAM_CORRECTIONS.md.", "",
          "## Phase constraints (optional)", "",
          table([[k, v["datasets"], v["points"], v["systems"]] for k, v in S["phase_constraints"]["folders"].items()],
                [("folder", 0), ("datasets", 1), ("ZPF points", 1), ("systems", 1)]), "",
          "Invariant brackets written: %s; skipped: %s." % (
              ", ".join("%s %d" % kv for kv in S["phase_constraints"]["brackets"]["made"].items()),
              ", ".join("%s %d" % kv for kv in S["phase_constraints"]["brackets"]["skipped"].items())),
          "Assessed points left out (phase not in the TDB, or ambiguous): %s." % (
              ", ".join("%s %d" % kv for kv in S["phase_constraints"]["assessed_points_left_out"].items()) or "none"),
          "Liquid-stability systems: %s." % ", ".join("%s %d" % kv for kv in
                                                     S["phase_constraints"]["liquid_stability_systems"].items())]
    L += ["", "## Inclusion gate (same rules as the published literature values)", "",
          "| reason | values removed |", "|---|---:|"]
    for k, v in S["values_dropped"].items():
        L.append("| %s: %s | %d |" % (k, DROP_REASONS[k], v))
    L += ["", "Datasets dropped (listed in `validation.json`, `dropped_datasets`):", ""]
    for k, v in S["datasets_dropped_by_gate"].items():
        L.append("- `%s`: %d%s" % (k, v, " (the corpus's DFT datasets are replaced by the endmember HM_FORM datasets "
                                         "in `dft/literature/`, written only where a TDB endmember is unambiguous)"
                                   if k == "dft_superseded" else ""))
    L += ["", "Datasets kept with some values removed: **%d** (`orchid.values_dropped`, index `n_values_dropped`). "
          "Observation match: %s (`id`: same ORCHID id as in the release; `row`: the only published value of the "
          "same table row, value check passed)." % (
              sum(1 for r in C["rows"].values() if r["n_values_dropped"]),
              ", ".join("%s %d" % kv for kv in S["observation_match"].items())), "",
          "## DFT endmember datasets", "",
          "Candidates (on-hull binary compounds of the published hull layers; accepted literature DFT HM_FORM values in "
          "kJ/mol-atom): %s. Values written (one per TDB model that has exactly one endmember at the compound's "
          "composition): %s. Exclusions (compound or compound x TDB model) by reason:" % (
              ", ".join("%s %d" % kv for kv in S["dft_candidates"].items()),
              ", ".join("%s %d" % kv for kv in S["dft_values_mapped"].items())), ""]
    L.append(table([[k.split("/")[0], k.split("/")[1], v] for k, v in S["dft_exclusions"].items()],
                   [("source", 0), ("reason", 0), ("count", 1)]))
    L += ["", "Every exclusion is listed in `validation.json` (`dft_excluded`).", "",
          "## Self-check against the corpus ESPEI evaluation", "",
          "%d invariant ZPF datasets carry the corpus driving-force evaluation against their assessment "
          "(`orchid.espei_evaluation`; only datasets kept unchanged)." % S["datasets_with_espei_evaluation"], ""]
    return "\n".join(L) + "\n"


def sigma_table(rows_iter):
    """rows of (folder, output, basis_counts, sigmas) -> summary table rows."""
    agg = defaultdict(lambda: [Counter(), []])
    for folder, output, bc, sig in rows_iter:
        a = agg[(folder, output)]
        a[0].update(bc)
        a[1].extend(sig)
    out = []
    for (folder, output) in sorted(agg):
        bc, sig = agg[(folder, output)]
        n = sum(bc.values())
        out.append([folder, output, n, bc.get("printed", 0), bc.get("dft_calibrated", 0),
                    bc.get("none_assumed_default", 0), bc.get("none", 0), bc.get("none_model", 0), bc.get("zpf", 0),
                    value_unit(output) if output != "ZPF" else "",
                    rnd(min(sig), 4) if sig else "", rnd(quantile(sig, .25), 4) if sig else "",
                    rnd(median(sig), 4) if sig else "", rnd(quantile(sig, .75), 4) if sig else "",
                    rnd(max(sig), 4) if sig else ""])
    return out


def uncertainty_md(C, D):
    rows = [(f, o, s["basis_counts"], s["sigmas"]) for f, o, s in C["sigma_stats"]]
    rows += [(r["folder"], r["output"], r["sigma_basis_counts"],
              [e["sigma"] for e in D["datasets"][r["path"]]["orchid_uncertainty"]["values"] if e["sigma"]])
             for r in D["rows"].values()]
    T = sigma_table(rows)
    return "\n".join([
        "# ORCHID ESPEI bundle: uncertainty summary", "",
        "Per folder and output: how many values carry a printed uncertainty, a calibrated DFT sigma, ESPEI's assumed "
        "experimental default (not printed), no sigma at all (no default exists, or model-calculated), or are ZPF "
        "(driving-force tolerance). Sigma distributions are over values with a printed or calibrated sigma, in the "
        "dataset's units (ESPEI units: J/mol-atom, J/mol-atom/K, activity). Rules: README.md, section "
        "\"Uncertainties\".", "",
        table(T, [("folder", 0), ("output", 0), ("values", 1), ("printed", 1), ("dft calibrated", 1),
                  ("assumed ESPEI default", 1), ("none", 1), ("none (model)", 1), ("ZPF", 1), ("sigma unit", 0),
                  ("min", 1), ("q25", 1), ("median", 1), ("q75", 1), ("max", 1)]), ""]) + "\n"


def calibration_md(cal, dp_release):
    L = ["# DFT formation-enthalpy uncertainty: calibration", "",
         "Calibrated from ORCHID's own DFT-vs-experiment pairs (`literature/dft_parity.json`, default subset: "
         "experimental, accepted, signed, not magnitude-suspect; residual = DFT - experiment, kJ/mol-atom). Materials "
         "Project and OQMD are calibrated separately (different energy references). Literature DFT: each accepted "
         "literature DFT HM_FORM value (kJ/mol-atom, binary) paired with the median experimental value of the same "
         "compound (same reduced composition) in that subset.", "",
         "## Method", "",
         "- bias: mean and median residual (**recorded, not applied** to any value);",
         "- robust sigma = 1.4826 x MAD (median absolute deviation about the median); 68 / 95 % quantiles of "
         "|residual|; RMSE;",
         "- 95%% confidence intervals: %d bootstrap resamples of the pairs (fixed seeds), percentile intervals;" % N_BOOT,
         "- strata: match quality (%s) x chemistry (%s); chemistry precedence: p-block / anion, then magnetic 3d, "
         "then all-metal; single-factor strata are listed too;" % (
             "; ".join("`%s` %s" % kv for kv in MATCH_TEXT.items()),
             "; ".join("`%s` %s" % kv for kv in CHEM_TEXT.items())),
         "- experimental-error-deconvolved sigma_DFT = sqrt(sigma_robust^2 - sigma_exp^2), sigma_exp = median printed "
         "experimental uncertainty of the matched records (reported, not used);",
         "- weights use the plain robust sigma (conservative), **leave-system-out**: for each system the sigma of its "
         "stratum is recomputed without that system's pairs; a stratum with fewer than %d pairs falls back to the "
         "source-level sigma (also leave-system-out);" % MIN_STRATUM,
         "- exported MP / OQMD values are on-hull vertices: stratum `on_hull_single` when the vertex is the only entry "
         "at its composition, else `polymorph_ambiguous`; literature DFT values use the chemistry stratum;",
         "- ESPEI weight = 500 J/mol / sigma_DFT (J/mol-atom). 500 J/mol is ESPEI's reference sigma for experimental "
         "HM data; the ratio only makes ESPEI use the calibrated DFT sigma.", "",
         "## Caveats", "",
         "- Heavy tails: RMSE exceeds the robust sigma (column `RMSE/sigma`) and the 95 % quantile of |residual| is "
         "several robust sigmas; a Gaussian likelihood with the robust sigma treats the tail compounds as far less "
         "likely than they are, so a few outliers can pull a fit hard. Consider a Student-t or a wider sigma if your "
         "fit is sensitive to a few compounds.",
         "- Pairs are treated as independent; several experimental records of one compound and correlated errors "
         "within a chemistry are not modelled.",
         "- The calibration is against literature experiments (machine-extracted, Version 1) and inherits their "
         "errors; sigma_deconvolved shows how much of the spread the printed experimental uncertainty explains.",
         "- 0 K DFT energies are used as 298.15 K enthalpies.", ""]
    for src, e in cal["sources"].items():
        L += ["## %s" % src, ""]
        rows = []
        for name, s in [("all", e["all"])] + sorted(e["strata"].items()):
            if not s:
                continue
            rows.append([name, s["n"], s["bias_mean"], s["bias_median"], "%s [%s, %s]" % (
                s["sigma_robust"], *s["ci95"]["sigma_robust"]), "%s [%s, %s]" % (s["q68_abs"], *s["ci95"]["q68_abs"]),
                "%s [%s, %s]" % (s["q95_abs"], *s["ci95"]["q95_abs"]), "%s [%s, %s]" % (s["rmse"], *s["ci95"]["rmse"]),
                s["rmse_over_sigma_robust"], s["sigma_exp_median"], s["sigma_deconvolved"],
                "" if s["n"] >= MIN_STRATUM else "< %d: falls back" % MIN_STRATUM])
        L.append(table(rows, [("stratum", 0), ("n", 1), ("bias mean", 1), ("bias median", 1),
                              ("robust sigma [95% CI]", 0), ("q68 |res|", 0), ("q95 |res|", 0), ("RMSE", 0),
                              ("RMSE/sigma", 1), ("sigma_exp", 1), ("sigma deconv.", 1), ("note", 0)]))
        L.append("")
    L.append("Values in kJ/mol-atom. Leave-system-out sigmas used for the weights: `CALIBRATION.json`, `lso` "
             "(`<source>|<stratum>|<system>`). Pairs file: %s." % dp_release)
    return "\n".join(L) + "\n"


def readme(S, C, cite):
    bf = S["by_folder"]
    pc = S["phase_constraints"]["folders"]
    fx = S["fixes"]
    pcontent = {"assessed": "assessed (handbook) phase-field points, Dataset assessed-phase-assemblages, in the phase "
                            "names of each TDB that has ZPF datasets here (tag `assessed_phase_diagram`)",
                "derived/invariant-brackets": "DERIVED: phase sets just above and below each printed invariant "
                                              "(tag `derived_bracket`)",
                "derived/liquid-stability": "DERIVED: single-phase liquid above the highest liquidus where the "
                                            "sources show complete liquid miscibility (tag `derived_liquid_stability`)"}
    content = {"experimental": "experimental thermochemical tables and phase-composition (ZPF) tables",
               "model": "values the papers calculated with a CALPHAD model, TDB (SHA-256) named where known",
               "invariants": "ZPF datasets from experimental invariant-reaction tables",
               "dft/literature": "literature DFT formation enthalpies, HM_FORM endmember datasets (`solver`)",
               "dft/mp": "Materials Project on-hull formation energies, HM_FORM endmember datasets",
               "dft/oqmd": "OQMD on-hull formation energies, HM_FORM endmember datasets"}
    return "\n".join([
        "# ORCHID ESPEI datasets from the literature", "",
        "ESPEI-format datasets built from the ORCHID literature values (Version 1, machine-extracted) and from DFT "
        "(literature, Materials Project, OQMD). Every literature value is traced to its published ORCHID observation "
        "id (`literature/`, field `id`) and the DOI of its source, passes the same inclusion rules as the published "
        "literature values, and carries its uncertainty in physical units next to the ESPEI weight derived from it.", "",
        table([[f, bf[f]["datasets"], bf[f]["values"], bf[f]["systems"], content[f]] for f in GROUPS],
              [("folder", 0), ("datasets", 1), ("values", 1), ("systems", 1), ("content", 0)]), "",
        "Optional phase constraints (separate folders; never mixed with measured data; see \"Phase constraints\"):", "",
        table([[f, pc[f]["datasets"], pc[f]["points"], pc[f]["systems"], pcontent[f]] for f in PGROUPS if f in pc],
              [("folder", 0), ("datasets", 1), ("ZPF points", 1), ("systems", 1), ("content", 0)]), "",
        "## Use with ESPEI", "",
        "Point ESPEI at one evidence class at a time; classes are never mixed in a folder, and for DFT use **one** of "
        "`dft/mp`, `dft/oqmd`, `dft/literature` per fit (MP uses corrected, OQMD fitted elemental references).", "",
        "```yaml", "# espei-in.yaml", "system:", "  phase_models: my-phases.json",
        "  datasets: orchid-espei/experimental      # or invariants/, dft/mp/, ... (separately, by choice)",
        "generate_parameters:", "  excess_model: linear", "  ref_state: SGTE91", "```", "",
        "```", "espei --input espei-in.yaml", "```", "",
        "or in Python (ESPEI 0.9.1):", "",
        "```python", "from espei.datasets import load_datasets, recursive_glob, apply_tags",
        "ds = load_datasets(recursive_glob('espei/experimental/AG-AL') + recursive_glob('espei/invariants/AG-AL'))",
        "# your choices, e.g. halve the weight of single-reader values, drop model values:",
        "apply_tags(ds, {'single_reader': {'weight': 0.5}, 'model': {'disabled': True}})",
        "```", "",
        "`apply_tags` replaces the dataset's `weight`; to scale the weight derived from the printed sigma, read "
        "`orchid_uncertainty` and set `weight` yourself. `index.csv` / `index.json` list every dataset (system, output, "
        "phases, values, temperature range, evidence class, DOIs, observation ids, reader agreement, sigma summary, "
        "weight).", "",
        "Tags: evidence class (`experimental`, `model`, `dft_literature`, `dft_mp`, `dft_oqmd`), quality "
        "(`single_reader`, `single_reader_text_layer`, `no_printed_uncertainty`, `polymorph_ambiguous`, `above_hull`, "
        "`orchid_caveat`), plus the corpus tags (`thermochemical_table`, `invariant_table`, `via_assessment`, ...).", "",
        "To restrict a fit to evidence certified against a CALPHAD prediction (`calphad_certified`: every value of "
        "the dataset is `espei_certified` in the literature release, i.e. espei_ready and its pycalphad calculation "
        "from the linked TDB agrees within 25 K / 0.05; per-value detail in `orchid_certification`):", "",
        "```yaml", "system:", "  phase_models: my-phases.json", "  datasets: orchid-espei/experimental",
        "  tags:", "    orchid: {disabled: true}             # every ORCHID dataset (literature and DFT) off ...",
        "    calphad_certified: {disabled: false}  # ... except the certified ones (tags apply in this order)",
        "```", "", "(ESPEI applies `tags` before it removes disabled datasets, so the second tag re-enables.)", "",
        "## Phase assignment", "",
        PR.RULE_TEXT, "",
        "Every value's phase comes from its published literature record (printed label -> phase registry -> the "
        "target TDB's phase name), never from a default. Values with different phases in one source table are written "
        "as separate datasets (`<TDB>.<file>-<PHASES>.json`). Formation enthalpies / entropies of a phase that is "
        "stoichiometric in the target TDB are written as ESPEI non-equilibrium datasets (`solver`, the compound's "
        "endmember, output `HM_FORM` / `SM_FORM`, tag `compound_formation_endmember`) when the printed composition "
        "matches the formula and the printed reference state is the TDB's SER; other properties of compounds are not "
        "written. `orchid_phase_assignment` gives the phases, the basis (printed labels or caption) and, flag only, "
        "values outside the phase's field in the calculated diagram of the target TDB (`outside_calculated_field`). "
        "Also excluded (VALIDATION.md, `drop_reasons`): activities outside (0, 1], negative absolute entropies, "
        "increments printed as absolute values, per-gram heat capacities, formation / mixing values outside "
        "physical ranges, values flagged unsigned_in_print / magnitude_suspect / system_mismatch / source_disputed, "
        "reference phases that the TDB does not have, binary three-phase ZPF regions away from every calculated "
        "invariant, ZPF regions with a gas and no printed pressure, and ZPF regions that carry no composition and no "
        "hyperplane (zero driving force in ESPEI 0.9.1).", "",
        "## Uncertainties", "",
        "ESPEI 0.9.1 has no per-value uncertainty field: the sigma it uses is a reference sigma divided by the "
        "dataset's `weight` (500 J/mol for HM, 0.2 J/K/mol for SM and CPM, 500 J/mol on the chemical potential for "
        "activities; ESPEI: \"estimated from NIST TRC uncertainties\", i.e. for experiments). Every value carries its "
        "sigma in physical units in `orchid_uncertainty.values[]` (`sigma`, `sigma_basis`, printed text and coverage "
        "or calibration provenance, `weight_used`). **Weights are derived from these sigmas; you can recompute your "
        "own weights from them.**", "",
        "Which sigma applies to which data:", "",
        table([["experimental, uncertainty printed", "printed (±x taken as 1 sigma: coverage is not printed with the "
                "values; a printed range gives sigma = half-width, bounds kept)", "ESPEI default / median sigma of the "
                "dataset (one weight per equilibrium dataset in ESPEI 0.9.1); per-value equivalents recorded"],
               ["experimental, nothing printed", "`sigma_printed` null; `sigma_assumed` = ESPEI's experimental default "
                "(500 J/mol, 0.2 J/K/mol), basis `espei_default_experimental (assumed, not printed)`; GM has no "
                "default (null)", "1 (tag `no_printed_uncertainty`)"],
               ["DFT (literature, MP, OQMD)", "calibrated against experiment (dft/CALIBRATION.md): robust sigma of "
                "the value's stratum, leave-system-out; bias recorded, not applied", "500 J/mol / sigma_DFT per value "
                "(a conversion, not a claim that DFT is accurate to 500 J/mol)"],
               ["model-calculated (model/)", "none (model-calculated; no measurement uncertainty)", "1 only because "
                "ESPEI needs a value: **choose your own weight** or exclude with the `model` tag"],
               ["ZPF (experimental/ ZPF, invariants/)", "none on the value; printed temperature uncertainty kept where "
                "printed; ESPEI weighs ZPF by a driving-force tolerance", "1 (not set)"]],
              [("data", 0), ("sigma", 0), ("ESPEI weight", 0)]), "",
        "Per folder and output: UNCERTAINTY_SUMMARY.md.", "",
        "## What each dataset carries", "",
        "ESPEI-native content (`components`, `phases`, `conditions`, `output`, `values`, `reference_states` / "
        "`reference_state`, `solver` for DFT endmembers, `tags`, `weight`), plus `reference` (source paper, DOI), "
        "`bibtex` (key in `references.bib`), a `comment` with the ORCHID observation ids, the evidence class and "
        "\"Version 1, machine-extracted\", an `orchid` block (ids, reader agreement, caveats, papers, DOIs, source) and "
        "`orchid_uncertainty`. Model datasets keep `orchid_model_evidence` (TDB name, SHA-256, link basis). DFT datasets "
        "name the TDB file (SHA-256), phase and endmember they map to (see SCHEMA.md).", "",
        "## Caveats", "",
        "- **Machine-extracted (Version 1).** Literature values were read by language-model readers and accepted on "
        "reader agreement plus deterministic checks; a hand-checked precision sample is in progress. Check a value "
        "against its paper (DOI) before relying on it.",
        "- **Single reader.** Values read by one visual reader carry `reader_agreement: single_reader` in "
        "`orchid.observations` and the tag `single_reader`. Reader agreement is taken from %s." % (
            "the literature release" if C["release_has_ra"] else
            "the record store (the published literature release predates the field)"),
        "- **Weights are your choice.** The weights written follow the rules above; tags and `orchid_uncertainty` let "
        "you set your own.",
        "- **DFT.** 0 K energies used as 298.15 K enthalpies; one DFT source per fit; a value is written only where the "
        "compound maps to exactly one TDB phase (phase registry: DFT entry id or paper label) and exactly one endmember "
        "of that TDB model at the compound's composition (one dataset per system, TDB file and phase; `solver` site "
        "ratios and configurations as in the TDB). Excluded compounds are listed in validation.json with reasons. ESPEI "
        "0.9.1 evaluates HM_FORM against each element's SER phase named on the TDB's ELEMENT line, which must exist as a "
        "phase in the TDB; in some binary TDB files it does not (e.g. GRAPHITE, ORTHORHOMBIC_GA, 1/2_MOLE_O2(G)), and "
        "ESPEI then stops with a KeyError. Such datasets are not written (likelihood evaluation gate, VALIDATION.md); "
        "to use those values, add the element phase to your working TDB and rebuild from the literature release.",
        "- **Energy basis.** ESPEI reads energies per mole of atoms. Literature values whose per-atom basis was not "
        "printed carry the caveat `energy_basis_unstated`.",
        "- **ZPF encoding.** Phase labels are mapped to the model phases of the assessment the source folder is named "
        "after (tag `via_assessment`); a printed compound composition that is its formula rounded is written at the "
        "exact formula value; a `__HYPERPLANE__` vertex is added only where every known composition comes from "
        "stoichiometry. Range temperatures are written at the midpoint (`orchid_value_bases`).",
        "- **Line compounds (E2).** A ZPF vertex of a phase that is stoichiometric in the target TDB is written with "
        "composition null: ESPEI 0.9.1 / pycalphad 0.11.2 return NaN for a line compound at a fixed composition (one "
        "NaN makes the ZPF likelihood -inf), and with null ESPEI evaluates the compound at its formula composition. "
        "Where that leaves a region with no known composition, a `__HYPERPLANE__` is added at the mean of the compound "
        "compositions. Known limitation of ESPEI 0.9.1: a single-phase vertex at a line-compound composition gives an "
        "undefined target hyperplane (zero driving force) once the compound is stable there, so a congruent melting "
        "point constrains the model in one direction only; the invariant brackets below constrain both.",
        "- **Duplicates (E3, E4).** Identical ZPF regions in several datasets (several assessment folders, or several "
        "papers quoting the same invariant) are written once, with every source in `orchid.observations[i]."
        "also_reported_as`; a printed row read twice (same row id, conditions and value) counts once (`also_read_as`). "
        "Records concerned are listed in UPSTREAM_CORRECTIONS.md.",
        "- **Activity reference states (E5).** Every activity table was checked against the paper: the reference state "
        "written is the one the paper prints (or, where nothing is printed, the pure liquid of an element that is "
        "liquid at the temperature, recorded as `not_printed_implied`); `orchid_reference_state` gives the basis and "
        "the printed text. Literature records whose reference state differs are listed in UPSTREAM_CORRECTIONS.md; "
        "the published literature values are not changed here.",
        "- Values removed by the inclusion gate and datasets ESPEI rejects are listed in VALIDATION.md.", "",
        "## Optional minimum-sigma policy (E6)", "",
        "Printed uncertainties are kept as printed. Some are much smaller than the scatter between independent "
        "measurements of the same quantity (e.g. formation enthalpies printed with ±80 J/mol-atom while other "
        "measurements of the same compound differ by 1-2 kJ/mol-atom); with the weight derived from such a sigma a few "
        "values can dominate a fit. Every experimental thermochemical dataset carries `orchid_weight_floor` "
        "(`suggested_weight` = min(weight, 1), i.e. sigma_eff = max(sigma, ESPEI's experimental default)). The floor "
        "changes the weight of %d of %d datasets. It is a suggestion, not applied: "
        "`d['weight'] = d['orchid_weight_floor']['suggested_weight']` before `load_datasets` applies it." % (
            fx["E6_datasets_where_floor_changes_weight"], fx["E6_datasets_with_floor_suggestion"]), "",
        "## Phase constraints (optional)", "",
        "ZPF data alone can let a refinement misfire: ESPEI 0.9.1 penalises a phase that should be present but is not, "
        "but gives no driving force for a liquid point at a line-compound composition once the compound is stable, so "
        "a compound can melt far above its printed melting point without cost. Three optional, tagged dataset groups "
        "add phase evidence. They share the phase names of the assessment folder (TDB) they are written for, so they "
        "load together with `experimental/` and `invariants/` of that TDB. Use them by choice, never as measured data.",
        "",
        "- `derived/invariant-brackets/<SYS>/<TDB>.<file>.json` (tags `derived`, `derived_bracket`): " + ec.BRACKET_RULE
        + " Each point records its source observation, T, u and offset (`orchid.points`).",
        "- `derived/liquid-stability/<SYS>/<TDB>.liquid-stability.json` (tags `derived`, `derived_liquid_stability`): "
        + ec.LIQUID_RULE + " `orchid.basis` gives the temperatures used.",
        "- `assessed/<SYS>/<TDB>.<DIAGRAM>.json` (tag `assessed_phase_diagram`): the ZPF datasets of "
        "assessed-phase-assemblages (one point inside each named field of an expert-drawn handbook diagram) with the "
        "phase names of that TDB (phase registry); points with a phase that the TDB does not model, or models twice, "
        "are left out; `orchid.confidence` per point. Evidence class assessed_phase_diagram: neither experimental nor "
        "CALPHAD.", "",
        "Counts (brackets written and skipped by reason, systems with liquid-stability points) are in validation.json "
        "(`phase_constraints`). Suggested use: weight each group so that it contributes comparably to the measured "
        "data (e.g. per-type weights proportional to 1/N), and check the fitted model's invariant temperatures "
        "against the printed ones after the fit.", "",
        "## Citation", "", cite["text"], "", "Part of: " + cite["part_of"]["text"], "",
        "Also cite " + "; ".join(cite["also_cite"]) + ".", "",
        "Licence: CC BY 4.0 (derived data; MP and OQMD data CC BY 4.0). The source papers are not redistributed; only "
        "their citations and DOIs.", ""])


SCHEMA = """# ORCHID ESPEI bundle: schema

## Layout

- `experimental/<SYS>/<source>.<file>.json`, `model/<SYS>/...`, `invariants/<SYS>/...`: corpus datasets.
  `<SYS>`: elements sorted, uppercase, joined by `-` (as in `literature/`); `<source>`: the corpus folder the dataset
  was written in (the assessment / TDB it was mapped through); `<file>`: the corpus file name.
- `dft/<literature|mp|oqmd>/<SYS>/<TDB stem>.<PHASE>.json`: DFT formation enthalpies of endmembers of one phase of one
  TDB file (`solver` datasets, output `HM_FORM`, T = 298.15 K, P = 101325 Pa).
- `dft/CALIBRATION.json|md`: DFT sigma calibration.
- Optional phase constraints (ZPF; tagged; never mixed with measured data): `derived/invariant-brackets/<SYS>/<TDB>.<file>.json`
  (one per source invariant dataset, same phase names), `derived/liquid-stability/<SYS>/<TDB>.liquid-stability.json`,
  `assessed/<SYS>/<TDB>.<DIAGRAM>.json` (assessed-phase-assemblages in that TDB's phase names).
- `UPSTREAM_CORRECTIONS.md`: literature records the export had to differ from (reference states, rows read twice).

## Dataset files (ESPEI 0.9.1 format)

ESPEI-native keys: `components`, `phases`, `conditions` (`P` Pa, `T` K, `X_<EL>` mole fractions), `output`, `values`
(J/mol-atom, J/mol-atom/K, activities; ZPF: `[phase, [components], [mole fractions]]` per vertex, `null` where not
printed), `reference_states` (or `reference_state` for activities), `solver` (DFT: `mode`, `sublattice_site_ratios`
and `sublattice_configurations` exactly as in the TDB phase model), `tags`, `weight` (scalar for equilibrium and
activity datasets; one per value for DFT `solver` datasets; absent = 1). Values removed by the inclusion gate are cut
from `values` and the matching `conditions` entries.

| key | meaning |
|---|---|
| `reference` | source paper(s) of the values (where measured / calculated): `Author Year, doi:...`; DFT databases by name and release |
| `bibtex` | BibTeX key(s) in `references.bib` (space-separated when several) |
| `comment` | literature: `ORCHID literature values, Version 1, machine-extracted. Evidence class: <class>. ORCHID observation ids (...): ORC-...`; DFT: `0 K DFT used as 298.15 K enthalpy. ...` with TDB, SHA-256, phase, ids and the eV conversion |
| `orchid.evidence_class` | `experimental`, `model` or `dft` |
| `orchid.observations[]` | per value, in order: `id` (ORCHID literature id), `reader_agreement` (`multi_reader`, `single_reader`, `single_reader_text_layer`; absent when not classified), `paper`, `measured_in`, `caveats` |
| `orchid.papers`, `orchid.dois` | paper ids (`literature/papers.json`) and DOIs involved: source, containing paper, printed cited works |
| `orchid.source_dataset`, `orchid.values_dropped`, `orchid.espei_evaluation` | corpus file; values removed by the gate; invariant ZPF driving forces against the assessment |
| `orchid.dft_source`, `orchid.tdb`, `orchid.registry_phase`, `orchid.mapping[]` | DFT: source; TDB file (name, path, SHA-256, phase); phase-registry id; per value: id (MP / OQMD entry or ORCHID id), compound, mapping basis, space group, match quality |
| `orchid_value_bases` | ZPF from ranges / central values: how a temperature was derived, keyed by ORCHID id |
| `orchid_model_evidence` | model datasets: `tdb` (name, path, SHA-256, catalog entry), `link_basis`, `model_citation`, `model_reproduction` |
| `orchid_uncertainty` | see below |
| `orchid_phase_assignment` | literature property datasets: `phases`, `basis` (counts: `printed_labels`, `caption_names_liquid`), `rule`, `outside_calculated_field[]` {id, phase, x, T_K} (flag only), `endmember` {tdb, path, sha256, phase, reference} for compound formation datasets |
| `orchid_certification` | literature datasets: `basis` (`pycalphad_invariant_agreement`), `all_values_certified`, `n_certified`, `values[]` {id, espei_certified, calphad_check}; tag `calphad_certified` when all values are certified |
| `orchid_reference_state` | activities (E5): `element`, `phase` written, `basis` (`printed` / `not_printed_implied`), `quote`, `where`, `record_reference_state_agrees`, `written_before_audit` |
| `orchid_weight_floor` | experimental thermochemical data (E6, optional): `suggested_weight` = min(weight, 1), `weight_written`, `changes_weight`, `sigma_floor` (+ unit), `rule` |
| `orchid_zpf_encoding` | ZPF (E2): number of line-compound vertices written with composition null, hyperplanes added, rule |
| `orchid.observations[i].also_reported_as` | E3: the same ZPF region from other datasets (`id`, `paper`, `measured_in`, `dataset`) |
| `orchid.observations[i].also_read_as` | E4: ids of the same printed row read again (not written as separate values) |
| `orchid.points` (derived brackets) | per point: `source_observation`, `side` (above/below), `T_invariant_K`, `offset_K`, `T_sigma_K`, `T_sigma_basis`, `x0`, `hyperplane_delta`, `reaction_printed`, `reaction_type`, `phases_above`, `phases_below` |
| `orchid.basis` (liquid stability) | `T_liq_max_K` and the temperatures it was taken from |
| `orchid.record_ids`, `orchid.confidence` (assessed) | assessed-phase-assemblages record id and confidence per point; `target_tdb`, `points_left_out` |

## Uncertainties (`orchid_uncertainty`)

`{unit, default_sigma {value, name}, weight_rule, values[]}`, one entry per value in value order:

| field | meaning |
|---|---|
| `sigma` | the value's standard uncertainty in the dataset's units (`unit`), or `null` when none is available (never invented) |
| `sigma_basis` | `printed`; `dft_calibrated`; `espei_default_experimental (assumed, not printed)`; `none (...)` with the reason (model-calculated; no ESPEI default); `zpf (...)` |
| `printed` | `{value, unit, interval, coverage, lower?, upper?, note?}`: printed uncertainty as published (`coverage`: `1sigma_assumed` for a bare ±x, `range_half_width` for a printed range) |
| `sigma_printed`, `sigma_assumed` | for experimental values without a printed uncertainty: `null` and ESPEI's experimental default (or `null` where ESPEI has none) |
| `sigma_mu_J_per_mol` | activities: R T sigma_a / a, the chemical-potential sigma ESPEI's activity residual uses |
| `calibration` | DFT: `{source, stratum, used (stratum or source_level), leave_system_out, system_left_out, n, sigma_robust_kJ, bias_median_kJ_not_applied, ci95_sigma_robust_kJ_full_stratum, n_full_stratum}` |
| `temperature_sigma_K`, `temperature_sigma_basis` | ZPF: printed temperature uncertainty |
| `weight_used` | the ESPEI weight that applies to this value; `weight_per_value_equivalent` = default sigma / this value's sigma |
| `printed_uncertainty_not_used` | model-calculated values: a printed uncertainty kept as printed, not used |

Caveat codes (published flags of the value): `magnitude_suspect`, `system_mismatch`, `source_disputed`,
`dft_outlier_unexplained`, `basis_inferred_from_magnitude`, `basis_corrected_from_formula`, `sign_from_heading`,
`unsigned_in_print`, `judge_adjudicated`, `page_unreadable`, `plausibility_fail`, `energy_basis_unstated`.

## index.json / index.csv

One row per dataset: `path`, `folder`, `system`, `components`, `output`, `phases`, `n_values`, `n_values_dropped`,
`T_min_K`, `T_max_K`, `evidence_class`, `reference`, `bibtex`, `dois`, `papers`, `observation_ids`, `reader_agreement`
(counts by code; `not_classified` where the record builder did not classify), `reader_agreement_source`, `caveats`,
`tags`, `weight` (`per_value` for DFT), `sigma_unit`, `sigma_min`, `sigma_median`, `sigma_max` (over printed or
calibrated sigmas), `sigma_basis_counts`, `dft_calibration_stratum`, `tdb`, `tdb_sha256`, `model_link_basis`,
`espei_median_abs_driving_force_J`, `source_dataset`. CSV: lists joined with `;`, objects as compact JSON.

## validation.json

`espei_version`; `folders` (files and datasets loaded per folder); `invalid_final` (empty); `excluded` (rejected by
ESPEI, with reason code); `dropped_datasets`; `dropped_values`; `dft_excluded` (compound or compound x TDB model, with
reason); `fixes` (E1-E6 counts, E3 merges, E5 audit per activity dataset); `upstream_corrections`; `phase_constraints`;
`summary`. index.json / index.csv also list the phase-constraint datasets (`evidence_class` `derived` or
`assessed_phase_diagram`); `totals.datasets` / `totals.values` count the literature and DFT datasets only.
"""


def fixes_summary(C, P):
    e5 = C["e5"]
    return {
        "E1_va_added_datasets": P["counts"].get("E1_va_added", 0),
        "E1_basis": {k[9:]: v for k, v in sorted(P["counts"].items()) if k.startswith("E1_basis_")},
        "E2_vertices_nulled": P["counts"].get("E2_vertices_nulled", 0),
        "E2_datasets": P["counts"].get("E2_datasets", 0),
        "E2_hyperplanes_added": P["counts"].get("E2_hyperplanes_added", 0),
        "E3_regions_removed": len(P["e3"]),
        "E3_datasets_emptied_and_dropped": P["counts"].get("E3_datasets_dropped", 0),
        "E4_values_removed": len(C["upstream"]["duplicate_reads"]),
        "E5_activity_datasets_audited": len(e5),
        "E5_reference_state_changed": sum(1 for x in e5 if x["changed"]),
        "E5_basis": dict(Counter(x["basis"] for x in e5)),
        "E5_records_disagreeing_with_paper": len(C["upstream"]["reference_state"]),
        "E6_datasets_with_floor_suggestion": sum(C["e6"].values()),
        "E6_datasets_where_floor_changes_weight": C["e6"].get(True, 0),
    }


def phase_constraint_summary(P):
    rows = P["rows"].values()
    by = defaultdict(lambda: {"datasets": 0, "points": 0, "systems": set()})
    for r in rows:
        b = by[r["folder"]]
        b["datasets"] += 1
        b["points"] += r["n_values"]
        b["systems"].add(r["system"])
    return {"folders": {k: {"datasets": v["datasets"], "points": v["points"], "systems": len(v["systems"])}
                        for k, v in sorted(by.items())},
            "brackets": P["brackets"], "assessed_points_left_out": P["assessed_skipped"],
            "liquid_stability_systems": P["liquid"]}


def upstream_md(C):
    L = ["# ORCHID ESPEI bundle: literature records to correct upstream", "",
         "The bundle never changes a published literature value. Where the export had to differ from a published "
         "literature record, the record is listed here for correction in the literature release.", "",
         "## Reference states that differ from what the paper states (E5)", ""]
    rs = C["upstream"]["reference_state"]
    if rs:
        L.append(table([[x["id"], x["table"], x["element"], x["record_reference_phase"], x["stated"], x["basis"]]
                        for x in rs], [("ORCHID id", 0), ("table", 0), ("element", 0), ("record reference phase", 0),
                                       ("stated in the paper", 0), ("basis", 0)]))
    else:
        L.append("None.")
    L += ["", "## The same printed row read twice (E4)", "",
          "Two published records with the same row id (`pid`), conditions and value; the bundle counts the value once.",
          ""]
    du = C["upstream"]["duplicate_reads"]
    L.append(table([[x["kept_id"], x["duplicate_id"], x["pid"]] for x in du],
                   [("kept id", 0), ("duplicate id", 0), ("row (pid)", 0)]) if du else "None.")
    L += ["", "## Corpus ESPEI export: phase of thermochemical datasets (root cause of LIQUID mislabels)", "",
          "The corpus record builder (`db-v2.1/bin-general/general_tables.py`, `espei_datasets`) writes "
          "`phases = sorted(mapped model phases | {\"LIQUID\"})`, or `[\"LIQUID\"]` when no printed label maps to a "
          "model phase, for every non-ZPF dataset: compound formation enthalpies at compound compositions and "
          "temperatures where no liquid exists were therefore labelled LIQUID, and LIQUID was added beside every "
          "mapped phase. The published literature records were correct (`property` HM_FORM / HM_MIX, printed phase "
          "labels). The bundle no longer takes the phase from the corpus dataset: each value's phase is derived from "
          "its published record (README.md, 'Phase assignment'). To fix upstream: drop the LIQUID default and the "
          "LIQUID union in `espei_datasets`.", "",
          "## Reference-state phase names", "",
          "The record builder writes the SER reference phase of O, N, H and F as `1` (truncated `1/2_MOLE_O2(G)`) "
          "and some metal reference phases as Strukturbericht symbols (`A1`, `A3`); the literature release completes "
          "the names (QUALITY_FIXES.md, section 5b) and the bundle drops values whose reference phase the target TDB "
          "does not define (`reference_phase_not_in_target_tdb`).", "",
          "## Compound values whose printed composition differs from the formula", "",
          "Composition copied from another table row, or a per-formula label on a solution composition; excluded "
          "from the bundle (`composition_not_at_compound`), to correct in the literature release:", ""]
    cc = sorted({(x["orchid_id"], x["source"]) for x in C["drops"] if x["reason"] == "composition_not_at_compound"})
    L.append(table([[a, b] for a, b in cc], [("ORCHID id", 0), ("corpus dataset", 0)]) if cc else "None.")
    return "\n".join(L) + "\n"


def build_files(C, D, val, data, release, dp_release, P):
    oc = cite_module()
    cite = oc.cite_as(NAME, data)
    S = summarize(C, D, val)
    S["fixes"] = fixes_summary(C, P)
    S["phase_constraints"] = phase_constraint_summary(P)
    files = {}
    alld = dict(C["datasets"])
    alld.update(D["datasets"])
    alld.update(P["datasets"])
    for rel in sorted(alld):
        files[rel] = dumps(alld[rel])
    allr = dict(C["rows"])
    allr.update(D["rows"])
    allr.update(P["rows"])
    rows = [allr[k] for k in sorted(allr)]
    files["index.json"] = dumps({"release": "orchid-espei %s (literature gate: %s)" % (oc.RELEASE_TAG, release),
                                 "cite_as": cite, "totals": S, "datasets": rows})
    cols = list(rows[0].keys())
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(cols)
    for r in rows:
        w.writerow([";".join(map(str, v)) if isinstance(v, list) else
                    (json.dumps(v, sort_keys=True) if isinstance(v, dict) else ("" if v is None else v))
                    for v in (r[c] for c in cols)])
    files["index.csv"] = buf.getvalue().encode("utf-8")
    papers = dict(C["papers"])
    papers.update(D["papers"])
    bib = [bib_entry(bib_key(pk, papers[pk]), pk, papers[pk]) for pk in sorted(papers)]
    bib += ["@article{jain2013_mp,\n  author   = {Jain, A. and others},\n  title    = {{Commentary: The Materials "
            "Project: A materials genome approach to accelerating materials innovation}},\n  journal  = {APL Materials},"
            "\n  volume   = {1},\n  pages    = {011002},\n  year     = {2013},\n  doi      = {10.1063/1.4812323}\n}",
            "@article{saal2013_oqmd,\n  author   = {Saal, J. E. and Kirklin, S. and Aykol, M. and Meredig, B. and "
            "Wolverton, C.},\n  title    = {{Materials design and discovery with high-throughput density functional "
            "theory: the Open Quantum Materials Database (OQMD)}},\n  journal  = {JOM},\n  volume   = {65},\n  pages    "
            "= {1501--1509},\n  year     = {2013},\n  doi      = {10.1007/s11837-013-0755-4}\n}",
            "@article{kirklin2015_oqmd,\n  author   = {Kirklin, S. and others},\n  title    = {{The Open Quantum "
            "Materials Database (OQMD): assessing the accuracy of DFT formation energies}},\n  journal  = {npj Comput. "
            "Mater.},\n  volume   = {1},\n  pages    = {15010},\n  year     = {2015},\n  doi      = "
            "{10.1038/npjcompumats.2015.10}\n}"]
    files["references.bib"] = ("% ORCHID ESPEI bundle: papers cited by the datasets (keys = the datasets' `bibtex`)\n\n"
                                + "\n\n".join(bib) + "\n").encode("utf-8")
    files["validation.json"] = dumps({"espei_version": val["espei_version"], "folders": val["folders"],
                                      "invalid_final": val["invalid_final"], "excluded": val["excluded"],
                                      "evaluation_excluded": val.get("evaluation_excluded"),
                                      "evaluation_regions_cut": val.get("evaluation_regions_cut"),
                                      "likelihood_evaluation": val.get("likelihood_evaluation"),
                                      "exclusion_reasons": ESPEI_EXCLUSION_TEXT,
                                      "dropped_datasets": C["dropped_datasets"], "dropped_values": C["drops"],
                                      "drop_reasons": DROP_REASONS,
                                      "dft_excluded": sorted(D["excluded"], key=lambda e: json.dumps(e, sort_keys=True)),
                                      "fixes": {"summary": S["fixes"], "e3_merged": P["e3"], "e5_activity": C["e5"]},
                                      "upstream_corrections": C["upstream"],
                                      "phase_constraints": S["phase_constraints"],
                                      "summary": S})
    files["UPSTREAM_CORRECTIONS.md"] = upstream_md(C).encode("utf-8")
    files["VALIDATION.md"] = validation_md(S, val, C, D, release).encode("utf-8")
    files["UNCERTAINTY_SUMMARY.md"] = uncertainty_md(C, D).encode("utf-8")
    files["dft/CALIBRATION.json"] = dumps({"method": {
        "pairs": "literature/dft_parity.json default subset; literature DFT vs median experiment of the same compound",
        "residual": "DFT - experiment, kJ/mol-atom", "sigma_robust": "1.4826 * MAD", "bootstrap": N_BOOT,
        "ci": "percentile 95 %", "min_stratum": MIN_STRATUM, "leave_system_out": True, "bias_applied": False,
        "weight": "500 J/mol / (sigma_robust_kJ * 1000)", "eV_to_J_per_mol": EV_TO_J_PER_MOL,
        "chemistry_classes": CHEM_TEXT, "match_quality": MATCH_TEXT}, "pairs_release": dp_release, **D["calibration"]})
    files["dft/CALIBRATION.md"] = calibration_md(D["calibration"], dp_release).encode("utf-8")
    files["README.md"] = readme(S, C, cite).encode("utf-8")
    files["SCHEMA.md"] = SCHEMA.encode("utf-8")
    files["CITATION.cff"] = oc.cff(NAME, data).replace(
        "keywords: [CALPHAD, thermodynamics, phase diagrams, convex hull, DFT, formation enthalpy, AI-ready data]",
        "keywords: [CALPHAD, ESPEI, thermodynamics, phase diagrams, ZPF, DFT, formation enthalpy, uncertainty, "
        "AI-ready data]").replace("original sources of the data.",
                                  "original sources of the data (references.bib).").encode("utf-8")
    return files, S


def prow(rel, d, cls):
    Ts = d["conditions"]["T"] if isinstance(d["conditions"]["T"], list) else [d["conditions"]["T"]]
    o = d.get("orchid") or {}
    return {"path": rel, "folder": "/".join(rel.split("/")[:-2]), "system": rel.split("/")[-2],
            "components": d["components"], "output": d["output"], "phases": d["phases"], "n_values": len(d["values"]),
            "n_values_dropped": 0, "T_min_K": min(Ts), "T_max_K": max(Ts), "evidence_class": cls,
            "reference": d.get("reference"), "bibtex": d.get("bibtex"), "dois": o.get("dois") or [],
            "papers": o.get("papers") or [], "observation_ids": [], "reader_agreement": {},
            "reader_agreement_source": None, "caveats": [], "tags": d["tags"], "weight": d.get("weight"),
            "sigma_unit": None, "sigma_min": None, "sigma_median": None, "sigma_max": None, "sigma_basis_counts": {},
            "dft_calibration_stratum": None, "tdb": o.get("target_tdb"), "tdb_sha256": None, "model_link_basis": None,
            "espei_median_abs_driving_force_J": None, "source_dataset": o.get("derived_from") or o.get("source")}


def refresh_row(C, rel):
    d, r = C["datasets"][rel], C["rows"][rel]
    Ts = d["conditions"]["T"] if isinstance(d["conditions"]["T"], list) else [d["conditions"]["T"]]
    r.update({"components": d["components"], "phases": d["phases"], "T_min_K": min(Ts), "T_max_K": max(Ts),
              "n_values": len(d["values"]) if d["output"] == "ZPF" else len(flat_values(d)),
              "n_values_dropped": d["orchid"].get("values_dropped") or 0, "dois": d["orchid"].get("dois") or [],
              "papers": d["orchid"].get("papers") or [],
              "observation_ids": [o["id"] for o in d["orchid"]["observations"]],
              "espei_median_abs_driving_force_J": (d["orchid"].get("espei_evaluation") or {}).get(
                  "median_abs_driving_force_J")})


def dedup_properties(C):
    """E3 for property datasets: a dataset identical to another of the same folder and system (same output, phases,
    conditions, values, solver and reference states; e.g. two assessment folders quoting the same table) is written
    once, the other's ids in `orchid.also_reported_as`. Returns the number dropped."""
    seen, removed = {}, []
    for rel in sorted(C["datasets"]):
        d = C["datasets"][rel]
        if d["output"] == "ZPF":
            continue
        sig = json.dumps([rel.split("/")[0], rel.split("/")[1], d["output"], d["phases"], d["conditions"],
                          d["values"], d.get("solver"), d.get("reference_states"), d.get("reference_state")],
                         sort_keys=True)
        if sig in seen:
            removed.append((rel, seen[sig]))
        else:
            seen[sig] = rel
    for rel, kept in removed:
        d = C["datasets"].pop(rel)
        r = C["rows"].pop(rel)
        C["datasets"][kept]["orchid"].setdefault("also_reported_as", []).append(
            {"dataset": rel, "ids": [o["id"] for o in d["orchid"]["observations"]]})
        C["dropped_datasets"].append({"source": r["source_dataset"].replace("db-v2.1/", "", 1),
                                      "reason": ["duplicate_identical_property_dataset"], "values": r["n_values"]})
    return len(removed)


def phase_constraints(C, D, data, pub):
    """E3 (dedup), the derived/assessed phase-constraint datasets, then E1/E2 on every dataset."""
    R = ec.Registry(data)
    counts = Counter()
    merged, changed = ec.dedup_zpf(C)
    counts["E3_property_datasets_dropped"] = dedup_properties(C)
    counts["E3_regions_removed"] = len(merged)
    counts["E3_datasets_dropped"] = sum(1 for x in C["dropped_datasets"] if x["reason"] == ["duplicate_identical_zpf"])
    brackets, b_skipped, b_made = ec.build_brackets(C, R)        # before E2: printed compositions
    for rel, d in list(C["datasets"].items()) + list(D["datasets"].items()):
        sysn = rel.split("/")[-2]
        stem = ec.stem_of(rel) if not rel.startswith("dft/") else os.path.basename(rel).rsplit(".", 2)[0]
        ec.apply_e1(d, sysn, stem, R, counts)
        if not rel.startswith("dft/") and ec.apply_e2(d, stem, R, counts):
            d["orchid"].pop("espei_evaluation", None)
            changed.add(rel)
    for rel in sorted(changed):
        if rel in C["rows"]:
            refresh_row(C, rel)
    for rel in C["rows"]:
        C["rows"][rel]["components"] = C["datasets"][rel]["components"]
    for rel in D["rows"]:
        D["rows"][rel]["components"] = D["datasets"][rel]["components"]
    for rel, d in brackets.items():
        ec.apply_e1(d, rel.split("/")[-2], ec.stem_of(rel), R, Counter())
    stems = defaultdict(set)
    for rel, d in C["datasets"].items():
        if rel.split("/")[0] in ("experimental", "invariants") and d["output"] == "ZPF":
            stems[rel.split("/")[1]].add(ec.stem_of(rel))
    assessed, a_skipped = ec.assessed_per_tdb(data, R, stems)
    liquid, l_why = ec.liquid_stability(data, R, pub, stems)
    P = {"datasets": {}, "rows": {}, "counts": counts, "e3": merged,
         "brackets": {"made": dict(b_made), "skipped": dict(sorted(b_skipped.items()))},
         "assessed_skipped": dict(sorted(a_skipped.items())), "liquid": dict(sorted(l_why.items()))}
    for group, cls in ((brackets, "derived"), (assessed, "assessed_phase_diagram"), (liquid, "derived")):
        for rel, d in group.items():
            P["datasets"][rel] = d
            P["rows"][rel] = prow(rel, d, cls)
    return P


# ------------------------------------------------------------------------------------------- likelihood gate --
EVAL_CACHE = os.path.join(SITE, "build-cache", "espei-eval")
EVAL_RULE = ("Every dataset is evaluated against its target TDB with the ESPEI 0.9.1 likelihood functions and the "
             "TDB's own parameters (tools/espei_evaluate.py: ZPF driving forces, equilibrium-property residuals, "
             "activity residuals and likelihood, non-equilibrium likelihood). A ZPF region with a non-finite driving "
             "force is cut (`espei_evaluation_nonfinite_region`); a dataset that raises an error (e.g. a KeyError for "
             "an SER phase such as GRAPHITE or 1/2_MOLE_O2(G) that the TDB does not define), gives a non-finite "
             "value, times out, or whose TDB pycalphad cannot load is not written (`espei_evaluation_<status>`).")


def _tdb_paths(data):
    reg = ec.Registry(data)
    idx = {}
    for stem, ms in reg.models.items():
        for m in ms.values():
            if m.get("path"):
                idx[stem] = os.path.join(lit.DEFAULT_OPEN_TDB, m["path"])
    for p in sorted(glob.glob(os.path.join(lit.DEFAULT_OPEN_TDB, "databases", "**", "*.[tT][dD][bB]"), recursive=True)):
        idx.setdefault(os.path.basename(p)[:-4], p)
    return idx


def _target_stem(rel, d):
    o = d.get("orchid") or {}
    if rel.startswith("dft/"):
        t = (o.get("tdb") or {}).get("name") or ""
    elif rel.startswith(("assessed/", "derived/")) and o.get("target_tdb"):
        t = o["target_tdb"]
    else:
        return ec.stem_of(rel)
    return t[:-4] if t.lower().endswith(".tdb") else t


def run_evaluation(datasets, data, espei_py):
    """{rel: result} for the given datasets; results cached by sha256(dataset bytes, TDB bytes)."""
    paths = _tdb_paths(data)
    os.makedirs(EVAL_CACHE, exist_ok=True)
    res, todo, keys = {}, [], {}
    tdb_sha = {}
    tmp = tempfile.mkdtemp(prefix="orchid-espei-eval-")
    try:
        for rel in sorted(datasets):
            b = dumps(datasets[rel])
            p = paths.get(_target_stem(rel, datasets[rel]))
            if not p or not os.path.exists(p):
                res[rel] = {"status": "tdb_not_found"}
                continue
            if p not in tdb_sha:
                tdb_sha[p] = hashlib.sha256(open(p, "rb").read()).hexdigest()
            k = hashlib.sha256(b + b"\0" + tdb_sha[p].encode() + b"\0espei-eval-v3").hexdigest()
            cp = os.path.join(EVAL_CACHE, k + ".json")
            if os.path.exists(cp):
                res[rel] = json.load(open(cp, encoding="utf-8"))
                continue
            f = os.path.join(tmp, "%06d.json" % len(todo))
            with open(f, "wb") as fh:
                fh.write(b)
            todo.append({"rel": rel, "file": f, "tdb": p})
            keys[rel] = cp
        if todo:
            man = os.path.join(tmp, "manifest.json")
            json.dump(todo, open(man, "w"))
            outp = os.path.join(tmp, "out.json")
            subprocess.run([espei_py, os.path.join(os.path.dirname(os.path.abspath(__file__)), "espei_evaluate.py"),
                            man, outp, str(EVAL_WORKERS)], check=True)
            new = json.load(open(outp, encoding="utf-8"))
            for rel, r in new.items():
                res[rel] = r
                with open(keys[rel] + ".tmp", "w", encoding="utf-8") as fh:
                    json.dump(r, fh, sort_keys=True)
                os.replace(keys[rel] + ".tmp", keys[rel])
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return res


def cut_regions(C, rel, idx, reason):
    """Remove ZPF regions `idx` of a literature dataset (values, T, observations, uncertainty, value bases, recs)."""
    d = C["datasets"][rel]
    keep = [i for i in range(len(d["values"])) if i not in idx]
    Ts = d["conditions"]["T"] if isinstance(d["conditions"]["T"], list) else [d["conditions"]["T"]] * len(d["values"])
    for i in sorted(idx):
        C["drops"].append({"source": C["rows"][rel]["source_dataset"].replace("db-v2.1/", ""), "index": i,
                           "reason": reason, "orchid_id": d["orchid"]["observations"][i].get("id")})
    d["values"] = [d["values"][i] for i in keep]
    d["conditions"]["T"] = [Ts[i] for i in keep]
    d["orchid"]["observations"] = [d["orchid"]["observations"][i] for i in keep]
    d["orchid_uncertainty"]["values"] = [d["orchid_uncertainty"]["values"][i] for i in keep]
    if "orchid_value_bases" in d:
        ids = {o["id"] for o in d["orchid"]["observations"]}
        d["orchid_value_bases"] = [b for b in d["orchid_value_bases"] if b.get("observation") in ids]
    d["phases"] = sorted({vx[0] for reg in d["values"] for vx in reg if vx[0] != "__HYPERPLANE__"})
    d["orchid"]["values_dropped"] = (d["orchid"].get("values_dropped") or 0) + len(idx)
    d["orchid"].pop("espei_evaluation", None)
    ids = [o["id"] for o in d["orchid"]["observations"]]
    d["comment"] = re.sub(r"(in value order\): ).*$", lambda mm: mm.group(1) + " ".join(ids), d["comment"])
    C["_recs"][rel] = [C["_recs"][rel][i] for i in keep]
    refresh_row(C, rel)


def evaluation_gate(C, D, P, data, espei_py):
    excluded, cut = [], []
    for it in range(8):
        alld = dict(C["datasets"])
        alld.update(D["datasets"])
        alld.update(P["datasets"])
        ev = run_evaluation(alld, data, espei_py)
        print("evaluation pass %d: %s" % (it + 1, dict(Counter(r["status"] for r in ev.values()))), file=sys.stderr)
        changed = False
        for rel in sorted(ev):
            r = ev[rel]
            if r["status"] == "ok":
                continue
            d = alld[rel]
            bad = r.get("bad_regions") or []
            if r["status"] == "nonfinite" and rel in C["datasets"] and d["output"] == "ZPF" and bad and \
                    len(bad) < len(d["values"]):
                cut_regions(C, rel, set(bad), "espei_evaluation_nonfinite_region")
                cut.append({"path": rel, "regions": len(bad)})
                changed = True
                continue
            reason = "espei_evaluation_" + r["status"]
            src = None
            for G in (C, D, P):
                if rel in G["datasets"]:
                    G["datasets"].pop(rel)
                    row_ = G["rows"].pop(rel, None) or {}
                    src = row_.get("source_dataset")
                    if G is C:
                        C["dropped_datasets"].append({"source": (src or rel).replace("db-v2.1/", "", 1),
                                                      "reason": [reason], "values": row_.get("n_values")})
            for prel in [x for x, pd_ in P["datasets"].items() if (pd_.get("orchid") or {}).get("derived_from") == rel]:
                P["datasets"].pop(prel)
                P["rows"].pop(prel)
                excluded.append({"path": prel, "source": rel, "reason": "derived_from_excluded_dataset",
                                 "error": None})
            excluded.append({"path": rel, "source": src, "reason": reason, "error": (r.get("error") or "")[:300] or None})
            changed = True
        if not changed:
            return excluded, cut, ev
    raise SystemExit("evaluation gate did not converge")


def evaluation_summary(ev, alld):
    out = defaultdict(lambda: {"datasets": 0, "values": 0, "median_abs": []})
    for rel in alld:
        r = ev.get(rel) or {}
        if r.get("status") != "ok":
            continue
        g = "/".join(rel.split("/")[:2]) if rel.startswith(("dft/", "derived/")) else rel.split("/")[0]
        k = "%s|%s" % (g, r["kind"])
        out[k]["datasets"] += 1
        out[k]["values"] += r["n"]
        if r.get("median_abs") is not None:
            out[k]["median_abs"].append(r["median_abs"])
    res = {}
    for k in sorted(out):
        v = out[k]
        m = sorted(v["median_abs"])
        res[k] = {"datasets": v["datasets"], "values": v["values"],
                  "median_of_dataset_median_abs": rnd(m[len(m) // 2], 4) if m else None,
                  "p90_of_dataset_median_abs": rnd(m[int(0.9 * (len(m) - 1))], 4) if m else None}
    return {"rule": EVAL_RULE, "datasets_evaluated": len(alld),
            "all_finite": all((ev.get(r) or {}).get("status") == "ok" for r in alld), "by_folder_kind": res}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--corpus", default=lit.DEFAULT_CORPUS, help="SGTE-corpus build directory (read only)")
    ap.add_argument("--pdfs", default=lit.DEFAULT_PDFS, help="PDF corpus root (read only; paper ids of no-DOI papers)")
    ap.add_argument("--data", default=orchid_cite.DEFAULT_DATA, help="orchid-data checkout (reads; writes espei/)")
    ap.add_argument("--espei-python", default=DEFAULT_ESPEI_PY)
    ap.add_argument("--check", action="store_true", help="verify the files on disk match a fresh build")
    a = ap.parse_args(argv)
    out = os.path.join(a.data, NAME)
    pub, alias, papers, lindex = published(os.path.join(a.data, "literature"))
    C = build_corpus(a.corpus, a.pdfs, a.data, pub, alias, papers)
    pairs, dp = calibration_pairs(a.data, pub)
    D = build_dft(a.data, pub, papers, pairs)
    P = phase_constraints(C, D, a.data, pub)
    dp_release = "dft_parity.json inputs literature sha256 %s" % ((dp.get("inputs") or {}).get("literature") or {}).get("sha256")
    alld = dict(C["datasets"])
    alld.update(D["datasets"])
    alld.update(P["datasets"])
    v1 = run_espei(a.espei_python, {rel: dumps(d) for rel, d in alld.items()}, "check")
    excluded, other = [], {}
    for rel, msg in sorted(v1["invalid"].items()):
        code = next((c for s, c in ESPEI_EXCLUSIONS.items() if s in msg), None)
        if code and rel in C["rows"]:
            excluded.append({"path": rel, "source": C["rows"][rel]["source_dataset"], "reason": code, "error": msg[:400]})
        else:
            other[rel] = msg
    if other:
        for rel, msg in list(other.items())[:20]:
            print("INVALID %s: %s" % (rel, msg[:300]), file=sys.stderr)
        raise SystemExit("build failed: %d dataset(s) invalid for ESPEI (not of an excluded class)" % len(other))
    for e in excluded:
        C["datasets"].pop(e["path"])
        C["rows"].pop(e["path"])
        for rel in [r for r, d in P["datasets"].items() if (d.get("orchid") or {}).get("derived_from") == e["path"]]:
            P["datasets"].pop(rel)
            P["rows"].pop(rel)
    ev_excluded, ev_cut, ev = evaluation_gate(C, D, P, a.data, a.espei_python)
    for rel, d in C["datasets"].items():                     # certification over the values finally written
        certification_block(d, [(None, r) for r in C["_recs"][rel]])
    alld = dict(C["datasets"])
    alld.update(D["datasets"])
    alld.update(P["datasets"])
    v2 = run_espei(a.espei_python, {rel: dumps(d) for rel, d in alld.items()}, "load")
    val = {"espei_version": v2["espei_version"], "folders": v2["folders"], "invalid_final": v2["invalid"],
           "excluded": excluded, "evaluation_excluded": ev_excluded, "evaluation_regions_cut": ev_cut,
           "likelihood_evaluation": evaluation_summary(ev, alld)}
    if v2["invalid"] or any("error" in x or x.get("loaded") != x.get("files") for x in v2["folders"].values()):
        raise SystemExit("build failed: final bundle does not validate: %s" % json.dumps(v2)[:2000])
    # sigma statistics only over datasets actually written
    C["sigma_stats"] = [(C["rows"][k]["folder"], C["rows"][k]["output"], C["_ustat"][k]) for k in sorted(C["rows"])]
    files, S = build_files(C, D, val, a.data, lindex.get("release"), dp_release, P)
    managed = set(files)
    existing = set()
    for g in ("experimental", "model", "invariants", "dft", "assessed", "derived"):
        for p in glob.glob(os.path.join(out, g, "**", "*.json"), recursive=True):
            existing.add(os.path.relpath(p, out).replace(os.sep, "/"))
    for p in glob.glob(os.path.join(out, "dft", "*.md")):
        existing.add(os.path.relpath(p, out).replace(os.sep, "/"))
    stale = sorted(existing - managed)
    if a.check:
        bad = [rel for rel, data in sorted(files.items())
               if not os.path.exists(os.path.join(out, rel)) or open(os.path.join(out, rel), "rb").read() != data]
        if bad or stale:
            print("out of date: %d changed/missing, %d stale" % (len(bad), len(stale)), file=sys.stderr)
            for x in (bad + stale)[:20]:
                print("  " + x, file=sys.stderr)
            return 1
        print("ok: %d files current" % len(files))
        return 0
    for rel in stale:
        os.remove(os.path.join(out, rel))
    for rel, data in files.items():
        p = os.path.join(out, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        if os.path.exists(p) and open(p, "rb").read() == data:
            continue
        with open(p + ".tmp", "wb") as fh:
            fh.write(data)
        os.replace(p + ".tmp", p)
    for g in GROUPS + PGROUPS:
        for dp_ in sorted(glob.glob(os.path.join(out, g, "*")), reverse=True):
            if os.path.isdir(dp_) and not os.listdir(dp_):
                os.rmdir(dp_)
    size = sum(len(v) for v in files.values())
    print(json.dumps({"summary": S, "excluded_by_espei": [(e["reason"], e["path"]) for e in excluded],
                      "files": len(files), "bytes": size}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
