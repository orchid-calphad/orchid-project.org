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

SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_ESPEI_PY = os.path.expanduser("~/.venvs/espei/bin/python")
NAME = "espei"
STATUS_LINE = "Version 1, machine-extracted"

SOURCES = [("general/espei", "experimental"), ("general/espei-model", "model"), ("invariants/espei", "invariants")]
FOLDER_CLASS = {"experimental": "experimental", "model": "model", "invariants": "experimental"}
GROUPS = ("experimental", "model", "invariants", "dft/literature", "dft/mp", "dft/oqmd")

CAVEAT_FLAGS = ("magnitude_suspect", "system_mismatch", "source_disputed", "dft_outlier_unexplained",
                "basis_inferred_from_magnitude", "basis_corrected_from_formula", "sign_from_heading",
                "unsigned_in_print", "judge_adjudicated", "page_unreadable")
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
                 "DFT against experiment) and the ESPEI weight derived from it.")
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
def build_corpus(corpus, pdfs, data, pub, alias, papers):
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
           "match_kind": Counter(), "sigma_stats": [], "_ustat": {}, "release_has_ra": release_has_ra}
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
                                                    ("caveats", cv or None)) if v is not None})
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
                ("values_dropped", (npts - len(keep)) or None), ("espei_evaluation", ev)) if v is not None}
            out["orchid_uncertainty"] = ublock
            syst = lit.system_name([c for c in d["components"] if c.upper() != "VA"])
            rel = "%s/%s/%s.%s" % (folder, syst, tdb_dir, fname)
            if rel in res["datasets"]:
                raise SystemExit("path collision: " + rel)
            res["datasets"][rel] = out
            Ts = out["conditions"]["T"] if isinstance(out["conditions"]["T"], list) else [out["conditions"]["T"]]
            nvals = len(out["values"]) if zpf else len(flat_values(out))
            res["rows"][rel] = row(rel, folder, syst, out, nvals, npts - len(keep), min(Ts), max(Ts), cls, dois, pids,
                                   ids, ra_counts, "release" if release_has_ra else "record_store", caveats,
                                   tdb=(me.get("tdb") or {}).get("name"), sha=(me.get("tdb") or {}).get("sha256"),
                                   link=me.get("link_basis"), ev=ev, src="db-v2.1/" + src_rel, weight=weight,
                                   ustat=ustat, unit=ublock.get("unit"))
            res["_ustat"][rel] = ustat
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
        r = subprocess.run([espei_py, "-I", "-c", VALIDATOR, tmp, mode, json.dumps(GROUPS)], capture_output=True,
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
                for g in GROUPS if g in val["folders"]],
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
        "ESPEI then stops with a KeyError: add the element phase to your working TDB or fit those datasets with another.",
        "- **Energy basis.** ESPEI reads energies per mole of atoms. Literature values whose per-atom basis was not "
        "printed carry the caveat `energy_basis_unstated`.",
        "- **ZPF encoding.** Phase labels are mapped to the model phases of the assessment the source folder is named "
        "after (tag `via_assessment`); a printed compound composition that is its formula rounded is written at the "
        "exact formula value; a `__HYPERPLANE__` vertex is added only where every known composition comes from "
        "stoichiometry. Range temperatures are written at the midpoint (`orchid_value_bases`).",
        "- Values removed by the inclusion gate and datasets ESPEI rejects are listed in VALIDATION.md.", "",
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
reason); `summary`.
"""


def build_files(C, D, val, data, release, dp_release):
    oc = cite_module()
    cite = oc.cite_as(NAME, data)
    S = summarize(C, D, val)
    files = {}
    alld = dict(C["datasets"])
    alld.update(D["datasets"])
    for rel in sorted(alld):
        files[rel] = dumps(alld[rel])
    allr = dict(C["rows"])
    allr.update(D["rows"])
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
                                      "exclusion_reasons": ESPEI_EXCLUSION_TEXT,
                                      "dropped_datasets": C["dropped_datasets"], "dropped_values": C["drops"],
                                      "drop_reasons": DROP_REASONS,
                                      "dft_excluded": sorted(D["excluded"], key=lambda e: json.dumps(e, sort_keys=True)),
                                      "summary": S})
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
    dp_release = "dft_parity.json inputs literature sha256 %s" % ((dp.get("inputs") or {}).get("literature") or {}).get("sha256")
    alld = dict(C["datasets"])
    alld.update(D["datasets"])
    v1 = run_espei(a.espei_python, {rel: dumps(d) for rel, d in alld.items()}, "check")
    excluded, other = [], {}
    for rel, msg in sorted(v1["invalid"].items()):
        code = next((c for s, c in ESPEI_EXCLUSIONS.items() if s in msg), None)
        if code and not rel.startswith("dft/"):
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
    alld = dict(C["datasets"])
    alld.update(D["datasets"])
    v2 = run_espei(a.espei_python, {rel: dumps(d) for rel, d in alld.items()}, "load")
    val = {"espei_version": v2["espei_version"], "folders": v2["folders"], "invalid_final": v2["invalid"],
           "excluded": excluded}
    if v2["invalid"] or any("error" in x or x.get("loaded") != x.get("files") for x in v2["folders"].values()):
        raise SystemExit("build failed: final bundle does not validate: %s" % json.dumps(v2)[:2000])
    # sigma statistics only over datasets actually written
    C["sigma_stats"] = [(C["rows"][k]["folder"], C["rows"][k]["output"], C["_ustat"][k]) for k in sorted(C["rows"])]
    files, S = build_files(C, D, val, a.data, lindex.get("release"), dp_release)
    managed = set(files)
    existing = set()
    for g in ("experimental", "model", "invariants", "dft"):
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
    for g in ("experimental", "model", "invariants", "dft/literature", "dft/mp", "dft/oqmd"):
        for dp_ in sorted(glob.glob(os.path.join(out, g, "*")), reverse=True):
            if os.path.isdir(dp_) and not os.listdir(dp_):
                os.rmdir(dp_)
    size = sum(len(v) for v in files.values())
    print(json.dumps({"summary": S, "excluded_by_espei": [(e["reason"], e["path"]) for e in excluded],
                      "files": len(files), "bytes": size}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
