#!/usr/bin/env python3
"""Assessed phase assemblages: one (T, x, P) point per named phase field of the expert-drawn (assessed) binary phase
diagrams annotated from Okamoto's Desk Handbook (run-20261021h-full), with ESPEI ZPF exports (2026-10-23).

Evidence class `assessed_phase_diagram` (an expert drew the diagram from experimental evidence, without a CALPHAD
model). The primary evidence of each record is the thermodynamic coordinate (T, x, P) and the assemblage present there;
the point is placed deep inside the field (area centroid, or the pole of inaccessibility when the centroid is not
safely inside) so that the claim does not depend on where the boundaries were drawn. Boundary positions carry low
credence: they appear only as a descriptive positional margin and, for two-phase points, as an optional `tie_line`
block marked drawing-based.

Licensing (PI decision 2026-10-21): data only -- phase names/assemblages, presence facts, temperatures, compositions,
phase metadata, citations (original evaluation first, the handbook as compilation), method, confidence, basis and
review state. Never images, crops, renders, paths to page images, or handbook text beyond phase names and citations.
`guard()` scans every output string and key and fails the build on anything of that kind.

Inputs (read only): the annotation run (calibrated.json, evidence.json, metadata.json, grid.npz per diagram;
summary.jsonl), the DFT-hull decisions of tier 4, the phase registry (<data>/phases/systems/<SYS>.json), the calculated
SGTE diagrams (<data>/phase-diagrams/sgte/*.json) and the SGTE TDB files (phase mapping, via the annotation pipeline's
tools/phase_mapper.py, version-controlled; its git blob id is recorded in the method). Output: <data>/assessed-phase-assemblages/. Deterministic; no network, no model calls.

  ~/orchid-py313/bin/python tools/build_assessed_assemblages.py            # write
  ~/orchid-py313/bin/python tools/build_assessed_assemblages.py --check    # exit 1 if the outputs differ from a rebuild
Needs numpy, scipy and pyarrow; ESPEI validation runs in --espei-python (espei.datasets.check_dataset on every
dataset and load_datasets on each directory; the build fails on any invalid dataset).
"""
from __future__ import annotations

import argparse
import collections
import copy
import csv
import glob
import gzip
import hashlib
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from scipy import ndimage

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dashed_boundaries  # noqa: E402  (solid vs dashed boundary lines; James Saal review 2026-10-24)
import doi_corrections  # noqa: E402
import ml_folds  # noqa: E402
import orchid_cite  # noqa: E402

DOI_CORR = doi_corrections.load()   # shared DOI correction table (tools/doi_corrections.json)
EXPERT_CORR = collections.defaultdict(list)   # tools/assessed_expert_corrections.json: expert decisions over machine names
for _c in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "assessed_expert_corrections.json"),
                         encoding="utf-8"))["corrections"]:
    EXPERT_CORR[_c["diagram"]].append(_c)

HOME = os.path.expanduser("~")
OPEN_TDB = os.path.join(HOME, "projects", "open-tdb")
DEF_RUN = os.path.join(OPEN_TDB, "build", "okamoto-annotation", "run-20261021h-full")
DEF_HULL = os.path.join(OPEN_TDB, "build", "okamoto-annotation", "tier4", "hull", "review_input", "decisions_hull.jsonl")
DEF_TOOLS = os.path.dirname(os.path.abspath(__file__))   # tools/phase_mapper.py (version-controlled since audit round 2)
DEF_TDB_DIR = os.path.join(OPEN_TDB, "databases", "sgte-binary")
DEF_ESPEI = os.path.join(HOME, ".venvs", "espei", "bin", "python")

DATASET = "assessed-phase-assemblages"
OUT_DIR = "assessed-phase-assemblages"
SCHEMA_VERSION = "assessed-phase-assemblage-v1"
METHOD_VERSION = "1.0.0"
EVIDENCE_CLASS = "assessed_phase_diagram"
MARGIN = 0.006            # minimum distance of a point to any drawn line or frame, in frame-normalised units (~1.4 pt)
EDGE_TOL = 0.006          # a field reaching within this (normalised) of a cut-off x frame edge is "partly drawn"
PURE_TOL = 0.005          # an x frame edge within this mole fraction of 0 or 1 is a pure-element edge, not a cut
SCHEMATIC_T_MIN_K = 298.15   # a point of a full-width field without a drawn lower boundary must lie above room temperature
LINE_COMPOUND_AT = 1.0    # a phase whose printed range is at most this wide (at.%) is treated as a line compound
GAP_MAX = 0.02            # tie-line: a drawn line wider than this (normalised) is not a clean single boundary
P_PA = 101325.0
MAX_FILE_BYTES = 45 * 1000 * 1000

# ---- dense field grid (PI decision 2026-10-24: "the dense assessed phase assemblages are facts"; staged change S1)
# Interior points on a regular (T, x) grid inside every named field, in addition to the primary point of each field.
# Method version 1.0.0 is the grid of the orchid-espei-demo Amendment 4 / PhaseForge+ Cu-Zr final fit, validated by
# reproducing its 177 points (tools/tests/test_dense_field_grid.py). Enabled by default only where validated.
DENSE_VERSION = "1.0.0"
DENSE_DIR = "dense-field-grid"
DENSE_ESPEI_DIR = "espei-zpf-dense"
DENSE_VALIDATED = {"CU-ZR-p371"}
DENSE_STEP_T_K, DENSE_STEP_X, DENSE_MAX_PER_FIELD, DENSE_LIQUID_BAND_K = 50.0, 0.05, 40, 150.0
# Phase labels the release leaves unmapped on its primary records but that the validated grid resolves (demo decisions,
# recorded on every dense record that uses them; the primary records keep phase_id null): label -> (registry id, reason)
DENSE_LABEL_OVERRIDES = {
    "CU-ZR-p371": {
        "CuZr": ("phase:cu-zr:cuzr-b2", "the printed structure cP2 (CsCl, B2) is the registry phase CuZr-B2; the release's "
                 "primary records leave the label unmapped because SGTE assessments disagree (BCC_B2 vs CUZR_B2)"),
        "Cu51Zr14": ("phase:cu-zr:cu51zr14", "the printed composition 21.5 at.% Zr and structure (hP65, Ag51Gd14 type) match "
                     "the only registry phase of that name; left unmapped on the primary records by the deterministic mapper"),
    }}
CAVEAT = ("drawing-based; boundary positions in assessed diagrams carry low credence -- use as weak/soft constraints")
COMPILATION = {"role": "compilation", "citation": "H. Okamoto, Desk Handbook: Phase Diagrams for Binary Alloys, 2nd ed., "
               "ASM International, Materials Park, OH (2010)", "doi": None}

# ---------------------------------------------------------------------------------------------------- guard --
FORBIDDEN_KEY = re.compile(r"image|png|jpe?g|crop|render|page_?text|base64|thumbnail|screenshot|svg|bitmap|raster|"
                           r"region_pdf|pdf_path|review_?html|packet", re.I)
FORBIDDEN_VAL = [re.compile(p, re.I) for p in (
    r"data:image", r"\.(png|jpe?g|gif|svg|webp|tiff?|bmp|pdf|html?)\b", r"base64", r"[A-Za-z0-9+/]{120,}={0,2}",
    r"/build/", r"okamoto-annotation", r"okamoto-pilot", r"/Users/", r"\bcrop", r"page\.png", r"review\.html")]
LONG_OK = {"caveat", "method", "comment", "reference", "note", "pressure_basis", "composition_basis", "description",
           "definition", "uncertainty_method", "margin_note", "mapping_basis", "sgte_results", "original_citation",
           "citation"}
MAX_STR, MAX_STR_DOC = 400, 1500


def guard(obj, where, path=""):
    """Raise SystemExit on any key or string that could carry handbook images, crops, page paths or page text."""
    bad = []

    def walk(o, p, key):
        if isinstance(o, dict):
            for k, v in o.items():
                if FORBIDDEN_KEY.search(str(k)):
                    bad.append("%s%s: forbidden key %r" % (where, p, k))
                walk(v, p + "." + str(k), str(k))
        elif isinstance(o, (list, tuple)):
            for i, v in enumerate(o):
                walk(v, "%s[%d]" % (p, i), key)
        elif isinstance(o, str):
            lim = MAX_STR_DOC if key in LONG_OK else MAX_STR
            if len(o) > lim:
                bad.append("%s%s: string of %d chars (limit %d)" % (where, p, len(o), lim))
            for rx in FORBIDDEN_VAL:
                if rx.search(o):
                    bad.append("%s%s: forbidden content %r" % (where, p, rx.pattern))
                    break
    walk(obj, path, "")
    if bad:
        raise SystemExit("LICENSING GUARD FAILED (%d):\n  " % len(bad) + "\n  ".join(bad[:30]))


def guard_text(text, where):
    """CSV and other flat text: every cell through the same value rules."""
    rows = list(csv.reader(io.StringIO(text)))
    head = rows[0] if rows else []
    for i, row in enumerate(rows):
        guard({(head[j] if j < len(head) else "c%d" % j): cell for j, cell in enumerate(row)}, "%s:row%d" % (where, i))


# -------------------------------------------------------------------------------------------------- helpers --
def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def git_blob(p):
    """{path, git_blob} of a version-controlled tool: the git object id is content-addressed (git hash-object), so the
    exact code is recoverable from the repository history; path relative to the repository root."""
    data = open(p, "rb").read()
    blob = hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest()
    return {"path": "orchid-project.org/tools/" + os.path.basename(p), "git_blob": blob,
            "repository": "https://github.com/orchid-calphad/orchid-project.org"}


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def dumps(o):
    return json.dumps(o, ensure_ascii=False, separators=(",", ":"), sort_keys=False)


def r6(v):
    return None if v is None else float("%.6g" % v)


def r10(v):
    """coordinates: full precision of the raster/calibration (10 significant digits strips float noise only)"""
    return None if v is None else float("%.10g" % v)


def cap(e):
    return e[0].upper() + e[1:].lower()


def sys_key(elements):
    return "-".join(sorted(e.upper() for e in elements))


def last_status(run):
    st = {}
    p = os.path.join(run, "summary.jsonl")
    if os.path.exists(p):
        for line in open(p, encoding="utf-8"):
            r = json.loads(line)
            st["%s-p%s" % (r["system"], r["page"])] = r.get("status")
    return st


# ------------------------------------------------------------------------------------------ basis / status --
def field_status(F):
    """named | open | undetermined_in_source | tracing_fault | awaiting_expert | cannot_tell (as the v1 field export)."""
    if F.get("name"):
        return "named"
    t5 = (F.get("tier5") or {}).get("decision")
    if t5 == "cannot_tell":
        return "cannot_tell"
    if F.get("source_undetermined") or t5 == "source_undetermined":
        return "undetermined_in_source"
    if F.get("tier4_verdict") == "tracing_fault" or t5 == "tracing_fault":
        return "tracing_fault"
    if F.get("tier5_pending"):
        return "awaiting_expert"
    return "open"


BASIS_LABEL = {
    "deterministic_printed_label": "deterministic: the field's printed label",
    "deterministic_phase_rules": "deterministic: phase rules (adjacency, vanishing, ZPF enclosure, invariant rule) solved globally",
    "deterministic_adjacency": "deterministic: adjacency to named neighbours (law of adjacent phase regions)",
    "reader_consensus": "model readers in agreement (tiers 1-3: Gemini 3.1 Pro + Claude Haiku; two Claude Opus; Fable + GPT-6 Astra)",
    "adjudicator": "adjudicator (tier 4, Claude Opus 5.5 with zoomed views and a written reason)",
    "dft_hull_corroborated": "adjudicator proposal corroborated by DFT convex hulls (Materials Project, OQMD; 0 K); not handbook evidence",
    "expert_reviewed": "decided or confirmed by the domain expert (tier 5)",
}


def basis_of(F, hull):
    """(basis, confidence, review_state, flags) of a named field."""
    b = str(F.get("basis") or "")
    t5 = F.get("tier5") or {}
    flags = []
    if hull is not None:
        sup = hull.get("hull_support")
        conf = {"both": "high", "disagree": "low"}.get(sup) or t5.get("confidence") or hull.get("confidence") or "medium"
        if t5 and t5.get("confidence") == "low":
            conf = "low"
        flags.append("hull_support:" + str(sup))
        return "dft_hull_corroborated", conf, ("expert-reviewed (tier 5)" if t5 else "not reviewed by the expert"), flags
    if t5 or b.startswith("tier 5"):
        return "expert_reviewed", (t5.get("confidence") or "high"), "expert decision (tier 5)", flags
    if "tier 4" in b:
        conf = F.get("confidence") or ("medium" if "medium" in b else "high")
        return "adjudicator", conf, "not reviewed by the expert", flags
    if "readers" in b:
        return "reader_consensus", F.get("confidence") or "high", "not reviewed by the expert", flags
    if b == "printed_label":
        return "deterministic_printed_label", "high", "not reviewed by the expert", flags
    if b == "all_solutions" or b.startswith("rules"):
        return "deterministic_phase_rules", "high", "not reviewed by the expert", flags
    if b == "local_conflict":
        flags.append("local_conflict")
        return "deterministic_adjacency", "low", "not reviewed by the expert", flags
    return "deterministic_adjacency", "medium", "not reviewed by the expert", flags


# elements with a single solid structure at ambient pressure (Pearson symbol): a terminal solution (X) printed with
# another structure is a misread phase-table row (checked: 11 of 1,871 rows in run h), so it is left unmapped
ELEMENT_ONLY_SOLID = {"Ni": "cF4", "Al": "cF4", "Cu": "cF4", "Ag": "cF4", "Au": "cF4", "Pt": "cF4", "Pd": "cF4", "Ir": "cF4",
                      "Rh": "cF4", "Pb": "cF4", "Cr": "cI2", "Mo": "cI2", "W": "cI2", "V": "cI2", "Nb": "cI2", "Ta": "cI2",
                      "K": "cI2", "Rb": "cI2", "Cs": "cI2", "Ba": "cI2", "Eu": "cI2", "Zn": "hP2", "Cd": "hP2", "Re": "hP2",
                      "Os": "hP2", "Ru": "hP2", "Mg": "hP2", "Ge": "cF8", "Si": "cF8", "Ga": "oC8", "In": "tI2",
                      "Bi": "hR2", "Sb": "hR2", "As": "hR2"}
# cF8 is not here (2026-10-24, James Saal review): it is shared by diamond (A4), rock salt (B1) and sphalerite (B3), so a
# Pearson symbol alone never names DIAMOND_A4; and the Pearson fallback applies only when no Strukturbericht is printed
STD_PEARSON = {"cF4": "A1", "cI2": "A2", "hP2": "A3", "tI4": "A5"}
ORDER_CONF = {"high": 2, "medium": 1, "low": 0}


# ------------------------------------------------------------------------------------------- phase mapping --
class Registry:
    """Handbook label -> phase-registry id, per binary system; TDB phase -> registry id. Never invents a mapping."""

    def __init__(self, data, tdb_dir, tools):
        self.data, self.tdb_dir = data, tdb_dir
        sys.path.insert(0, tools)
        import phase_mapper  # noqa: E402  (annotation pipeline, deterministic)
        self.PM = phase_mapper
        self.pm_sha = sha256_file(os.path.join(tools, "phase_mapper.py"))
        self.pm_git = git_blob(os.path.join(tools, "phase_mapper.py"))
        self.sys_cache, self.tdb_cache, self.od, self.partners = {}, {}, {}, {}
        self.grid_index = collections.defaultdict(list)   # system -> calculated diagram files (loaded per system)
        for p in sorted(glob.glob(os.path.join(data, "phase-diagrams", "sgte", "*.json"))):
            els = re.findall(r"[A-Z][a-z]?", os.path.basename(p).split("-")[0])
            if len(els) == 2:
                self.grid_index[sys_key(els)].append(p)
        self._grids = (None, [])

    def grids_of(self, key):
        if self._grids[0] != key:
            out = []
            for p in self.grid_index.get(key, []):
                d = json.load(open(p, encoding="utf-8"))
                if d.get("status") == "ok" and d.get("grid") and sys_key(d["elements"]) == key:
                    out.append(d)
            self._grids = (key, out)
        return self._grids[1]

    def system(self, key):
        if key not in self.sys_cache:
            p = os.path.join(self.data, "phases", "systems", key + ".json")
            phases = json.load(open(p, encoding="utf-8"))["phases"] if os.path.exists(p) else []
            by_tdb, by_paper, info = {}, collections.defaultdict(set), {}
            od = set()
            for e in phases:
                for mo in e.get("models", []):
                    subs = [tuple(sorted(c for c in s_.get("constituents", []) if c != "VA")) for s_ in mo.get("sublattices", [])]
                    multi = [t for t in subs if len(t) >= 2]
                    if mo.get("disordered_part") or len(multi) != len(set(multi)):
                        od.add((mo.get("tdb"), str(mo.get("phase_name")).upper()))
                tdbn = []
                for n in e["names"]:
                    nm = n["name"]
                    if nm.startswith("tdb:"):
                        _, f, ph = nm.split(":", 2)
                        by_tdb[(f, ph.upper())] = e["id"]
                        tdbn.append(ph.upper())
                    elif nm.startswith("paper:"):
                        rest = nm[len("paper:"):]
                        doi, _, lab = rest.rpartition(":")
                        by_paper[(doi.lower(), lab)].add(e["id"])
                va = any("VA" in s.get("constituents", []) for m in e.get("models", []) for s in m.get("sublattices", []))
                indep = {str(mo.get("phase_name")).upper() for mo in e.get("models", []) if not mo.get("disordered_part_of")}
                info[e["id"]] = {"preferred_name": e["preferred_name"], "tdb_names": sorted(set(tdbn)), "va": va,
                                 "independent_names": sorted(indep)}
            self.sys_cache[key] = (by_tdb, by_paper, info)
            self.od[key] = od
            # order/disorder partners: the ordered model and its DIS_PART disordered part are one Gibbs-energy model
            # (audit round 2, A2 class 6); canonical id = the disordered part's registry id
            part = {}
            for e in phases:
                for mo in e.get("models", []):
                    dp = mo.get("disordered_part")
                    if dp and (mo.get("tdb"), str(dp).upper()) in by_tdb:
                        part[e["id"]] = by_tdb[(mo.get("tdb"), str(dp).upper())]
            self.partners[key] = part
        return self.sys_cache[key]

    def tdb(self, fname, axis_order):
        k = (fname, axis_order)
        if k not in self.tdb_cache:
            p = os.path.join(self.tdb_dir, fname)
            if not os.path.exists(p):
                self.tdb_cache[k] = None
            else:
                self.tdb_cache[k] = (self.PM.tdb_phases(p), self.PM.tdb_compositions(p, list(axis_order)), sha256_file(p))
        return self.tdb_cache[k]

    def map_label(self, key, label, info, axis_order, source_dois):
        """(registry id or None, mapping_basis). axis_order = (other element, axis element), upper case."""
        by_tdb, by_paper, _ = self.system(key)
        el = label[1:-1] if label.startswith("(") and label.endswith(")") else None
        pear = (info.get("pearson_symbol") or "").rstrip("*")
        if el in ELEMENT_ONLY_SOLID and pear and pear != ELEMENT_ONLY_SOLID[el]:
            return None, ("unmapped: printed structure %s contradicts the only solid structure of %s (%s); phase-table "
                          "reading suspect" % (pear, el, ELEMENT_ONLY_SOLID[el]))
        base = re.sub(r"#\d+$", "", label) if not re.fullmatch(r"L#\d+", label) else label
        # A compound can share hP2/cF4 with a terminal solution. That is not identity evidence.
        formula_named = len(self.PM.formula_dict(base) or {}) >= 2 or bool(re.search(r"[+-]\s*x", base))
        ids, bases = collections.defaultdict(list), []
        for g in self.grids_of(key):
            t = self.tdb(g["tdb"], axis_order)
            if not t:
                continue
            tdbp, comps, sha = t
            if sha != g["sha256"]:
                continue
            m, basis, _ = self.PM.map_phase(base, info, tdbp, comps, list(axis_order))
            if not m and not formula_named and not (base.startswith("(") and base.endswith(")")):
                # a letter-named phase with a standard solution structure (e.g. Ag-Al delta, A3) is that TDB solution
                # phase, as in the pipeline's compare_handbook_tdb.py: Strukturbericht first, else the A-type Pearson symbol
                sb = re.sub(r"[_\s]", "", info.get("strukturbericht") or "").upper()
                sb = sb if sb in self.PM.STD else (None if sb else STD_PEARSON.get((info.get("pearson_symbol") or "").rstrip("*")))
                if sb and self.PM.STD[sb] in tdbp:
                    m, basis = self.PM.STD[sb], "letter_name_with_standard_solution_structure"
            if not m and not (base.startswith("(") and base.endswith(")")):
                # an ordered compound whose Strukturbericht token (B2, L12, D03 ...) names exactly one TDB phase
                # (e.g. AlNi B2 -> B2_BCC) whose model composition range overlaps the printed range
                tok = re.sub(r"[_\s]", "", info.get("strukturbericht") or "").upper()
                rng = info.get("composition_at_pct")
                if re.fullmatch(r"[B-L]\d+[A-Z]?", tok):
                    cand = [p for p in tdbp if tok in p.upper().split("_")]
                    if rng:
                        cand = [p for p in cand if p in comps and comps[p][0] <= rng[1] / 100 + 0.03 and comps[p][1] >= rng[0] / 100 - 0.03]
                    if len(cand) == 1:
                        m, basis = cand[0], "ordered_compound_by_structure_token_in_tdb_name"
            if m and (g["tdb"], m.upper()) in by_tdb:
                ids[by_tdb[(g["tdb"], m.upper())]].append("%s:%s (%s)" % (g["tdb"], m, basis))
        paper = set()
        for doi in source_dois:
            paper |= by_paper.get((doi.lower(), base), set())
        if len(ids) > 1:
            return None, "unmapped: assessments map it to different registry phases (%s)" % ", ".join(sorted(ids))
        if len(paper) > 1:
            return None, "unmapped: the original evaluation's label is linked to several registry phases"
        if ids and paper and set(ids) != paper:
            return None, "unmapped: phase_mapper and the original evaluation's label disagree"
        if ids:
            rid = next(iter(ids))
            ev = sorted(ids[rid])
            return rid, "phase_mapper via %d SGTE file(s), e.g. %s" % (len(ev), ev[0])
        if paper:
            return next(iter(paper)), "registry paper label of this diagram's original evaluation (same DOI, same label)"
        return None, "unmapped: no deterministic mapping"


# ------------------------------------------------------------------------------------------------ geometry --
def field_point(mask, xs, Ts):
    """One interior point: area centroid when inside with >= MARGIN clearance, else the pole of inaccessibility.
    Distances in frame-normalised units (each axis scaled to the frame span). Returns dict or None (too thin)."""
    nT, nx = mask.shape
    dxn, dTn = 1.0 / (nx - 1), 1.0 / (nT - 1)
    rows, cols = np.nonzero(mask)
    r0, r1, c0, c1 = rows.min(), rows.max(), cols.min(), cols.max()
    sub = np.zeros((r1 - r0 + 3, c1 - c0 + 3), dtype=bool)
    sub[1:-1, 1:-1] = mask[r0:r1 + 1, c0:c1 + 1]
    dist = ndimage.distance_transform_edt(sub, sampling=(dTn, dxn))
    # centroid in the physical axes (uniform grid: the mean of the cell coordinates)
    xc, Tc = float(xs[cols].mean()), float(Ts[rows].mean())
    ci, ri = int(np.argmin(np.abs(xs - xc))), int(np.argmin(np.abs(Ts - Tc)))
    inside = bool(mask[ri, ci])
    dc = float(dist[ri - r0 + 1, ci - c0 + 1]) if (r0 <= ri <= r1 and c0 <= ci <= c1) else 0.0
    centroid = {"x": r6(xc), "T_C": r6(Tc), "inside": inside, "clearance_normalized": r6(dc if inside else 0.0)}
    if inside and dc >= MARGIN:
        return {"method": "centroid", "row": ri, "col": ci, "clearance": dc, "centroid": centroid}
    k = int(np.argmax(dist))   # first maximum in row-major order: deterministic
    pr, pc = np.unravel_index(k, dist.shape)
    dp = float(dist[pr, pc])
    if dp < MARGIN:
        return {"method": None, "clearance": dp, "centroid": centroid}
    return {"method": "pole_of_inaccessibility", "row": int(pr + r0 - 1), "col": int(pc + c0 - 1), "clearance": dp,
            "centroid": centroid}


def dense_grid_cells(mask, xs, Ts, liquid_only):
    """Regular (T, x) grid inside one field (method DENSE_VERSION): spacing dT = min(50 K, height / 4), dx = min(0.05,
    width / 4) anchored at the lower-left corner of the field plus half a step; a node is kept when it lies in the field
    at least MARGIN (frame-normalised) from every other field or drawn line; at most DENSE_MAX_PER_FIELD points (even
    subsample); the single-phase liquid only within DENSE_LIQUID_BAND_K above its lower boundary at that x. Returns
    [(row, col, clearance)]."""
    nT, nx = mask.shape
    dxn, dTn = 1.0 / (nx - 1), 1.0 / (nT - 1)
    rows, cols = np.nonzero(mask)
    r0, r1, c0, c1 = rows.min(), rows.max(), cols.min(), cols.max()
    # distances on the bounding box plus a one-cell ring of outside cells: identical to the full raster (the nearest
    # outside cell of an inside cell is never farther than that ring), much faster
    dsub = ndimage.distance_transform_edt(np.pad(mask[r0:r1 + 1, c0:c1 + 1], 1), sampling=(dTn, dxn))[1:-1, 1:-1]
    dist = np.zeros(mask.shape)
    dist[r0:r1 + 1, c0:c1 + 1] = dsub
    T0, T1, x0, x1 = Ts[r0], Ts[r1], xs[c0], xs[c1]
    dT, dx = min(DENSE_STEP_T_K, (T1 - T0) / 4), min(DENSE_STEP_X, (x1 - x0) / 4)
    cand = []
    for T in np.arange(T0 + dT / 2, T1, dT):
        i = int(np.argmin(abs(Ts - T)))
        for x in np.arange(x0 + dx / 2, x1, dx):
            j = int(np.argmin(abs(xs - x)))
            if not mask[i, j] or dist[i, j] < MARGIN:
                continue
            if liquid_only and Ts[i] - Ts[np.nonzero(mask[:, j])[0].min()] > DENSE_LIQUID_BAND_K:
                continue
            cand.append((i, j, float(dist[i, j])))
    if len(cand) > DENSE_MAX_PER_FIELD:
        cand = [cand[k] for k in np.linspace(0, len(cand) - 1, DENSE_MAX_PER_FIELD).astype(int)]
    return cand


def dense_enabled(a, diag):
    mode = getattr(a, "dense", "default")
    return mode in ("all", "report") or (mode == "default" and diag in DENSE_VALIDATED)


def espei_names(asm, rinfo):
    """(TDB phase names, labels without a name, any vacancy sublattice) of an assemblage."""
    tnames, miss = [], []
    for e in asm:
        ri = rinfo.get(e["phase_id"]) if e["phase_id"] else None
        if ri and ri["preferred_name"].upper() in ri["tdb_names"] and ri["preferred_name"].upper() in ri.get("independent_names", []):
            tnames.append(ri["preferred_name"].upper())
        else:
            miss.append(e["label"])
    return tnames, miss, (any(rinfo[e["phase_id"]]["va"] for e in asm if e["phase_id"]) if not miss else False)


def run_extent(lab_row, c, fid, step):
    """From column c (inside field fid) walk along the row: (first column outside the field, first nonzero label
    beyond the drawn line or None at the frame, column of that label)."""
    n = len(lab_row)
    j = c
    while 0 <= j < n and lab_row[j] == fid:
        j += step
    j0 = j
    while 0 <= j < n and lab_row[j] == 0:
        j += step
    if not (0 <= j < n):
        return j0, None, None
    return j0, int(lab_row[j]), j


def col_extent(lab_col, r, fid):
    out = []
    for step in (-1, 1):
        j = r
        while 0 <= j < len(lab_col) and lab_col[j] == fid:
            j += step
        out.append(j)
    return out


# ------------------------------------------------------------------------------------------------- the build --
def build(a, log=print):
    run = a.run
    status_by_diag = last_status(run)
    hull = {}
    if a.hull and os.path.exists(a.hull):
        for line in open(a.hull, encoding="utf-8"):
            d = json.loads(line)
            hull[(d["diagram"], int(d["field"]))] = d
    reg = Registry(a.data, a.tdb_dir, a.tools)
    manifest = json.load(open(os.path.join(run, "run_manifest.json"), encoding="utf-8"))
    run_id = os.path.basename(os.path.normpath(run))
    method = {"name": "assessed phase-assemblage points (ORCHID)", "version": METHOD_VERSION, "annotation_run": run_id,
              "annotation_tools_sha256": manifest.get("tools_sha256", {}).get("pd_fields.py"),
              "phase_mapper_sha256": reg.pm_sha, "phase_mapper": reg.pm_git, "builder": "orchid-project.org/tools/build_assessed_assemblages.py"}
    records, exclusions, espei_excl = [], [], []
    a._dense = {"records": [], "espei_excl": [], "espei": {}, "espei_pts": collections.defaultdict(list),
                "counts": collections.Counter()}
    espei_default, espei_tl = {}, {}
    cnt = collections.Counter()
    for sd in sorted(glob.glob(os.path.join(run, "systems", getattr(a, "only", None) or "*"))):
        diag = os.path.basename(sd)
        need = [os.path.join(sd, f) for f in ("metadata.json", "evidence.json", "calibrated.json", "grid.npz")]
        dstat = status_by_diag.get(diag)
        if not all(os.path.exists(p) for p in need):
            exclusions.append({"diagram": diag, "field": None, "reason": "diagram_not_annotated",
                               "detail": dstat or "no annotation outputs"})
            cnt["diagram_not_annotated"] += 1
            continue
        m = json.load(open(need[0], encoding="utf-8"))
        ev = json.load(open(need[1], encoding="utf-8"))
        cal = json.load(open(need[2], encoding="utf-8"))
        system, elements = m["system"], [cap(e) for e in m["elements"]]
        key = sys_key(elements)
        comp = cap(cal.get("component") or elements[1])
        other = [e for e in elements if e != comp]
        other = other[0] if other else elements[0]
        axis_order = (other.upper(), comp.upper())
        fields = ev.get("fields", {})
        for c in EXPERT_CORR.get(diag, []):
            F = fields.get(str(c["field"]))
            if F is None or F.get("name") != c["machine_name"]:
                raise SystemExit("expert correction does not match the annotation: %s field %s" % (diag, c["field"]))
            F = fields[str(c["field"])] = dict(F)
            F["name_history"] = [{"name": c["machine_name"], "basis": F.get("basis"), "superseded_by": c["basis"],
                                  "reviewed_at": c["reviewed_at"]}]
            F["name"] = c["expert_name"]
            if c["basis"] == "page_check":   # read on the rendered page (audit round 2); not an expert decision
                F["page_check"] = {k: c[k] for k in ("reviewer", "reviewed_at", "source", "confidence") if c.get(k)}
                ev = c.get("evidence") or {}
                F["page_check"].update({"page_pdf": ev.get("page_pdf"), "pdf_sha256_prefix": ev.get("pdf_sha256_prefix"),
                                        "page_view_sha256": ev.get("image_sha256"), "seen": ev.get("seen")})
                cnt["page_check_corrections_applied"] += 1
                continue
            F["tier5"] = {"decision": c["expert_name"], "confidence": c.get("confidence") or "high",
                          "reviewer": c["reviewer"], "reviewed_at": c["reviewed_at"], "source": c["source"]}
            cnt["expert_corrections_applied"] += 1
        cnt["diagrams_seen"] += 1
        if comp not in elements:   # the composition axis was read as an element outside the system
            for fid in sorted(fields, key=int):
                exclusions.append({"diagram": diag, "field": int(fid), "reason": "axis_component_not_in_system",
                                   "detail": "composition axis read as %s, not an element of %s" % (comp, key)})
                cnt["excl:axis_component_not_in_system"] += 1
            continue
        if dstat == "tracing_suspect":
            for fid in sorted(fields, key=int):
                exclusions.append({"diagram": diag, "field": int(fid), "reason": "diagram_tracing_suspect",
                                   "detail": "the diagram's traced lines failed the self-check as a whole"})
                cnt["excl:diagram_tracing_suspect"] += 1
            continue
        g = np.load(need[3])
        lab, xs, Ts = g["lab"], g["x"], g["T"]
        nT, nx = lab.shape
        fpath = os.path.join(sd, "fields.json")
        fj = json.load(open(fpath, encoding="utf-8")) if os.path.exists(fpath) else {}
        lines = dashed_boundaries.DiagramLines(cal, fj, (float(xs[0]), float(xs[-1])), (float(Ts[0]), float(Ts[-1])))
        if fj and fj.get("dash_joins") is not None and fj["dash_joins"] != lines.n_dash_joins:
            raise SystemExit("%s: re-derived dash joins (%d) differ from the tracer's count (%s)"
                             % (diag, lines.n_dash_joins, fj["dash_joins"]))
        cnt["diagrams_with_dashed_lines"] += bool(lines.n_dash_joins)
        x_span, T_span = float(xs[-1] - xs[0]), float(Ts[-1] - Ts[0])
        dx_cell, dT_cell = x_span / (nx - 1), T_span / (nT - 1)
        res_x = float(cal["calibration"]["x"].get("max_residual") or 0.0) / 100.0   # at.% -> mole fraction
        res_T = float(cal["calibration"]["y"].get("max_residual_C") or 0.0)
        u_x, u_T = res_x + dx_cell / 2, res_T + dT_cell / 2
        cut_left, cut_right = xs[0] > PURE_TOL, xs[-1] < 1 - PURE_TOL
        edge_cols = max(1, int(round(EDGE_TOL * (nx - 1))))
        table = {p["name"]: p for p in m.get("phases") or []}
        srcs = []
        for s in m["sources"].get("diagram_sources") or []:
            s = dict(s)
            e = DOI_CORR.handbook(key, s["doi"]) if s.get("doi") else None
            if e:   # the DOI matched to the printed source line describes a different work: keep the printed citation
                DOI_CORR.mark(s, e)
            srcs.append(s)
        src_dois = [s["doi"] for s in srcs if s.get("doi")]
        hb = m["sources"]["handbook"]
        citations = [dict({"role": "original_evaluation", "citation": s.get("citation"), "doi": s.get("doi"),
                           "doi_matched_by": s.get("matched_by"), "kind": s.get("kind")},
                          **({"doi_status": s["doi_status"], "doi_correction": s["doi_correction"]} if s.get("doi_status") else {}))
                     for s in srcs]
        citations.append(dict(COMPILATION, page_pdf=hb.get("page"), locator=hb.get("locator")))
        by_tdb, by_paper, rinfo = reg.system(key)
        names = {int(r): (F["name"].split("+") if F.get("name") else None) for r, F in fields.items()}
        label_map = {}

        def phase_entry(label):
            if label not in label_map:
                base = re.sub(r"#\d+$", "", label) if not label.startswith("L#") else label
                t = table.get(base) or {}
                info = {"pearson_symbol": t.get("pearson_symbol"), "strukturbericht": t.get("strukturbericht"),
                        "prototype": t.get("prototype"), "composition_at_pct": t.get("composition_at_pct")}
                rid, why = reg.map_label(key, label, info, axis_order, src_dois)
                e = {"label": label, "canonical_name": rinfo[rid]["preferred_name"] if rid else label,
                     "phase_id": rid, "mapping_basis": why}
                for k_ in ("pearson_symbol", "strukturbericht", "prototype", "space_group_printed"):
                    if t.get(k_):
                        e[k_] = t[k_]
                if t.get("composition_at_pct"):
                    e["composition_range_at_pct"] = t["composition_at_pct"]
                label_map[label] = e
            return label_map[label]

        sgte = reg.grids_of(key)
        for k_ in [k_ for k_ in reg.sys_cache if k_ != key]:
            reg.sys_cache.pop(k_, None)
        espei_pts, espei_pts_tl = [], []
        for fid_s in sorted(fields, key=int):
            F = fields[fid_s]
            fid = int(fid_s)
            st = field_status(F)
            base_ex = {"diagram": diag, "field": fid}
            if st != "named":
                exclusions.append(dict(base_ex, reason=st, detail=None))
                cnt["excl:" + st] += 1
                continue
            h = hull.get((diag, fid)) if F.get("source_undetermined") else None
            mask = lab == fid
            if not mask.any():
                exclusions.append(dict(base_ex, reason="field_not_in_raster", detail=None))
                cnt["excl:field_not_in_raster"] += 1
                continue
            touches = (cut_left and mask[:, :edge_cols + 1].any()) or (cut_right and mask[:, nx - edge_cols - 1:].any())
            if touches and h is None:
                exclusions.append(dict(base_ex, reason="partly_drawn_truncated_edge",
                                       detail="the field reaches a cut-off composition edge of the drawing"))
                cnt["excl:partly_drawn_truncated_edge"] += 1
                continue
            phases = names[fid]
            if len(phases) > 2:   # a binary field holds one or two phases; three or more is a naming/tracing fault
                exclusions.append(dict(base_ex, reason="more_than_two_phases", detail=None, assemblage=phases))
                cnt["excl:more_than_two_phases"] += 1
                continue
            pt = field_point(mask, xs, Ts)
            basis, conf, review, flags = basis_of(F, h)
            if F.get("carried_against_rules"):
                flags.append("carried_against_rules")
                if ORDER_CONF.get(conf, 0) > 1:
                    conf = "medium"
            if F.get("phase_not_in_table"):
                flags.append("phase_not_in_table")
            if touches:
                flags.append("partly_drawn_dft_corroborated")
            if pt["method"] is None:
                exclusions.append(dict(base_ex, reason="too_thin_for_margin",
                                       detail="max interior clearance %.4f < margin %.4f (normalised)" % (pt["clearance"], MARGIN),
                                       assemblage=phases))
                cnt["excl:too_thin_for_margin"] += 1
                continue
            r, c = pt["row"], pt["col"]
            x, T_C = float(xs[c]), float(Ts[r])
            T_K = T_C + 273.15
            # schematic fields (audit round 2, A1 C8): a full-width two-phase field that reaches the bottom of the frame
            # has no drawn lower boundary; its centroid below room temperature (down to -347 degC) is an artefact of
            # where the frame ends, not a drawn state
            bottom = int(np.argmin(Ts))
            unbounded_below = bool(mask[bottom, :].any())
            full_width = bool(mask[:, :edge_cols + 1].any() and mask[:, nx - edge_cols - 1:].any())
            if T_K <= 0 or (unbounded_below and full_width and T_K < SCHEMATIC_T_MIN_K):
                exclusions.append(dict(base_ex, reason="schematic_field_unbounded_below",
                                       detail="point at %.1f K in a full-width field with no drawn lower boundary" % T_K,
                                       assemblage=phases))
                cnt["excl:schematic_field_unbounded_below"] += 1
                continue
            if unbounded_below and full_width:
                flags.append("schematic_field_unbounded_below")
            row = lab[r]
            jl, nl, kl = run_extent(row, c, fid, -1)
            jr, nr, kr = run_extent(row, c, fid, +1)
            cu, cd = col_extent(lab[:, c], r, fid)
            margin = {"normalized": r6(pt["clearance"]),
                      "isotherm_mole_fraction": r6(min(x - float(xs[max(jl, 0)]), float(xs[min(jr, nx - 1)]) - x)),
                      "isopleth_K": r6(min(T_C - float(Ts[max(cu, 0)]), float(Ts[min(cd, nT - 1)]) - T_C)),
                      "margin_note": "distance from the point to the nearest drawn line or frame edge; drawing-based, "
                                     "low credence; describes placement, not an uncertainty of any phase composition"}
            # boundary lines that bound the field along the point's isotherm and isopleth (the lines that set its
            # positional margin): printed solid, dashed (estimated by the evaluator) or both
            cross = [(float(xs[j_]), T_C) for j_ in (jl, jr) if 0 <= j_ < nx] + \
                    [(x, float(Ts[i_])) for i_ in (cu, cd) if 0 <= i_ < nT]
            bbasis, n_solid, n_dashed = lines.basis_at(cross)
            asm = [phase_entry(p) for p in phases]
            if any("reading suspect" in e["mapping_basis"] for e in asm):
                flags.append("phase_table_structure_suspect")
            rec = {"id": "aspa-" + sha256_bytes(("%s|%s|%d|0" % (system, diag, fid)).encode())[:16],
                   "record_type": "assessed_phase_assemblage", "evidence_class": EVIDENCE_CLASS,
                   "system": key, "elements": sorted(elements), "diagram": diag, "field": fid,
                   "assemblage": asm, "assemblage_label": "+".join(phases), "n_phases": len(phases),
                   "T_K": r10(T_K), "T_C": r10(T_C), "temperature_axis": "degC (as drawn); T_K = T_C + 273.15",
                   "x": r10(x), "x_component": comp, "composition_axis": "at.%% %s" % comp,
                   "composition_basis": "diagram drawn in at.%; x = at.%/100 (mole fraction); no wt.% conversion",
                   "P_Pa": P_PA, "pressure_basis": "1 atm assumed (binary diagram at ambient pressure; the diagram states no other pressure)",
                   "coordinate_uncertainty": {"T_K": r6(u_T), "x": r6(u_x), "calibration_max_residual_T_K": r6(res_T),
                                              "calibration_max_residual_x": r6(res_x),
                                              "calibration_pairs_T": cal["calibration"]["y"].get("pairs"),
                                              "calibration_pairs_x": cal["calibration"]["x"].get("pairs"),
                                              "uncertainty_method": "bound = max axis-calibration residual (printed tick labels "
                                              "vs fitted linear axis) + half a raster cell; applies to the point's coordinate"},
                   "point_method": pt["method"], "centroid": pt["centroid"] if pt["method"] != "centroid" else None,
                   "positional_margin": margin, "phase_compositions": None,
                   "boundary_basis": bbasis,
                   "boundary_lines": {"solid": n_solid, "dashed": n_dashed, "diagram_has_dashed_lines": bool(lines.n_dash_joins),
                                      "note": (dashed_boundaries.BASIS_TEXT.get(bbasis) or "no traced boundary line along the "
                                               "point's isotherm and isopleth (frame edges only)")},
                   "basis": basis, "basis_description": BASIS_LABEL[basis], "confidence": conf, "review_state": review,
                   "flags": flags, "citations": citations, "method": method}
            if F.get("name_history"):
                rec["name_history"] = F["name_history"]
                if F.get("tier5"):
                    rec["expert_review"] = {k: F["tier5"][k] for k in ("reviewer", "reviewed_at", "source", "confidence")}
                if F.get("page_check"):
                    rec["page_check"] = F["page_check"]
                    flags.append("name_corrected_on_page_check")
            # ---- tie line (two-phase points only; drawing-based, low credence)
            tl, why = None, None
            if len(phases) == 2:
                tl, why = tie_line(phases, (jl, nl, kl), (jr, nr, kr), xs, names, table, u_x, x, T_K, cut_left, cut_right)
                cnt["tie_line" if tl else "tie_line_omitted:" + why] += 1
            rec["tie_line"] = tl
            rec["uncertainty_budget"] = assessed_budget(rec)
            rec["tie_line_omitted"] = why
            # ---- SGTE cross-check
            rec["sgte_crosscheck"] = [crosscheck(gd, x, T_K, comp, asm, by_tdb, reg.od.get(key, frozenset()),
                                                 reg.partners.get(key, {})) for gd in sgte]
            records.append(rec)
            cnt["points"] += 1
            cnt["boundary_basis:%s" % bbasis] += 1
            cnt["method:" + pt["method"]] += 1
            # ---- ESPEI
            tnames, miss = [], []
            for e in asm:
                ri = rinfo.get(e["phase_id"]) if e["phase_id"] else None
                # an ESPEI ZPF phase must be an independent PHASE (never only a DIS_PART disordered part; audit
                # round 2, A2 class 7)
                if ri and ri["preferred_name"].upper() in ri["tdb_names"] and \
                        ri["preferred_name"].upper() in ri.get("independent_names", []):
                    tnames.append(ri["preferred_name"].upper())
                else:
                    miss.append(e["label"])
            if miss:
                espei_excl.append({"id": rec["id"], "diagram": diag, "field": fid, "reason": "phase without TDB name in the registry",
                                   "phases": miss})
                cnt["espei_excluded"] += 1
            else:
                va = any(rinfo[e["phase_id"]]["va"] for e in asm)
                espei_pts.append((rec, tnames, va, None))
                espei_pts_tl.append((rec, tnames, va, tl))
            # ---- dense field grid (additional records; the primary point above is unchanged)
            if dense_enabled(a, diag):
                cells = dense_grid_cells(mask, xs, Ts, phases == ["L"])
                dd = a._dense
                dd["counts"]["fields_with_grid"] += 1
                dd["counts"]["points"] += len(cells)
                dd["counts"]["diagram:" + diag] += len(cells)
                if getattr(a, "dense", "default") == "report":
                    continue
                asm_d = []
                for e in asm:
                    e2 = dict(e)
                    ov = DENSE_LABEL_OVERRIDES.get(diag, {}).get(e["label"])
                    if e2["phase_id"] is None and ov and ov[0] in rinfo:
                        e2.update(phase_id=ov[0], canonical_name=rinfo[ov[0]]["preferred_name"],
                                  mapping_basis="dense_field_grid override: " + ov[1])
                    asm_d.append(e2)
                for k_, (r_, c_, clr) in enumerate(cells, 1):
                    xg, T_Cg = round(float(xs[c_]), 6), float(Ts[r_])
                    row_ = lab[r_]
                    jl_, _, _ = run_extent(row_, c_, fid, -1)
                    jr_, _, _ = run_extent(row_, c_, fid, +1)
                    cu_, cd_ = col_extent(lab[:, c_], r_, fid)
                    mg = {"normalized": r6(clr),
                          "isotherm_mole_fraction": r6(min(xg - float(xs[max(jl_, 0)]), float(xs[min(jr_, nx - 1)]) - xg)),
                          "isopleth_K": r6(min(T_Cg - float(Ts[max(cu_, 0)]), float(Ts[min(cd_, nT - 1)]) - T_Cg)),
                          "margin_note": margin["margin_note"]}
                    cross_ = [(float(xs[j_]), T_Cg) for j_ in (jl_, jr_) if 0 <= j_ < nx] + \
                             [(xg, float(Ts[i_])) for i_ in (cu_, cd_) if 0 <= i_ < nT]
                    bb_, ns_, nd_ = lines.basis_at(cross_)
                    dr = copy.deepcopy(rec)
                    dr.update({
                        "id": "aspa-" + sha256_bytes(("%s|%s|%d|%d" % (system, diag, fid, k_)).encode())[:16],
                        "assemblage": asm_d, "T_K": T_Cg + 273.15, "T_C": T_Cg, "x": xg,   # raster values unrounded (grid of the validated fit)
                        "point_method": "dense_field_grid", "centroid": None, "positional_margin": mg, "boundary_basis": bb_,
                        "boundary_lines": {"solid": ns_, "dashed": nd_, "diagram_has_dashed_lines": bool(lines.n_dash_joins),
                                           "note": (dashed_boundaries.BASIS_TEXT.get(bb_) or "no traced boundary line along the "
                                                    "point's isotherm and isopleth (frame edges only)")},
                        "tie_line": None, "tie_line_omitted": "dense_field_grid_point",
                        "flags": list(flags) + ["dense_field_grid", "interior_grid"],
                        "dense_grid": {"method_version": DENSE_VERSION, "primary_record_id": rec["id"], "index_in_field": k_,
                                       "n_points_in_field": len(cells), "field_weight_share": r6(1.0 / len(cells)),
                                       "interior_margin_normalized_min": MARGIN, "clearance_normalized": r6(clr),
                                       "spacing": {"T_K_max": DENSE_STEP_T_K, "x_max": DENSE_STEP_X},
                                       "weight_note": "points of one diagram are not independent: the shares of one field sum "
                                                      "to one observation (a fit may instead weight every point equally, as the "
                                                      "PhaseForge+ Cu-Zr example does)"}})
                    dr["uncertainty_budget"] = assessed_budget(dr)
                    dr["sgte_crosscheck"] = [crosscheck(gd, xg, T_Cg + 273.15, comp, asm_d, by_tdb, reg.od.get(key, frozenset()),
                                                        reg.partners.get(key, {})) for gd in sgte]
                    dd["records"].append(dr)
                    tn_, miss_, va_ = espei_names(asm_d, rinfo)
                    if miss_:
                        dd["espei_excl"].append({"id": dr["id"], "diagram": diag, "field": fid,
                                                 "reason": "phase without TDB name in the registry", "phases": miss_})
                    else:
                        dd["espei_pts"][diag].append((dr, tn_, va_, None))
        if a._dense["espei_pts"].get(diag):
            a._dense["espei"][diag] = espei_dataset(diag, a._dense["espei_pts"][diag], comp, elements, citations, False, dense=True)
        if espei_pts:
            espei_default[diag] = espei_dataset(diag, espei_pts, comp, elements, citations, False)
            espei_tl[diag] = espei_dataset(diag, espei_pts_tl, comp, elements, citations, True)
    return records, exclusions, espei_excl, espei_default, espei_tl, cnt, method


def tie_line(phases, left, right, xs, names, table, u_x, x, T_K, cut_left, cut_right):
    """Isotherm through the point intersected with the field's drawn boundaries; each endpoint assigned to the phase
    whose field (or line compound) lies on that side. (block, None) or (None, reason)."""
    ends = []
    nx = len(xs)
    dxn = 1.0 / (nx - 1)
    for side, (j0, nb, k) in (("left", left), ("right", right)):
        if nb is None:   # the isotherm reaches the frame
            at_edge = float(xs[0] if side == "left" else xs[-1])
            if (side == "left" and cut_left) or (side == "right" and cut_right):
                return None, "isotherm_reaches_cut_off_frame_edge"
            pure = 0.0 if side == "left" else 1.0
            cands = [p for p in phases if p in table and table[p].get("composition_at_pct")
                     and min(abs(v - 100 * pure) for v in table[p]["composition_at_pct"]) <= 0.5]
            if len(cands) != 1:
                return None, "frame_edge_phase_ambiguous"
            ends.append({"phase": cands[0], "x": pure, "x_uncertainty": r6(u_x), "side": side,
                         "kind": "pure-element edge (no solubility drawn)"})
            continue
        gap = abs(k - j0) * dxn
        if gap > GAP_MAX:
            return None, "boundary_not_a_single_clean_line"
        xb = float((xs[j0] + xs[k - (1 if side == "right" else -1)]) / 2)
        nn = names.get(nb)
        if not nn:
            return None, "neighbour_field_unnamed"
        shared = [p for p in nn if p in phases]
        if len(nn) == 1 and shared:
            ends.append({"phase": shared[0], "x": r10(xb), "x_uncertainty": r6(u_x + gap / 2 * (xs[-1] - xs[0])),
                         "side": side, "kind": "drawn boundary to the single-phase field"})
            continue
        if len(shared) == 1:
            t = table.get(shared[0]) or {}
            rng = t.get("composition_at_pct")
            if rng and rng[1] - rng[0] <= LINE_COMPOUND_AT and min(abs(100 * xb - rng[0]), abs(100 * xb - rng[1])) <= 1.5:
                ends.append({"phase": shared[0], "x": r10((rng[0] + rng[1]) / 200.0), "x_drawn": r10(xb),
                             "x_uncertainty": r6(max(u_x, (rng[1] - rng[0]) / 200.0)), "side": side,
                             "kind": "line compound (composition from the printed phase table)", "line_compound": True})
                continue
        return None, "side_assignment_ambiguous"
    if len({e["phase"] for e in ends}) != 2 or set(e["phase"] for e in ends) != set(phases):
        return None, "endpoints_do_not_cover_the_two_phases"
    if not (ends[0]["x"] <= x <= ends[1]["x"]):
        return None, "point_outside_endpoints"
    return {"T_K": r10(T_K), "endpoints": ends, "method": "isotherm through field point, intersected with drawn boundaries",
            "caveat": CAVEAT}, None


def grid_cells(gd):
    if "_dec" not in gd:
        g = gd["grid"]
        dec = np.zeros((len(g["rows"]), g["x_points"]), dtype=np.int32)
        for i, row in enumerate(g["rows"]):
            j = 0
            for idx, n in row:
                dec[i, j:j + n] = idx
                j += n
        gd["_dec"] = dec
    return gd["_dec"]


def crosscheck(gd, x, T_K, comp, asm, by_tdb, od=frozenset(), partners=None):
    g = gd["grid"]
    base = {"tdb": gd["tdb"], "tdb_sha256": gd["sha256"]}
    xg = x if gd["x_component"].upper() == comp.upper() else 1.0 - x
    dec = grid_cells(gd)
    i = int(round((T_K - g["T0"]) / g["dT"]))
    k = int(round(xg / g["x_step"]))
    if i < 0 or i >= len(dec) or not 0 <= k < g["x_points"]:
        return dict(base, result="not_comparable", calculated=None, reason="outside the calculated T-x range")
    labl = g["labels"][int(dec[i, k])]
    base.update(cell={"T_K": r6(g["T0"] + i * g["dT"]), "x": r6(k * g["x_step"])}, calculated=labl or None)
    if not labl:
        return dict(base, result="not_comparable", reason="no calculated equilibrium at this cell")
    ids = [e["phase_id"] for e in asm]
    if any(i_ is None for i_ in ids):
        return dict(base, result="not_comparable", reason="assemblage phase without registry mapping")
    if len(set(ids)) < len(ids):
        return dict(base, result="not_comparable", reason="the same registry phase twice (the grid collapses composition sets)")
    calc = [(by_tdb.get((gd["tdb"], p.upper()), "tdb:" + p.upper()), (gd["tdb"], p.upper()) in od) for p in labl.split("+")]
    if {c for c, _ in calc} == set(ids):
        return dict(base, result="agrees", reason=None)
    canon = (lambda i_: (partners or {}).get(i_, i_))
    if len(calc) == len(ids) and {canon(c) for c, _ in calc} == {canon(i_) for i_ in ids}:
        # the ordered and the disordered state of one order/disorder model (B2_BCC / A2_BCC, L12 / FCC_A1): the grid
        # label does not resolve the order parameter, the handbook names the state; treated as agreement
        return dict(base, result="agrees", reason="order/disorder partners (one Gibbs-energy model)")
    if any(o for _, o in calc) and {c for c, o in calc if not o} <= set(ids) and len(calc) == len(ids):
        return dict(base, result="not_comparable",
                    reason="order-disorder model: the grid label does not tell the ordered from the disordered state")
    return dict(base, result="disagrees", reason=None)


# ---- uncertainty budget (PI decision 2026-10-24: components kept separate, as in the literature release)
BUDGET_CALIBRATION = "assessed-phase-assemblages handbook-class prior 2026-10-24"
BUDGET_RULE = ("max(sqrt(max(printed, reading, method)^2 + composition^2 + conversion^2), scatter); boundary margin "
               "reported, not combined")
HANDBOOK_PRIOR = {"T": 10.0, "x": 0.02}   # K, mole fraction (ORCHID assumption; see SCHEMA.md)
HANDBOOK_PRIOR_BASIS = (
    "ORCHID assumption (PI decision 2026-10-24), not fitted: an expert-drawn assessed diagram reproduces its "
    "experimental evidence to about 10 K and 0.02 mole fraction, the same priors the literature release uses for "
    "model / assessment values (tools/uncertainty_budget.py MODEL_PRIORS T_INVARIANT, X_PHASE)")


# dashed (estimated) boundaries (James Saal review 2026-10-24): the evaluator marks such a line as estimated, i.e. not
# fixed by measurements, so the reading component of a point bounded by one is enlarged to a multiple of the handbook
# class prior (ORCHID assumption, not fitted): 2x for a field bounded only by dashed lines along the point's isotherm
# and isopleth, 1.5x when solid and dashed lines both bound it; solid-only fields keep the digitisation bound.
DASHED_READING_FACTOR = {"dashed_estimated": 2.0, "mixed": 1.5}


def _bq(unit, vt, reading, prior, margin, boundary_basis=None):
    f = DASHED_READING_FACTOR.get(boundary_basis)
    rd = {"value": reading, "basis": "digitisation (max axis-calibration residual + half a raster cell)"} \
        if reading is not None else None
    if f:
        rd = {"value": r6(max(reading or 0.0, f * prior)), "digitisation": reading,
              "basis": "boundary_basis %s: %.1f x the handbook class prior (dashed lines are estimated by the evaluator; "
                       "ORCHID assumption), or the digitisation bound if larger" % (boundary_basis, f)}
    reading = rd["value"] if rd else None
    eff = max(x for x in (reading, prior) if x is not None)
    return {"unit": unit, "value_type": vt, "sigma_printed": None,
            "sigma_reading": rd,
            "sigma_method": {"class": "assessed_phase_diagram", "value": prior, "floor": prior, "floor_unit": unit,
                             "basis": "prior", "n_pairs": 0, "ci95": None, "entry": "assessed handbook class prior",
                             "prior_source": "orchid_assumed"},
            "sigma_composition": None, "sigma_scatter": None, "sigma_conversion": None,
            "sigma_boundary_margin": {"value": margin, "basis": "distance of the point to the nearest drawn boundary "
                                      "(drawing-based, low credence; describes placement, NOT added to "
                                      "sigma_effective)"},
            "sigma_effective": r6(eff)}


def assessed_budget(rec):
    cu, pm = rec["coordinate_uncertainty"], rec["positional_margin"]
    return {"calibration": BUDGET_CALIBRATION, "sigma_rule": BUDGET_RULE,
            "boundary_basis": rec.get("boundary_basis"),
            "quantities": {"temperature_K": _bq("K", "T_ASSESSED", cu["T_K"], HANDBOOK_PRIOR["T"], pm["isopleth_K"],
                                                rec.get("boundary_basis")),
                           "composition": _bq("mole fraction", "X_ASSESSED", cu["x"], HANDBOOK_PRIOR["x"],
                                              pm["isotherm_mole_fraction"], rec.get("boundary_basis"))}}


def espei_dataset(diag, pts, comp, elements, citations, with_tl, dense=False):
    comps = sorted(e.upper() for e in elements)
    if any(va for _, _, va, _ in pts):
        comps.append("VA")
    C = comp.upper()
    temps, values, phases = [], [], set()
    for rec, tn, _, tl in pts:
        phases.update(tn)
        temps.append(rec["T_K"])
        if len(tn) == 1:
            values.append([[tn[0], [C], [rec["x"]]]])
        elif tl:
            lab2t = dict(zip([e["label"] for e in rec["assemblage"]], tn))
            # a line-compound endpoint is written with composition null (ESPEI 0.9.1 / pycalphad 0.11.2 return NaN
            # for a line compound at a fixed composition); hyperplane at the point when no endpoint is left known
            reg = [[lab2t[e["phase"]], [C], [None if e.get("line_compound") else e["x"]]] for e in tl["endpoints"]]
            if all(v[2][0] is None for v in reg):
                reg.append(["__HYPERPLANE__", [C], [rec["x"]]])
            values.append(reg)
        else:
            values.append([[t, [C], [None]] for t in tn] + [["__HYPERPLANE__", [C], [rec["x"]]]])
    orig = "; ".join(c["citation"] + (" doi:%s" % c["doi"] if c.get("doi") else "") for c in citations
                     if c["role"] == "original_evaluation" and c.get("citation"))
    hb = [c for c in citations if c["role"] == "compilation"][0]
    ref = "%s%scompiled in %s, p. %s" % (orig, "; " if orig else "", hb["citation"], hb.get("page_pdf"))
    comment = ("evidence_class=assessed_phase_diagram (expert-drawn diagram, not experimental, not CALPHAD); ORCHID "
               "assessed-phase-assemblages v%s, diagram %s; one point per phase field (centroid or pole of inaccessibility); "
               "single-phase: the phase at the point composition; two-phase: phases present, compositions unknown (null) with "
               "__HYPERPLANE__ at the point's overall composition." % (METHOD_VERSION, diag))
    if dense:
        comment = ("evidence_class=assessed_phase_diagram (expert-drawn diagram, not experimental, not CALPHAD); ORCHID "
                   "assessed-phase-assemblages v%s, dense field grid (method %s), diagram %s; regular interior (T, x) grid "
                   "inside every named field, at least %.3f (frame-normalised) from every drawn line; single-phase: the phase at "
                   "the point composition; two-phase: phases present, compositions unknown (null) with __HYPERPLANE__ at the "
                   "point's overall composition. Points of one diagram are not independent: orchid.field_weight_share gives "
                   "each point's share of one observation per field." % (METHOD_VERSION, DENSE_VERSION, diag, MARGIN))
    if with_tl:
        comment += (" Variant with tie-lines: two-phase points with a tie_line block carry its endpoint compositions "
                    "(isotherm through the point intersected with the drawn boundaries); caveat: " + CAVEAT + ".")
    return {"components": comps, "phases": sorted(phases), "conditions": {"P": P_PA, "T": temps}, "output": "ZPF",
            "values": values, "reference": ref, "comment": comment,
            "tags": ["assessed_phase_diagram", "orchid"] + (["tie_line_variant"] if with_tl else [])
            + (["dense_field_grid"] if dense else []),
            "orchid": dict({"evidence_class": EVIDENCE_CLASS, "diagram": diag, "record_ids": [r["id"] for r, _, _, _ in pts],
                            "boundary_basis": [r.get("boundary_basis") for r, _, _, _ in pts]},
                           **({"point_method": "dense_field_grid", "dense_method_version": DENSE_VERSION,
                               "field": [r["field"] for r, _, _, _ in pts],
                               "field_weight_share": [r["dense_grid"]["field_weight_share"] for r, _, _, _ in pts]}
                              if dense else {})),
            "orchid_uncertainty_budget": budget_points([r for r, _, _, _ in pts])}


def budget_points(recs):
    """ESPEI metadata: per point (value order) the component values of the record's uncertainty budget."""
    def comp(q):
        return {k: (q.get(k) or {}).get("value") if isinstance(q.get(k), dict) else q.get(k)
                for k in ("sigma_printed", "sigma_reading", "sigma_method", "sigma_boundary_margin", "sigma_effective")}
    return {"sigma_rule": BUDGET_RULE, "calibration": BUDGET_CALIBRATION,
            "points": [{"record_id": r["id"], "boundary_basis": r.get("boundary_basis"),
                        "temperature_K": comp(r["uncertainty_budget"]["quantities"]["temperature_K"]),
                        "composition_mole_fraction": comp(r["uncertainty_budget"]["quantities"]["composition"])}
                       for r in recs]}


# ------------------------------------------------------------------------------------------------ outputs --
CSV_COLS = ["id", "system", "diagram", "field", "assemblage_label", "n_phases", "phase_labels", "canonical_names",
            "phase_ids", "T_K", "T_C", "x", "x_component", "P_Pa", "u_T_K", "u_x", "sigma_T_reading_K", "sigma_x_reading",
            "sigma_T_method_K", "sigma_x_method", "sigma_T_effective_K", "sigma_x_effective", "point_method",
            "margin_normalized",
            "margin_isotherm_x", "margin_isopleth_K", "boundary_basis", "basis", "chemical_system", "paper_id", "fold_lso", "fold_lpo", "confidence", "review_state", "flags",
            "tie_line_x_left", "tie_line_phase_left", "tie_line_x_right", "tie_line_phase_right", "tie_line_omitted",
            "sgte_n", "sgte_agrees", "sgte_disagrees", "sgte_not_comparable", "sgte_results",
            "original_citation", "original_doi", "handbook_page_pdf", "evidence_class"]
INT_COLS = {"fold_lso", "fold_lpo", "field", "n_phases", "sgte_n", "sgte_agrees", "sgte_disagrees", "sgte_not_comparable", "handbook_page_pdf"}
FLT_COLS = {"T_K", "T_C", "x", "P_Pa", "u_T_K", "u_x", "sigma_T_reading_K", "sigma_x_reading", "sigma_T_method_K",
            "sigma_x_method", "sigma_T_effective_K", "sigma_x_effective", "margin_normalized", "margin_isotherm_x", "margin_isopleth_K",
            "tie_line_x_left", "tie_line_x_right"}


def flat(rec):
    tl = rec["tie_line"] or {}
    e = tl.get("endpoints") or [{}, {}]
    sg = rec["sgte_crosscheck"]
    oc = [c for c in rec["citations"] if c["role"] == "original_evaluation"]
    hb = [c for c in rec["citations"] if c["role"] == "compilation"][0]
    bT = rec["uncertainty_budget"]["quantities"]["temperature_K"]
    bX = rec["uncertainty_budget"]["quantities"]["composition"]
    return {"id": rec["id"], "system": rec["system"], "diagram": rec["diagram"], "field": rec["field"],
            "assemblage_label": rec["assemblage_label"], "n_phases": rec["n_phases"],
            "phase_labels": ";".join(p["label"] for p in rec["assemblage"]),
            "canonical_names": ";".join(p["canonical_name"] for p in rec["assemblage"]),
            "phase_ids": ";".join(p["phase_id"] or "" for p in rec["assemblage"]),
            "T_K": rec["T_K"], "T_C": rec["T_C"], "x": rec["x"], "x_component": rec["x_component"], "P_Pa": rec["P_Pa"],
            "u_T_K": rec["coordinate_uncertainty"]["T_K"], "u_x": rec["coordinate_uncertainty"]["x"],
            "sigma_T_reading_K": (bT["sigma_reading"] or {}).get("value"),
            "sigma_x_reading": (bX["sigma_reading"] or {}).get("value"),
            "sigma_T_method_K": bT["sigma_method"]["value"], "sigma_x_method": bX["sigma_method"]["value"],
            "sigma_T_effective_K": bT["sigma_effective"], "sigma_x_effective": bX["sigma_effective"],
            "point_method": rec["point_method"], "margin_normalized": rec["positional_margin"]["normalized"],
            "margin_isotherm_x": rec["positional_margin"]["isotherm_mole_fraction"],
            "margin_isopleth_K": rec["positional_margin"]["isopleth_K"], "boundary_basis": rec.get("boundary_basis"),
            "basis": rec["basis"], "chemical_system": rec["chemical_system"], "paper_id": rec["paper_id"],
            "fold_lso": rec["fold_lso"], "fold_lpo": rec["fold_lpo"],
            "confidence": rec["confidence"], "review_state": rec["review_state"], "flags": ";".join(rec["flags"]),
            "tie_line_x_left": e[0].get("x"), "tie_line_phase_left": e[0].get("phase"),
            "tie_line_x_right": e[1].get("x"), "tie_line_phase_right": e[1].get("phase"),
            "tie_line_omitted": rec["tie_line_omitted"],
            "sgte_n": len(sg), "sgte_agrees": sum(s["result"] == "agrees" for s in sg),
            "sgte_disagrees": sum(s["result"] == "disagrees" for s in sg),
            "sgte_not_comparable": sum(s["result"] == "not_comparable" for s in sg),
            "sgte_results": ";".join("%s=%s" % (s["tdb"], s["result"]) for s in sg),
            "original_citation": " | ".join(c["citation"] or "" for c in oc),
            "original_doi": ";".join(c["doi"] for c in oc if c.get("doi")), "handbook_page_pdf": hb.get("page_pdf"),
            "evidence_class": rec["evidence_class"]}


def to_csv(rows, cols):
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols, lineterminator="\n")
    w.writeheader()
    for r in rows:
        w.writerow({k: ("" if r.get(k) is None else r.get(k)) for k in cols})
    return buf.getvalue()


def parquet_bytes(rows, cite):
    import pyarrow as pa
    import pyarrow.parquet as pq
    fields = [(c, pa.int32() if c in INT_COLS else pa.float64() if c in FLT_COLS else pa.string()) for c in CSV_COLS]
    meta = {"orchid.dataset": DATASET, "orchid.table": "assessed_phase_assemblages", "orchid.license": "CC-BY-4.0",
            "orchid.generator": "orchid-project.org/tools/build_assessed_assemblages.py", "orchid.doi": cite["doi"],
            "orchid.citation": cite["text"], "orchid.evidence_class": EVIDENCE_CLASS,
            "orchid.cite_as": json.dumps(cite, sort_keys=True, ensure_ascii=False)}
    data = [{c: r.get(c) for c in CSV_COLS} for r in rows]
    t = pa.Table.from_pylist(data, schema=pa.schema(fields).with_metadata(meta))
    sink = pa.BufferOutputStream()
    pq.write_table(t, sink, compression="zstd", compression_level=9, row_group_size=65536, use_dictionary=True,
                   data_page_version="2.0", write_statistics=True)
    return sink.getvalue().to_pybytes()


def parquet_content(b):
    import pyarrow as pa
    import pyarrow.parquet as pq
    t = pq.read_table(pa.BufferReader(b))
    h = hashlib.sha256()
    for r in t.to_pylist():
        h.update(json.dumps(r, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode())
        h.update(b"\n")
    kv = {k.decode(): v.decode() for k, v in (t.schema.metadata or {}).items() if k.startswith(b"orchid.")}
    return h.hexdigest(), sha256_bytes(json.dumps(kv, sort_keys=True).encode())


ESPEI_CHECK = r'''
import json, os, sys, glob
from espei.datasets import check_dataset, load_datasets, recursive_glob, DatasetError
out = {}
for d in sys.argv[1:]:
    files = sorted(glob.glob(os.path.join(d, "*.json")))
    bad = []
    for f in files:
        try:
            check_dataset(json.load(open(f)))
        except Exception as e:
            bad.append("%s: %s" % (os.path.basename(f), str(e)[:200]))
    try:
        db = load_datasets(sorted(recursive_glob(d, "*.json")))
        n = len(db.all())
    except Exception as e:
        n = -1; bad.append("load_datasets: %s" % str(e)[:300])
    import espei
    out[os.path.basename(d)] = {"files": len(files), "valid": len(files) - len(bad), "loaded": n, "errors": bad[:20],
                                "espei": espei.__version__}
print(json.dumps(out))
'''


def validate_espei(py, dirs):
    if not py or not os.path.exists(py):
        raise SystemExit("ESPEI interpreter not found (%s): the ESPEI exports must be validated" % py)
    r = subprocess.run([py, "-I", "-c", ESPEI_CHECK] + dirs, capture_output=True, text=True)
    if r.returncode:
        raise SystemExit("ESPEI validation crashed: " + r.stderr[-2000:])
    res = json.loads(r.stdout.strip().splitlines()[-1])
    for d, v in res.items():
        if v["errors"] or v["valid"] != v["files"] or v["loaded"] != v["files"]:
            raise SystemExit("ESPEI validation failed in %s: %s" % (d, json.dumps(v)[:3000]))
    return res


def write_outputs(a, records, exclusions, espei_excl, espei_default, espei_tl, cnt, method):
    """{relative path: bytes}"""
    cite = orchid_cite.cite_as(DATASET, a.data)
    out = {}
    for r in records:
        # machine-learning folds (tools/ml_folds.py); the source is the diagram's original evaluation (its DOI), else
        # the handbook diagram itself
        oc = [c["doi"] for c in r["citations"] if c["role"] == "original_evaluation" and c.get("doi")]
        r.update(ml_folds.fields(r["system"], ("doi:" + oc[0].lower()) if oc else "handbook:" + r["diagram"]))
        guard(r, "records[%s]" % r["id"])
    for e in exclusions + espei_excl:
        guard(e, "exclusions")
    # gzip (mtime 0: deterministic bytes); the uncompressed JSONL exceeds the 45 MB per-file limit
    out["records.jsonl.gz"] = gzip.compress("".join(dumps(r) + "\n" for r in records).encode(), compresslevel=9, mtime=0)
    rows = [flat(r) for r in records]
    csv_text = to_csv(rows, CSV_COLS)
    guard_text(csv_text, "records.csv")
    out["records.csv"] = csv_text.encode()
    pqb = parquet_bytes(rows, cite)
    out["records.parquet"] = pqb
    ex_cols = ["diagram", "field", "reason", "detail", "assemblage"]
    out["exclusions.jsonl"] = "".join(dumps(e) + "\n" for e in exclusions).encode()
    out["exclusions.csv"] = to_csv([dict(e, assemblage="+".join(e["assemblage"]) if e.get("assemblage") else None)
                                    for e in exclusions], ex_cols).encode()
    out["espei_exclusions.jsonl"] = "".join(dumps(e) + "\n" for e in espei_excl).encode()
    out["sgte_not_comparable.csv"] = to_csv(
        [{"id": r["id"], "system": r["system"], "diagram": r["diagram"], "field": r["field"],
          "assemblage_label": r["assemblage_label"], "tdb": x["tdb"], "calculated": x.get("calculated"), "reason": x["reason"],
          "unmapped_labels": ";".join(p["label"] for p in r["assemblage"] if not p["phase_id"])}
         for r in records for x in r["sgte_crosscheck"] if x["result"] == "not_comparable"],
        ["id", "system", "diagram", "field", "assemblage_label", "tdb", "calculated", "reason", "unmapped_labels"]).encode()
    for sub, dsets in (("espei-zpf", espei_default), ("espei-zpf-with-tie-lines", espei_tl)):
        for diag, ds in sorted(dsets.items()):
            guard(ds, "%s/%s" % (sub, diag))
            out["%s/%s.json" % (sub, diag)] = (json.dumps(ds, ensure_ascii=False, indent=1) + "\n").encode()
    dd = getattr(a, "_dense", None)
    dense_recs = dd["records"] if dd else []
    if dense_recs:
        for r in dense_recs:
            oc = [c["doi"] for c in r["citations"] if c["role"] == "original_evaluation" and c.get("doi")]
            r.update(ml_folds.fields(r["system"], ("doi:" + oc[0].lower()) if oc else "handbook:" + r["diagram"]))
            guard(r, "dense[%s]" % r["id"])
        out[DENSE_DIR + "/records.jsonl.gz"] = gzip.compress("".join(dumps(r) + "\n" for r in dense_recs).encode(),
                                                             compresslevel=9, mtime=0)
        drows = [flat(r) for r in dense_recs]
        dcsv = to_csv(drows, CSV_COLS)
        guard_text(dcsv, DENSE_DIR + "/records.csv")
        out[DENSE_DIR + "/records.csv"] = dcsv.encode()
        out[DENSE_DIR + "/records.parquet"] = parquet_bytes(drows, cite)
        out[DENSE_DIR + "/espei_exclusions.jsonl"] = "".join(dumps(e) + "\n" for e in dd["espei_excl"]).encode()
        for diag, ds in sorted(dd["espei"].items()):
            guard(ds, "%s/%s" % (DENSE_ESPEI_DIR, diag))
            out["%s/%s.json" % (DENSE_ESPEI_DIR, diag)] = (json.dumps(ds, ensure_ascii=False, indent=1) + "\n").encode()
    # ---- summary
    sg = collections.Counter()
    sg_any = collections.Counter()
    for r in records:
        res = {s["result"] for s in r["sgte_crosscheck"]}
        for s in r["sgte_crosscheck"]:
            sg[s["result"]] += 1
            if s["result"] == "not_comparable":
                sg["nc:" + s["reason"]] += 1
        sg_any["no SGTE assessment" if not r["sgte_crosscheck"] else "agrees with >=1" if "agrees" in res
               else "disagrees with all comparable" if "disagrees" in res else "not comparable"] += 1
    comparable = sg["agrees"] + sg["disagrees"]
    # metric definitions (James Saal review 2026-10-24): agreement is reported over COMPARABLE comparisons only, and
    # per point as the share of points with >= 1 comparable assessment that agree with at least one of them;
    # non-comparable comparisons are counted separately, with their reasons (sgte_not_comparable.csv)
    pts_cmp = [r for r in records if any(x["result"] in ("agrees", "disagrees") for x in r["sgte_crosscheck"])]
    pts_any = sum(1 for r in pts_cmp if any(x["result"] == "agrees" for x in r["sgte_crosscheck"]))
    pts_all = sum(1 for r in pts_cmp if all(x["result"] != "disagrees" for x in r["sgte_crosscheck"]))
    by_sys = collections.defaultdict(lambda: [0, 0])
    for r in records:
        for x in r["sgte_crosscheck"]:
            if x["result"] in ("agrees", "disagrees"):
                by_sys[r["system"]][x["result"] == "disagrees"] += 1
    rates = sorted(a_ / (a_ + d_) for a_, d_ in by_sys.values() if a_ + d_ >= 5)
    med = (rates[len(rates) // 2] if len(rates) % 2 else (rates[len(rates) // 2 - 1] + rates[len(rates) // 2]) / 2) if rates else None
    two = [r for r in records if r["n_phases"] == 2]
    summ = {
        "dataset": DATASET, "schema": SCHEMA_VERSION, "method": method,
        "diagrams_with_points": len({r["diagram"] for r in records}), "diagrams_seen": cnt["diagrams_seen"],
        "diagrams_not_annotated": cnt["diagram_not_annotated"],
        "fields_seen": len(records) + sum(1 for e in exclusions if e.get("field") is not None),
        "points": len(records),
        "point_method": {"centroid": cnt["method:centroid"], "pole_of_inaccessibility": cnt["method:pole_of_inaccessibility"]},
        "points_by_n_phases": dict(sorted(collections.Counter(r["n_phases"] for r in records).items())),
        "points_by_basis": dict(collections.Counter(r["basis"] for r in records).most_common()),
        "points_by_confidence": dict(collections.Counter(r["confidence"] for r in records).most_common()),
        "points_by_review_state": dict(collections.Counter(r["review_state"] for r in records).most_common()),
        "excluded_fields_by_reason": dict(collections.Counter(e["reason"] for e in exclusions if e.get("field") is not None).most_common()),
        "phase_labels_mapped_to_registry": {
            "mapped": sum(1 for r in records for p in r["assemblage"] if p["phase_id"]),
            "unmapped": sum(1 for r in records for p in r["assemblage"] if not p["phase_id"])},
        "two_phase_tie_lines": {"with_tie_line": sum(1 for r in two if r["tie_line"]),
                                "without": sum(1 for r in two if not r["tie_line"]),
                                "omitted_by_reason": dict(collections.Counter(r["tie_line_omitted"] for r in two if not r["tie_line"]).most_common()),
                                "line_compound_endpoints": sum(1 for r in two if r["tie_line"] for e in r["tie_line"]["endpoints"] if e.get("line_compound"))},
        "sgte_crosscheck": {"comparisons": sg["agrees"] + sg["disagrees"] + sg["not_comparable"], "agrees": sg["agrees"],
                            "disagrees": sg["disagrees"], "not_comparable": sg["not_comparable"],
                            "agreement_rate_of_comparable": round(sg["agrees"] / comparable, 4) if comparable else None,
                            "not_comparable_by_reason": {k[3:]: v for k, v in sorted(sg.items()) if k.startswith("nc:")},
                            "points": dict(sg_any.most_common()),
                            "metric_definitions": {
                                "comparison": "one point x one calculated SGTE diagram of its system (nearest grid cell)",
                                "comparable": "both phase sets can be expressed as registry phases and the grid cell has an "
                                              "equilibrium: result agrees or disagrees",
                                "agreement_rate_of_comparable": "agrees / (agrees + disagrees), over all comparable "
                                                                "comparisons (aggregate, not a per-system mean)",
                                "points_agreeing_with_at_least_one_rate": "points that agree with >= 1 SGTE assessment / "
                                                                          "points with >= 1 comparable comparison",
                                "median_system_rate": "median over systems with >= 5 comparable comparisons of the "
                                                      "system's agreement rate",
                                "not_comparable": "counted separately, never in a rate; reasons in not_comparable_by_reason "
                                                  "and, per comparison, in sgte_not_comparable.csv"},
                            "comparable": comparable,
                            "points_with_comparable_assessment": len(pts_cmp),
                            "points_agreeing_with_at_least_one": pts_any,
                            "points_agreeing_with_at_least_one_rate": round(pts_any / len(pts_cmp), 4) if pts_cmp else None,
                            "points_agreeing_with_all_comparable": pts_all,
                            "points_disagreeing_with_all_comparable": len(pts_cmp) - pts_any,
                            "systems_with_5_comparable": len(rates),
                            "median_system_rate": round(med, 4) if med is not None else None,
                            "mean_system_rate": round(sum(rates) / len(rates), 4) if rates else None},
        "points_by_boundary_basis": {str(k): v for k, v in collections.Counter(
            r.get("boundary_basis") or "none" for r in records).most_common()},
        "diagrams_with_dashed_lines": cnt["diagrams_with_dashed_lines"],
        "ml_folds": dict(ml_folds.balance(records), method=ml_folds.DOC),
        "espei": {"datasets_default": len(espei_default), "datasets_with_tie_lines": len(espei_tl),
                  "points_in_espei": len(records) - len(espei_excl), "points_excluded_from_espei": len(espei_excl)},
        "doi_corrections": dict(DOI_CORR.provenance(), records_with_dropped_original_doi=sum(
            1 for r in records if any(c.get("doi_status") for c in r["citations"])),
            systems=sorted({r["system"] for r in records if any(c.get("doi_status") for c in r["citations"])})),
        "dense_field_grid": None if not dense_recs else {
            "method_version": DENSE_VERSION, "diagrams": sorted({r["diagram"] for r in dense_recs}),
            "points": len(dense_recs), "fields": len({(r["diagram"], r["field"]) for r in dense_recs}),
            "points_by_n_phases": dict(sorted(collections.Counter(r["n_phases"] for r in dense_recs).items())),
            "points_in_espei": len(dense_recs) - len(dd["espei_excl"]), "points_excluded_from_espei": len(dd["espei_excl"]),
            "default_enabled_for": sorted(DENSE_VALIDATED), "files": DENSE_DIR + "/ and " + DENSE_ESPEI_DIR + "/",
            "label_overrides": {d: {l: v[0] for l, v in m.items()} for d, m in DENSE_LABEL_OVERRIDES.items()
                                if d in {r["diagram"] for r in dense_recs}}},
        "cite_as": cite,
    }
    summ["files"] = {}
    out["summary.json"] = b""
    out["CITATION.cff"] = orchid_cite.cff(DATASET, a.data, " and, for each diagram, its original evaluation and the "
                                          "Okamoto handbook as the compilation").encode()
    return out, summ, cite, pqb


def finalize(a, out, summ, cite, pqb, espei_res):
    summ["espei"]["validation"] = espei_res
    summ["parquet"] = {"table": "records.parquet", "rows": summ["points"], "csv_sha256": sha256_bytes(out["records.csv"]),
                       "content_sha256": parquet_content(pqb)[0],
                       "note": "content_sha256 = SHA-256 of the canonical JSON rows (sorted keys); independent of the Parquet encoder"}
    top = collections.Counter()
    for rel, b in out.items():
        k = rel.split("/")[0] if "/" in rel else rel
        top[k] += len(b)
    summ["files"] = {"bytes_by_file_or_dir": dict(sorted(top.items()))}
    guard({k: v for k, v in summ.items() if k != "cite_as"}, "summary")
    out["summary.json"] = (json.dumps(summ, ensure_ascii=False, indent=1) + "\n").encode()
    out["parquet.json"] = (json.dumps(summ["parquet"], indent=1) + "\n").encode()
    dr_, ds_ = dense_docs(summ)
    out["README.md"] = (readme(summ, cite) + dr_).encode()
    out["SCHEMA.md"] = (schema_md(summ) + ds_).encode()
    big = [rel for rel, b in out.items() if len(b) > MAX_FILE_BYTES]
    if big:
        raise SystemExit("files over 45 MB: %s" % ", ".join(big))
    return out


def fmt(n):
    return format(n, ",")


def dense_docs(s):
    """README and SCHEMA sections of the dense field grid ('' when the release has none)."""
    d = s.get("dense_field_grid")
    if not d:
        return "", ""
    ov = "; ".join("%s: %s" % (dg, ", ".join("%s -> %s" % kv for kv in m.items())) for dg, m in d["label_overrides"].items())
    readme_ = """

## Dense field grid (derived points of the assessed-diagram class)

Besides the primary point of each field, a regular interior grid of (T, x) points inside every named field is released
for the diagrams where it was validated (%s; method %s, %s points in %s fields, %s of them usable in ESPEI
datasets). The grid is the one used by the PhaseForge+ Cu-Zr example (`cuzr-phaseforgeplus` in the ORCHID examples):
spacing min(50 K, field height / 4) by min(0.05, field width / 4) anchored at the lower-left corner of the field plus half
a step, kept when the node lies at least %.3f (frame-normalised) from every other field or drawn line, at most %d points
per field (even subsample), the single-phase liquid only within %d K above its lower boundary. Like the primary points, they
are derived coordinates and phase assemblages only (never page images, crops or handbook text), with the same record schema,
uncertainty budget (handbook-class prior), basis, review state and citations (original evaluation first, then the
handbook); they are distinguished by `point_method` = `dense_field_grid`, the flags `dense_field_grid` and `interior_grid`
and the `dense_grid` block (`primary_record_id`, `n_points_in_field`, `field_weight_share`, margin used). The points of one
diagram are not independent observations: `field_weight_share` = 1 / (points in the field), so that each field counts as
one observation when a fit weights by it. Files: `%s/` (records.jsonl.gz, records.csv, records.parquet,
espei_exclusions.jsonl) and `%s/` (ESPEI ZPF datasets, dense points only; the bundle carries them as
`assessed/<SYSTEM>/<TDB>.<DIAGRAM>.dense.json`, tag `dense_field_grid`). The counts above (`points`, `espei`) are the
primary points only. Phase labels the primary records leave unmapped but the validated grid resolves (recorded in each
record's `mapping_basis`): %s. Enable for other diagrams only after validation (`--dense all` in the builder).
""" % (", ".join(d["diagrams"]), d["method_version"], fmt(d["points"]), fmt(d["fields"]), fmt(d["points_in_espei"]), MARGIN,
       DENSE_MAX_PER_FIELD, int(DENSE_LIQUID_BAND_K), DENSE_DIR, DENSE_ESPEI_DIR, ov)
    schema_ = """

## Dense field grid records (`%s/`)

Same fields as above. Differences: `point_method` = `dense_field_grid`; `centroid` null; `id` = `aspa-` + hash of
(system, diagram, field, grid index >= 1) (the primary record of a field has index 0); `tie_line` null with
`tie_line_omitted` = `dense_field_grid_point`; `flags` add `dense_field_grid` and `interior_grid`; `phase_id` of a label
the primary record leaves unmapped may be filled by a documented override (`mapping_basis` starts with `dense_field_grid
override`); extra block `dense_grid` = `{method_version, primary_record_id, index_in_field, n_points_in_field,
field_weight_share (1 / n_points_in_field), interior_margin_normalized_min, clearance_normalized, spacing, weight_note}`.
`x` is rounded to 1e-6 (raster resolution 5e-4). The ESPEI datasets (`%s/`) carry `orchid.point_method`,
`orchid.field` and `orchid.field_weight_share` per point and the `dense_field_grid` tag.
""" % (DENSE_DIR, DENSE_ESPEI_DIR)
    return readme_, schema_


def readme(s, cite):
    pm, sg, tl, es = s["point_method"], s["sgte_crosscheck"], s["two_phase_tie_lines"], s["espei"]
    ex = "\n".join("| %s | %s |" % (k, fmt(v)) for k, v in s["excluded_fields_by_reason"].items())
    bas = "\n".join("| %s | %s |" % (k, fmt(v)) for k, v in s["points_by_basis"].items())
    return """# ORCHID assessed phase assemblages (binary systems)

Status: **Version 1, machine-extracted.** Release %s. Licence CC BY 4.0 (derived data and citations).

One record per named phase field of the expert-drawn (assessed) binary phase diagrams compiled in H. Okamoto,
*Desk Handbook: Phase Diagrams for Binary Alloys*, 2nd ed. (ASM International, 2010): the thermodynamic coordinate
(T, x, P) of a point deep inside the field and the phase assemblage present there. Evidence class
`assessed_phase_diagram`: an expert drew the diagram from experimental evidence without a CALPHAD model; these are
neither raw measurements nor model values, and they are kept separable at fitting time.

## Counts (this build)

| | |
|---|---|
| diagrams with points | %s |
| fields seen | %s |
| points | %s (centroid %s, pole of inaccessibility %s) |
| points by number of phases | %s |
| two-phase points with / without a tie_line block | %s / %s |
| SGTE cross-check (point x assessment) | %s comparisons: %s comparable (%s agree, %s disagree; agreement %s of comparable), %s not comparable (counted separately, reasons in `sgte_not_comparable.csv`) |
| SGTE cross-check (points) | %s points with a comparable assessment; %s (%s) agree with at least one SGTE assessment; median system agreement %s (%s systems with >= 5 comparable comparisons) |
| boundary basis (lines bounding the point's field along its isotherm and isopleth) | %s |
| ESPEI datasets (one per diagram) | %s default, %s with tie-lines; %s points exported, %s excluded (phase without TDB name) |

Points by basis:

| basis | points |
|---|---|
%s

Fields without a point, by reason:

| reason | fields |
|---|---|
%s

## Why one point deep inside each field

The PI's guidance (2026-10-23): the positions of drawn boundaries and tie-lines deserve little credence; the
thermodynamic coordinates of the assemblage are what matters. Each field therefore gives exactly one point: the area
centroid of the field (in the diagram's T and composition axes) when it lies inside the field at least %.3f
(frame-normalised) from every drawn line, otherwise the pole of inaccessibility (the interior point farthest from all
lines). Fields too thin for that margin have no point (listed in `exclusions.jsonl`). The claim "this assemblage is
stable at (T, x, 1 atm)" then does not depend on where exactly the expert drew the boundaries. Coordinate uncertainty
comes from the axis calibration (printed tick labels vs the fitted axis) plus half a raster cell.

Each record carries an `uncertainty_budget` with the components kept separate, in the same structure as the ORCHID
literature release (PI decision 2026-10-24): digitisation (`sigma_reading`, the coordinate uncertainty above), a
handbook-class prior (`sigma_method`: 10 K and 0.02 mole fraction for an expert-drawn assessed diagram; an ORCHID
assumption, the same as the release's model / assessment priors, not fitted) and the distance to the nearest drawn
boundary (`sigma_boundary_margin`, reported but not combined, since boundary positions carry low credence);
`sigma_effective` = max(reading, handbook prior). The ESPEI datasets carry the per-point components and the
`boundary_basis` in `orchid_uncertainty_budget` (ZPF data are not weighted by a sigma in ESPEI).

## Solid and dashed (estimated) boundaries

Handbook diagrams draw boundaries that the evaluator could not fix from measurements as dashed lines. The tracer joins
the dashes of such a line into a chain and uses it as a field boundary like a solid line (dashed Curie-temperature
lines, cooling curves and the dashed second assessment of an overlay are not boundaries and are dropped). Each point
records `boundary_basis` for the lines that bound its field along its isotherm and isopleth, i.e. the lines that set
its positional margin: `printed_solid`, `dashed_estimated` or `mixed` (null when only frame edges bound it), with the
counts in `boundary_lines`. The dash chains are re-derived from the traced lines with the tracer's own rule
(`tools/dashed_boundaries.py`; the build fails if the number of dash joins differs from the tracer's count in any
diagram). For points bounded by dashed lines the reading component of the budget is enlarged (ORCHID assumption, not
fitted): 2x the handbook-class prior (20 K, 0.04) for `dashed_estimated` and 1.5x (15 K, 0.03) for `mixed`.

## SGTE cross-check: metric definitions

A comparison is one point against one SGTE assessment of its system (the phase set at the nearest cell of the
calculated 10 K x 0.01 grid). It is comparable when every phase of the point maps to a registry phase and the cell
has an equilibrium; the result is then agrees or disagrees. Rates are computed over comparable comparisons only:
`agreement_rate_of_comparable` = agrees / (agrees + disagrees) over all comparable comparisons (an aggregate, not a
per-system mean); `points_agreeing_with_at_least_one_rate` = share of the points with at least one comparable
comparison that agree with at least one SGTE assessment of their system (systems have up to 16 assessments, and an
older assessment that differs still counts); `median_system_rate` over systems with at least 5 comparable
comparisons. Non-comparable comparisons are never part of a rate; they are counted by reason in `summary.json` and
listed one per row in `sgte_not_comparable.csv`.

Two-phase points carry the two phases with **unknown compositions** (`phase_compositions: null`; ESPEI null vertices
with a `__HYPERPLANE__` vertex at the point's overall composition). A separate `tie_line` block gives the isotherm
through the point intersected with the drawn boundaries, with the caveat that it is drawing-based and of low credence
(weak/soft constraint); it is left out with a reason when the isotherm does not meet two clean, unambiguous boundaries.

## Files

`records.jsonl.gz` (full records, gzip JSON Lines), `records.csv` / `records.parquet` (flat), `exclusions.jsonl` / `.csv` (fields and
diagrams without points, with reasons), `espei-zpf/<diagram>.json` (ESPEI ZPF datasets, null endpoints),
`espei-zpf-with-tie-lines/<diagram>.json` (variant using tie_line endpoints where present), `espei_exclusions.jsonl`,
`sgte_not_comparable.csv` (each non-comparable SGTE comparison with its reason),
`summary.json`, `parquet.json`, `SCHEMA.md`, `CITATION.cff`. ESPEI datasets were validated with
`espei.datasets.check_dataset` and `load_datasets` (ESPEI %s).

## Not included (licence)

No page images, crops, renders, paths to page images, or handbook text beyond phase names and citations. The build
fails if any output key or string looks like such content (forbidden keys, image/PDF paths, base64, long free text).
Traced boundary geometry, reader/adjudicator reasons and review pages are not part of this dataset.

## Citation

Cite this dataset, the umbrella ORCHID data record, and, for every diagram you use, its **original evaluation**
(`citations[role=original_evaluation]`, DOI where matched) followed by the Okamoto handbook as the compilation.

%s

Part of: %s
""" % (orchid_cite.RELEASE_TAG, fmt(s["diagrams_with_points"]), fmt(s["fields_seen"]), fmt(s["points"]),
       fmt(pm["centroid"]), fmt(pm["pole_of_inaccessibility"]),
       ", ".join("%s-phase %s" % (k, fmt(v)) for k, v in s["points_by_n_phases"].items()),
       fmt(tl["with_tie_line"]), fmt(tl["without"]), fmt(sg["comparisons"]), fmt(sg["comparable"]), fmt(sg["agrees"]),
       fmt(sg["disagrees"]), sg["agreement_rate_of_comparable"], fmt(sg["not_comparable"]),
       fmt(sg["points_with_comparable_assessment"]), fmt(sg["points_agreeing_with_at_least_one"]),
       sg["points_agreeing_with_at_least_one_rate"], sg["median_system_rate"], fmt(sg["systems_with_5_comparable"]),
       ", ".join("%s %s" % (k, fmt(v)) for k, v in s["points_by_boundary_basis"].items()), fmt(es["datasets_default"]),
       fmt(es["datasets_with_tie_lines"]), fmt(es["points_in_espei"]), fmt(es["points_excluded_from_espei"]), bas, ex, MARGIN,
       next(iter(s["espei"]["validation"].values()))["espei"], cite["text"], cite["part_of"]["text"])


def schema_md(s):
    return """# Assessed phase assemblages: schema (%s)

Builder: `tools/build_assessed_assemblages.py` (orchid-project.org; deterministic, `--check`). Source annotation run:
`%s`. Evidence class of every record: `assessed_phase_diagram`.

## Record (`records.jsonl.gz`, gzip-compressed JSON Lines, one object per line)

| field | meaning |
|---|---|
| `id` | `aspa-` + first 16 hex of SHA-256 of `system|diagram|field|0` (stable across rebuilds) |
| `record_type`, `evidence_class` | `assessed_phase_assemblage`, `assessed_phase_diagram` (never `experimental`, never a model value) |
| `system`, `elements` | binary system (upper-case elements, sorted, `-`-joined; the phase-registry key) |
| `diagram`, `field` | diagram id `<SYSTEM>-p<PDF page>` and field id of the annotation run |
| `assemblage[]` | phases present: `label` (handbook label as printed, with the annotation conventions: `L#2` second liquid, `δ#2` second composition of one phase, `G` gas), `canonical_name` (registry preferred name when mapped, else the label), `phase_id` (phase-registry id or `null`; never invented), `mapping_basis`, and phase metadata from the diagram's phase table (`pearson_symbol`, `strukturbericht`, `prototype`, `space_group_printed`, `composition_range_at_pct`) |
| `assemblage_label`, `n_phases` | `+`-joined labels; 1 or 2 (three-phase invariants are other records) |
| `T_K`, `T_C`, `temperature_axis` | temperature of the point; the diagram axis is in °C |
| `x`, `x_component`, `composition_axis`, `composition_basis` | mole fraction of `x_component`; every diagram here is drawn in at.%% (no wt.%% conversion) |
| `P_Pa`, `pressure_basis` | 101325; 1 atm assumed |
| `coordinate_uncertainty` | `{T_K, x}` bounds = max axis-calibration residual + half a raster cell; also the residuals and number of calibration pairs per axis |
| `uncertainty_budget` | `{calibration, sigma_rule, quantities{temperature_K, composition}}`, the same separate-component structure as the literature release (PI decision 2026-10-24). Per quantity (`unit` K / mole fraction, `value_type` `T_ASSESSED` / `X_ASSESSED`): `sigma_printed` null (the handbook prints no uncertainty); `sigma_reading` = digitisation, i.e. the `coordinate_uncertainty` bound (max axis-calibration residual + half a raster cell), enlarged for points bounded by dashed (estimated) lines to 2x the handbook prior (`dashed_estimated`: 20 K, 0.04) or 1.5x (`mixed`: 15 K, 0.03), ORCHID assumption, with the digitisation value kept as `digitisation`; `boundary_basis` repeated at budget level; `sigma_method` = handbook-class prior `{class: assessed_phase_diagram, value, floor, basis: prior, n_pairs: 0, entry, prior_source: orchid_assumed}`, 10 K and 0.02 mole fraction: ORCHID assumption, not fitted: an expert-drawn assessed diagram reproduces its experimental evidence to about 10 K and 0.02 mole fraction, the priors the literature release uses for model / assessment values; `sigma_composition`, `sigma_scatter`, `sigma_conversion` null (do not apply); `sigma_boundary_margin` = `positional_margin` isopleth_K / isotherm_mole_fraction (drawing-based, reported, **not** combined); `sigma_effective` = max(sigma_reading, sigma_method). Rule: `max(sqrt(max(printed, reading, method)^2 + composition^2 + conversion^2), scatter); boundary margin reported, not combined`. Flat columns: `sigma_T_reading_K`, `sigma_x_reading`, `sigma_T_method_K`, `sigma_x_method`, `sigma_T_effective_K`, `sigma_x_effective`. Tie-line endpoints keep their own `x_uncertainty` (drawing-based, not part of the budget) |
| `point_method` | `centroid` or `pole_of_inaccessibility`; for the fallback `centroid` gives the centroid, whether it was inside, and its clearance |
| `positional_margin` | `{normalized, isotherm_mole_fraction, isopleth_K}`: distance to the nearest drawn line or frame edge (frame-normalised; along the isotherm; along the isopleth). Drawing-based, low credence: describes placement only, not an uncertainty of any phase composition |
| `chemical_system`, `paper_id`, `fold_lso`, `fold_lpo` | machine-learning folds (`tools/ml_folds.py`): leave-chemical-system-out and leave-paper-out, int(sha256(key), 16) mod 5; here `paper_id` is the original evaluation (`doi:<doi>`), else `handbook:<diagram>`; subsystems are systems of their own |
| `boundary_basis` | `printed_solid` / `dashed_estimated` / `mixed` / null: the traced boundary lines that bound the point's field along its isotherm and isopleth (the lines that set `positional_margin`) are all solid, all dashed (estimated by the evaluator) or both; null when only frame edges bound it. Dash chains re-derived with the tracer's rule (`tools/dashed_boundaries.py`) |
| `boundary_lines` | `{solid, dashed, diagram_has_dashed_lines, note}`: how many of those crossings lie on solid and on dashed lines, and whether the diagram has any dashed boundary |
| `phase_compositions` | always `null`: the record does not claim phase compositions |
| `tie_line` | two-phase points only, or `null`: `{T_K, endpoints[{phase, x, x_uncertainty, side, kind, line_compound?, x_drawn?}], method, caveat}`. Isotherm through the point intersected with the field's drawn boundaries; each endpoint belongs to the phase whose single-phase field lies on that side, a line compound (composition from the printed phase table, `x_drawn` the drawn line), or a pure-element frame edge. Caveat: drawing-based, weak/soft constraint |
| `tie_line_omitted` | reason when a two-phase point has no tie_line: `neighbour_field_unnamed`, `side_assignment_ambiguous`, `boundary_not_a_single_clean_line`, `isotherm_reaches_cut_off_frame_edge`, `frame_edge_phase_ambiguous`, `endpoints_do_not_cover_the_two_phases`, `point_outside_endpoints` |
| `basis`, `basis_description` | how the field was named: `deterministic_printed_label`, `deterministic_phase_rules`, `deterministic_adjacency`, `reader_consensus`, `adjudicator`, `dft_hull_corroborated` (adjudicator + Materials Project/OQMD hulls; not handbook evidence), `expert_reviewed` |
| `confidence` | `high` / `medium` / `low` (blind-audit defaults per basis; hull rule: both hulls agree high, disagree low; `local_conflict` low; `carried_against_rules` capped at medium) |
| `review_state` | `expert decision (tier 5)`, `expert-reviewed (tier 5)`, or `not reviewed by the expert` |
| `name_history`, `expert_review` | only where the expert corrected the machine name (tools/assessed_expert_corrections.json): the superseded machine name with its basis, and the reviewer, review time and source of the decision |
| `flags` | `local_conflict`, `carried_against_rules`, `phase_not_in_table`, `hull_support:<both|disagree|mp_only|oqmd_only|none>`, `partly_drawn_dft_corroborated` |
| `sgte_crosscheck[]` | per calculated SGTE diagram of the system: `{tdb, tdb_sha256, cell{T_K, x}, calculated (TDB phase set at the nearest grid cell), result: agrees / disagrees / not_comparable, reason}`. Phase sets compared as registry ids. Never changes the assessed record |
| `citations[]` | original evaluation(s) first (`citation`, `doi`, `doi_matched_by`, `kind`; `doi` null with `doi_status: dropped_mismatch` and `doi_correction` when the DOI matched to the printed source line describes a different work, `tools/doi_corrections.json`), then the handbook (`role: compilation`, `page_pdf`, `locator`) |
| `method` | builder, version, annotation run, tool hashes |

## Point placement

Exactly one point per named field. The field is the set of raster cells of the traced field (drawn lines are not part
of any field). Distances use frame-normalised axes (each axis scaled to the drawn frame). The area centroid (mean T
and x of the field's cells) is used when it lies inside the field with clearance >= %.3f; otherwise the pole of
inaccessibility (the cell of maximum clearance; first in row-major order on ties). If even that is < %.3f the field
gets no point (`too_thin_for_margin`). Points are never placed on invariant lines or line compounds (both are lines,
not fields). Rationale: an interior point makes the assemblage claim independent of the exact boundary positions,
which in an assessed diagram carry little credence.

## Exclusions (`exclusions.jsonl`)

`{diagram, field, reason, detail, assemblage?}`. Reasons: `tracing_fault`, `open`, `awaiting_expert` (tier-5
pending), `undetermined_in_source` (source_undetermined, e.g. printed '?', drawing ends), `cannot_tell`,
`partly_drawn_truncated_edge` (the field reaches a cut-off composition edge of the drawing; included instead when
DFT-hull corroborated), `too_thin_for_margin`, `diagram_tracing_suspect` (whole diagram), `more_than_two_phases` (a binary field named with three or more phases), `axis_component_not_in_system` (the composition axis was read as an element outside the system; whole diagram), `diagram_not_annotated`
(field null: calibration failed, overlay of several assessments, no axis title).

## ESPEI (`espei-zpf/`, `espei-zpf-with-tie-lines/`)

One ZPF dataset per diagram: `components` (upper case, `VA` when a model of an included phase has vacancies),
`phases` (TDB phase names: the registry preferred name, i.e. the SGTE name used by most assessments of the system),
`conditions {P: 101325, T: [K...]}`, `values` aligned with `T`. Single-phase point: `[[PHASE, [X], [x]]]`. Two-phase
point (default): `[[P1, [X], [null]], [P2, [X], [null]], ["__HYPERPLANE__", [X], [x]]]` -- phases present,
compositions unknown, the target hyperplane at the point's overall composition. The tie-line variant uses
`[[P1, [X], [x1]], [P2, [X], [x2]]]` where a `tie_line` block exists, else the default form; a line-compound endpoint
is written with composition null (ESPEI 0.9.1 / pycalphad 0.11.2 return NaN for a line compound at a fixed composition),
and where both endpoints are line compounds the hyperplane is set at the point. `tags`: `assessed_phase_diagram`,
`orchid` (+ `tie_line_variant`). `reference`: original evaluation(s) and the handbook; `comment`: evidence class, point
method and (variant) the caveat; `orchid`: record ids; `orchid_uncertainty_budget`: `{sigma_rule, calibration,
points[{record_id, temperature_K, composition_mole_fraction}]}` with the budget component values of each point (value
order); ESPEI weighs ZPF by its driving-force tolerance, so the budget is information, not a weight. Points whose phases do not all map to a registry phase with a TDB
name are left out (`espei_exclusions.jsonl`) but stay in `records.*`. The ESPEI bundle (`espei/assessed/`) carries the
same points in the phase names of each assessment (TDB) that has ESPEI ZPF datasets, for use together with them.

## CSV / Parquet columns

%s. Lists are `;`-joined. `sgte_results`: `tdb=result` per assessment.
""" % (SCHEMA_VERSION, s["method"]["annotation_run"], MARGIN, MARGIN, ", ".join("`%s`" % c for c in CSV_COLS))


def zenodo_metadata(s, data):
    import project_credits
    e = orchid_cite.entry(DATASET, data)
    u = orchid_cite.entry(orchid_cite.UMBRELLA, data)
    desc = ("<p><b>%s</b>: one record per named phase field of the expert-drawn (assessed) binary phase diagrams compiled "
            "in H. Okamoto's Desk Handbook (ASM International, 2010): the thermodynamic coordinate (T, x, P) of a point deep "
            "inside the field and the phase assemblage present there, with phase-registry ids, coordinate uncertainty from "
            "the axis calibration, naming basis, confidence and review state, a cross-check against the calculated SGTE "
            "diagrams, and citations (original evaluation first, the handbook as compilation). Evidence class "
            "<code>assessed_phase_diagram</code>. Part of the ORCHID data release %s.</p>"
            "<h3>What</h3><ul><li>%s points from %s diagrams (centroid %s, pole of inaccessibility %s); %s two-phase points "
            "with a drawing-based tie_line block (low credence).</li><li>ESPEI ZPF datasets, one per diagram: default "
            "(two-phase compositions null with a __HYPERPLANE__ at the point) and a variant with tie-line endpoints. "
            "The files pass ESPEI dataset-format checks; this establishes format compatibility and does not validate the "
            "scientific correctness of all assemblages or their locations.</li><li>JSONL, CSV, Parquet, exclusions with reasons, SCHEMA.md.</li></ul>"
            "<h3>Limitations</h3><ul><li>Version 1, machine-extracted: fields named by deterministic rules, model readers, "
            "an adjudicator and partial expert review (basis and review state per record).</li><li>Drawing-based tie "
            "lines and interior points are soft constraints with recorded naming basis and review state.</li><li>No page "
            "images, crops or handbook text are distributed.</li></ul><h3>Licence</h3><p>CC BY 4.0. %s</p><h3>How to cite</h3><p>%s</p><p>Part of: %s</p><p>Cite each diagram's "
            "original evaluation and the Okamoto handbook as the compilation.</p>"
            % (e["title"], orchid_cite.RELEASE_TAG, fmt(s["points"]), fmt(s["diagrams_with_points"]),
               fmt(s["point_method"]["centroid"]), fmt(s["point_method"]["pole_of_inaccessibility"]),
               fmt(s["two_phase_tie_lines"]["with_tie_line"]), orchid_cite.LICENCE_SCOPE, e["text"], u["text"]))
    return {"_comment": "Zenodo deposition metadata (legacy REST format). Written by tools/build_assessed_assemblages.py. "
                        "DOI placeholders 10.5281/zenodo.PENDING-<name> are replaced by tools/zenodo_release.py writeback.",
            "record": DATASET, "doi": e["doi"],
            "metadata": {"title": e["title"], "upload_type": "dataset", "version": orchid_cite.VERSION, "description": desc + project_credits.html_block() + "<p>" + orchid_cite.FUNDING_NOTE + "</p>",
                         "creators": [{"name": "%s, %s" % (f, g), "affiliation": af, "orcid": o} for f, g, af, o in orchid_cite.AUTHORS],
                         "keywords": ["phase diagrams", "assessed phase diagram", "phase equilibria", "ZPF", "ESPEI", "CALPHAD",
                                      "binary alloys", "AI-ready data"],
                         "license": "cc-by-4.0", "access_right": "open", "language": "eng",
                         "communities": [{"identifier": orchid_cite.ZENODO_COMMUNITY}],
                         "related_identifiers": [
                             {"identifier": u["doi"], "relation": "isPartOf", "resource_type": "dataset", "scheme": "doi"},
                             {"identifier": orchid_cite.entry("phase-registry", data)["doi"], "relation": "references", "resource_type": "dataset", "scheme": "doi"},
                             {"identifier": orchid_cite.entry("phase-diagrams", data)["doi"], "relation": "references", "resource_type": "dataset", "scheme": "doi"},
                             {"identifier": "https://orchid-project.org", "relation": "isDocumentedBy", "scheme": "url"},
                             {"identifier": "https://data.orchid-project.org", "relation": "isDocumentedBy", "scheme": "url"}]
                         + ([{"identifier": orchid_cite.PAPER["identifier"], "relation": "isSupplementTo",
                              "scheme": orchid_cite.PAPER["scheme"]}] if orchid_cite.PAPER else []) + project_credits.related_identifiers(),
                         "notes": "Status: Version 1, machine-extracted. " + orchid_cite.FUNDING_NOTE}}


# --------------------------------------------------------------------------------------------------- main --
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=orchid_cite.DEFAULT_DATA, help="orchid-data checkout (writes %s/)" % OUT_DIR)
    ap.add_argument("--run", default=DEF_RUN, help="annotation run directory (read only)")
    ap.add_argument("--hull", default=DEF_HULL, help="tier-4 DFT-hull decisions (JSONL, read only)")
    ap.add_argument("--tools", default=DEF_TOOLS, help="directory of phase_mapper.py (default: this tools/ directory)")
    ap.add_argument("--tdb-dir", default=DEF_TDB_DIR, help="SGTE binary TDB files (read only)")
    ap.add_argument("--espei-python", default=DEF_ESPEI, help="interpreter with ESPEI for validation")
    ap.add_argument("--dense", choices=("default", "none", "all", "report"), default="default",
                    help="dense field grid: default = validated diagrams only (%s); none; all = every diagram (not for release "
                         "without validation); report = count the points of every diagram and exit, writing nothing"
                         % ", ".join(sorted(DENSE_VALIDATED)))
    ap.add_argument("--only", default=None, help="development: build only this diagram directory (glob); use with a scratch --data")
    ap.add_argument("--check", action="store_true", help="exit 1 if the outputs differ from a fresh build")
    ap.add_argument("--zenodo-only", action="store_true",
                    help="only (re)write (or with --check, compare) zenodo/%s.metadata.json from the released "
                         "%s/summary.json, without rebuilding the data" % (DATASET, OUT_DIR))
    a = ap.parse_args(argv)
    if a.zenodo_only:
        zp = os.path.join(a.data, "zenodo", DATASET + ".metadata.json")
        with open(os.path.join(a.data, OUT_DIR, "summary.json"), encoding="utf-8") as fh:
            zen = (json.dumps(zenodo_metadata(json.load(fh), a.data), ensure_ascii=False, indent=1) + "\n").encode()
        cur = open(zp, "rb").read() if os.path.exists(zp) else None
        if a.check:
            print(("%s is current" if cur == zen else "stale: %s") % os.path.relpath(zp, a.data),
                  file=sys.stdout if cur == zen else sys.stderr)
            return 0 if cur == zen else 1
        if cur != zen:
            with open(zp + ".tmp", "wb") as fh:
                fh.write(zen)
            os.replace(zp + ".tmp", zp)
            print("wrote " + os.path.relpath(zp, a.data))
        return 0
    cache = os.environ.get("ASPA_DEV_CACHE")   # development only: reuse build() results (pickle) between runs
    if cache and os.path.exists(cache):
        import pickle
        records, exclusions, espei_excl, ed, et, cnt, method = pickle.load(open(cache, "rb"))
    else:
        records, exclusions, espei_excl, ed, et, cnt, method = build(a)
        if cache:
            import pickle
            pickle.dump((records, exclusions, espei_excl, ed, et, cnt, method), open(cache, "wb"))
    if a.dense == "report":
        c = a._dense["counts"]
        per = {k[8:]: v for k, v in c.items() if k.startswith("diagram:")}
        print(json.dumps({"fields_with_grid": c["fields_with_grid"], "points": c["points"], "diagrams": len(per),
                          "primary_points": len(records), "points_by_diagram": dict(sorted(per.items()))}, indent=1))
        return 0
    out, summ, cite, pqb = write_outputs(a, records, exclusions, espei_excl, ed, et, cnt, method)
    tmp = tempfile.mkdtemp(prefix="aspa-espei-")
    try:
        for rel, b in out.items():
            if rel.startswith("espei-zpf"):
                p = os.path.join(tmp, rel)
                os.makedirs(os.path.dirname(p), exist_ok=True)
                with open(p, "wb") as fh:
                    fh.write(b)
        vdirs = [os.path.join(tmp, "espei-zpf"), os.path.join(tmp, "espei-zpf-with-tie-lines")]
        if os.path.isdir(os.path.join(tmp, DENSE_ESPEI_DIR)):
            vdirs.append(os.path.join(tmp, DENSE_ESPEI_DIR))
        espei_res = validate_espei(a.espei_python, vdirs)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    out = finalize(a, out, summ, cite, pqb, espei_res)
    zen_rel = os.path.join("zenodo", DATASET + ".metadata.json")
    zen = (json.dumps(zenodo_metadata(summ, a.data), ensure_ascii=False, indent=1) + "\n").encode()
    out_dir = os.path.join(a.data, OUT_DIR)
    if a.check:
        bad = []
        for rel, b in sorted(out.items()):
            p = os.path.join(out_dir, rel)
            if not os.path.exists(p):
                bad.append("missing " + rel)
                continue
            cur = open(p, "rb").read()
            if rel.endswith("records.parquet"):
                if parquet_content(cur) != parquet_content(b):
                    bad.append(rel)
            elif cur != b:
                bad.append(rel)
        disk = {os.path.relpath(os.path.join(d, f), out_dir) for d, _, fs in os.walk(out_dir) for f in fs}
        bad += ["orphan " + x for x in sorted(disk - set(out))]
        zp = os.path.join(a.data, zen_rel)
        if not os.path.exists(zp) or open(zp, "rb").read() != zen:
            bad.append(zen_rel)
        if bad:
            print("stale: " + ", ".join(bad[:20]) + (" ... (%d)" % len(bad) if len(bad) > 20 else ""), file=sys.stderr)
            return 1
        print("%s/ is up to date (%d files, %d points; ESPEI valid)" % (OUT_DIR, len(out), summ["points"]))
        return 0
    os.makedirs(out_dir, exist_ok=True)
    disk = {os.path.relpath(os.path.join(d, f), out_dir) for d, _, fs in os.walk(out_dir) for f in fs}
    for rel in sorted(disk - set(out)):
        os.remove(os.path.join(out_dir, rel))
    for rel, b in sorted(out.items()):
        p = os.path.join(out_dir, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        if rel.endswith("records.parquet") and os.path.exists(p) and parquet_content(open(p, "rb").read()) == parquet_content(b):
            continue   # same content: keep the existing bytes
        if not os.path.exists(p) or open(p, "rb").read() != b:
            with open(p + ".tmp", "wb") as fh:
                fh.write(b)
            os.replace(p + ".tmp", p)
    zp = os.path.join(a.data, zen_rel)
    if not os.path.exists(zp) or open(zp, "rb").read() != zen:
        with open(zp, "wb") as fh:
            fh.write(zen)
    print(json.dumps({k: v for k, v in summ.items() if k not in ("cite_as", "method")}, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
