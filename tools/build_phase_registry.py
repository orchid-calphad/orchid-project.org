#!/usr/bin/env python3
"""ORCHID phase registry v0 (binary systems): one entry per (binary system, phase).

Reads (never writes) local sources and writes <orchid-data>/phases/:

  structure_types.json  system-independent structure types (prototype, Strukturbericht, Pearson, space group), cited
  registry.json         compact index: one row per phase of each binary system (id, system, names, structure summary
                        with tier and prototype, counts); no TDB definitions or reference details
  systems/<SYS>.json    full records of one system: names, models with tdb_definition, stability, literature counts,
                        example references, structure assignment or review candidates
  registry.jsonld       the same terms as a schema.org DefinedTermSet (context: context.jsonld)
  registry.csv          one flat row per phase
  structure_types.csv   one flat row per structure type (also structure_types.parquet)
  review_queue.json     phases without an assigned structure that have candidates (with the evidence)
  CONFLICTS.md          same name -> different structures, polymorph families, models that disagree
  index.json            compact search index for the website
  registry.parquet      flat table (only when the running interpreter has pyarrow; --parquet-python otherwise)
  SCHEMA.md, CITATION.cff

Sources (read only):
  open-tdb catalog/catalog.json            licence tier (bundled + confirmed, two-element records only)
  open-tdb databases/**                    PHASE / CONSTITUENT / TYPE_DEFINITION declarations, header comments
  open-tdb catalog/model-audit/bundled-summary.json   SHA-256 cross-check of every TDB read
  open-tdb hulls/dft-mp, hulls/dft-oqmd    on-hull vertices (formula, space group, degenerate alternatives)
  sgte-corpus db-v2.1/canonical_phase_aliases.json and phase_mapping records (corpus ids, paper label -> TDB phase)
  sgte-corpus tdb_pdf_links.json           DOI of the assessment behind each TDB (fallback)
  orchid-data literature/systems/*.json, literature/papers.json   paper labels, value counts, citations
  orchid-data phase-diagrams/sgte/*.json   where each phase is stable in the calculated diagrams
  Hallstedt 2025 supplementary table (mmc2.xlsx)  crystallographic table of the SGTE collection (CC BY 4.0)

Structure assignment (PI rule: assign only when certain; otherwise list candidates and mark needs_review):
  A  standard_name     the phase name encodes one structure (FCC_A1, BCC_B2, L12_FCC, AL3NI_D011, SIGMA, ...) and the
                       model's stoichiometry / sublattices do not contradict it, and DFT does not contradict it
  B  dft_unique_match  a stoichiometric compound whose composition matches exactly one on-hull entry in BOTH Materials
                       Project and OQMD, with the same space group, a unique structure type for that space group and
                       stoichiometry, no other phase of the system at that composition, and no structure in the name
  C  assessment_declared  the Hallstedt 2025 crystallographic table (mmc2) lists exactly one structure for this phase
                       name in this system, nothing in A/B contradicts it, and the TDB model matches the listed compound
                       (confidence "assessor-stated"; source cited)
  liquid / gas         structure_type "liquid" / "gas" (no prototype)
  everything else      no structure type; `candidates` with their evidence; status needs_review

TDB definitions: every model carries `tdb_definition` (verbatim PHASE / CONSTITUENT / TYPE_DEFINITION statements with
line numbers, a parsed view, parameter counts and provenance links), using open-tdb's statement and parameter parsing.

Deterministic: sorted output, no timestamps; input hashes recorded. `--check` rebuilds in memory and compares.

  python3 tools/build_phase_registry.py [--data ~/projects/orchid-data] [--check]
"""
from __future__ import annotations

import argparse
import collections
import csv
import glob
import hashlib
import io
import json
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
import zipfile
from fractions import Fraction

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import doi_corrections  # noqa: E402
import orchid_cite  # noqa: E402

DOI_CORR = doi_corrections.load()   # shared DOI correction table

HOME = os.path.expanduser("~")
DEF_OPEN_TDB = os.path.join(HOME, "projects", "open-tdb")
DEF_CORPUS = os.path.join(DEF_OPEN_TDB, "build", "sgte-corpus")
DEF_MMC2 = os.path.join(HOME, "data", "PD-extraction", "sources", "hallstedt2025_mmc2.xlsx")
DATASET = "phase-registry"
REGISTRY_VERSION = "0.1.0"
SITE_DATA = "https://data.orchid-project.org/phases/"
SITE = "https://orchid-project.org"
PIPELINE_GITHUB = "https://github.com/orchid-calphad/orchid-pipeline"
MAX_REFS = 10
TOL = 1e-6

HALLSTEDT = {"citation": "B. Hallstedt, The SGTE collection of binary datasets, Calphad 89 (2025) 102833, "
                         "supplementary crystallographic table (mmc2)", "doi": "10.1016/j.calphad.2025.102833",
             "license": "CC BY 4.0"}
REF_STD = [
    {"citation": "P. Villars, L. D. Calvert, Pearson's Handbook of Crystallographic Data for Intermetallic Phases, "
                 "2nd ed., ASM International (1991)", "doi": None},
    {"citation": "M. J. Mehl et al., The AFLOW Library of Crystallographic Prototypes: Part 1, Comput. Mater. Sci. "
                 "136 (2017) S1-S828", "doi": "10.1016/j.commatsci.2017.01.017"},
]

# Built-in standard structure types: (Strukturbericht, prototype, Pearson, space-group number, symbol).
# Written from memory against the standard tables above; rows the Hallstedt table confirms (same Strukturbericht,
# Pearson and space-group number) are marked verified, the others carry "%TODO verify".
BUILTIN = [
    ("A1", "Cu", "cF4", 225, "Fm-3m"), ("A2", "W", "cI2", 229, "Im-3m"), ("A3", "Mg", "hP2", 194, "P6_3/mmc"),
    ("A3'", "alpha-La", "hP4", 194, "P6_3/mmc"), ("A4", "C (diamond)", "cF8", 227, "Fd-3m"),
    ("A5", "beta-Sn", "tI4", 141, "I4_1/amd"), ("A6", "In", "tI2", 139, "I4/mmm"),
    ("A7", "alpha-As", "hR2", 166, "R-3m"), ("A8", "gamma-Se", "hP3", 152, "P3_121"),
    ("A9", "C (graphite)", "hP4", 194, "P6_3/mmc"), ("A10", "alpha-Hg", "hR1", 166, "R-3m"),
    ("A11", "alpha-Ga", "oC8", 64, "Cmce"), ("A12", "alpha-Mn", "cI58", 217, "I-43m"),
    ("A13", "beta-Mn", "cP20", 213, "P4_132"), ("A15", "Cr3Si", "cP8", 223, "Pm-3n"),
    ("A20", "alpha-U", "oC4", 63, "Cmcm"),
    ("B1", "NaCl", "cF8", 225, "Fm-3m"), ("B2", "CsCl", "cP2", 221, "Pm-3m"), ("B3", "ZnS (sphalerite)", "cF8", 216, "F-43m"),
    ("B4", "ZnS (wurtzite)", "hP4", 186, "P6_3mc"), ("B8_1", "NiAs", "hP4", 194, "P6_3/mmc"),
    ("B19", "AuCd", "oP4", 51, "Pmma"), ("B20", "FeSi", "cP8", 198, "P2_13"), ("B27", "FeB", "oP8", 62, "Pnma"),
    ("B32", "NaTl", "cF16", 227, "Fd-3m"), ("B33", "CrB", "oC8", 63, "Cmcm"),
    ("C1", "CaF2", "cF12", 225, "Fm-3m"), ("C11_b", "MoSi2", "tI6", 139, "I4/mmm"),
    ("C14", "MgZn2", "hP12", 194, "P6_3/mmc"), ("C15", "Cu2Mg", "cF24", 227, "Fd-3m"),
    ("C16", "Al2Cu", "tI12", 140, "I4/mcm"), ("C32", "AlB2", "hP3", 191, "P6/mmm"),
    ("C36", "MgNi2", "hP24", 194, "P6_3/mmc"), ("C40", "CrSi2", "hP9", 180, "P6_222"),
    ("C49", "ZrSi2", "oC12", 63, "Cmcm"), ("C54", "TiSi2", "oF24", 70, "Fddd"),
    ("D0_3", "BiF3", "cF16", 225, "Fm-3m"), ("D0_11", "Fe3C", "oP16", 62, "Pnma"),
    ("D0_19", "Ni3Sn", "hP8", 194, "P6_3/mmc"), ("D0_22", "Al3Ti", "tI8", 139, "I4/mmm"),
    ("D0_23", "Al3Zr", "tI16", 139, "I4/mmm"), ("D0_24", "Ni3Ti", "hP16", 194, "P6_3/mmc"),
    ("D0_a", "beta-Cu3Ti", "oP8", 59, "Pmmn"), ("D5_13", "Ni2Al3", "hP5", 164, "P-3m1"),
    ("D8_2", "Cu5Zn8", "cI52", 217, "I-43m"), ("D8_5", "Fe7W6", "hR13", 166, "R-3m"),
    ("D8_8", "Mn5Si3", "hP16", 193, "P6_3/mcm"), ("D8_b", "sigma-CrFe", "tP30", 136, "P4_2/mnm"),
    ("D8_m", "W5Si3", "tI32", 140, "I4/mcm"), ("L1_0", "AuCu", "tP4", 123, "P4/mmm"),
    ("L1_2", "AuCu3", "cP4", 221, "Pm-3m"), ("L2_1", "Cu2MnAl", "cF16", 225, "Fm-3m"),
]
# words that name one structure in SGTE-style phase names (only used when no Strukturbericht token is present)
WORD_SB = {"SIGMA": "D8B", "CHI": "A12", "MU": "D85", "DIAMOND": "A4", "GRAPHITE": "A9", "ZINCBLENDE": "B3",
           "HALITE": "B1", "CEMENTITE": "D011", "DHCP": "A3'", "FCC": "A1", "BCC": "A2", "HCP": "A3"}
LATTICE_FAMILY = {"FCC": {"A1", "L12", "L10", "D022", "D023", "D1A"}, "BCC": {"A2", "B2", "D03", "L21", "B32"},
                  "HCP": {"A3", "D019", "B19"}, "LAVES": {"C14", "C15", "C36"}}
DISORDERED = {"A1", "A2", "A3"}
# crystal-system words: allowed beside a Strukturbericht token when they agree with its Pearson symbol
CRYSTAL_WORDS = {"CUB": "c", "CUBIC": "c", "CBCC": "c", "TETRAGONAL": "t", "TETRAG": "t", "TET": "t", "BCT": "t",
                 "ORTHORHOMBIC": "o", "ORTHO": "o", "HEXAGONAL": "hP", "HEX": "hP", "RHOMBOHEDRAL": "hR",
                 "RHOMBO": "hR", "RHOMB": "hR", "MONOCLINIC": "m", "TRICLINIC": "a"}
SB_TOKEN = re.compile(r"^[A-L](\d+[A-Z]?\d*|[A-Z])$")
BENIGN = {"OLD", "NEW", "PHASE", "OXIDE", "HT", "LT", "RT", "HT1", "HT2", "LT1", "LT2"}
LIQ_TOKENS = {"LIQUID", "LIQ", "L", "IONIC", "ION", "IONIC_LIQ"}
GREEK_WORDS = ("ALPHA", "BETA", "GAMMA", "DELTA", "EPSILON", "ZETA", "ETA", "THETA", "KAPPA", "LAMBDA", "MU", "NU",
               "XI", "PI", "RHO", "SIGMA", "TAU", "PHI", "CHI", "PSI", "OMEGA")


# ----------------------------------------------------------------------------------------------- helpers --
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def dump(obj, indent=1):
    return json.dumps(obj, indent=indent, ensure_ascii=False, sort_keys=False) + "\n"


def sb_norm(sb):
    if not sb:
        return None
    s = str(sb).strip().replace("_", "").replace(" ", "").upper()
    return None if s in ("", "-", "?") else s


def parse_sg(text):
    m = re.match(r"\s*0*(\d+)\s*,\s*(\S+)", str(text or ""))
    return (int(m.group(1)), m.group(2)) if m else (None, None)


def formula_counts(text):
    """'Ni2Al3' -> {'Al': 3, 'Ni': 2}; None unless the whole text is a formula of element symbols."""
    t = re.sub(r"^\(|\)$", "", str(text or "").strip())
    t = re.sub(r"^(alpha|beta|gamma|delta|epsilon|sigma|eta|zeta|theta|kappa|lambda|mu|omega)-", "", t, flags=re.I)
    parts = re.findall(r"([A-Z][a-z]?)(\d+(?:\.\d+)?)?", t)
    if not parts or "".join(a + (b or "") for a, b in parts) != t:
        return None
    out = collections.Counter()
    for el, n in parts:
        out[el] += Fraction(n) if n else 1
    return dict(out)


def ratio_key(counts):
    """Sorted normalised amounts (stoichiometry pattern, element-independent): {'Al':3,'Ni':2} -> (2/5, 3/5)."""
    if not counts:
        return None
    tot = sum(Fraction(v) for v in counts.values())
    return tuple(sorted(Fraction(v) / tot for v in counts.values()))


def proto_key(p):
    c = formula_counts(p)
    if c and len(c) > 1:
        pre = re.match(r"^\(?((?:alpha|beta|gamma|delta|sigma|eta|zeta|theta|omega)-)", p.strip(), flags=re.I)
        return (pre.group(1).lower() if pre else "") + "".join("%s%s" % (k, c[k]) for k in sorted(c))
    return re.sub(r"\s+", "", p.strip().strip("()")).lower()


def slug(t):
    return re.sub(r"[^a-z0-9]+", "-", str(t).lower()).strip("-")


def pearson_count(p):
    m = re.match(r"^([amothc])([PCFIRSAB])(\d+)$", str(p or ""))
    return (m.group(2), int(m.group(3))) if m else (None, None)


def nsites_ok(pearson, nsites):
    cen, n = pearson_count(pearson)
    if n is None or nsites is None:
        return None
    div = {"P": 1, "C": 2, "A": 2, "B": 2, "S": 2, "I": 2, "F": 4, "R": 1}[cen]
    ok = {n, n // div if n % div == 0 else n}
    if cen == "R":
        ok.add(3 * n)
    return nsites in ok


# ------------------------------------------------------------------------------------------- xlsx (mmc2) --
def read_xlsx_rows(path):
    ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).iter(ns + "si"):
            ss.append("".join(t.text or "" for t in si.iter(ns + "t")))
    sheets = sorted(n for n in z.namelist() if re.match(r"xl/worksheets/sheet\d+\.xml$", n))
    root = ET.fromstring(z.read(sheets[0]))
    rows = []
    for row in root.iter(ns + "row"):
        cells = {}
        for c in row.iter(ns + "c"):
            ref = re.match(r"([A-Z]+)", c.get("r")).group(1)
            col = 0
            for ch in ref:
                col = col * 26 + ord(ch) - 64
            t, v = c.get("t"), c.find(ns + "v")
            if t == "s" and v is not None:
                val = ss[int(v.text)]
            elif t == "inlineStr":
                val = "".join(x.text or "" for x in c.iter(ns + "t"))
            else:
                val = v.text if v is not None else None
            cells[col - 1] = val
        if cells:
            rows.append([cells.get(i) for i in range(max(cells) + 1)])
    return rows


def load_mmc2(path):
    rows = read_xlsx_rows(path)
    out = []
    for r in rows[1:]:
        r = (list(r) + [None] * 7)[:7]
        proto, pearson, sb, sg, name, model, actual = [(str(x).strip() if x is not None else None) for x in r]
        num, sym = parse_sg(sg)
        out.append({"prototype": proto if proto and proto not in ("?", "-") else None,
                    "pearson": pearson if pearson and "?" not in pearson and "*" not in pearson else None,
                    "sb": sb if sb and sb not in ("-", "?") else None, "sg_number": num, "sg_symbol": sym,
                    "phase_name": name.upper() if name and name != "-" else None, "model": model, "actual": actual})
    return out


# ----------------------------------------------------------------------------------------- structure types --
def build_structure_types(mmc2):
    types = {}       # key -> type

    def key_of(sbn, proto, pearson, sgn):
        return ("sb", sbn, pearson, sgn) if sbn else ("proto", proto_key(proto or "?"), pearson, sgn)

    for r in mmc2:
        if not r["prototype"] and not r["sb"]:
            continue
        sbn = sb_norm(r["sb"])
        k = key_of(sbn, r["prototype"], r["pearson"], r["sg_number"])
        t = types.setdefault(k, {"sb": r["sb"], "sbn": sbn, "prototypes": collections.Counter(), "pearson": r["pearson"],
                                 "sg_number": r["sg_number"], "sg_symbol": r["sg_symbol"], "from_mmc2": 0,
                                 "builtin": False, "verified": False})
        if r["prototype"]:
            t["prototypes"][r["prototype"].replace("\t", "").strip()] += 1
        t["from_mmc2"] += 1
    for sb, proto, pearson, sgn, sgs in BUILTIN:
        sbn = sb_norm(sb)
        k = key_of(sbn, proto, pearson, sgn)
        if k in types:
            types[k].update(builtin=True, verified=True, builtin_prototype=proto, sb=sb, sg_symbol=types[k]["sg_symbol"] or sgs)
        else:
            types[k] = {"sb": sb, "sbn": sbn, "prototypes": collections.Counter({proto: 0}), "pearson": pearson,
                        "sg_number": sgn, "sg_symbol": sgs, "from_mmc2": 0, "builtin": True, "verified": False,
                        "builtin_prototype": proto}
    out, used = [], set()
    for k in sorted(types, key=lambda k: tuple("" if x is None else str(x) for x in k)):
        t = types[k]
        protos = sorted(t["prototypes"], key=lambda p: (-t["prototypes"][p], p))
        proto = t.get("builtin_prototype") or (protos[0] if protos else None)
        proto = re.sub(r"\s*\?$", "", proto or "").strip() or None
        counts = formula_counts(proto) if proto else None
        if counts and len(counts) == 1:
            counts = None
        base = "-".join(x for x in [slug(t["sb"]) if t["sb"] else None, slug(proto or "x"), slug(t["pearson"] or "")] if x)
        sid, i = "st:" + base, 2
        while sid in used:
            sid, i = "st:%s-%d" % (base, i), i + 1
        used.add(sid)
        src = []
        if t["from_mmc2"]:
            src.append(dict(HALLSTEDT, rows=t["from_mmc2"]))
        if t["builtin"]:
            src += [dict(x) for x in REF_STD]
        out.append({
            "id": sid, "prototype": proto,
            "prototype_spellings": [p.replace("\t", "").strip() for p in protos if p != proto and "?" not in p],
            "strukturbericht": t["sb"], "pearson": t["pearson"],
            "space_group": {"number": t["sg_number"], "symbol": t["sg_symbol"]} if t["sg_number"] else None,
            "stoichiometry": "".join("%s%s" % (k2, counts[k2]) for k2 in sorted(counts)) if counts else None,
            "sources": src,
            "verification": ("in Hallstedt 2025 mmc2" if t["from_mmc2"] and not t["builtin"] else
                             "built-in row confirmed by Hallstedt 2025 mmc2 (Strukturbericht, Pearson, space group)"
                             if t["verified"] else "%TODO verify (built-in row written from memory; not in mmc2)"),
            "_sbn": t["sbn"], "_ratio": ratio_key(counts) if counts else None})
    out.append({"id": "st:liquid", "prototype": None, "prototype_spellings": [], "strukturbericht": None,
                "pearson": None, "space_group": None, "stoichiometry": None, "sources": [],
                "verification": "definition (non-crystalline)", "_sbn": None, "_ratio": None})
    out.append({"id": "st:gas", "prototype": None, "prototype_spellings": [], "strukturbericht": None,
                "pearson": None, "space_group": None, "stoichiometry": None, "sources": [],
                "verification": "definition (non-crystalline)", "_sbn": None, "_ratio": None})
    return out


# ------------------------------------------------------------------------------------------------- TDBs --
def parse_tdb(path):
    text = open(path, "rb").read().decode("latin-1")
    phases, consts, ordis, comments_for = {}, {}, {}, {}
    buf, recent = [], []
    header = []
    seen_stmt = False
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("$"):
            c = line.lstrip("$").strip()
            if c:
                recent.append(c)
                if not seen_stmt:
                    header.append(c)
            continue
        if not line:
            continue
        buf.append(line)
        if "!" in line:
            stmt = " ".join(buf).split("!")[0].strip()
            buf = []
            seen_stmt = True
            toks = stmt.split()
            if not toks:
                recent = []
                continue
            kw = toks[0].upper()
            if kw in ("PHASE", "PHAS") and len(toks) >= 4:
                name = toks[1].split(":")[0].upper()
                try:
                    n = int(toks[3])
                    ratios = [float(x) for x in toks[4:4 + n]]
                except ValueError:
                    n, ratios = None, []
                phases[name] = {"type_code": toks[2], "n_sublattices": n, "site_ratios": ratios,
                                "comments": list(recent[-6:])}
            elif kw.startswith("CONST") and len(toks) >= 2:
                rest = stmt.split(None, 1)[1]
                nm, _, subl = rest.partition(":")
                nm = nm.strip().split()[0].upper() if nm.strip() else ""
                parts = subl.split(":")
                if parts and re.match(r"^[A-Z]$", parts[0].strip().upper()):
                    parts = parts[1:]  # NAME:L : ...
                sl = []
                for p in parts:
                    cs = [c.strip().rstrip("%").upper() for c in re.split(r"[,\s]+", p.strip()) if c.strip()]
                    if cs:
                        sl.append(sorted(set(cs)))
                consts[nm] = sl
            elif kw.startswith("TYPE"):
                up = [t.upper() for t in toks]
                for marker in ("DIS_PART", "DISORDERED_PART", "DIS_PA", "DISORD"):
                    if marker in up and "AMEND_PHASE_DESCRIPTION" in up:
                        i, j = up.index("AMEND_PHASE_DESCRIPTION"), up.index(marker)
                        if i + 1 < len(up) and j + 1 < len(up):
                            ordis[up[i + 1].split(":")[0]] = up[j + 1].split(":")[0]
                        break
            recent = []
    out = {}
    for name, p in phases.items():
        p = dict(p)
        p["constituents"] = consts.get(name, [])
        out[name] = p
    return out, ordis, header


# ------------------------------------------------------------------------------- TDB phase definitions --
# Verbatim PHASE / CONSTITUENT / TYPE_DEFINITION statements with line numbers, a parsed view and a parameter summary
# per phase model. Statement boundaries, regular expressions and parameter keys are open-tdb's (scripts/validate_tdb.py,
# scripts/check_subsystem_compat.py, scripts/audit_tdb_models.py); TDB text is read as text and never executed.
_OT_MODULES = {}


def open_tdb_modules(open_tdb):
    """validate_tdb, check_subsystem_compat, audit_tdb_models from <open-tdb>/scripts (standard library only)."""
    if not _OT_MODULES:
        sd = os.path.join(open_tdb, "scripts")
        sys.path.insert(0, sd)
        try:
            import validate_tdb  # noqa: E402
            import check_subsystem_compat  # noqa: E402
            import audit_tdb_models  # noqa: E402
        finally:
            sys.path.remove(sd)
        _OT_MODULES.update(vt=validate_tdb, compat=check_subsystem_compat, audit=audit_tdb_models)
    return _OT_MODULES


def statement_spans(text, vt):
    """Character spans of the statements validate_tdb.statements() yields, in the same order and with the same
    rules ('$' comments, several statements per line, LIST_OF_REFERENCE blocks, '!' terminators). Each item:
    (start_offset, end_offset_after_terminator, start_line, end_line, boundary)."""
    out = []
    pos = 0
    pending, start, start_line, last_end, last_line = False, 0, 0, 0, 0
    reference_block = False
    for ln, raw in enumerate(text.splitlines(keepends=True), 1):
        base = pos
        pos += len(raw)
        content = raw.splitlines()[0] if raw.splitlines() else ""
        cut = content.split("$", 1)[0]
        if not cut.strip():
            continue
        if re.match(r"^LIST[-_]OF[-_]REFERENCE", cut.lstrip(), re.IGNORECASE):
            reference_block = True
            continue
        if reference_block:
            d = cut.lstrip()
            if not (re.match(r"^(?:ELEMENT|SPECIES)\s+\S+\s+", d, re.I) or
                    re.match(r"^FUNCTION\s+\S+\s+[+-]?(?:\d|\.)", d, re.I) or
                    re.match(r"^PHASE\s+\S+\s+%", d, re.I) or re.match(r"^TYPE_DEFINITION\s+\S+", d, re.I)):
                continue
            reference_block = False
        if pending and vt.COMMAND_START_RE.match(cut.lstrip()):
            out.append((start, last_end, start_line, last_line, "missing"))
            pending = False
        off = 0
        broke = False
        while "!" in cut[off:]:
            i = cut.index("!", off)
            before = cut[off:i]
            if before.strip() or pending:
                if not pending:
                    start, start_line = base + off + (len(before) - len(before.lstrip())), ln
                out.append((start, base + i + 1, start_line, ln, "complete"))
                pending = False
            off = i + 1
            rem = cut[off:]
            if rem.strip() and not vt.COMMAND_START_RE.match(rem.lstrip()):
                broke = True
                break
        rem = "" if broke else cut[off:]
        if rem.strip():
            if not pending:
                start, start_line = base + off + (len(rem) - len(rem.lstrip())), ln
                pending = True
            last_end, last_line = base + off + len(rem.rstrip()), ln
    if pending:
        out.append((start, last_end, start_line, last_line, "trailing"))
    return out


def verbatim(text, s, e):
    """The statement exactly as in the file; only its line breaks are normalised to LF."""
    return re.sub(r"\r\n|\r|\n|\x0b|\x0c|\x1c|\x1d|\x1e|\x85| | ", "\n", text[s:e])


def _const_lists(body):
    """CONSTITUENT body ' : AL NI% : VA : ' -> [['AL', 'NI%'], ['VA']] (spelling and '%' markers as printed)."""
    pieces = body.strip().split(":")
    pieces = pieces[1:-1] if not pieces[-1].strip() else pieces[1:]
    return [[t for t in re.split(r"[\s,]+", p.strip()) if t] for p in pieces]


def _num(t):
    try:
        v = float(t)
    except ValueError:
        return t
    return int(v) if v == int(v) and abs(v) < 1e15 else v


def tdb_definitions(path, open_tdb, rel_path, sha, links):
    """{PHASE_NAME: tdb_definition} for every PHASE statement of one TDB file, plus per-file statistics."""
    M = open_tdb_modules(open_tdb)
    vt, compat, audit = M["vt"], M["compat"], M["audit"]
    data = open(path, "rb").read()
    text, issue = vt.decode(data)
    stmts = list(vt.statements(text))
    spans = statement_spans(text, vt)
    if len(spans) != len(stmts):
        raise SystemExit("statement spans disagree with validate_tdb.statements: " + rel_path)
    phases, consts, typedefs = {}, {}, []
    for (st, bnd), (s, e, l0, l1, b2) in zip(stmts, spans):
        vb = verbatim(text, s, e)
        # the verbatim span must reduce to validate_tdb's statement text (comments dropped, parts joined)
        red = " ".join(x.split("$", 1)[0].strip() for x in vb.split("\n") if x.split("$", 1)[0].strip())
        red = red[:-1].rstrip() if red.endswith("!") else red
        if st.line != l0 or bnd != b2 or red != st.text:
            raise SystemExit("statement span mismatch in %s line %d" % (rel_path, st.line))
        rec = {"line_start": l0, "line_end": l1, "text": vb}
        value = st.text
        tm = compat.TYPE_DEFINITION_RE.match(value)
        if tm:
            typedefs.append((tm.group(1), value, rec))
            continue
        pm = compat.PHASE_RE.match(value)
        if pm:
            phases[vt.normalize_name(pm.group(1))] = (pm, rec)
            continue
        cm = vt.CONSTITUENT_RE.match(value)
        if cm:
            consts[vt.normalize_name(cm.group(1))] = (cm, rec)
    parsed_src = compat.parse_source(path, {})
    params = collections.defaultdict(list)
    unmatched = 0
    for key, recs in parsed_src["parameters"].items():
        parts = audit.parameter_parts(key)
        if not parts:
            unmatched += len(recs)
            continue
        kind, ph, conf, order = parts
        if ph not in phases:
            unmatched += len(recs)
            continue
        inter = any(len(a) > 1 for a in audit.configuration_arrays(conf))
        for _ in recs:
            params[ph].append((kind, 0 if order is None else order, inter))
    out = {}
    for name, (pm, prec) in phases.items():
        code = pm.group(2)
        codes = list(dict.fromkeys(code))
        tds, magnetic = [], None
        for ident, value, rec in typedefs:
            if ident not in codes:
                continue
            dm = compat.PHASE_DIRECTIVE_RE.match(value) or re.match(
                r"^TYPE[-_](?:DEFINITION|DEF)\s+(\S+)\s+GES\s+A_P_D\s+(\S+)\s+(.*)$", value, re.I)
            target = vt.normalize_name(dm.group(2)) if dm else None
            tds.append(dict(rec, type_code=ident, amends=target))
            if dm and target in ("@", name):
                toks = dm.group(3).split()
                if toks and toks[0].upper().startswith("MAG") and len(toks) >= 3:
                    magnetic = {"type_code": ident, "afm_factor": _num(toks[1]), "structure_fraction": _num(toks[2]),
                                "line": rec["line_start"]}
        n = int(pm.group(3))
        ratios = [_num(t) for t in pm.group(4).split()[:n]]
        cm, crec = consts.get(name, (None, None))
        cl = _const_lists(cm.group(2)) if cm else []
        ps = params.get(name, [])
        kinds = collections.Counter(k for k, _, _ in ps)
        lorders = collections.Counter(o for k, o, _ in ps if k == "L")
        out[name] = {
            "verbatim": {"phase": prec, "constituent": crec, "type_definitions": tds},
            "parsed": {"phase_name": name, "type_code": code, "type_codes": codes, "n_sublattices": n,
                       "site_ratios": ratios, "constituents": cl,
                       "major_constituents": [[c.rstrip("%") for c in sl if c.endswith("%")] for sl in cl],
                       "magnetic": magnetic},
            "parameters": {"total": len(ps), "by_kind": dict(sorted(kinds.items())),
                           "L_by_order": {str(k): v for k, v in sorted(lorders.items())},
                           "interaction_orders": sorted({o for _, o, i in ps if i}),
                           "endmember_G": sum(1 for k, _, i in ps if k == "G" and not i)},
            "provenance": dict({"path": rel_path, "sha256": sha, "decoded_as": "latin-1" if issue else "utf-8"},
                               **links)}
    stats = {"unmatched_parameters": unmatched, "type_codes_without_definition":
             sorted({c for _, (pm, _) in phases.items() for c in pm.group(2)} - {t[0] for t in typedefs})}
    return out, stats


def model_composition(p, elements):
    """Fixed composition of a stoichiometric model ({EL: amount}) or None (a sublattice with >1 constituent, a
    non-element species, or missing data)."""
    els = {e.upper() for e in elements}
    if not p["constituents"] or len(p["constituents"]) != len(p["site_ratios"]):
        return None
    comp = collections.Counter()
    for ratio, cs in zip(p["site_ratios"], p["constituents"]):
        real = [c for c in cs if c != "VA"]
        if len(cs) != 1:
            return None
        if real:
            if real[0] not in els:
                return None
            comp[real[0]] += Fraction(ratio).limit_denominator(10000)
    if len(comp) < 1:
        return None
    return dict(comp)


def partitioned(p):
    """>=2 non-vacancy sublattices with identical constituent sets of >=2 species (order/disorder capable)."""
    sets = [tuple(c for c in cs if c != "VA") for cs in p["constituents"]]
    sets = [s for s in sets if len(s) >= 2]
    return len(sets) >= 2 and len(set(sets)) == 1


def site_ratio_key(p):
    rs = [r for r, cs in zip(p["site_ratios"], p["constituents"]) if [c for c in cs if c != "VA"]]
    if len(rs) != 2:
        return None
    return ratio_key({"a": Fraction(rs[0]).limit_denominator(1000), "b": Fraction(rs[1]).limit_denominator(1000)})


# --------------------------------------------------------------------------- canonical ids (corpus convention) --
# Same rules as the SGTE corpus (tools-release-v2.1.18/phase_aliases.canonical_phase_id), so ids carry over.
_LATTICE = ("LIQUID", "GAS", "FCC", "BCC", "HCP", "DHCP", "CBCC", "CUB", "TET", "ORTHO", "RHOMB", "ION", "IONIC")
_SB = re.compile(r"^[A-L]\d{1,2}[A-Z]?$|^[A-L]\d_\d{1,2}$")
_TEMP_TAG = re.compile(r"^(HT\d?|LT\d?|RT)$")
_LIQ = {"LIQUID", "LIQ", "L", "LIQUID_PHASE", "LIQ_PHASE"}


def _is_formula(tok, elements):
    syms = sorted((e.upper() for e in elements), key=len, reverse=True)
    if not syms:
        return False
    rest = tok.upper()
    pat = re.compile(r"(" + "|".join(map(re.escape, syms)) + r")(\d+(?:\.\d+)?)?")
    n = 0
    while rest:
        m = pat.match(rest)
        if not m:
            return False
        rest, n = rest[m.end():], n + 1
    return n > 0


def canonical_phase_id(elements, model_phase, tdb_phases=(), mmc2_sb=None):
    system = "-".join(sorted(e.lower() for e in elements))
    name, _, inst = str(model_phase).partition("#")
    up = name.upper().strip()
    if up in _LIQ:
        return "phase:%s:liquid" % system
    toks = [t for t in re.split(r"[_\-\s]+", up) if t]
    tagged = [t for t in toks if _TEMP_TAG.match(t)]
    core = [t for t in toks if not _TEMP_TAG.match(t)]
    if tagged:
        siblings = {"_".join(t for t in re.split(r"[_\-\s]+", str(p).partition("#")[0].upper()) if t and not _TEMP_TAG.match(t))
                    for p in tdb_phases if p and str(p).partition("#")[0].upper() != up}
        if "_".join(core) not in siblings:
            tagged = []
    formula = [t for t in core if _is_formula(t, elements)]
    lattice = [t for t in core if t not in formula and t in _LATTICE]
    sb = [t for t in core if t not in formula and t not in lattice and _SB.match(t)]
    other = [t for t in core if t not in formula and t not in lattice and t not in sb]
    if not formula and not sb and not other and len(lattice) == 1:
        known = sorted((mmc2_sb or {}).get(up, set()))
        guess = known[0] if len(known) == 1 else {"FCC": "A1", "BCC": "A2", "HCP": "A3"}.get(lattice[0])
        if guess:
            sb = [guess]
    ordered = formula + sorted(lattice, key=_LATTICE.index) + sorted(sb) + sorted(other) + sorted(tagged)
    s = re.sub(r"[^a-z0-9]+", "-", "-".join(ordered).lower()).strip("-")[:80] or "x"
    return "phase:%s:%s" % (system, s)


# ------------------------------------------------------------------------------------------ name analysis --
def classify_name(name, elements, sb_types, model, st_by_id=None):
    """Tier-A analysis of a TDB phase name. Returns dict(kind, sbn, type_id, ok, reasons, caveats)."""
    up = name.upper()
    toks = [t for t in re.split(r"[_\-\s]+", up) if t]
    if up in _LIQ or set(toks) <= LIQ_TOKENS:
        return {"kind": "liquid"}
    if up == "GAS" or re.match(r"^[A-Z0-9]*GAS$", up):
        return {"kind": "gas"}
    sbs, words, unknown, formula_sb, crystal = [], [], [], [], []
    for t in toks:
        tn = sb_norm(t)
        if tn in sb_types and SB_TOKEN.match(t):
            if _is_formula(t, elements):
                formula_sb.append(tn)
            else:
                sbs.append(tn)
        elif t in WORD_SB or t == "LAVES":
            words.append(t)
        elif t in CRYSTAL_WORDS:
            crystal.append(t)
        elif t in GREEK_WORDS and t not in WORD_SB:
            continue
        elif _is_formula(t, elements) or t in BENIGN or _TEMP_TAG.match(t):
            continue
        else:
            unknown.append(t)
    res = {"kind": "crystal", "sbn": None, "type_id": None, "ok": False, "reasons": [], "caveats": [],
           "tokens": toks, "unknown": unknown}
    lattice_words = [w for w in words if w in LATTICE_FAMILY]
    if formula_sb and not sbs:
        # a token that is both a Strukturbericht code and a formula of this system (B2 in a boron system): only with a
        # confirming lattice word
        if any(fs in LATTICE_FAMILY.get(w, ()) for fs in formula_sb for w in lattice_words):
            sbs = formula_sb
        else:
            res["reasons"].append("structure_token_is_also_a_formula")
    sset = sorted(set(sbs))
    if crystal and not sset:
        unknown += crystal
    if len(sset) > 1:
        res["reasons"].append("several_structure_tokens:" + ",".join(sset))
        res["candidate_sbn"] = sset
        return res
    if len(sset) == 1:
        sbn = sset[0]
        for w in words:
            if w in LATTICE_FAMILY and sbn not in LATTICE_FAMILY[w]:
                res["reasons"].append("lattice_word_%s_contradicts_%s" % (w, sbn))
            elif w not in LATTICE_FAMILY and WORD_SB.get(w) != sbn:
                res["reasons"].append("word_%s_contradicts_%s" % (w, sbn))
    else:
        implied = {WORD_SB[w] for w in words if w in WORD_SB}
        if len(implied) != 1 or "LAVES" in words and not implied:
            if "LAVES" in words:
                res["candidate_sbn"] = sorted(LATTICE_FAMILY["LAVES"])
            res["reasons"].append("no_structure_token" if not implied else "words_imply_several_structures")
            return res
        sbn = implied.pop()
    res["sbn"] = sbn
    res["candidate_sbn"] = [sbn]
    if unknown:
        res["reasons"].append("unrecognised_name_tokens:" + ",".join(unknown))
    tids = sb_types.get(sbn, [])
    if len(tids) != 1:
        res["reasons"].append("strukturbericht_%s_maps_to_%d_structure_types" % (sbn, len(tids)))
        return res
    res["type_id"] = tids[0]
    pear = (st_by_id or {}).get(tids[0], {}).get("pearson") or ""
    for w in crystal:
        if not pear.startswith(CRYSTAL_WORDS[w]):
            res["reasons"].append("crystal_word_%s_contradicts_%s" % (w, pear or "unknown_pearson"))
    return res


# ----------------------------------------------------------------------------------------------- build --
class Inputs:
    def __init__(self):
        self.hashes = {}

    def add(self, label, path):
        self.hashes[label] = sha256_file(path)
        return path


def load_hull_vertices(open_tdb, inputs):
    out = {}  # (SYS, source) -> {"x_element":..., "vertices": [...]}
    for src, layer in (("mp", "dft-mp"), ("oqmd", "dft-oqmd")):
        cat_path = inputs.add("open-tdb/hulls/%s/catalog.json" % layer, os.path.join(open_tdb, "hulls", layer, "catalog.json"))
        cat = json.load(open(cat_path))
        for s in cat["systems"]:
            h = json.load(open(os.path.join(open_tdb, "hulls", s["hull_path"])))
            vs = []
            for iso in h["isotherms"]:
                vs += iso["vertices"]
            out[(h["system_id"], src)] = {"x_element": h["x_element"], "vertices": vs,
                                          "release": (h.get("provenance") or {}).get("source_release", (h.get("provenance") or {}).get("mp_release"))}
    return out


def build(a):
    inputs = Inputs()
    open_tdb, corpus_db, data = a.open_tdb, os.path.join(a.corpus, "db-v2.1"), a.data
    mmc2 = load_mmc2(inputs.add("Hallstedt 2025 mmc2.xlsx", a.mmc2))
    stypes = build_structure_types(mmc2)
    st_by_id = {t["id"]: t for t in stypes}
    sb_types = collections.defaultdict(list)
    for t in stypes:
        if t["_sbn"]:
            sb_types[t["_sbn"]].append(t["id"])
    mmc2_sb_by_name = collections.defaultdict(set)
    mmc2_by_name = collections.defaultdict(list)
    for r in mmc2:
        if r["phase_name"]:
            mmc2_by_name[r["phase_name"]].append(r)
            if r["sb"]:
                mmc2_sb_by_name[r["phase_name"]].add(str(r["sb"]).replace("_", "").upper())

    def mmc2_rows_for(name, elements):
        """Hallstedt rows for this phase name that concern this system: a listed compound with exactly the system's
        elements, or (generic names) every row when they all agree on Pearson symbol and space group."""
        rows = mmc2_by_name.get(name, [])
        els = set(elements)
        hit = []
        for r in rows:
            for comp in re.split(r",\s*", r["actual"] or ""):
                c = formula_counts(comp.split("(")[0].strip().replace(" ", ""))
                if c and {e.upper() for e in c} == els:
                    hit.append(r)
                    break
        if hit:
            return hit
        if rows and len({(r["pearson"], r["sg_number"]) for r in rows}) == 1:
            return rows
        return []

    def type_for_row(r):
        sbn = sb_norm(r["sb"])
        for t in stypes:
            if sbn and t["_sbn"] == sbn and t["pearson"] == r["pearson"] and (t["space_group"] or {}).get("number") == r["sg_number"]:
                return t["id"]
        pk = proto_key(r["prototype"] or "?")
        for t in stypes:
            if not t["_sbn"] and proto_key(t["prototype"] or "?") == pk and t["pearson"] == r["pearson"] and \
                    (t["space_group"] or {}).get("number") == r["sg_number"]:
                return t["id"]
        return None

    # ---- catalog + TDB files
    catalog = json.load(open(inputs.add("open-tdb/catalog/catalog.json", os.path.join(open_tdb, "catalog", "catalog.json"))))
    audit = json.load(open(inputs.add("open-tdb/catalog/model-audit/bundled-summary.json",
                                      os.path.join(open_tdb, "catalog", "model-audit", "bundled-summary.json"))))
    audit_sha = {e["path"]: e["sha256"] for e in audit["entries"]}
    links = json.load(open(inputs.add("sgte-corpus/tdb_pdf_links.json", os.path.join(a.corpus, "tdb_pdf_links.json"))))
    inputs.add("orchid-project.org/tools/doi_corrections.json", DOI_CORR.path)
    link_by_tdb = {}
    for l in links.get("links", []):
        if l.get("doi"):
            link_by_tdb.setdefault(l["tdb"], l)
    tdbs = []
    for e in catalog:
        if e["tier"] == "bundled" and e["license_status"] == "confirmed" and e["n_elements"] == 2 and e["local_path"]:
            tdbs.append(e)
    tdbs.sort(key=lambda e: e["local_path"])
    tdb_digest = hashlib.sha256()
    # where each TDB file is published: the website copy (same bytes) and the orchid-pipeline GitHub blob (same git
    # blob id on the recorded ref, read from the local clone; no network)
    pipe_blobs, pipe_tree = {}, None
    if a.pipeline and os.path.isdir(os.path.join(a.pipeline, ".git")):
        r = subprocess.run(["git", "-C", a.pipeline, "ls-tree", "-r", a.pipeline_ref, "databases"],
                           capture_output=True, text=True)
        for line in r.stdout.splitlines():
            meta, _, rel = line.partition("\t")
            if meta.split()[1:2] == ["blob"]:
                pipe_blobs[rel] = meta.split()[2]
        pipe_tree = subprocess.run(["git", "-C", a.pipeline, "rev-parse", a.pipeline_ref + ":databases"],
                                   capture_output=True, text=True).stdout.strip() or None
    if pipe_tree:
        inputs.hashes["orchid-pipeline %s:databases (git tree id)" % a.pipeline_ref] = pipe_tree
    tdb_stats = collections.Counter()
    tdb_codes_undefined = collections.Counter()

    def file_links(rel, data_bytes, sha):
        out = {"website_url": None, "github_url": None}
        site_copy = os.path.join(orchid_cite.SITE, "data", rel)
        if os.path.isfile(site_copy) and sha256_file(site_copy) == sha:
            out["website_url"] = SITE + "/data/" + rel
        blob = hashlib.sha1(b"blob %d\0" % len(data_bytes) + data_bytes).hexdigest()
        if pipe_blobs.get(rel) == blob:
            out["github_url"] = "%s/blob/%s/%s" % (PIPELINE_GITHUB, a.pipeline_ref.split("/")[-1], rel)
        return out

    # ---- calculated diagrams
    pd_dir = os.path.join(data, "phase-diagrams", "sgte")
    diagrams = {}
    pd_digest = hashlib.sha256()
    for fp in sorted(glob.glob(os.path.join(pd_dir, "*.json"))):
        pd_digest.update(("%s  %s\n" % (sha256_file(fp), os.path.basename(fp))).encode())
        d = json.load(open(fp))
        diagrams[d["tdb"]] = (os.path.basename(fp), d)
    inputs.hashes["orchid-data/phase-diagrams/sgte/*.json (sha256 of 'sha256  name' lines)"] = pd_digest.hexdigest()

    papers = json.load(open(inputs.add("orchid-data/literature/papers.json", os.path.join(data, "literature", "papers.json"))))

    entries = {}

    def entry(pid, system, elements):
        if pid not in entries:
            entries[pid] = {"id": pid, "system": system, "elements": elements, "models": [], "names": {},
                            "stability": [], "lit": collections.Counter(), "lit_acc": collections.Counter(),
                            "refs": {}, "corpus": {"in_corpus": False, "aliases": set(), "mappings": collections.Counter()},
                            "lit_struct": collections.Counter(), "dft_names": []}
        return entries[pid]

    def add_name(e, name, ns, **kw):
        e["names"].setdefault(name, dict({"name": name, "namespace": ns}, **kw))

    system_phase_names = collections.defaultdict(set)   # SYS -> TDB phase names (for literature linking)
    tdb_models_by_file = {}
    for e in tdbs:
        path = os.path.join(open_tdb, e["local_path"])
        sha = sha256_file(path)
        tdb_digest.update(("%s  %s\n" % (sha, e["local_path"])).encode())
        if e["local_path"] in audit_sha and audit_sha[e["local_path"]] != sha:
            raise SystemExit("TDB changed since the model audit: " + e["local_path"])
        elements = sorted(x.upper() for x in e["elements"].split("-"))
        system = "-".join(elements)
        phases, ordis, header = parse_tdb(path)
        fname = os.path.basename(e["local_path"])
        defs, dstats = tdb_definitions(path, open_tdb, e["local_path"], sha,
                                       file_links(e["local_path"], open(path, "rb").read(), sha))
        if set(defs) != set(phases):
            raise SystemExit("PHASE statements disagree between parsers in " + e["local_path"])
        tdb_stats["unmatched_parameters"] += dstats["unmatched_parameters"]
        for c in dstats["type_codes_without_definition"]:
            tdb_codes_undefined[c] += 1
        tdb_models_by_file[fname] = (system, phases)
        dis_of = {v: k for k, v in ordis.items()}
        for name in sorted(phases):
            p = phases[name]
            pid = canonical_phase_id(elements, name, phases.keys(), mmc2_sb_by_name)
            en = entry(pid, system, elements)
            system_phase_names[system].add(name)
            comp = model_composition(p, elements)
            m = {"tdb": fname, "path": e["local_path"], "sha256": sha, "phase_name": name,
                 "type_code": p["type_code"],
                 "sublattices": [{"site_ratio": r, "constituents": cs}
                                 for r, cs in zip(p["site_ratios"], p["constituents"] or [[]] * len(p["site_ratios"]))],
                 "stoichiometric": comp is not None,
                 "composition": {k: float(v / sum(comp.values())) for k, v in sorted(comp.items())} if comp else None}
            if name in ordis:
                m["disordered_part"] = ordis[name]
            if name in dis_of:
                m["disordered_part_of"] = dis_of[name]
            td = defs[name]
            td["parsed"]["disordered_part"] = ordis.get(name)
            td["parsed"]["ordered_partner"] = dis_of.get(name)
            m["tdb_definition"] = td
            protos = [c for c in p["comments"] if re.match(r"(?i)^prototypes?\b", c)]
            if protos:
                m["tdb_comment_structure"] = protos
            m["_raw"] = p
            m["_comp"] = comp
            en["models"].append(m)
            add_name(en, "tdb:%s:%s" % (fname, name), "tdb", label=name, tdb_sha256=sha)
            # assessment reference
            dg = diagrams.get(fname)
            osrc = (dg[1].get("original_source") if dg else None) or {}
            doi = osrc.get("doi") or (None if osrc.get("doi_status") else (link_by_tdb.get(fname) or {}).get("doi"))
            if doi and DOI_CORR.tdb(fname, doi):   # TDB header citation -> DOI mismatch (tools/doi_corrections.json)
                doi = None
            if doi:
                cit = (papers.get("doi:" + doi) or {}).get("citation") or osrc.get("citation") or \
                      (link_by_tdb.get(fname) or {}).get("header_citation")
                r = en["refs"].setdefault(doi.lower(), {"doi": doi, "citation": cit, "roles": set(), "labels": set(),
                                                        "values": 0, "tdb_files": set()})
                r["roles"].add("assessment that models this phase")
                r["labels"].add(name)
                r["tdb_files"].add(fname)
    inputs.hashes["open-tdb/databases (sha256 of 'sha256  path' lines, bundled confirmed binaries)"] = tdb_digest.hexdigest()

    # ---- stability from calculated diagrams
    for fname, (system, phases) in sorted(tdb_models_by_file.items()):
        if fname not in diagrams:
            continue
        pfile, d = diagrams[fname]
        elements = sorted(x.upper() for x in d["elements"])
        pts = collections.defaultdict(list)
        for seg in d.get("regions") or []:
            for r in seg:
                for x, T in zip(r.get("x") or [], r.get("T") or []):
                    pts[r["phase"].upper()].append((x, T))
        for inv in d.get("invariants") or []:
            for p in inv.get("points") or []:
                pts[p["phase"].upper()].append((p["x"], inv["T"]))
        g = d.get("grid") or {}
        trow = collections.defaultdict(list)
        if g.get("rows"):
            labels = g["labels"]
            for i, row in enumerate(g["rows"]):
                T = g["T0"] + i * g["dT"]
                xi = 0
                for li, cnt in row:
                    lab = labels[li]
                    ph = [q for q in lab.split("+") if q]
                    for q in ph:
                        trow[q.upper()].append(T)
                    if len(ph) == 1:
                        for k in range(cnt):
                            pts[ph[0].upper()].append((round((xi + k) * g["x_step"], 6), T))
                    xi += cnt
        tr = d.get("T_range") or [None, None]
        for ph in sorted(set(phases) | {q.upper() for q in d.get("phases") or []}):
            pid = canonical_phase_id(elements, ph, phases.keys(), mmc2_sb_by_name)
            if pid not in entries:
                continue
            P, Ts = pts.get(ph, []), trow.get(ph, [])
            st = {"tdb": fname, "tdb_sha256": d.get("sha256"), "diagram": "phase-diagrams/sgte/" + pfile,
                  "status": d.get("status"), "x_element": d.get("x_component")}
            if tr != [300, 2500]:
                st["calc_T_range_K"] = tr
            if d.get("status") != "ok":
                st["stable"] = None
            elif not P and not Ts:
                st["stable"] = False
            else:
                allT = [t for _, t in P] + Ts
                xs = [x for x, _ in P]
                st.update(stable=True, T_min_K=round(min(allT), 2), T_max_K=round(max(allT), 2),
                          x_min=round(min(xs), 5) if xs else None, x_max=round(max(xs), 5) if xs else None,
                          stable_at_top_of_range=bool(tr[1] and max(allT) >= tr[1] - g.get("dT", 10)))
            entries[pid]["stability"].append(st)

    # ---- corpus canonical ids + mappings
    cal = json.load(open(inputs.add("sgte-corpus/db-v2.1/canonical_phase_aliases.json",
                                    os.path.join(corpus_db, "canonical_phase_aliases.json"))))
    corpus_only = 0
    tok_index = {}
    for pid0 in sorted(entries):
        sys0, sl0 = pid0.split(":")[1], pid0.split(":")[2]
        tok_index.setdefault((sys0, tuple(sorted(sl0.split("-")))), pid0)
    redirect = {}
    for old, cid in sorted(cal["aliases"].items()):
        system_l = cid.split(":")[1]
        els = system_l.split("-")
        if len(els) != 2:
            continue
        if cid not in entries:
            # same name tokens in another order as a registry phase (B2-BCC = BCC_B2): an alias of that phase
            tgt = tok_index.get((system_l, tuple(sorted(cid.split(":")[2].split("-")))))
            if tgt:
                redirect[cid] = tgt
                entries[tgt]["corpus"]["aliases"].add(cid)
                if old != cid:
                    entries[tgt]["corpus"]["aliases"].add(old)
                entries[tgt]["corpus"]["in_corpus"] = True
                continue
            entry(cid, "-".join(e.upper() for e in els), sorted(e.upper() for e in els))
            corpus_only += 1
        en = entries[cid]
        en["corpus"]["in_corpus"] = True
        if old != cid:
            en["corpus"]["aliases"].add(old)
    map_digest = hashlib.sha256()
    n_map = 0
    for f in sorted(glob.glob(os.path.join(corpus_db, "invariants", "master", "*.jsonl")) +
                    glob.glob(os.path.join(corpus_db, "general", "records", "*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                if '"record_type": "phase_mapping"' not in line:
                    continue
                r = json.loads(line)
                cid = r.get("canonical_phase")
                cid = cal["aliases"].get(cid, cid)
                cid = redirect.get(cid, cid)
                if cid in entries:
                    entries[cid]["corpus"]["mappings"][(r.get("basis"), r.get("confidence"))] += 1
                    entries[cid]["corpus"].setdefault("labels", collections.Counter())[r.get("printed_label")] += 1
                    n_map += 1
                    map_digest.update(line.strip().encode() + b"\n")
    inputs.hashes["sgte-corpus phase_mapping records (sha256 of the linked lines)"] = map_digest.hexdigest()

    # ---- literature
    lit_digest = hashlib.sha256()
    lit_stats = collections.Counter()
    for fp in sorted(glob.glob(os.path.join(data, "literature", "systems", "*.json"))):
        lit_digest.update(("%s  %s\n" % (sha256_file(fp), os.path.basename(fp))).encode())
        d = json.load(open(fp))
        els = [x.upper() for x in d.get("elements") or d["system"].split("-")]
        if len(els) != 2:
            continue
        system = "-".join(sorted(els))
        names = system_phase_names.get(system, set())
        for r in d["records"]:
            seen = {}
            for ph in r.get("phases") or []:
                lab, mp = ph.get("label"), ph.get("model_phase")
                pid = None
                if mp:
                    pid = canonical_phase_id(els, mp, names, mmc2_sb_by_name)
                    pid = cal["aliases"].get(pid, pid)
                    pid = redirect.get(pid, pid)
                    basis = "corpus_mapping:" + (ph.get("mapping_confidence") or "unstated")
                elif lab and str(lab).strip().upper() in names:
                    pid = canonical_phase_id(els, str(lab).strip().upper(), names, mmc2_sb_by_name)
                    basis = "exact_tdb_name"
                if not pid or pid not in entries:
                    lit_stats["phase_mentions_unlinked"] += 1
                    continue
                lit_stats["phase_mentions_linked"] += 1
                seen[pid] = (lab, basis)
                if ph.get("pearson_symbol") or ph.get("prototype"):
                    entries[pid]["lit_struct"][(ph.get("prototype"), ph.get("pearson_symbol"))] += 1
            acc = bool((r.get("flags") or {}).get("extraction_accepted"))
            paper = r.get("paper") or ""
            for pid, (lab, basis) in seen.items():
                en = entries[pid]
                en["lit"]["values"] += 1
                if acc:
                    en["lit"]["accepted"] += 1
                en["lit_acc"][paper] += 1
                if paper.startswith("doi:") and paper in papers:
                    doi = papers[paper]["doi"]
                    add_name(en, "paper:%s:%s" % (doi, lab), "paper", label=lab, basis=basis)
                    rr = en["refs"].setdefault(doi.lower(), {"doi": doi, "citation": papers[paper].get("citation"),
                                                             "roles": set(), "labels": set(), "values": 0,
                                                             "tdb_files": set()})
                    rr["roles"].add("reports values for this phase")
                    rr["labels"].add(lab)
                    rr["values"] += 1
                elif paper:
                    add_name(en, "paper:%s:%s" % (paper.replace("sha256:", "sha256-"), lab), "paper", label=lab, basis=basis)
                for c in r.get("cited_sources") or []:
                    cd = c.get("doi")
                    if cd and ("doi:" + cd) in papers:
                        rr = en["refs"].setdefault(cd.lower(), {"doi": cd, "citation": papers["doi:" + cd].get("citation"),
                                                                "roles": set(), "labels": set(), "values": 0,
                                                                "tdb_files": set()})
                        rr["roles"].add("cited as the source of values for this phase")
                        rr["values"] += 1
    inputs.hashes["orchid-data/literature/systems/*.json (sha256 of 'sha256  name' lines)"] = lit_digest.hexdigest()

    # ---- DFT
    hulls = load_hull_vertices(open_tdb, inputs)

    def dft_at(system, comp, src):
        h = hulls.get((system, src))
        if not h or not comp:
            return None, []
        tot = sum(comp.values())
        x = float(Fraction(comp.get(h["x_element"], 0)) / tot)
        hits = [v for v in h["vertices"] if abs(v["x"] - x) < TOL]
        return x, hits

    def dft_in_range(system, lo, hi, src):
        h = hulls.get((system, src))
        if not h:
            return []
        return [v for v in h["vertices"] if lo - TOL <= v["x"] <= hi + TOL and 0 < v["x"] < 1]

    # ---- per-system composition index (polymorph families)
    comp_index = collections.defaultdict(set)
    for pid, en in entries.items():
        for m in en["models"]:
            if m["_comp"] and len(m["_comp"]) == 2:
                tot = sum(m["_comp"].values())
                key = tuple(sorted((k, v / tot) for k, v in m["_comp"].items()))
                comp_index[(en["system"], key)].add(pid)

    # ---- assignment
    conflicts = {"name_dft_disagreement": [], "polymorph_families": [], "model_composition_disagreement": [],
                 "same_name_different_structures": [], "mp_oqmd_disagreement": [], "assessment_declared_blocked": []}
    for key, pids in sorted(comp_index.items()):
        if len(pids) > 1:
            conflicts["polymorph_families"].append({"system": key[0], "composition": {k: float(v) for k, v in key[1]},
                                                    "phases": sorted(pids)})

    def candidate_from_type(tid, source, detail, **kw):
        t = st_by_id.get(tid) if tid else None
        c = {"source": source, "structure_type": tid, "prototype": t["prototype"] if t else kw.get("prototype"),
             "strukturbericht": t["strukturbericht"] if t else kw.get("strukturbericht"),
             "pearson": t["pearson"] if t else kw.get("pearson"),
             "space_group": (t["space_group"] if t else kw.get("space_group")), "detail": detail}
        for k in ("material_id", "count", "doi"):
            if kw.get(k) is not None:
                c[k] = kw[k]
        return c

    def match_types(sgn, ratio, nsites_list):
        hits = []
        for t in stypes:
            if not t["space_group"] or t["space_group"]["number"] != sgn or t["_ratio"] != ratio:
                continue
            oks = [nsites_ok(t["pearson"], n) for n in nsites_list]
            if any(o is False for o in oks):
                continue
            hits.append(t["id"])
        return hits

    C_CONFLICT = ("stoichiometry_contradicts_", "site_ratios_contradict_", "order_disorder_ambiguity", "lattice_word_",
                  "word_", "crystal_word_", "several_structure_tokens", "names_imply_different_structures",
                  "sgte_mmc2_row_disagrees", "name_suggests_", "models_disagree_on_composition", "words_imply_several")

    def mmc2_system_rows(name, elements):
        """Hallstedt rows for this phase name whose listed compound has exactly this system's elements."""
        els, hit = set(elements), []
        for r in mmc2_by_name.get(name, []):
            for c in re.split(r",\s*", r["actual"] or ""):
                cc = formula_counts(c.split("(")[0].strip().replace(" ", ""))
                if cc and {x.upper() for x in cc} == els:
                    hit.append((r, cc))
                    break
        return hit

    def tier_c(en, tdb_names, analyses, reasons, dft, comp, polymorph):
        """Tier C: the SGTE crystallographic table (mmc2) states exactly one structure for this phase in this system,
        nothing in tiers A/B contradicts it, and the TDB model is consistent with the listed compound."""
        rows = [(n, r, cc) for n in tdb_names for r, cc in mmc2_system_rows(n, en["elements"])]
        if not rows:
            return {}
        brief = [{"phase_name": n, "actual": r["actual"], "model": r["model"], "prototype": r["prototype"],
                  "strukturbericht": r["sb"], "pearson": r["pearson"],
                  "space_group": {"number": r["sg_number"], "symbol": r["sg_symbol"]} if r["sg_number"] else None}
                 for n, r, cc in rows]
        brief = [json.loads(x) for x in sorted({json.dumps(b, sort_keys=True) for b in brief})]
        res = {"rows_brief": brief}
        tids = {type_for_row(r) for _, r, _ in rows}
        keys = {(r["pearson"], r["sg_number"]) for _, r, _ in rows}
        if len(keys) > 1 or len(tids) > 1:
            return dict(res, block="sgte_mmc2_rows_disagree")
        if None in tids or not all(r["sg_number"] and r["pearson"] for _, r, _ in rows):
            return dict(res, block="sgte_mmc2_row_incomplete")
        tid = tids.pop()
        t = st_by_id[tid]
        bad = [x for x in reasons if x.startswith(C_CONFLICT)]
        if bad:
            return dict(res, block="assessment_declared_blocked_by:" + bad[0].split(":")[0])
        sbns = {s for x in analyses for s in (x.get("candidate_sbn") or [])}
        if sbns and t["_sbn"] not in sbns:
            return dict(res, block="assessment_declared_conflicts_with_name")
        if comp and not polymorph and dft.get("mp") and dft.get("oqmd"):
            sgs = {v["spacegroup_number"] for sx in ("mp", "oqmd") for v in dft[sx]}
            if len(dft["mp"]) == 1 and len(dft["oqmd"]) == 1 and len(sgs) == 1 and not dft["mp"][0]["degenerate"] \
                    and not dft["oqmd"][0]["degenerate"] and t["space_group"] and sgs.pop() != t["space_group"]["number"]:
                return dict(res, block="assessment_declared_conflicts_with_dft")
        cks = {ratio_key(cc) for _, _, cc in rows}
        for m in en["models"]:
            p = m["_raw"]
            if t["_sbn"] in DISORDERED and partitioned(p):
                return dict(res, block="assessment_declared_order_disorder_ambiguity")
            if m["_comp"] and len(m["_comp"]) == 2:
                if ratio_key(m["_comp"]) not in cks:
                    return dict(res, block="assessment_declared_compound_differs_from_model_composition")
            elif not partitioned(p):
                sk = site_ratio_key(p)
                if sk is not None and sk not in cks:
                    return dict(res, block="assessment_declared_compound_differs_from_site_ratios")
        res["struct"] = {"structure_type": tid, "basis": "assessment_declared", "tier": "C",
                         "confidence": "assessor-stated", "label": "assessment-declared (Hallstedt 2025, mmc2)",
                         "source": {"citation": HALLSTEDT["citation"], "doi": HALLSTEDT["doi"],
                                    "license": HALLSTEDT["license"], "table": "mmc2", "rows": brief}}
        return res

    name_struct = collections.defaultdict(set)   # TDB phase name -> structure types assigned/declared (conflicts)
    for pid in sorted(entries):
        en = entries[pid]
        reasons, cands, caveats = [], [], []
        struct = None
        tdb_names = sorted({m["phase_name"] for m in en["models"]})
        comps = {}
        for m in en["models"]:
            if m["_comp"]:
                tot = sum(m["_comp"].values())
                comps[tuple(sorted((k, v / tot) for k, v in m["_comp"].items()))] = m["_comp"]
        all_stoich = bool(en["models"]) and all(m["_comp"] for m in en["models"])
        if len(comps) > 1:
            reasons.append("models_disagree_on_composition")
            conflicts["model_composition_disagreement"].append(
                {"phase": pid, "compositions": [{k: float(v) for k, v in c} for c in sorted(comps)],
                 "tdb_files": sorted({m["tdb"] for m in en["models"]})})
        comp = next(iter(comps.values())) if len(comps) == 1 and all_stoich else None
        compkey = next(iter(comps)) if comp else None
        polymorph = bool(comp) and len(comp_index.get((en["system"], compkey), ())) > 1
        # name analysis over every TDB name of this phase (they share a canonical id, i.e. the same name tokens)
        analyses = [classify_name(n, en["elements"], sb_types, None, st_by_id) for n in tdb_names] or \
                   [classify_name(pid.split(":")[2].upper().replace("-", "_"), en["elements"], sb_types, None, st_by_id)]
        kinds = {x["kind"] for x in analyses}
        # DFT at the fixed composition
        dft = {}
        if comp and len(comp) == 2:
            for src in ("mp", "oqmd"):
                x, hits = dft_at(en["system"], comp, src)
                dft[src] = hits
        if kinds == {"liquid"}:
            struct = {"structure_type": "st:liquid", "basis": "liquid", "confidence": "high", "label": "liquid"}
        elif kinds == {"gas"}:
            struct = {"structure_type": "st:gas", "basis": "gas", "confidence": "high", "label": "gas"}
        else:
            an = analyses[0]
            tids = {x.get("type_id") for x in analyses}
            if len(tids) > 1:
                reasons.append("names_imply_different_structures")
            for x in analyses:
                reasons += x.get("reasons") or []
                for sbn in x.get("candidate_sbn") or []:
                    for tid in sb_types.get(sbn, []):
                        cands.append(candidate_from_type(tid, "name_token", "Strukturbericht/word '%s' in the TDB name" % sbn))
            tid = an.get("type_id") if len(tids) == 1 else None
            if tid:
                t = st_by_id[tid]
                ok = not reasons
                # stoichiometry / sublattices vs. structure type
                for m in en["models"]:
                    p = m["_raw"]
                    if t["_sbn"] in DISORDERED and partitioned(p):
                        reasons.append("order_disorder_ambiguity:%s_model_is_partitioned" % m["tdb"])
                    elif t["_ratio"]:
                        if m["_comp"] and len(m["_comp"]) == 2:
                            if ratio_key(m["_comp"]) != t["_ratio"]:
                                reasons.append("stoichiometry_contradicts_%s:%s" % (t["strukturbericht"], m["tdb"]))
                        elif partitioned(p):
                            caveats.append("partitioned order/disorder model in %s also describes the disordered "
                                           "state" % m["tdb"])
                        else:
                            sk = site_ratio_key(p)
                            if sk is not None and sk != t["_ratio"]:
                                reasons.append("site_ratios_contradict_%s:%s" % (t["strukturbericht"], m["tdb"]))
                    if m.get("disordered_part") and t["_sbn"] not in DISORDERED:
                        caveats.append("%s: ordered phase with disordered part %s (order/disorder model)" %
                                       (m["tdb"], m["disordered_part"]))
                # the SGTE crystallographic table must not contradict the name
                for n in tdb_names:
                    for r in mmc2_rows_for(n, en["elements"]):
                        if r["sg_number"] and t["space_group"] and r["sg_number"] != t["space_group"]["number"]:
                            reasons.append("sgte_mmc2_row_disagrees_with_name:%s" % n)
                # DFT must not contradict a named compound structure
                if comp and not polymorph and dft.get("mp") and dft.get("oqmd"):
                    sgs = {v["spacegroup_number"] for s in ("mp", "oqmd") for v in dft[s]}
                    if len(dft["mp"]) == 1 and len(dft["oqmd"]) == 1 and len(sgs) == 1 and \
                            not dft["mp"][0]["degenerate"] and not dft["oqmd"][0]["degenerate"]:
                        sg = sgs.pop()
                        if t["space_group"] and sg != t["space_group"]["number"]:
                            reasons.append("name_suggests_%s_but_dft_space_group_%s" % (t["strukturbericht"], sg))
                            conflicts["name_dft_disagreement"].append(
                                {"phase": pid, "names": tdb_names, "name_structure": tid,
                                 "name_space_group": t["space_group"]["number"], "dft_space_group": sg,
                                 "mp": dft["mp"][0]["phase"], "oqmd": dft["oqmd"][0]["phase"]})
                        else:
                            caveats.append("DFT agrees: on-hull %s and %s, space group %s" %
                                           (dft["mp"][0]["phase"], dft["oqmd"][0]["phase"], sg))
                            en["dft_names"] = [dft["mp"][0]["phase"], dft["oqmd"][0]["phase"]]
                if not reasons:
                    struct = {"structure_type": tid, "basis": "standard_name", "tier": "A", "confidence": "high",
                              "label": "from the standard phase name"}
            elif not an.get("candidate_sbn") and not any(x.get("candidate_sbn") for x in analyses):
                # tier B: no structure in the name
                desc = sorted({u for x in analyses for u in x.get("unknown") or []})
                if desc:
                    reasons.append("name_has_descriptive_tokens:" + ",".join(desc))
                elif not en["models"]:
                    reasons.append("no_bundled_tdb_model")
                elif not all_stoich or not comp:
                    if not all_stoich:
                        reasons.append("off_stoichiometric_model")
                elif len(comp) != 2:
                    reasons.append("unary_composition")
                elif polymorph:
                    reasons.append("polymorph_family:other phases of this system have the same composition")
                else:
                    mp, oq = dft.get("mp") or [], dft.get("oqmd") or []
                    if not mp or not oq:
                        reasons.append("no_on_hull_entry_in_%s" % ("both" if not mp and not oq else ("mp" if not mp else "oqmd")))
                    elif len(mp) > 1 or len(oq) > 1 or mp[0]["degenerate"] or oq[0]["degenerate"]:
                        reasons.append("dft_degenerate_or_multiple_on_hull_entries")
                    elif mp[0]["spacegroup_number"] != oq[0]["spacegroup_number"]:
                        reasons.append("mp_oqmd_space_groups_differ")
                        conflicts["mp_oqmd_disagreement"].append(
                            {"phase": pid, "names": tdb_names, "mp": mp[0]["phase"], "mp_space_group": mp[0]["spacegroup_number"],
                             "oqmd": oq[0]["phase"], "oqmd_space_group": oq[0]["spacegroup_number"]})
                    else:
                        sgn = mp[0]["spacegroup_number"]
                        hits = match_types(sgn, ratio_key(comp), [mp[0].get("nsites"), oq[0].get("nsites")])
                        mm = [r for n in tdb_names for r in mmc2_rows_for(n, en["elements"]) if r["sg_number"]]
                        if any(r["sg_number"] != sgn for r in mm):
                            reasons.append("sgte_mmc2_row_disagrees_with_dft")
                        if len(hits) == 1 and "sgte_mmc2_row_disagrees_with_dft" not in reasons:
                            struct = {"structure_type": hits[0], "basis": "dft_unique_match", "tier": "B",
                                      "confidence": "high", "label": "DFT-supported",
                                      "dft": {"mp": mp[0]["phase"], "oqmd": oq[0]["phase"], "space_group_number": sgn,
                                              "space_group_symbol": mp[0]["spacegroup_symbol"]}}
                            en["dft_names"] = [mp[0]["phase"], oq[0]["phase"]]
                        elif len(hits) != 1:
                            reasons.append("dft_space_group_%s_matches_%d_structure_types" % (sgn, len(hits)))
                            for h in hits[:6]:
                                cands.append(candidate_from_type(h, "dft_space_group_and_stoichiometry",
                                                                 "space group %s and stoichiometry match" % sgn))
            else:
                if not reasons:
                    reasons.append("name_structure_not_unique")
        # ---- tier C: assessment-declared (Hallstedt 2025 mmc2 lists one structure for this phase in this system)
        if not struct and kinds == {"crystal"} and en["models"]:
            tc = tier_c(en, tdb_names, analyses, reasons, dft, comp, polymorph)
            if tc.get("struct"):
                struct = tc["struct"]
                caveats.append("not certain from the name or DFT alone (%s); structure as stated in the SGTE "
                               "crystallographic table for this phase name and system" % (", ".join(sorted(set(reasons))) or "no name/DFT basis"))
            elif tc.get("block"):
                reasons.append(tc["block"])
                conflicts["assessment_declared_blocked"].append(
                    {"phase": pid, "names": tdb_names, "reason": tc["block"], "mmc2": tc.get("rows_brief", [])})
        # ---- candidates (evidence for review, recorded for every crystalline phase without an assignment)
        if not struct:
            for n in tdb_names:
                for r in mmc2_rows_for(n, en["elements"]):
                    cands.append(candidate_from_type(type_for_row(r), "sgte_mmc2",
                                                     "Hallstedt 2025 mmc2 row for %s: %s" % (n, r["actual"] or r["model"] or ""),
                                                     prototype=r["prototype"], strukturbericht=r["sb"], pearson=r["pearson"],
                                                     space_group={"number": r["sg_number"], "symbol": r["sg_symbol"]}
                                                     if r["sg_number"] else None))
            for m in en["models"]:
                for c in m.get("tdb_comment_structure") or []:
                    mm = re.match(r"(?i)^prototypes?\s+(.+?)\s*\((\w+),\s*([^)]+)\)", c)
                    cands.append({"source": "tdb_comment", "structure_type": None,
                                  "prototype": mm.group(1) if mm else None, "strukturbericht": None,
                                  "pearson": mm.group(2) if mm else None,
                                  "space_group": {"number": None, "symbol": mm.group(3).strip()} if mm else None,
                                  "detail": "%s: '%s'" % (m["tdb"], c)})
            for src in ("mp", "oqmd"):
                vs = dft.get(src)
                if vs is None and en["models"]:
                    lo_hi = None
                    for st in en["stability"]:
                        if st.get("stable") and st.get("x_min") is not None:
                            h = hulls.get((en["system"], src))
                            if h and h["x_element"] == st["x_element"]:
                                lo_hi = (min(st["x_min"], (lo_hi or (9, 0))[0]), max(st["x_max"], (lo_hi or (9, -1))[1]))
                    vs = dft_in_range(en["system"], lo_hi[0], lo_hi[1], src) if lo_hi else []
                for v in (vs or [])[:6]:
                    if not (0 < v["x"] < 1):
                        continue
                    cands.append({"source": "dft_" + src, "structure_type": None, "prototype": None,
                                  "strukturbericht": None, "pearson": None,
                                  "space_group": {"number": v["spacegroup_number"], "symbol": v["spacegroup_symbol"]}
                                  if v.get("spacegroup_number") else None,
                                  "detail": "on-hull %s %s at x(%s)=%.4f%s" % (
                                      src.upper(), v["formula"], hulls[(en["system"], src)]["x_element"], v["x"],
                                      "; degenerate with " + ",".join(dd["phase"] for dd in v["degenerate"]) if v["degenerate"] else ""),
                                  "material_id": v["phase"]})
            for (proto, pear), n in sorted(en["lit_struct"].items(), key=lambda kv: (-kv[1], str(kv[0]))):
                cands.append({"source": "literature_printed", "structure_type": None, "prototype": proto,
                              "strukturbericht": None, "pearson": pear, "space_group": None,
                              "detail": "printed with a label mapped to this phase", "count": n})
        # dedup candidates
        seen, uc = set(), []
        for c in cands:
            k = json.dumps(c, sort_keys=True)
            if k not in seen:
                seen.add(k)
                uc.append(c)
        cands = sorted(uc, key=lambda c: (c["source"], str(c.get("structure_type")), c["detail"]))
        if struct:
            t = st_by_id[struct["structure_type"]]
            struct.update(prototype=t["prototype"], strukturbericht=t["strukturbericht"], pearson=t["pearson"],
                          space_group=t["space_group"])
            if caveats:
                struct["caveats"] = sorted(set(caveats))
            status = "assigned"
            for n in tdb_names:
                name_struct[n].add(struct["structure_type"])
        else:
            status = "needs_review" if cands else "unresolved_no_candidates"
        en["structure"] = struct
        en["status"] = status
        en["reasons"] = sorted(set(reasons)) if not struct else []
        en["candidates"] = cands if not struct else []
        en["comp"] = comp

    # mmc2 view of "same name -> different structures"
    mm_name = collections.defaultdict(set)
    for r in mmc2:
        if r["phase_name"] and (r["pearson"] or r["sg_number"]):
            mm_name[r["phase_name"]].add((r["prototype"], r["pearson"], r["sg_number"]))
    used_names = set()
    for en in entries.values():
        used_names |= {m["phase_name"] for m in en["models"]}
    for n in sorted(used_names):
        protos = mm_name.get(n, set())
        pk = {(p[1], p[2]) for p in protos}
        if len(pk) > 1:
            conflicts["same_name_different_structures"].append(
                {"name": n, "source": "sgte_mmc2", "structures": [{"prototype": p, "pearson": q, "space_group_number": s}
                                                                  for p, q, s in sorted(protos, key=str)],
                 "registry_phases": sorted(pid for pid, en in entries.items() if n in {m["phase_name"] for m in en["models"]})})

    return finalize(entries, stypes, inputs, conflicts, lit_stats, corpus_only, n_map, data,
                    {"unmatched_parameters": tdb_stats["unmatched_parameters"],
                     "type_codes_without_definition": dict(sorted(tdb_codes_undefined.items()))})


def finalize(entries, stypes, inputs, conflicts, lit_stats, corpus_only, n_map, data, tdb_extra):
    cite = orchid_cite.cite_as(DATASET, data)
    reg = []
    for pid in sorted(entries):
        en = entries[pid]
        names = [en["names"][k] for k in sorted(en["names"])]
        names.append({"name": "corpus:" + pid, "namespace": "corpus", "label": pid,
                      "in_sgte_corpus": en["corpus"]["in_corpus"]})
        for al in sorted(en["corpus"]["aliases"]):
            names.append({"name": "corpus:" + al, "namespace": "corpus", "label": al, "alias_of": pid})
        if en["structure"]:
            for mid in en["dft_names"]:
                ns = "mp" if mid.startswith("mp-") else "oqmd"
                names.append({"name": "%s:%s" % (ns, mid), "namespace": ns, "label": mid,
                              "basis": en["structure"]["basis"]})
        tdbn = collections.Counter(m["phase_name"] for m in en["models"])
        pref = sorted(tdbn, key=lambda n: (-tdbn[n], n))[0] if tdbn else pid.split(":")[2].upper()
        refs = []
        for k, r in en["refs"].items():
            refs.append({"doi": r["doi"], "citation": r["citation"], "roles": sorted(r["roles"]),
                         "labels": sorted(str(x) for x in r["labels"]), "values": r["values"],
                         **({"tdb_files": sorted(r["tdb_files"])} if r["tdb_files"] else {})})
        refs.sort(key=lambda r: (0 if "assessment that models this phase" in r["roles"] else 1, -r["values"], r["doi"].lower()))
        models = []
        for m in sorted(en["models"], key=lambda m: (m["tdb"], m["phase_name"])):
            models.append({k: v for k, v in m.items() if not k.startswith("_")})
        lit_papers = sorted(p for p in en["lit_acc"] if p)
        maps = en["corpus"]["mappings"]
        e = {"id": pid, "system": en["system"], "elements": en["elements"], "preferred_name": pref,
             "status": en["status"],
             "structure": en["structure"],
             "candidates": en["candidates"], "review_reasons": en["reasons"],
             "names": names, "models": models,
             "stoichiometry": ({k: float(v / sum(en["comp"].values())) for k, v in sorted(en["comp"].items())}
                               if en["comp"] else None),
             "stability": sorted(en["stability"], key=lambda s: s["tdb"]),
             "literature": {"values": en["lit"]["values"], "accepted_values": en["lit"]["accepted"],
                            "papers": len(lit_papers)},
             "corpus": {"in_sgte_corpus": en["corpus"]["in_corpus"],
                        "phase_mappings": sum(maps.values()),
                        "mapping_basis": {"%s/%s" % k: v for k, v in sorted(maps.items(), key=lambda kv: str(kv[0]))},
                        "printed_labels": [l for l, _ in sorted((en["corpus"].get("labels") or {}).items(),
                                                                key=lambda kv: (-kv[1], str(kv[0])))][:15]},
             "references": refs[:MAX_REFS], "total_papers": len(refs)}
        reg.append(e)

    st_public = [{k: v for k, v in t.items() if not k.startswith("_")} for t in stypes]
    used = collections.Counter((e["structure"] or {}).get("structure_type") for e in reg if e["structure"])
    for t in st_public:
        t["registry_phases"] = used.get(t["id"], 0)
    systems = sorted({e["system"] for e in reg})
    by_basis = collections.Counter((e["structure"] or {}).get("basis") for e in reg)
    status = collections.Counter(e["status"] for e in reg)
    ns = collections.Counter(n["namespace"] for e in reg for n in e["names"])
    reasons = collections.Counter(r.split(":")[0] for e in reg for r in e["review_reasons"])
    counts = {"systems": len(systems), "phases": len(reg), "names": sum(ns.values()), "names_by_namespace": dict(sorted(ns.items())),
              "structure_types": len(st_public), "structure_types_used": len([t for t in st_public if t["registry_phases"]]),
              "assigned_tier_A_standard_name": by_basis.get("standard_name", 0),
              "assigned_tier_B_dft_unique_match": by_basis.get("dft_unique_match", 0),
              "assigned_tier_C_assessment_declared": by_basis.get("assessment_declared", 0),
              "liquid": by_basis.get("liquid", 0), "gas": by_basis.get("gas", 0),
              "no_prototype": sum(1 for e in reg if not e["structure"]),
              "status": dict(sorted(status.items())), "review_queue": status.get("needs_review", 0),
              "review_reasons": dict(sorted(reasons.items())),
              "conflicts": {k: len(v) for k, v in conflicts.items()},
              "tdb_models": sum(len(e["models"]) for e in reg),
              "phases_with_models": sum(1 for e in reg if e["models"]),
              "corpus_only_phases": corpus_only, "corpus_phase_mappings_linked": n_map,
              "phases_with_literature_values": sum(1 for e in reg if e["literature"]["values"]),
              "literature_phase_mentions_linked": lit_stats["phase_mentions_linked"],
              "literature_phase_mentions_unlinked": lit_stats["phase_mentions_unlinked"],
              "phases_with_references": sum(1 for e in reg if e["references"]),
              "references_listed": sum(len(e["references"]) for e in reg),
              "distinct_reference_dois": len({r["doi"].lower() for e in reg for r in e["references"]}),
              "phases_with_stability": sum(1 for e in reg if any(s.get("stable") for s in e["stability"])),
              "tdb_definitions": tdb_def_counts(reg, tdb_extra)}
    head = {"dataset": "ORCHID phase registry (binary systems)", "registry_version": REGISTRY_VERSION,
            "release": orchid_cite.RELEASE_TAG, "license": "CC-BY-4.0",
            "description": "One entry per (binary system, phase). Structure types are assigned only when certain "
                           "(tier A: standard phase name; tier B: unique MP+OQMD on-hull match; tier C: one structure stated by "
                           "the assessment's crystallographic table, assessor-stated); otherwise candidates "
                           "are listed and the phase is marked needs_review. Field reference: SCHEMA.md.",
            "id_rule": "phase:<sorted lower-case elements>:<slug>, the SGTE-corpus canonical id convention",
            "stability_rule": "calculated diagrams: pycalphad 0.11.2, 300-2500 K unless calc_T_range_K is given, "
                              "1 atm, 0.01 composition grid (orchid-data phase-diagrams/sgte)",
            "inputs": dict(sorted(inputs.hashes.items())), "counts": counts, "cite_as": cite}
    out = {}
    out["structure_types.json"] = dump(dict(head, counts={"structure_types": len(st_public)}, structure_types=st_public))
    out["structure_types.csv"] = structure_types_csv(st_public)
    # registry.json is a compact index (kept well under GitHub's 50 MB file limit); the full records, including
    # models[].tdb_definition and references, are in systems/<SYS>.json
    out["registry.json"] = json.dumps(dict(head, record_format="compact", detail="systems/<SYS>.json",
                                           phases=[compact_entry(e) for e in reg]),
                                      ensure_ascii=False, separators=(",", ":")) + "\n"
    queue = [{"id": e["id"], "system": e["system"], "preferred_name": e["preferred_name"],
              "tdb_names": sorted({m["phase_name"] for m in e["models"]}), "stoichiometry": e["stoichiometry"],
              "review_reasons": e["review_reasons"], "candidates": e["candidates"],
              "literature_values": e["literature"]["values"]}
             for e in reg if e["status"] == "needs_review"]
    out["review_queue.json"] = dump({"dataset": head["dataset"], "registry_version": REGISTRY_VERSION,
                                     "description": "Phases without an assigned structure type that have candidates. "
                                                    "A decision records structure_type, reviewer identity and basis.",
                                     "n": len(queue), "reasons": counts["review_reasons"], "phases": queue,
                                     "cite_as": cite})
    out["CONFLICTS.md"] = conflicts_md(conflicts, counts)
    out["registry.csv"] = registry_csv(reg)
    ctx, ld = jsonld(reg, cite)
    out["context.jsonld"] = dump(ctx)
    out["registry.jsonld"] = json.dumps(ld, ensure_ascii=False, separators=(",", ":")) + "\n"
    # website: compact index + per-system files
    idx = []
    spell = {t["id"]: t["prototype_spellings"] for t in st_public}
    for e in reg:
        s = e["structure"] or {}
        search = sorted({n["label"] for n in e["names"] if n["namespace"] in ("tdb", "paper") and n.get("label")}, key=str)
        alt = []
        for pr in [s.get("prototype")] + spell.get(s.get("structure_type"), []):
            parts = re.findall(r"([A-Z][a-z]?)(\d*)", re.sub(r"^\(|\)$", "", pr or ""))
            if len(parts) == 2 and "".join(a + b for a, b in parts) == re.sub(r"^\(|\)$", "", pr):
                alt.append("".join(a + b for a, b in reversed(parts)))
        search = spell.get(s.get("structure_type"), []) + alt + search
        idx.append([e["id"], e["system"], e["preferred_name"], s.get("structure_type"), s.get("prototype"),
                    s.get("strukturbericht"), s.get("pearson"), s.get("basis"), e["status"],
                    e["literature"]["values"], len(e["models"]), " ".join(str(x) for x in search)[:400]])
    out["index.json"] = json.dumps({"release": orchid_cite.RELEASE_TAG, "registry_version": REGISTRY_VERSION,
                                    "columns": ["id", "system", "name", "structure_type", "prototype",
                                                "strukturbericht", "pearson", "basis", "status", "values", "models",
                                                "search"],
                                    "counts": counts, "phases": idx, "cite_as": cite},
                                   ensure_ascii=False, separators=(",", ":")) + "\n"
    stl = {t["id"]: t for t in st_public}
    for sysn in systems:
        ph = [e for e in reg if e["system"] == sysn]
        out["systems/%s.json" % sysn] = json.dumps(
            {"system": sysn, "release": orchid_cite.RELEASE_TAG, "phases": ph,
             "structure_types": {k: stl[k] for k in sorted({(e["structure"] or {}).get("structure_type") for e in ph} - {None})},
             "cite_as": cite}, ensure_ascii=False, separators=(",", ":")) + "\n"
    out["SCHEMA.md"] = schema_md(counts)
    out["CITATION.cff"] = orchid_cite.cff(DATASET, data)
    return out, reg, counts


MAX_FILE_BYTES = 45 * 1000 * 1000  # every committed file stays below this (GitHub warns at 50 MB)
COMPACT_STRUCTURE = ("structure_type", "basis", "tier", "confidence", "prototype", "strukturbericht", "pearson",
                     "space_group")


def compact_entry(e):
    """One registry.json row: identity, names, structure summary and counts; full record in systems/<SYS>.json."""
    s = e["structure"]
    return {"id": e["id"], "system": e["system"], "elements": e["elements"], "preferred_name": e["preferred_name"],
            "status": e["status"],
            "structure": {k: s[k] for k in COMPACT_STRUCTURE if k in s} if s else None,
            "names": [n["name"] for n in e["names"]],
            "tdb_names": sorted({m["phase_name"] for m in e["models"]}),
            "stoichiometry": e["stoichiometry"],
            "counts": {"models": len(e["models"]),
                       "models_with_tdb_definition": sum(1 for m in e["models"] if m.get("tdb_definition")),
                       "candidates": len(e["candidates"]),
                       "calculated_diagrams": len(e["stability"]),
                       "stable_in_diagrams": sum(1 for x in e["stability"] if x.get("stable")),
                       "literature_values": e["literature"]["values"],
                       "literature_accepted_values": e["literature"]["accepted_values"],
                       "literature_papers": e["literature"]["papers"],
                       "references_listed": len(e["references"]), "total_papers": e["total_papers"]},
            "detail": "systems/%s.json" % e["system"]}


def tdb_def_counts(reg, extra):
    tds = [m["tdb_definition"] for e in reg for m in e["models"] if m.get("tdb_definition")]
    stm = {(t["provenance"]["path"], d["line_start"]) for t in tds for d in t["verbatim"]["type_definitions"]}
    return {"models_with_definition": len(tds),
            "phases_with_definitions": sum(1 for e in reg if any(m.get("tdb_definition") for m in e["models"])),
            "constituent_statements": sum(1 for t in tds if t["verbatim"]["constituent"]),
            "type_definition_references": sum(len(t["verbatim"]["type_definitions"]) for t in tds),
            "distinct_type_definition_statements": len(stm),
            "magnetic_models": sum(1 for t in tds if t["parsed"]["magnetic"]),
            "ordered_with_disordered_part": sum(1 for t in tds if t["parsed"]["disordered_part"]),
            "parameters_counted": sum(t["parameters"]["total"] for t in tds),
            "parameters_not_matched_to_a_declared_phase": extra["unmatched_parameters"],
            "files_with_type_codes_without_definition": extra["type_codes_without_definition"],
            "with_website_url": sum(1 for t in tds if t["provenance"]["website_url"]),
            "with_github_url": sum(1 for t in tds if t["provenance"]["github_url"])}


def structure_types_csv(st_public):
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(["id", "prototype", "prototype_spellings", "strukturbericht", "pearson", "space_group_number",
                "space_group_symbol", "stoichiometry", "verification", "registry_phases", "n_sources", "source_citations",
                "source_dois"])
    for t in st_public:
        sg = t.get("space_group") or {}
        src = t.get("sources") or []
        st = t.get("stoichiometry")
        w.writerow([t["id"], t.get("prototype") or "", ";".join(t.get("prototype_spellings") or []),
                    t.get("strukturbericht") or "", t.get("pearson") or "", sg.get("number") or "", sg.get("symbol") or "",
                    json.dumps(st, sort_keys=True) if isinstance(st, (dict, list)) else (st or ""),
                    t.get("verification") or "", t.get("registry_phases", 0), len(src),
                    " | ".join(x.get("citation") or "" for x in src), ";".join(x["doi"] for x in src if x.get("doi"))])
    return buf.getvalue()


def registry_csv(reg):
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(["id", "system", "preferred_name", "status", "structure_basis", "confidence", "structure_type",
                "prototype", "strukturbericht", "pearson", "space_group_number", "space_group_symbol", "label",
                "tdb_names", "n_models", "n_names", "stoichiometry", "stable_T_min_K", "stable_T_max_K",
                "literature_values", "literature_accepted_values", "n_candidates", "review_reasons",
                "references_first5", "total_papers", "tdb_first_model", "tdb_phase_line", "tdb_constituent_line"])
    for e in reg:
        s = e["structure"] or {}
        sg = s.get("space_group") or {}
        stab = [x for x in e["stability"] if x.get("stable")]
        w.writerow([e["id"], e["system"], e["preferred_name"], e["status"], s.get("basis") or "",
                    s.get("confidence") or "", s.get("structure_type") or "", s.get("prototype") or "",
                    s.get("strukturbericht") or "", s.get("pearson") or "", sg.get("number") or "", sg.get("symbol") or "",
                    s.get("label") or "", ";".join(sorted({m["phase_name"] for m in e["models"]})), len(e["models"]),
                    len(e["names"]),
                    ";".join("%s:%.6g" % kv for kv in (e["stoichiometry"] or {}).items()),
                    min(x["T_min_K"] for x in stab) if stab else "", max(x["T_max_K"] for x in stab) if stab else "",
                    e["literature"]["values"], e["literature"]["accepted_values"], len(e["candidates"]),
                    ";".join(e["review_reasons"]), ";".join(r["doi"] for r in e["references"][:5]), e["total_papers"]]
                   + first_def(e))
    return buf.getvalue()


def first_def(e):
    """CSV: file, verbatim PHASE and CONSTITUENT statements of the first model (models sorted by file, name)."""
    for m in e["models"]:
        td = m.get("tdb_definition")
        if td:
            v = td["verbatim"]
            return [m["path"], v["phase"]["text"], (v["constituent"] or {}).get("text", "")]
    return ["", "", ""]


def jsonld(reg, cite):
    ctx = {"@context": {"@vocab": "https://schema.org/", "schema": "https://schema.org/",
                        "orchid": "https://orchid-project.org/ns/phase#",
                        "prototype": "orchid:prototype", "strukturbericht": "orchid:strukturbericht",
                        "pearsonSymbol": "orchid:pearsonSymbol", "spaceGroupNumber": "orchid:spaceGroupNumber",
                        "structureBasis": "orchid:structureBasis", "reviewStatus": "orchid:reviewStatus",
                        "chemicalSystem": "orchid:chemicalSystem", "structureType": {"@id": "orchid:structureType", "@type": "@id"}}}
    terms = []
    for e in reg:
        s = e["structure"] or {}
        t = {"@id": SITE_DATA + "#" + e["id"], "@type": "DefinedTerm", "termCode": e["id"], "name": e["preferred_name"],
             "inDefinedTermSet": SITE_DATA, "chemicalSystem": e["system"], "reviewStatus": e["status"],
             "alternateName": sorted({n["name"] for n in e["names"]})}
        if s:
            t["structureBasis"] = s.get("basis")
            t["structureType"] = SITE_DATA + "structure_types.json#" + s["structure_type"]
            for k, kk in (("prototype", "prototype"), ("strukturbericht", "strukturbericht"), ("pearson", "pearsonSymbol")):
                if s.get(k):
                    t[kk] = s[k]
            if s.get("space_group"):
                t["spaceGroupNumber"] = s["space_group"]["number"]
        props = []
        for m in e["models"]:
            td = m.get("tdb_definition")
            if not td:
                continue
            v = td["verbatim"]
            txt = "\n".join([v["phase"]["text"]] + ([v["constituent"]["text"]] if v["constituent"] else []) +
                             [x["text"] for x in v["type_definitions"]])
            props.append({"@type": "PropertyValue", "propertyID": "orchid:tdbPhaseDefinition",
                          "name": "TDB definition of %s in %s (line %d)" % (m["phase_name"], m["tdb"], v["phase"]["line_start"]),
                          "value": txt, "url": td["provenance"]["website_url"] or td["provenance"]["github_url"],
                          "description": "verbatim PHASE, CONSTITUENT and TYPE_DEFINITION statements; sha256 " + m["sha256"]})
        if props:
            t["additionalProperty"] = props
        if e["references"]:
            t["citation"] = [{"@id": "https://doi.org/" + r["doi"], "@type": "ScholarlyArticle", "name": r["citation"]}
                             for r in e["references"]]
        based = ["https://doi.org/" + r["doi"] for r in e["references"] if "assessment that models this phase" in r["roles"]]
        based += sorted({"https://github.com/orchid-calphad/orchid-project.org/blob/main/data/" + m["path"] for m in e["models"]})
        if based:
            t["isBasedOn"] = based
        terms.append(t)
    ld = {"@context": ctx["@context"], "@id": SITE_DATA, "@type": "DefinedTermSet",
          "name": "ORCHID phase registry (binary systems)", "version": orchid_cite.VERSION,
          "license": "https://creativecommons.org/licenses/by/4.0/", "identifier": "https://doi.org/" + cite["doi"],
          "citation": cite["text"], "hasDefinedTerm": terms}
    return ctx, ld


def conflicts_md(conflicts, counts):
    o = ["# ORCHID phase registry: conflicts", "",
         "Generated by `tools/build_phase_registry.py` (orchid-project.org). Every item here blocks an automatic "
         "structure assignment for the phases it names; they are in `review_queue.json` with their candidates.", "",
         "| kind | items |", "|---|---|"]
    for k, v in conflicts.items():
        o.append("| %s | %d |" % (k, len(v)))
    o += ["", "## Same TDB name, different structures (Hallstedt 2025 mmc2)", "",
          "A phase name the SGTE collection uses for more than one structure (Pearson symbol / space group). The name "
          "alone therefore cannot fix the structure.", ""]
    for c in conflicts["same_name_different_structures"]:
        o.append("- `%s`: %s (%d registry phases)" % (c["name"], "; ".join(
            "%s %s sg %s" % (s["prototype"], s["pearson"], s["space_group_number"]) for s in c["structures"]),
            len(c["registry_phases"])))
    o += ["", "## Name suggests a structure, DFT disagrees", ""]
    for c in conflicts["name_dft_disagreement"]:
        o.append("- `%s` (%s): name -> %s (sg %s); MP %s / OQMD %s on hull with sg %s" % (
            c["phase"], ", ".join(c["names"]), c["name_structure"], c["name_space_group"], c["mp"], c["oqmd"],
            c["dft_space_group"]))
    o += ["", "## Materials Project and OQMD disagree (single on-hull entries, different space groups)", ""]
    for c in conflicts["mp_oqmd_disagreement"]:
        o.append("- `%s` (%s): MP %s sg %s; OQMD %s sg %s" % (c["phase"], ", ".join(c["names"]), c["mp"],
                                                            c["mp_space_group"], c["oqmd"], c["oqmd_space_group"]))
    o += ["", "## Models disagree on composition (same canonical id, different fixed stoichiometry)", ""]
    for c in conflicts["model_composition_disagreement"]:
        o.append("- `%s`: %s in %s" % (c["phase"], " | ".join(
            ",".join("%s %.4g" % kv for kv in sorted(x.items())) for x in c["compositions"]), ", ".join(c["tdb_files"])))
    o += ["", "## Tier C (assessment-declared) not applied", "",
          "The SGTE crystallographic table (Hallstedt 2025, mmc2) has rows for this phase name and system, but they "
          "list more than one structure, are incomplete, or something in tiers A/B or the TDB model contradicts them. "
          "The phase stays in the review queue.", ""]
    for c in conflicts["assessment_declared_blocked"]:
        o.append("- `%s` (%s): %s; mmc2: %s" % (c["phase"], ", ".join(c["names"]), c["reason"], "; ".join(
            "%s %s %s sg %s (%s)" % (r["phase_name"], r["prototype"] or "?", r["pearson"] or "?",
                                    (r["space_group"] or {}).get("number") or "?", r["actual"]) for r in c["mmc2"])))
    o += ["", "## Polymorph families (several phases of one system at the same composition)", "",
          "Either polymorphs (alpha/beta, LT/HT) or the same phase under different names in different assessments; "
          "a unique DFT match cannot tell which phase it is.", ""]
    for c in conflicts["polymorph_families"]:
        o.append("- %s %s: %s" % (c["system"], ",".join("%s %.4g" % kv for kv in sorted(c["composition"].items())),
                                  ", ".join("`%s`" % p for p in c["phases"])))
    return "\n".join(o) + "\n"


def schema_md(counts):
    return """# ORCHID phase registry (binary systems): schema

Registry version %s, release %s. Generated by `tools/build_phase_registry.py` in the orchid-project.org repository
(deterministic; `--check` verifies the files). Licence: CC BY 4.0 (facts and citations; cite ORCHID and the sources
named in `cite_as`).

**Status.** Version 0 of the registry. Structure types are assigned only when certain (tiers A, B) or when the
assessment's crystallographic table states exactly one consistent structure (tier C, assessor-stated); everything else is
a review candidate. Paper labels come from the machine-extracted literature values (Version 1, machine-extracted).

## Files

| file | content |
|---|---|
| `registry.json` | compact index of every phase: `{dataset, registry_version, inputs, counts, cite_as, record_format: "compact", detail, phases[]}`; each row is a compact entry (below), without `tdb_definition` text or reference details |
| `systems/<SYS>.json` | full records of one binary system: `{system, release, phases[], structure_types{}, cite_as}`; `phases[]` are full phase entries (below, with `models[].tdb_definition`, `stability`, `references`) |
| `structure_types.json` | system-independent structure types, each with its sources |
| `registry.jsonld`, `context.jsonld` | the phases as a schema.org `DefinedTermSet` of `DefinedTerm`s (`citation`, `isBasedOn`; `additionalProperty`: one `PropertyValue` `orchid:tdbPhaseDefinition` per model with the verbatim statements) |
| `registry.csv` | one flat row per phase (`references_first5`: first five reference DOIs, `;`-joined; `tdb_first_model`, `tdb_phase_line`, `tdb_constituent_line`: file and verbatim PHASE / CONSTITUENT statements of the first model, models sorted by file and phase name) |
| `registry.parquet` | the CSV columns as a typed Parquet table (when built with pyarrow) |
| `structure_types.csv`, `structure_types.parquet` | one row per structure type: `id, prototype, prototype_spellings (;), strukturbericht, pearson, space_group_number, space_group_symbol, stoichiometry, verification, registry_phases, n_sources, source_citations (" \\| "), source_dois (;)` |
| `parquet.json` | row counts, CSV SHA-256 and content SHA-256 (canonical JSON rows) of the Parquet tables |
| `pages/<SYS>.html`, `pages/<SYS>.md`, `pages/index.*` | static, JavaScript-free page per binary system (written by `tools/build_static_pages.py`): names, structure, stability, verbatim PHASE/CONSTITUENT statements, example DOIs, DefinedTerm JSON-LD |
| `review_queue.json` | phases with candidates but no assigned structure (`status: needs_review`) |
| `CONFLICTS.md` | same name -> different structures, name vs DFT, MP vs OQMD, model compositions, polymorph families |
| `index.json` | website search index (column arrays) |

Every file stays below 45 MB: the full records are split by system, and `registry.json` is the index over them.
To get all full records, read `registry.json` and then each distinct `detail` file (1 file per system).

## Compact entry (`registry.json`)

`{id, system, elements, preferred_name, status, structure, names[], tdb_names[], stoichiometry, counts, detail}`.
`structure` is `null` or `{structure_type, basis, tier, confidence, prototype, strukturbericht, pearson, space_group}`
(the full entry adds `label`, `caveats`, `dft`, `source`). `names` are the `name` strings of the full entry's `names`.
`counts`: `{models, models_with_tdb_definition, candidates, calculated_diagrams, stable_in_diagrams, literature_values,
literature_accepted_values, literature_papers, references_listed, total_papers}`. `detail`: the per-system file
(relative to `phases/`) holding the full entry with the same `id`.

## Phase entry (full, `systems/<SYS>.json`)

| field | meaning |
|---|---|
| `id` | `phase:<sorted lower-case elements>:<slug>`; the SGTE-corpus canonical id convention (token order folded, liquid synonyms, `#n` instance and non-distinguishing HT/LT tags dropped), so corpus ids carry over |
| `system`, `elements` | chemical system (upper-case elements joined by `-`) |
| `preferred_name` | the TDB phase name used by most assessments |
| `status` | `assigned` (structure below), `needs_review` (no structure; `candidates` listed), `unresolved_no_candidates` |
| `structure` | `{structure_type, basis, tier, confidence, label, prototype, strukturbericht, pearson, space_group, caveats[], dft?}` or `null`. `basis`: `standard_name` (tier A: the TDB name encodes one structure and the model stoichiometry / sublattices and DFT do not contradict it), `dft_unique_match` (tier B, label "DFT-supported": a stoichiometric compound matching exactly one on-hull entry in both Materials Project and OQMD with the same space group, one structure type for that space group and stoichiometry, and no other phase of the system at that composition), `assessment_declared` (tier C, confidence `assessor-stated`: tiers A and B do not assign it, the Hallstedt 2025 crystallographic table (mmc2) rows for this TDB phase name whose listed compound has exactly this system's elements all give one structure type with Pearson symbol and space group, nothing in tiers A/B contradicts it (name structure token, DFT ground state, model stoichiometry/sublattices, order/disorder), and every model's fixed composition or two-sublattice site ratio equals the listed compound; `source` cites mmc2 with its rows), `liquid`, `gas`. Phases whose mmc2 rows disagree stay `needs_review` (`sgte_mmc2_rows_disagree`) |
| `candidates` | for unassigned phases: `{source, structure_type, prototype, strukturbericht, pearson, space_group, detail, material_id?, count?}`; `source` is `name_token`, `sgte_mmc2` (Hallstedt 2025 crystallographic table row for that name), `tdb_comment` (the compiler's comment above the PHASE line), `dft_mp` / `dft_oqmd` (on-hull entries at the composition or within the calculated homogeneity range), `dft_space_group_and_stoichiometry`, `literature_printed` (Pearson / prototype printed beside a label mapped to the phase) |
| `review_reasons` | why no structure was assigned (`off_stoichiometric_model`, `polymorph_family`, `no_on_hull_entry_in_*`, `mp_oqmd_space_groups_differ`, `dft_space_group_*_matches_n_structure_types`, `order_disorder_ambiguity`, `stoichiometry_contradicts_*`, `name_suggests_*_but_dft_space_group_*`, `unrecognised_name_tokens`, `no_structure_token`, ...) |
| `names` | every name: `{name, namespace, label, ...}`; namespaces `tdb` (`tdb:<file>:<NAME>`), `corpus` (`corpus:<canonical id>` and its older aliases), `paper` (`paper:<doi>:<label>`, or `paper:sha256-<prefix>:<label>` for papers without DOI; `basis` says how the label was linked), `mp` / `oqmd` (only for assigned phases whose DFT entry was matched) |
| `models` | `{tdb, path, sha256, phase_name, type_code, sublattices[{site_ratio, constituents}], stoichiometric, composition, disordered_part?, disordered_part_of?, tdb_comment_structure?, tdb_definition}` |
| `models[].tdb_definition` | the model as written in the TDB (bundled, licence-confirmed files only). `verbatim`: `{phase, constituent, type_definitions[]}`, each `{line_start, line_end, text}` with the exact statement from its first character to the terminating `!` (only line breaks normalised to LF; bytes decoded as UTF-8, Latin-1 fallback, see `provenance.decoded_as`); `type_definitions` are every TYPE_DEFINITION whose code is in the phase's type code (`%%` SEQ, magnetic `GES AMEND_PHASE_DESCRIPTION ... MAGNETIC`, `DIS_PART`, ...), each with `type_code` and `amends` (`@` = every phase with that code). `parsed`: `{phase_name, type_code, type_codes[], n_sublattices, site_ratios[], constituents[[...]] (spelling and `%%` major-constituent markers as printed), major_constituents[[...]], magnetic {type_code, afm_factor, structure_fraction, line} or null, disordered_part, ordered_partner}`. `parameters`: counts only, never expressions: `{total, by_kind {G, L, TC, BM, BMAG, BMAGN, THETA, ...as spelled}, L_by_order, interaction_orders[] (orders of parameters with an interacting sublattice, any kind; a missing order counts as 0), endmember_G}`. `provenance`: `{path, sha256, decoded_as, website_url (https://orchid-project.org/data/<path>, when the site copy has the same SHA-256), github_url (orchid-calphad/orchid-pipeline blob, when that repository holds the same git blob)}`. Statement boundaries follow open-tdb `scripts/validate_tdb.py` (each span is checked to reduce to its statement); parameter keys follow `check_subsystem_compat.py` / `audit_tdb_models.py`. TDB text is never executed |
| `stoichiometry` | mole fractions when every model of the phase is stoichiometric and they agree |
| `stability` | per calculated diagram (pycalphad, SGTE TDB): `{tdb, tdb_sha256, diagram, status, stable, T_min_K, T_max_K, x_element, x_min, x_max, stable_at_top_of_range, calc_T_range_K?}` (calculated range 300-2500 K unless `calc_T_range_K` says otherwise); `x` is the mole fraction of `x_element` of the phase itself (boundary, invariant and single-phase grid points; 0.01 grid); `stable: false` = not stable anywhere in the calculated range |
| `literature` | `{values, accepted_values, papers}`: literature values whose phase list contains this phase |
| `corpus` | `{in_sgte_corpus, phase_mappings, mapping_basis, printed_labels}` from the SGTE-corpus phase-mapping records |
| `references`, `total_papers` | up to %d example references, deterministic: assessments first, then papers with most values, ties by DOI: `{doi, citation, roles[], labels[], values, tdb_files?}`. Roles: `assessment that models this phase` (source of a TDB declaring it), `reports values for this phase`, `cited as the source of values for this phase`. Only DOIs present in ORCHID records |

## Structure type

`{id, prototype, prototype_spellings[], strukturbericht, pearson, space_group{number, symbol}, stoichiometry, sources[],
verification, registry_phases}`. Rows come from the Hallstedt 2025 crystallographic table (mmc2; CC BY 4.0) and a
small built-in table of standard types (Villars & Calvert; AFLOW prototype library). `verification` says which: rows
not confirmed by mmc2 carry `%%TODO verify`. `st:liquid` and `st:gas` are definitions.

## Counts (this build)

```
%s
```
""" % (REGISTRY_VERSION, orchid_cite.RELEASE_TAG, MAX_REFS, json.dumps(counts, indent=1))


# ------------------------------------------------------------------------------------------------ parquet --
PARQUET_HELPER = r'''
import csv, io, json, sys, hashlib
import pyarrow as pa, pyarrow.parquet as pq
src, dst, meta = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
rows = list(csv.DictReader(open(src, encoding="utf-8")))
ints = {"space_group_number", "n_models", "n_names", "literature_values", "literature_accepted_values",
        "n_candidates", "total_papers", "registry_phases", "n_sources"}
flts = {"stable_T_min_K", "stable_T_max_K"}
cols = list(rows[0].keys()) if rows else []
fields = [(c, pa.int32() if c in ints else pa.float64() if c in flts else pa.string()) for c in cols]
def conv(c, v):
    if v == "":
        return None
    return int(v) if c in ints else float(v) if c in flts else v
data = [{c: conv(c, r[c]) for c in cols} for r in rows]
t = pa.Table.from_pylist(data, schema=pa.schema(fields).with_metadata(meta))
pq.write_table(t, dst, compression="zstd", compression_level=9, row_group_size=65536, use_dictionary=True,
               data_page_version="2.0", write_statistics=True)
t = pq.read_table(dst)
h = hashlib.sha256()
for r in t.to_pylist():
    h.update(json.dumps(r, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()); h.update(b"\n")
kv = {k.decode(): v.decode() for k, v in (t.schema.metadata or {}).items() if k.startswith(b"orchid.")}
print(h.hexdigest() + " " + hashlib.sha256(json.dumps(kv, sort_keys=True).encode()).hexdigest())
'''


def parquet_python(a):
    try:
        import pyarrow  # noqa: F401
        return sys.executable
    except ImportError:
        return a.parquet_python if a.parquet_python and os.path.exists(a.parquet_python) else None


PARQUET_TABLES = (("registry.parquet", "registry.csv", "phase_registry"),
                  ("structure_types.parquet", "structure_types.csv", "structure_types"))
PAGES_DIR = "pages"  # phases/pages/ is written by tools/build_static_pages.py, not by this builder


def parquet_meta(cite, table="phase_registry"):
    return {"orchid.dataset": DATASET, "orchid.table": table, "orchid.license": "CC-BY-4.0",
            "orchid.generator": "orchid-project.org/tools/build_phase_registry.py", "orchid.doi": cite["doi"],
            "orchid.citation": cite["text"], "orchid.cite_as": json.dumps(cite, sort_keys=True, ensure_ascii=False)}


def meta_hash(cite, table="phase_registry"):
    return sha256_bytes(json.dumps(parquet_meta(cite, table), sort_keys=True).encode())


def write_parquet(a, out_dir, csv_path, cite, name="registry.parquet", table="phase_registry"):
    py = parquet_python(a)
    if not py:
        print("no pyarrow: %s not written" % name)
        return None
    meta = parquet_meta(cite, table)
    dst = os.path.join(out_dir, name)
    r = subprocess.run([py, "-c", PARQUET_HELPER, csv_path, dst + ".tmp", json.dumps(meta)], capture_output=True, text=True)
    if r.returncode:
        raise SystemExit("parquet failed: " + r.stderr)
    return dst, r.stdout.strip()


def parquet_content_hash(a, path):
    py = parquet_python(a)
    if not py or not os.path.exists(path):
        return None
    code = ("import sys,json,hashlib,pyarrow.parquet as pq\nt=pq.read_table(sys.argv[1])\nh=hashlib.sha256()\n"
            "for r in t.to_pylist():\n"
            "    h.update(json.dumps(r,sort_keys=True,ensure_ascii=False,separators=(',',':')).encode());h.update(b'\\n')\n"
            "kv={k.decode():v.decode() for k,v in (t.schema.metadata or {}).items() if k.startswith(b'orchid.')}\n"
            "print(h.hexdigest()+' '+hashlib.sha256(json.dumps(kv,sort_keys=True).encode()).hexdigest())")
    return subprocess.run([py, "-c", code, path], capture_output=True, text=True).stdout.strip()


# --------------------------------------------------------------------------------------------------- main --
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=orchid_cite.DEFAULT_DATA, help="orchid-data checkout (writes phases/)")
    ap.add_argument("--open-tdb", default=DEF_OPEN_TDB, help="open-tdb checkout (read only)")
    ap.add_argument("--corpus", default=DEF_CORPUS, help="SGTE corpus build directory (read only)")
    ap.add_argument("--mmc2", default=DEF_MMC2, help="Hallstedt 2025 supplementary table mmc2.xlsx")
    ap.add_argument("--pipeline", default=os.path.join(HOME, "projects", "orchid-pipeline"),
                    help="local orchid-pipeline clone (git blob ids for the GitHub links; read only)")
    ap.add_argument("--pipeline-ref", default="origin/main", help="ref of the orchid-pipeline clone to link to")
    default_ppy = os.path.expanduser("~/orchid-py313/bin/python")
    ap.add_argument("--parquet-python", default=default_ppy, help="interpreter with pyarrow for registry.parquet")
    ap.add_argument("--check", action="store_true", help="exit 1 if phases/ differs from a fresh build")
    a = ap.parse_args(argv)
    out, reg, counts = build(a)
    out_dir = os.path.join(a.data, "phases")
    big = sorted((len(t.encode()), rel) for rel, t in out.items() if len(t.encode()) > MAX_FILE_BYTES)
    if big:
        raise SystemExit("files over %d bytes (GitHub limit 50 MB): %s" % (MAX_FILE_BYTES, ", ".join("%s %d" % (r, n) for n, r in big)))
    cite = orchid_cite.cite_as(DATASET, a.data)
    man_rel = "parquet.json"
    if a.check:
        bad = []
        for rel, text in sorted(out.items()):
            p = os.path.join(out_dir, rel)
            if not os.path.exists(p) or open(p, encoding="utf-8").read() != text:
                bad.append(rel)
        disk = set()
        for d, _, fs in os.walk(out_dir):
            for f in fs:
                disk.add(os.path.relpath(os.path.join(d, f), out_dir))
        extra = sorted(x for x in disk - set(out) - {pq_ for pq_, _, _ in PARQUET_TABLES} - {man_rel}
                       if not x.startswith(PAGES_DIR + os.sep))
        bad += ["stale file: " + x for x in extra]
        pm = os.path.join(out_dir, man_rel)
        if parquet_python(a):
            man = json.load(open(pm)) if os.path.exists(pm) else {}
            for pq_, csv_, table in PARQUET_TABLES:
                mt = man if table == "phase_registry" else man.get(table, {})
                want = sha256_bytes(out[csv_].encode())
                got = (parquet_content_hash(a, os.path.join(out_dir, pq_)) or " ").split(" ")
                if mt.get("csv_sha256") != want or got[1] != meta_hash(cite, table) or got[0] != mt.get("content_sha256"):
                    bad.append(pq_)
        if bad:
            print("stale: " + ", ".join(bad[:20]) + (" ... (%d)" % len(bad) if len(bad) > 20 else ""), file=sys.stderr)
            return 1
        print("phases/ is up to date (%d files, %d phases)" % (len(out), len(reg)))
        return 0
    os.makedirs(os.path.join(out_dir, "systems"), exist_ok=True)
    keep = set(out)
    for d, _, fs in os.walk(out_dir):
        for f in fs:
            rel = os.path.relpath(os.path.join(d, f), out_dir)
            if rel not in keep and rel not in {pq_ for pq_, _, _ in PARQUET_TABLES} | {man_rel} \
                    and not rel.startswith(PAGES_DIR + os.sep):
                os.remove(os.path.join(d, f))
    for rel, text in sorted(out.items()):
        p = os.path.join(out_dir, rel)
        cur = open(p, encoding="utf-8").read() if os.path.exists(p) else None
        if cur != text:
            with open(p + ".tmp", "w", encoding="utf-8") as fh:
                fh.write(text)
            os.replace(p + ".tmp", p)
    pm = os.path.join(out_dir, man_rel)
    manifest = {}
    for pq_, csv_, table in PARQUET_TABLES:
        res = write_parquet(a, out_dir, os.path.join(out_dir, csv_), cite, pq_, table)
        if not res:
            continue
        dst, ch = res
        ch, mh = ch.split(" ")
        old = (parquet_content_hash(a, dst) or " ").split(" ") if os.path.exists(dst) else ["", ""]
        if old == [ch, mh]:
            os.remove(dst + ".tmp")  # unchanged content: keep the existing bytes
        else:
            os.replace(dst + ".tmp", dst)
        rows = out[csv_].count("\n") - 1 if table != "phase_registry" else len(reg)
        entry = {"table": pq_, "rows": rows, "csv_sha256": sha256_bytes(out[csv_].encode()), "content_sha256": ch}
        if table == "phase_registry":
            manifest.update(entry)
            manifest["note"] = ("content_sha256 = SHA-256 of the canonical JSON rows (sorted keys); independent of "
                                "the Parquet encoder")
        else:
            manifest[table] = entry
    if manifest:
        with open(pm, "w", encoding="utf-8") as fh:
            fh.write(dump(manifest))
    print(json.dumps(counts, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
