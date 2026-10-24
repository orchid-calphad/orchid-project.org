#!/usr/bin/env python3
"""Build the DFT-vs-experiment parity data for formation enthalpies (data/literature/dft_parity.json).

Reads (never writes):

  data/literature/orchid-literature.jsonl.gz      literature records (this site)
  data/literature/papers.json                     bibliography (citation, DOI, PDF hash prefix)
  <open-tdb>/hulls/snapshot/snapshot.json         Materials Project 2026.04.13 (GGA/GGA+U, MP corrections)
  <open-tdb>/hulls/snapshot-oqmd/snapshot.json    OQMD qmdb v1.8 (fitted elemental references)
  <open-tdb>/build/sgte-corpus/dft-check/sign-check/sign_check.json   page-image sign verdicts (optional)

Selection: property HM_FORM in the normalized unit kJ/mol-atom, evidence class experimental / model /
unspecified (literature DFT values are not compared with DFT). Values printed with basis unstated
("kJ/mol (basis unstated)") are excluded: per-atom vs per-formula cannot be settled here.

Matching: the record's phase label is parsed into a formula (Greek / polymorph prefixes, sublattice brackets,
subscript markers (NiSi_2) and formula groups handled; tokens outside the record's system dropped); if no label
parses, a phase composition derived from stoichiometry is used. The reduced composition is matched to the
lowest-energy non-deprecated entry at exactly that composition (4-decimal mole fractions) in each DFT snapshot,
unless the record names a crystal structure (Strukturbericht designator or prototype in the printed label, or in
a model-phase name whose formula agrees): then the lowest entry of that structure's space group is used
(`dft_structure_match` exact), and when the snapshot has no entry of that space group the point is paired with the
ground state but marked structure_unmatched and left out of the default (calibration) subset (audit A4 E1,
2026-10-24). Energy above hull is derived from ORCHID's own released hull of the same snapshot (hulls/dft-<src>/),
the upstream value is kept as ehull_upstream (A4 E2). MP and
OQMD are separate sources with different references and are never averaged. DFT eV/atom x 96.485 =
kJ/mol-atom. Residual = DFT - literature.

One point per (literature record x DFT source). Aggregate statistics per source and subset; exclusions are
counted by reason. Deterministic (sorted keys, no timestamps); standard library only; no network.

  python3 tools/build_dft_parity.py            # write data/literature/dft_parity.json (+ croissant entry)
  python3 tools/build_dft_parity.py --check    # exit 1 if the files on disk are stale
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
import re
import sys
import unicodedata
from collections import Counter, defaultdict
from fractions import Fraction

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orchid_cite  # noqa: E402

SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPEN_TDB = os.path.expanduser("~/projects/open-tdb")
LIT = os.path.join(SITE, "data", "literature")
OUT = os.path.join(LIT, "dft_parity.json")
HULLS = os.path.join(SITE, "data", "hulls")   # copy of <orchid-data>/hulls (ORCHID's released DFT hulls)
CROISSANT = os.path.join(SITE, "data", "croissant.json")
KJ_PER_EV = 96.485
UNIT = "kJ/mol-atom"
CLASSES = ("experimental", "model", "unspecified")
BAND = 10.0  # kJ/mol-atom, the band drawn on the page and used for "within band" counts

SOURCES = {
    "mp": {"label": "Materials Project", "rel": os.path.join("hulls", "snapshot", "snapshot.json"),
           "citation": "A. Jain et al., Commentary: The Materials Project: A materials genome approach to "
                       "accelerating materials innovation, APL Materials 1 (2013) 011002.",
           "doi": "10.1063/1.4812323", "url": "https://next-gen.materialsproject.org/",
           "license": "CC BY 4.0",
           "reference_note": "GGA/GGA+U mixing scheme with MP fitted energy corrections (notably anions and "
                             "some transition-metal oxides); elemental ground states at 0 K."},
    "oqmd": {"label": "OQMD", "rel": os.path.join("hulls", "snapshot-oqmd", "snapshot.json"),
             "citation": "J. E. Saal et al., JOM 65 (2013) 1501-1509; S. Kirklin et al., npj Comput. Mater. 1 "
                         "(2015) 15010.",
             "doi": "10.1007/s11837-013-0755-4", "doi2": "10.1038/npjcompumats.2015.10",
             "url": "https://oqmd.org/", "license": "CC BY 4.0",
             "reference_note": "'standard' fit: elemental reference energies fitted to experimental formation "
                               "enthalpies (fitted elemental reference, 0 K)."},
}

ELEMENTS = set("""H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn
Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm
Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu
Am Cm Bk Cf Es Fm Md No Lr""".split())
OXIDE_HALIDE = {"O", "F", "Cl", "Br", "I"}
OTHER_NONMETAL = {"H", "B", "C", "N", "P", "S", "Se", "Te", "He", "Ne", "Ar", "Kr", "Xe"}
SUBSCRIPTS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")

METHOD_RULES = [  # (class, regex) on the lower-cased printed method; first match wins
    ("emf_vapour", r"\be\.?m\.?f\b|emf|galvanic|potentiomet|electrochem|coulometr|vapou?r|knudsen|effusion|"
                   r"\bkem[ls]\b|keml|kems|torsion|isopiestic|mass[- ]spectr|p[-–]c[-–]t|\bpct\b|equilibri|"
                   r"second[- ]law|third[- ]law|dissociation pressure|transpiration"),
    ("calorimetry", r"calori|\bcal\b|\bcal\.|sol\. ?cal|solution|\bdrop\b|combustion|\bbomb\b|direct (reaction|"
                    r"synthesis)|\bdsc\b|htdsc|\bdta\b|reaction cal|react\. ?cal|aneroid|enthalp"),
]
METHOD_LABELS = {"calorimetry": "calorimetry (solution, drop, combustion, direct reaction, DSC)",
                 "emf_vapour": "EMF / vapour pressure / Knudsen effusion / equilibrium",
                 "other_unknown": "other or not printed"}
CHEM_LABELS = {"oxide_halide": "oxides and halides (contains O, F, Cl, Br or I)",
               "other_nonmetal": "hydrides, borides, carbides, nitrides, phosphides, chalcogenides",
               "intermetallic": "intermetallics (metals and metalloids only)"}


# ----------------------------------------------------------------- formula parsing (from open-tdb dft_check) --
TOKEN_RE = re.compile(r"^(?:[A-Z][a-z]?\d*(?:\.\d+)?)+$")
LIQUID_AMORPHOUS = re.compile(r"\((?:[^()]*[-‐ ,/])?(?:liq\w*|am|amorph\w*|glass\w*|melt|vitreous)(?:[-‐ ,/][^()]*)?\)|"
                              r"\b(?:liquid|amorphous|glass(?:y)?|melt|vitreous)\b", re.I)
PART_RE = re.compile(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)")


def _counts(run):
    out = {}
    for el, n in PART_RE.findall(run):
        if el not in ELEMENTS:
            return None
        out[el] = out.get(el, 0) + (Fraction(n) if n else Fraction(1))
    return out


def _runs(text):
    """Maximal runs of formula tokens; runs join across whitespace only."""
    out = []
    for seg in re.split(r"[-‐-―,;:/'′’\"`+=~_]", text):
        cur = []
        for tok in seg.split():
            t = re.sub(r"^[^A-Z0-9]+(?=[A-Z])", "", tok)  # glued prefix: aGeSe, αGeSe
            if re.fullmatch(r"[A-Z]{3,}", t):  # all-caps words (UNKNOWN, HCP) are not formulas
                pass
            elif TOKEN_RE.match(t) and _counts(t) is not None:
                cur.append(t)
            elif cur and re.fullmatch(r"\d+(?:\.\d+)?", tok) and re.search(r"[A-Za-z]$", cur[-1]):
                cur.append(tok)
            else:
                if cur:
                    out.append(cur)
                cur = []
        if cur:
            out.append(cur)
    return out


def _pick(runs, sysel):
    cands = []
    for run in runs:
        whole = _counts("".join(run))
        if sysel is None or set(whole) <= sysel:
            if len(whole) >= 2:
                cands.append(("".join(run), whole, ""))
            continue
        best = None
        for i in range(len(run)):
            for j in range(len(run), i, -1):
                c = _counts("".join(run[i:j]))
                if c and len(c) >= 2 and set(c) <= sysel and (best is None or j - i > best[0]):
                    best = (j - i, "".join(run[i:j]), c)
        if best:
            cands.append((best[1], best[2], "tokens outside declared system dropped"))
        elif len(whole) >= 2:
            cands.append(("".join(run), whole, "elements outside declared system"))
    return cands


def parse_phase(phase, system):
    """Return (counts, note) or (None, reason)."""
    if phase is None or not str(phase).strip():
        return None, "no phase label"
    s = unicodedata.normalize("NFKC", str(phase)).translate(SUBSCRIPTS)
    # subscript markers between an element (or a closing group bracket) and its count: NiSi_2, NiSi_{2}, Al^3
    # (audit A4 E4: "αNiSi_2" was read as NiSi because "_" separates formula tokens)
    s = re.sub(r"(?<=[A-Za-z)\]])[_^]\{?(\d+(?:\.\d+)?)\}?", r"\1", s)
    # comma decimal subscripts: PuO1,52 = PuO1.52 (audit round 3, POSTFIX A3-N6; it was read as PuO, so a value per
    # formula was divided by 2 instead of 2.52). Only a digit-comma-digits run that directly follows an element symbol.
    s = re.sub(r"(?<=[A-Za-z])(\d+),(\d+)(?![\d,])", r"\1.\2", s)
    notes = []
    if re.search(r"\((g|gas|l|liq)\)", s, re.I):
        return None, "gas or liquid state marker in phase label"
    if LIQUID_AMORPHOUS.search(s):
        # CaO(liq-am), glassy / amorphous / melt states are not the crystalline compound of a 0 K DFT entry (A4 E7)
        return None, "liquid, amorphous or glass state in phase label"
    had_group = False
    while True:  # expand formula groups: (Re)24 -> Re24, Mg(AlH4)2 -> MgAl2H8
        m = re.search(r"\(([A-Z][A-Za-z0-9.]*)\)(\d+(?:\.\d+)?)", s)
        c = _counts(m.group(1)) if m and TOKEN_RE.match(m.group(1)) else None
        if not c:
            break
        k = Fraction(m.group(2))
        s = s[:m.start()] + "".join(el + (str((c[el] * k).numerator) if (c[el] * k).denominator == 1
                                          else "%g" % float(c[el] * k)) for el in c) + s[m.end():]
        had_group = True
    pairs = re.findall(r"\(([^()]*)\)|\[([^\[\]]*)\]", s)
    inner = " , ".join(a or b for a, b in pairs)
    outer = re.sub(r"\([^()]*\)|\[[^\[\]]*\]", " , ", s)
    sysel = set(system) if system else None
    cands = _pick(_runs(outer), sysel)
    if not cands and inner:
        cands = _pick(_runs(inner), sysel)
        if cands:
            notes.append("formula taken from parenthetical")
    if not cands:
        return None, "no multi-element formula in phase label"
    if len({tuple(sorted(c.items())) for _, c, _ in cands}) > 1:
        return None, "ambiguous: several formulas in phase label"
    run, counts, n = cands[0]
    if n == "elements outside declared system":
        return None, "formula elements outside the record's system"
    if n:
        notes.append(n)
    if len(PART_RE.findall(run)) != len(counts):
        if not had_group:
            return None, "element repeated in formula without grouping"
        notes.append("element repeated in printed formula; counts summed")
    return counts, "; ".join(notes)


def reduced_key(counts):
    tot = sum(counts.values())
    return tuple(sorted((el, round(float(c) / float(tot), 4)) for el, c in counts.items()))


def _integerize(counts):
    """Non-integer printed ratios (U0.5C0.5, Si0.3333Ta0.6667) -> the small-integer formula they denote, when one
    with denominators <= 24 reproduces every mole fraction within 5e-4 (audit A4 E7); else unchanged."""
    if all(Fraction(c).denominator == 1 for c in counts.values()):
        return counts
    tot = sum(Fraction(c) for c in counts.values())
    fr = {el: Fraction(c) / tot for el, c in counts.items()}
    approx = {el: f.limit_denominator(24) for el, f in fr.items()}
    if any(abs(float(approx[el] - fr[el])) > 5e-4 for el in fr) or sum(approx.values()) != 1:
        return counts
    lcm = 1
    for f in approx.values():
        lcm = lcm * f.denominator // math.gcd(lcm, f.denominator)
    return {el: f * lcm for el, f in approx.items()}


def formula_str(counts):
    """Reduced formula: divide by the GCD of integer counts (non-integer counts that denote a small-integer formula
    are written as that formula; others kept as printed ratios)."""
    counts = _integerize(counts)
    if all(Fraction(c).denominator == 1 for c in counts.values()):
        g = 0
        for c in counts.values():
            g = math.gcd(g, int(c))
        counts = {el: Fraction(c) / g for el, c in counts.items()}

    def fmt(c):
        c = Fraction(c)
        if c == 1:
            return ""
        return str(c.numerator) if c.denominator == 1 else "%g" % float(c)
    return "".join(el + fmt(counts[el]) for el in sorted(counts, key=lambda e: (e in OXIDE_HALIDE | OTHER_NONMETAL, e)))


# --------------------------------------------------------------------- crystal structure designators (A4 E1) --
# Strukturbericht designator -> space-group numbers (standard structure types; B19' is the NiTi martensite).
STRUKTUR_SG = {
    "A1": {225}, "A2": {229}, "A3": {194}, "A3'": {194}, "A4": {227}, "A5": {141}, "A6": {139}, "A7": {166},
    "A9": {194}, "A12": {217}, "A13": {213}, "A15": {223}, "A20": {63},
    "B1": {225}, "B2": {221}, "B3": {216}, "B4": {186}, "B81": {194}, "B82": {194}, "B10": {129}, "B11": {129},
    "B19": {51}, "B19'": {11}, "B20": {198}, "B27": {62}, "B31": {62}, "B32": {227}, "B33": {63}, "B35": {191},
    "BH": {187}, "C1": {225}, "C2": {205}, "C6": {164}, "C11B": {139}, "C14": {194}, "C15": {227}, "C16": {140},
    "C19": {166}, "C23": {62}, "C32": {191}, "C36": {194}, "C37": {62}, "C38": {129}, "C40": {180}, "C49": {63},
    "C54": {70}, "D03": {225}, "D011": {62}, "D019": {194}, "D022": {139}, "D023": {139}, "D024": {194},
    "D0A": {59}, "D1A": {87}, "D2B": {139}, "D2D": {191}, "D513": {164}, "D52": {164}, "D58": {62}, "D73": {220},
    "D82": {217}, "D85": {166}, "D88": {193}, "D8B": {136}, "D8M": {140}, "E93": {227},
    "L10": {123}, "L12": {221}, "L21": {225}, "L60": {123},
}
# prototype formula (as printed, "MgZn2-type") -> designator
PROTOTYPE_SB = {"CU": "A1", "W": "A2", "MG": "A3", "CSCL": "B2", "NACL": "B1", "NIAS": "B81", "AUCD": "B19",
                "FESI": "B20", "FEB": "B27", "MNP": "B31", "NATL": "B32", "CRB": "B33", "CAF2": "C1",
                "MOSI2": "C11B", "MGZN2": "C14", "MGCU2": "C15", "CU2MG": "C15", "AL2CU": "C16", "ALB2": "C32",
                "MGNI2": "C36", "CRSI2": "C40", "ZRSI2": "C49", "TISI2": "C54", "BIF3": "D03", "FE3C": "D011",
                "NI3SN": "D019", "AL3TI": "D022", "AL3ZR": "D023", "NI3TI": "D024", "CU3TI": "D0A",
                "MN5SI3": "D88", "W5SI3": "D8M", "AUCU": "L10", "AUCU3": "L12", "CU3AU": "L12", "CR3SI": "A15",
                "CU2MNAL": "L21", "CACU5": "D2D", "THMN12": "D2B", "CU2SB": "C38", "WC": "BH", "COSN": "B35",
                "CUTI": "B11", "NI2AL3": "D513", "TH3P4": "D73"}
_SB_TOKEN = re.compile(r"^([A-EL])(\d{1,2})(?:_?\{?(\d{1,2}|[a-z])\}?)?(['′’]?)$")
_DISORDERED_SB = {"A1", "A2", "A3"}   # solution designators: never a compound's structure


def _sb_norm(tok):
    """'D0_3' -> 'D03', 'C11b' -> 'C11B', "B19′" -> "B19'"; None if not a designator."""
    m = _SB_TOKEN.match(tok)
    if not m:
        return None
    k = (m.group(1) + m.group(2) + (m.group(3) or "")).upper() + ("'" if m.group(4) else "")
    return k if k in STRUKTUR_SG else None


def structure_designators(labels, model_phases, counts):
    """{designator: provenance} for a record's compound: Strukturbericht tokens or 'X-type' prototypes printed in
    the phase label; Strukturbericht tokens of a model-phase name (CUZR_B2, C14_LAVES) when the name's formula part,
    if any, is the record's compound. Solution designators A1/A2/A3 are ignored. A token that is also a formula of
    the record's elements (B2 in a B-containing system) is not read as a designator."""
    els = {e.upper() for e in counts}
    want = reduced_key(counts)
    out = {}
    for lab in labels:
        t = unicodedata.normalize("NFKC", str(lab))
        for tok in re.split(r"[\s()\[\],;:/+=~-]+", t):
            k = _sb_norm(tok.strip())
            if k and k not in _DISORDERED_SB and not (tok[0].upper() in els and re.fullmatch(r"[A-Z]\d+", tok)):
                out.setdefault(k, "printed_label")
        for m in re.finditer(r"([A-Z][A-Za-z0-9]*)\s*[-‐ ]?\s*type\b", t):
            k = PROTOTYPE_SB.get(m.group(1).upper())
            if k and k not in _DISORDERED_SB:
                out.setdefault(k, "printed_prototype")
    for mp in model_phases:
        parts = [p for p in str(mp).upper().split("_") if p]
        # SGTE writes a prime as a trailing P (NITI_B19P = B19')
        sbs = [_sb_norm(p) or (_sb_norm(p[:-1] + "'") if p.endswith("P") and len(p) > 2 else None) for p in parts]
        if not any(sbs):
            continue
        formula_ok = True
        for p, sb in zip(parts, sbs):
            if sb:
                continue
            # SGTE names are upper case: a formula-looking part (NISI, CUZR, FE2TI) must be the record's compound
            fe = _upper_formula(p, {e.capitalize() for e in els})
            if fe is not None:
                formula_ok = formula_ok and reduced_key(fe) == want
        if not formula_ok:
            continue
        for sb in sbs:
            if sb and sb not in _DISORDERED_SB:
                out.setdefault(sb, "model_phase_name")
    return out


def _upper_formula(part, alphabet):
    """Counts of an upper-case formula word (NISI2 -> {Ni:1, Si:2}) when it splits unambiguously into elements of
    `alphabet` (the record's elements); None for words that are not formulas (LAVES, BETA)."""
    sols = []

    def rec(i, acc):
        if len(sols) > 1:
            return
        if i == len(part):
            if acc:
                sols.append(dict(acc))
            return
        for L in (2, 1):
            el = part[i:i + L].capitalize()
            if len(el) == L and el in alphabet:
                j = i + L
                m = re.match(r"\d+(?:\.\d+)?", part[j:])
                n = Fraction(m.group(0)) if m else Fraction(1)
                j += len(m.group(0)) if m else 0
                if el in acc:
                    continue
                acc[el] = n
                rec(j, acc)
                del acc[el]
    rec(0, {})
    if len(sols) != 1 or len(sols[0]) < 2:
        return None
    return sols[0]


# ------------------------------------------------------------------------------------------------- DFT --
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class DFT:
    def __init__(self, key, root):
        cfg = SOURCES[key]
        path = os.path.join(root, cfg["rel"])
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
        self.by_comp = defaultdict(list)
        for e in d["entries"]:
            if e.get("deprecated"):
                continue
            comp = {k: float(v) for k, v in e["composition_reduced"].items()}
            tot = sum(comp.values())
            k = tuple(sorted((el, round(v / tot, 4)) for el, v in comp.items()))
            self.by_comp[k].append(e)
        self.meta = {"label": cfg["label"], "snapshot": "open-tdb/" + cfg["rel"].replace(os.sep, "/"),
                     "snapshot_sha256": sha256_file(path), "release": d.get("release"),
                     "thermo_type": d.get("thermo_type"),
                     "reference": d.get("reference") or {"kind": "dft_elemental_ground_state_0K", "units": "eV/atom",
                                                         "per": "atom", "temperature_K": 0.0},
                     "reference_note": cfg["reference_note"], "citation": cfg["citation"], "doi": cfg["doi"],
                     "doi2": cfg.get("doi2"), "url": cfg["url"], "license": cfg["license"]}

    def lookup(self, key, designators=None, hulls=None):
        cands = self.by_comp.get(key)
        if not cands:
            return None
        order = lambda e: (float(e["formation_energy_per_atom"]), bool(e.get("theoretical")), e["material_id"])
        ground = min(cands, key=order)
        best, smatch, sgs = ground, "no_structure_printed", set()
        if designators:
            for d in designators:
                sgs |= STRUKTUR_SG[d]
            same = [e for e in cands if e.get("spacegroup_number") in sgs]
            if same:
                best, smatch = min(same, key=order), "exact"
            else:
                smatch = "structure_unmatched"
        eh_up = best.get("energy_above_hull")
        ef = float(best["formation_energy_per_atom"])
        eh, slope = None, None
        if hulls is not None:
            cr = {k: float(v) for k, v in best["composition_reduced"].items()}
            tot = sum(cr.values())
            eh, slope = hulls.distance(tuple(sorted((k, v / tot) for k, v in cr.items())), ef)
        return {"id": best["material_id"], "ef_eV": round(ef, 6), "dft": round(ef * KJ_PER_EV, 3),
                "sg": best.get("spacegroup_symbol"), "sgn": best.get("spacegroup_number"),
                "ehull": eh, "ehull_upstream": None if eh_up is None else round(float(eh_up) * KJ_PER_EV, 3),
                "npoly": len(cands), "theoretical": best.get("theoretical"), "smatch": smatch,
                "designators": sorted(designators) if designators else None,
                "ground_id": ground["material_id"] if best is not ground else None, "hull_slopes": slope}


class Hulls:
    """ORCHID's released binary hulls of one source (hulls/dft-<src>/<A>/<A-B>.json): the hull energy at a
    composition and the slopes of the two adjacent hull segments (eV/atom per unit x of the system's x element)."""

    def __init__(self, root, src, snapshot_sha):
        self.dir = os.path.join(root, "dft-" + src)
        self.sha = snapshot_sha
        self.cache = {}
        self.meta = {"path": "hulls/dft-%s/" % src}

    def _doc(self, els):
        sid = "-".join(sorted(e.upper() for e in els))
        if sid not in self.cache:
            p = os.path.join(self.dir, sid.split("-")[0], sid + ".json")
            d = None
            if os.path.exists(p):
                with open(p, encoding="utf-8") as fh:
                    d = json.load(fh)
                psha = (d.get("provenance") or {}).get("snapshot_sha256")
                if psha != self.sha:
                    raise SystemExit("hull %s was built from snapshot %s, parity uses %s: rebuild the hulls"
                                     % (p, psha, self.sha))
            self.cache[sid] = d
        return self.cache[sid]

    def distance(self, key, ef):
        """(energy above ORCHID's hull in kJ/mol-atom, {left, right} hull slopes in kJ/mol-atom per unit x of the
        x element, x element) for formation energy ef (eV/atom) at reduced composition key."""
        comp = dict(key)
        d = self._doc(comp)
        if d is None:
            return None, None
        xel = d["x_element"].capitalize()
        x = comp.get(xel, 0.0)
        iso = d["isotherms"][0]
        hull, left, right = None, None, None
        for it in iso["intervals"]:
            if it["x_lo"] - 1e-4 <= x <= it["x_hi"] + 1e-4:
                e = it["intercept"] + it["slope"] * x
                hull = e if hull is None else min(hull, e)
                if abs(x - it["x_hi"]) <= 1e-4:
                    left = it["slope"]
                elif abs(x - it["x_lo"]) <= 1e-4:
                    right = it["slope"]
                else:
                    left = right = it["slope"]
        if hull is None:
            return None, None
        eh = (ef - hull) * KJ_PER_EV
        return round(eh, 4), {"x_element": xel.upper(), "x": x,
                              "left": None if left is None else round(left * KJ_PER_EV, 3),
                              "right": None if right is None else round(right * KJ_PER_EV, 3)}


# ------------------------------------------------------------------------------------------ classifiers --
def method_class(method):
    if not method or not str(method).strip():
        return "other_unknown"
    s = str(method).lower()
    for cls, rx in METHOD_RULES:
        if re.search(rx, s):
            return cls
    return "other_unknown"


def chem_class(elements):
    els = set(elements)
    if els & OXIDE_HALIDE:
        return "oxide_halide"
    if els & OTHER_NONMETAL:
        return "other_nonmetal"
    return "intermetallic"


def _ranks(xs):
    """Average ranks (1-based) with ties sharing their mean rank."""
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        for k in range(i, j + 1):
            r[order[k]] = (i + j) / 2.0 + 1
        i = j + 1
    return r


def stats(pairs):
    """pairs: list of (lit, dft). Residual = dft - lit."""
    n = len(pairs)
    if not n:
        return {"n": 0}
    res = [d - l for l, d in pairs]
    out = {"n": n, "mae": sum(abs(r) for r in res) / n, "rmse": math.sqrt(sum(r * r for r in res) / n),
           "mse": sum(res) / n, "median_abs": sorted(abs(r) for r in res)[n // 2] if n % 2 else
           sum(sorted(abs(r) for r in res)[n // 2 - 1:n // 2 + 1]) / 2,
           "same_sign": sum(1 for l, d in pairs if (l > 0) == (d > 0) and l != 0 and d != 0) / n,
           "within_band": sum(1 for r in res if abs(r) <= BAND) / n, "r": None, "spearman": None}
    if n >= 3:
        mx = sum(l for l, _ in pairs) / n
        my = sum(d for _, d in pairs) / n
        sxx = sum((l - mx) ** 2 for l, _ in pairs)
        syy = sum((d - my) ** 2 for _, d in pairs)
        sxy = sum((l - mx) * (d - my) for l, d in pairs)
        out["r"] = sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else None
        rl, rd = _ranks([l for l, _ in pairs]), _ranks([d for _, d in pairs])
        m = (n + 1) / 2.0
        a = sum((x - m) * (y - m) for x, y in zip(rl, rd))
        b = math.sqrt(sum((x - m) ** 2 for x in rl) * sum((y - m) ** 2 for y in rd))
        out["spearman"] = a / b if b > 0 else None
    return {k: (round(v, 4) if isinstance(v, float) else v) for k, v in out.items()}


def median(xs):
    xs = sorted(xs)
    n = len(xs)
    return xs[n // 2] if n % 2 else (xs[n // 2 - 1] + xs[n // 2]) / 2


# ------------------------------------------------------------------------------------------------ build --
def load_sign_check(path, papers):
    """Map (paper id, page, |printed value|) -> verdict for HM_FORM-like sign-check entries."""
    if not path or not os.path.exists(path):
        return {}, None
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    prefix = {p["pdf_sha256_prefix"]: pid for pid, p in papers.items() if p.get("pdf_sha256_prefix")}
    out = {}
    for s in d.get("records", []):
        pid = prefix.get((s.get("pdf_sha256") or "")[:16])
        v = s.get("printed_text_layer_value")
        if pid is None or not isinstance(v, (int, float)):
            continue
        out[(pid, s.get("page"), round(abs(float(v)), 6))] = s.get("verdict")
    return out, {"path": "open-tdb/build/sgte-corpus/dft-check/sign-check/sign_check.json",
                 "sha256": sha256_file(path), "schema": d.get("meta", {}).get("schema"),
                 "entries": len(d.get("records", [])), "entries_mapped_to_corpus_papers": len(out)}


def record_page(r):
    loc = r.get("locator") or {}
    return loc.get("page") if loc.get("page") is not None else r.get("quote_page")


OFF_STOICH_MAX = 0.02   # mole fraction; beyond it a measured composition is not the DFT compound (A4 E5)


def measured_offset(r, key):
    """(|dx| max over elements, measured composition) between a printed (not formula-derived) overall composition
    and the matched compound; (None, None) when the record carries no printed composition."""
    oc = r.get("overall_composition") or {}
    if not oc or r.get("overall_composition_basis") == "from_formula":
        return None, None
    comp = dict(key)
    meas = {k.capitalize(): float(v) for k, v in oc.items() if isinstance(v, (int, float))}
    if not meas:
        return None, None
    if len(meas) == len(comp) - 1:   # balance element implied
        rest = [e for e in comp if e not in meas]
        if len(rest) == 1:
            meas[rest[0]] = 1.0 - sum(meas.values())
    if set(meas) - set(comp):
        return None, None
    dx = max(abs(meas.get(e, 0.0) - x) for e, x in comp.items())
    return round(dx, 4), {e: round(meas.get(e, 0.0), 4) for e in sorted(comp)}


def sigma_composition(dx, meas, slopes):
    """|dx| x local hull slope (kJ/mol-atom), the slope of the hull segment on the side of the offset (Cu-Zr
    historical-formula rule, CHANGELOG 2026-10-24)."""
    if not dx or not slopes or not meas:
        return None
    xm = meas.get(slopes["x_element"].capitalize())
    if xm is None:
        return None
    sl = slopes["right"] if xm > slopes["x"] else slopes["left"]
    if sl is None:
        sl = slopes["left"] if sl is slopes["right"] else slopes["right"]
    return None if sl is None else round(abs(dx * sl), 3)


_LIQ_GAS_REF = re.compile(r"LIQ|GAS|\(G\)|\(L\)|^L$|VAPO|MELT", re.I)


def reference_state_mismatch(r, compound_elements):
    """Why a record's printed reference state is not the solid standard element reference that 0 K DFT formation
    energies use (audit round 3, POSTFIX A4-N2), or None: a liquid / gas element reference (Pt3Sn referred to liquid
    Pt and Sn; 'Na(l) + Sb(s)'), or a reference element that is not in the compound (a garbled 'IN' in Ca-Sb; the
    stored key, kept in reference_state_stored, counts)."""
    from espei_constraints import ELEMENT_TM_TB
    els = {e.upper() for e in compound_elements}
    why = []
    src_els = set(str(r.get("source_system") or "").upper().split("-")) | \
        set(str(r.get("system") or "").upper().split("-"))
    keys = [k.upper() for rs in (r.get("reference_state"), r.get("reference_state_stored")) if isinstance(rs, dict)
            for k in rs if re.fullmatch(r"[A-Z][a-z]?", k.capitalize()) and k.lower() != "printed"]
    if any(k not in els | src_els for k in keys):
        # keys of the paper's other elements (a state narrowed to its subsystem) are not a mismatch; an element of
        # neither the compound nor the paper's system is a garbled reference (Ca-Sb 'IN')
        why.append("reference element not in the compound")
    rs = r.get("reference_state") or {}
    T = r.get("temperature_K")
    for el, st in rs.items() if isinstance(rs, dict) else []:
        ph = str((st or {}).get("phase") if isinstance(st, dict) else st)
        if el.upper() in els and _LIQ_GAS_REF.search(ph) and not el.upper() in ("H", "N", "O", "F", "CL"):
            tm = (ELEMENT_TM_TB.get(el.upper()) or (None,))[0]
            why.append("liquid reference of %s%s" % (el.upper(), " below its melting point" if tm and isinstance(T, (int, float)) and T < tm else ""))
    return "; ".join(why) or None


OUTLIER_ABS = 50.0      # kJ/mol-atom: |DFT - lit| beyond it from every matched source
OUTLIER_SIGN = 10.0     # kJ/mol-atom: sign opposite to DFT with |DFT - lit| beyond it
BASIS_RULES = ("table_basis", "per_formula", "basis_correct")
SIGN_RULES = ("sign", "unsigned")


def as_extracted_value(r):
    """The record's value before the page-read basis / sign corrections (the first `before` of a value correction)."""
    for c in r.get("corrections") or []:
        if c.get("field") == "value" and isinstance(c.get("before"), (int, float)):
            return c["before"], c.get("rule") or ""
    return r.get("value"), None


def is_outlier(v, dfts):
    if not dfts or not isinstance(v, (int, float)):
        return False
    big = all(abs(d - v) > OUTLIER_ABS for d in dfts)
    sign = v != 0 and any((d < 0) != (v < 0) and abs(d - v) > OUTLIER_SIGN for d in dfts)
    return big or sign

def build(open_tdb, sign_check_path, hulls_root=None):
    with open(os.path.join(LIT, "papers.json"), encoding="utf-8") as fh:
        papers = json.load(fh)
    recs = []
    with gzip.open(os.path.join(LIT, "orchid-literature.jsonl.gz"), "rt", encoding="utf-8") as fh:
        for line in fh:
            r = json.loads(line)
            if r.get("property") == "HM_FORM":
                recs.append(r)
    lit_sha = sha256_file(os.path.join(LIT, "orchid-literature.jsonl.gz"))
    sign, sign_meta = load_sign_check(sign_check_path, papers)
    dft = {k: DFT(k, open_tdb) for k in SOURCES}
    hulls_root = hulls_root or HULLS
    hulls = {k: Hulls(hulls_root, k, dft[k].meta["snapshot_sha256"]) for k in SOURCES}

    excl = Counter()
    unit_counts = Counter()
    rows = []
    for r in sorted(recs, key=lambda x: x["id"]):
        unit_counts[r.get("unit")] += 1
        if r.get("unit") != UNIT:
            excl["unit is not kJ/mol-atom" + (" (basis unstated)" if r.get("unit") == "kJ/mol (basis unstated)"
                                                else "")] += 1
            continue
        ev = r.get("evidence_class")
        if ev not in CLASSES:
            excl["evidence class %s (literature DFT values are not compared with DFT)" % ev] += 1
            continue
        v = r.get("value")
        if isinstance(v, bool) or not isinstance(v, (int, float)):
            excl["no numeric value"] += 1
            continue
        system = [e.capitalize() for e in (r.get("elements") or r.get("system", "").split("-")) if e]
        counts, rule, note = None, None, ""
        reasons = []
        labels = [p.get("label") for p in (r.get("phases") or []) if p.get("label")]
        found = {}
        for lab in labels:
            c, n = parse_phase(lab, system)
            if c:
                found[reduced_key(c)] = (c, n, lab)
            else:
                reasons.append(n)
        variants = {}
        for lab in [lv for p in (r.get("phases") or []) for lv in (p.get("label_variants") or [])]:
            c, n = parse_phase(lab, system)
            if c:
                variants[reduced_key(c)] = (c, n, lab)
        if not found:
            # spellings of the same printed label by other readers (literature_gates G7 `label_variants`)
            found = dict(variants)
        elif variants:
            # a reader's spelling that agrees with the printed composition wins over a parse that does not
            # (A4 E4: "αNiSi_2" vs "αNiSi2" at x_Si 0.667)
            agree = {}
            for k, cand in {**found, **variants}.items():
                d = measured_offset(r, k)[0]
                if d is not None and d <= OFF_STOICH_MAX:
                    agree[k] = cand
            if len(agree) == 1 and (len(found) > 1 or next(iter(agree)) not in found):
                found = agree
        if len(found) == 1:
            key = next(iter(found))
            counts, note, label = found[key]
            rule = "phase_label"
        elif len(found) > 1:
            excl["unparseable phase: several distinct formulas among the record's phases"] += 1
            continue
        else:
            stoich = [pc for pc in (r.get("phase_compositions") or [])
                      if pc.get("basis") == "mole_fraction" and pc.get("composition_source") in ("stoichiometry",
                                                                                                    "formula_nominal")]
            comps = {tuple(sorted((k.capitalize(), round(float(x), 4)) for k, x in pc["composition"].items()
                                  if float(x) > 0)) for pc in stoich}
            if len(comps) == 1:
                key = next(iter(comps))
                if len(key) >= 2 and abs(sum(x for _, x in key) - 1) < 2e-3:
                    counts = {el: Fraction(x).limit_denominator(60) for el, x in key}
                    rule, label = "phase_composition_stoichiometry", (labels[0] if labels else "")
            if counts is None:
                why = sorted(set(reasons)) or ["no phase label"]
                excl["unparseable phase: " + why[0]] += 1
                continue
        if len(counts) > 2:
            excl["ternary or higher compound (DFT snapshots hold unary/binary entries only)"] += 1
            continue
        key = reduced_key(counts)
        P = papers.get(r.get("paper"), {})
        sv = sign.get((r.get("paper"), record_page(r), round(abs(float(r.get("value_printed", v))), 6))) \
            if isinstance(r.get("value_printed", v), (int, float)) else None
        flags = r.get("flags") or {}
        ref_mismatch = reference_state_mismatch(r, counts)
        gate_codes = sorted({g.get("code") for g in (r.get("gate_failures") or []) if g.get("code")})
        v0, v0_rule = as_extracted_value(r)
        sign_fixed = any(c.get("field") == "value" and any(k in (c.get("rule") or "") for k in SIGN_RULES)
                         for c in r.get("corrections") or []) or bool(flags.get("sign_pdf_page_verified")) \
            or bool(flags.get("sign_from_heading"))
        basis_fixed = any(c.get("field") == "value" and any(k in (c.get("rule") or "") for k in BASIS_RULES)
                          for c in r.get("corrections") or [])
        unsigned = bool(flags.get("unsigned_in_print")) or sv == "no_minus_printed"
        sign_lost = sv == "minus_printed" and v > 0 and not flags.get("sign_pdf_page_verified")
        cited = [{"ref": c.get("reference") or "", "doi": c.get("doi") or ""} for c in (r.get("cited_sources") or [])][:3]
        loc = r.get("locator") or {}
        labels_all = labels + [lv for p in (r.get("phases") or []) for lv in (p.get("label_variants") or [])]
        desig = structure_designators(labels_all, [p.get("model_phase") for p in (r.get("phases") or [])
                                                   if p.get("model_phase")], counts)
        dx, meas = measured_offset(r, key)
        state = r.get("chemical_system") or r.get("system")
        base = {
            "rid": r["id"], "sys": r.get("system"), "formula": formula_str(counts),
            # machine-learning keys exactly as the literature record carries them (A1 C5: one fold per record)
            "chemical_system": state, "fold_lso": r.get("fold_lso"), "fold_lpo": r.get("fold_lpo"),
            "designators": desig or None, "dx_measured": dx, "comp_measured": meas,
            "elements": sorted(counts), "comp": {el: x for el, x in key}, "phase": label, "match": rule,
            "match_note": note or None, "lit": round(float(v), 4), "ev": ev,
            "acc": bool(flags.get("extraction_accepted")), "unsigned": unsigned,
            "sign_check": sv, "sign_lost": sign_lost,
            "sign_verified": bool(flags.get("sign_pdf_page_verified")),
            "suspect": [k for k in ("magnitude_suspect", "system_mismatch", "source_disputed") if flags.get(k)],
            "method": r.get("method"), "mcls": method_class(r.get("method")), "ccls": chem_class(counts),
            "T": r.get("temperature_K"), "paper": r.get("paper"), "cited": cited,
            "quote": r.get("quote"), "page": record_page(r), "kind": r.get("record_kind"),
            "table": ("Table %s" % loc["label"]) if loc.get("artifact") == "table" and loc.get("label") else None,
            "own": (r.get("origin") or {}).get("own_work_of_paper"),
            "qfix": next((k for k in ("basis_corrected_from_formula", "sign_from_heading", "dft_outlier_unexplained")
                          if flags.get(k)), None),
            "basis": r.get("unit_basis"), "basis_inferred": bool(flags.get("basis_inferred_from_magnitude")),
            "ref_mismatch": ref_mismatch, "gate_failures": gate_codes or None,
            "lit_as_extracted": round(float(v0), 4) if isinstance(v0, (int, float)) and v0 != v else None,
            "basis_fixed": basis_fixed, "sign_fixed": sign_fixed,
        }
        matched = {}
        for s, src in dft.items():
            m = src.lookup(key, set(desig), hulls[s])
            if m:
                matched[s] = m
        if not matched:
            excl["no DFT entry at this composition in either source"] += 1
            continue
        for s in SOURCES:
            if s not in matched:
                excl["no %s entry at this composition (other source matched)" % SOURCES[s]["label"]] += 1
        rows.append((base, matched))

    points, used_papers, records_out = [], set(), {}
    DFT_KEYS = ("src", "dft", "dft_eV", "mid", "sg", "ehull", "ehull_upstream", "npoly", "theoretical", "res",
                "smatch", "ground_id", "sigma_comp", "default", "excluded_by", "outlier_class")

    def record_default(p, classes=("experimental",)):
        return (p["ev"] in classes and p["acc"] and not p["unsigned"] and not p["sign_lost"]
                and "magnitude_suspect" not in p["suspect"])

    def point_exclusions(p):
        """Point-level reasons that keep an otherwise default record out of the default (calibration) subset."""
        why = []
        if p["smatch"] == "structure_unmatched":
            why.append("structure_unmatched")
        if p["dx_measured"] is not None and p["dx_measured"] > OFF_STOICH_MAX:
            why.append("off_stoichiometry" if p["dx_measured"] <= 0.1 else "composition_conflicts_label")
        if p["ev"] == "experimental" and p["T"] == 0:
            why.append("experimental_value_at_0K")
        if p.get("ref_mismatch"):
            why.append("reference_state")
        if p.get("gate_failures"):
            why.append("gate_flagged")
        return why

    for base, matched in rows:
        used_papers.add(base["paper"])
        records_out[base["rid"]] = {k: v for k, v in base.items() if k != "rid"}
        for s in SOURCES:
            if s in matched:
                m = matched[s]
                p = dict(base)
                p.update({"src": s, "dft": m["dft"], "dft_eV": m["ef_eV"], "mid": m["id"], "sg": m["sg"],
                          "ehull": m["ehull"], "ehull_upstream": m["ehull_upstream"], "npoly": m["npoly"],
                          "theoretical": m["theoretical"], "res": round(m["dft"] - base["lit"], 3),
                          "smatch": m["smatch"], "ground_id": m["ground_id"],
                          "sigma_comp": sigma_composition(base["dx_measured"], base["comp_measured"], m["hull_slopes"])})
                p["excluded_by"] = point_exclusions(p) or None
                p["default"] = record_default(p) and not p["excluded_by"]
                points.append(p)
    points.sort(key=lambda p: (p["src"], p["rid"]))

    # ---- aggregate statistics
    def default(p, classes=("experimental",)):
        return record_default(p, classes) and not p["excluded_by"]

    # ---- one outlier classification (audit round 3): population = accepted experimental records not flagged
    # magnitude_suspect (the default subset before its point exclusions and before the unsigned / sign-lost exclusion,
    # so that the classes partition one set); a record is examined when its value as extracted (before the page-read
    # basis / sign corrections) differs from every matched DFT source by more than OUTLIER_ABS, or has the opposite
    # sign with |DFT - lit| > OUTLIER_SIGN for a source. Classes, in this order: per_formula (a basis correction read on
    # the page), unsigned (printed without a sign, a lost sign, or a sign read on the page), unexplained (the rest).
    by_rec = defaultdict(list)
    for p in points:
        by_rec[p["rid"]].append(p)
    outl = {}
    for rid, ps in by_rec.items():
        p0 = ps[0]
        if not (p0["ev"] == "experimental" and p0["acc"] and "magnitude_suspect" not in p0["suspect"]):
            continue
        v0 = p0["lit_as_extracted"] if p0["lit_as_extracted"] is not None else p0["lit"]
        dfts = [q["dft"] for q in ps]
        if not is_outlier(v0, dfts):
            continue
        if p0["basis_fixed"]:
            c = "per_formula"
        elif p0["unsigned"] or p0["sign_lost"] or p0["sign_fixed"]:
            c = "unsigned"
        else:
            c = "unexplained"
        outl[rid] = c
        for q in ps:
            q["outlier_class"] = c
    for p in points:
        p.setdefault("outlier_class", None)
    oc = Counter(outl.values())
    unexpl_default = {rid for rid, c in outl.items() if c == "unexplained" and any(default(q) for q in by_rec[rid])}
    outlier_classification = {
        "population": "accepted experimental HM_FORM records matched to MP and/or OQMD, not flagged magnitude_suspect "
                      "(the default subset before its point exclusions and before the unsigned / sign-lost exclusion)",
        "criterion": "value as extracted (before the page-read basis / sign corrections) differs from every matched DFT "
                     "source by more than %g kJ/mol-atom, or has the opposite sign of a source with |DFT - lit| > %g"
                     % (OUTLIER_ABS, OUTLIER_SIGN),
        "classes": "per_formula: a basis correction read on the page; unsigned: printed without a sign, a lost sign or "
                   "a sign read on the page; unexplained: no deterministic explanation (kept as printed)",
        "examined": len(outl), "per_formula": oc.get("per_formula", 0), "unsigned": oc.get("unsigned", 0),
        "unexplained": oc.get("unexplained", 0),
        "unexplained_in_default_subset": len(unexpl_default),
        "unexplained_excluded_by_point_rules": oc.get("unexplained", 0) - len(unexpl_default),
    }
    assert outlier_classification["per_formula"] + outlier_classification["unsigned"] + \
        outlier_classification["unexplained"] == outlier_classification["examined"]

    subsets = {
        "experimental": ("Experimental, accepted, signed, not magnitude-suspect (default)", lambda p: default(p)),
        "experimental_model": ("+ model values", lambda p: default(p, ("experimental", "model"))),
        "experimental_model_unspecified": ("+ model + unspecified", lambda p: default(p, CLASSES)),
        "experimental_with_unsigned": ("Experimental incl. unsigned-in-print",
                                       lambda p: p["ev"] == "experimental" and p["acc"]
                                       and "magnitude_suspect" not in p["suspect"] and not p["excluded_by"]),
        "experimental_incl_structure_unmatched_offstoich": (
            "Default incl. structure-unmatched, off-stoichiometric and 0 K experimental points (diagnostic only)",
            lambda p: record_default(p)),
        "experimental_with_magnitude_suspect": ("Experimental incl. magnitude-suspect extractions",
                                                lambda p: p["ev"] == "experimental" and p["acc"]
                                                and not p["unsigned"] and not p["sign_lost"]
                                                and not p["excluded_by"]),
        "experimental_explicit_basis": ("Default, per-atom basis printed (not inferred from magnitude)",
                                        lambda p: default(p) and not p["basis_inferred"]),
        "experimental_excl_unexplained_outliers": ("Default without the unexplained outliers of the outlier "
                                                   "classification (diagnostic only)",
                                                   lambda p: default(p) and p["outlier_class"] != "unexplained"),
        "experimental_incl_not_accepted": ("Experimental incl. not-accepted extractions",
                                           lambda p: p["ev"] == "experimental" and not p["unsigned"]
                                           and not p["sign_lost"] and "magnitude_suspect" not in p["suspect"]
                                           and not p["excluded_by"]),
    }
    for c in METHOD_LABELS:
        subsets["experimental_method_" + c] = ("Experimental, method: " + METHOD_LABELS[c],
                                               (lambda c: lambda p: default(p) and p["mcls"] == c)(c))
    for c in CHEM_LABELS:
        subsets["experimental_chem_" + c] = ("Experimental, " + CHEM_LABELS[c],
                                             (lambda c: lambda p: default(p) and p["ccls"] == c)(c))
    st = {}
    for s in SOURCES:
        st[s] = {}
        for k, (label, f) in subsets.items():
            sel = [p for p in points if p["src"] == s and f(p)]
            st[s][k] = dict(stats([(p["lit"], p["dft"]) for p in sel]), label=label)
        # per-compound median of the default subset
        g = defaultdict(list)
        for p in points:
            if p["src"] == s and default(p):
                # one compound = one reduced composition (UC and U0.5C0.5 are one compound; A4 E7)
                g[tuple(sorted(p["comp"].items()))].append(p)
        med = [(median([q["lit"] for q in v]), v[0]["dft"]) for v in g.values()]
        st[s]["experimental_per_compound_median"] = dict(stats(med), label="Experimental, median per compound")
        # largest outliers of the default subset
        sel = sorted((p for p in points if p["src"] == s and default(p)), key=lambda p: (-abs(p["res"]), p["rid"]))
        st[s]["largest_outliers"] = [{k: p[k] for k in ("rid", "formula", "sys", "lit", "dft", "res", "mid", "sg",
                                                        "ehull", "npoly", "method", "paper")} for p in sel[:10]]

    paper_out = {}
    for pid in sorted(used_papers):
        P = papers.get(pid, {})
        paper_out[pid] = {"citation": P.get("citation") or pid, "doi": P.get("doi"), "year": P.get("year")}

    n_matched_rec = len(rows)
    out = {
        "schema": "orchid-dft-parity-v1",
        "cite_as": orchid_cite.cite_as("dft-benchmark"),
        "description": "Literature formation enthalpies (HM_FORM, kJ/mol-atom) against DFT formation energies at "
                       "the same composition: the lowest entry of the space group the record names (Strukturbericht "
                       "designator or prototype), else the ground state (lowest-energy entry). Materials Project and "
                       "OQMD are separate sources with different energy references and are never averaged. Residual "
                       "= DFT - literature.",
        "unit": UNIT, "kJ_per_mol_per_eV": KJ_PER_EV, "band_kJ_per_mol_atom": BAND,
        "inputs": {"literature": {"path": "data/literature/orchid-literature.jsonl.gz", "sha256": lit_sha},
                   "sign_check": sign_meta},
        "dft_sources": {s: dft[s].meta for s in SOURCES},
        "rules": {
            "selection": "property HM_FORM; unit kJ/mol-atom; evidence class experimental, model or unspecified. "
                         "'kJ/mol (basis unstated)' is excluded (per-atom vs per-formula cannot be settled).",
            "formula": "Parsed from the record's printed phase label(s): Greek and polymorph prefixes stripped, "
                       "formula groups expanded, sublattice brackets read, tokens outside the record's system "
                       "dropped; several distinct formulas in one record are refused. Fallback: a phase "
                       "composition derived from stoichiometry.",
            "match": "Reduced composition (mole fractions rounded to 4 decimals) matched exactly to non-deprecated "
                     "DFT entries; the lowest entry of the named space group (see structure), else the lowest "
                     "formation energy at that composition, is used; npoly = number of "
                     "entries at that composition. Binary compounds only (snapshots hold unary/binary entries).",
            "units": "DFT eV/atom x 96.485 = kJ/mol-atom; ehull = energy above ORCHID's released hull of the same "
                     "snapshot (hulls/dft-<src>/, kJ/mol-atom; one source of truth for MP and OQMD); ehull_upstream = "
                     "the snapshot's upstream value (MP energy_above_hull; OQMD max(stability, 0), computed upstream "
                     "against a different phase set).",
            "structure": "dft_structure_match (`smatch`): exact = the record names a structure (Strukturbericht "
                         "designator or 'X-type' prototype in a printed label, or in a model-phase name whose formula "
                         "is the compound; A1/A2/A3 ignored) and the lowest entry of that space group is used "
                         "(`ground_id` = the ground state when it differs); structure_unmatched = named structure "
                         "absent from the snapshot, paired with the ground state, not in the default subset; "
                         "no_structure_printed = ground state.",
            "composition_offset": "dx_measured = max |x_printed - x_formula| for a printed (not formula-derived) "
                                  "overall composition; sigma_comp = |dx| x local hull slope (kJ/mol-atom, the "
                                  "Cu-Zr historical-formula rule); dx > %.2f is not in the default subset "
                                  "(off_stoichiometry, or composition_conflicts_label beyond 0.1)." % OFF_STOICH_MAX,
            "default_subset": "evidence class experimental, extraction accepted, not unsigned-in-print, no "
                              "uncorrected sign loss found by the page-image sign check, not flagged "
                              "magnitude_suspect by the extraction checks (a flag set before and independently "
                              "of this DFT comparison); per point (`default`): not structure_unmatched, measured "
                              "composition within 0.02 of the compound, not an experimental value stored at 0 K, "
                              "the record's reference state is the solid standard element reference that 0 K DFT "
                              "uses (no liquid / gas element reference, no reference element outside the compound; "
                              "reference_state), and the record carries no literature gate failure (gate_flagged) "
                              "(`excluded_by` lists the reason; audit round 3).",
            "layout": "records{rid: literature fields} + points[] (one per record x DFT source, joined by rid) "
                      "+ papers{paper id: citation, DOI}.",
            "quality_fixes": "data/literature/QUALITY_FIXES.md: the per-formula / per-atom basis of a column is decided "
                             "on the printed page only (gate R2 table_basis; the earlier rule that converted a value "
                             "when the converted value agreed with DFT was withdrawn in audit round 2 because it used "
                             "DFT to correct the values that calibrate DFT); unsigned magnitudes are flagged "
                             "unsigned_in_print (or signed from a printed negated heading); see "
                             "counts.outlier_classification for the classes of the large residuals.",
            "unsigned": "flags.unsigned_in_print, or a sign-check verdict no_minus_printed matched by paper "
                        "(PDF hash prefix), page and |printed value|: such values are magnitudes by the paper's "
                        "convention and are kept out of the default set.",
            "method_class": {k: v for k, v in METHOD_LABELS.items()},
            "chem_class": {k: v for k, v in CHEM_LABELS.items()},
        },
        "counts": {
            "hm_form_records": len(recs), "hm_form_by_unit": dict(sorted((str(k), v) for k, v in unit_counts.items())),
            "records_matched": n_matched_rec,
            "records_matched_by_class": dict(sorted(Counter(b["ev"] for b, _ in rows).items())),
            "points": len(points),
            "points_by_source": {s: sum(1 for p in points if p["src"] == s) for s in SOURCES},
            "compounds": len({b["formula"] for b, _ in rows}),
            "unsigned_points": sum(1 for p in points if p["unsigned"]),
            "sign_lost_points": sum(1 for p in points if p["sign_lost"]),
            "default_points_by_source": {s: sum(1 for p in points if p["src"] == s and p["default"]) for s in SOURCES},
            "structure_match_by_source": {s: dict(sorted(Counter(p["smatch"] for p in points if p["src"] == s).items()))
                                          for s in SOURCES},
            "excluded_from_default_by_source": {s: dict(sorted(Counter(x for p in points if p["src"] == s
                                                                       and record_default(p)
                                                                       for x in (p["excluded_by"] or [])).items()))
                                                for s in SOURCES},
            "exclusions": dict(sorted(excl.items(), key=lambda kv: (-kv[1], kv[0]))),
            "outlier_classification": outlier_classification,
        },
        "stats": st,
        "papers": paper_out,
        "records": records_out,
        "points": [dict({k: p[k] for k in DFT_KEYS}, rid=p["rid"]) for p in points],
    }
    return (json.dumps(out, sort_keys=True, ensure_ascii=False, separators=(",", ":")) + "\n").encode()


def croissant_update(path, data):
    with open(path, encoding="utf-8") as fh:
        d = json.load(fh)
    dist = [x for x in d.get("distribution", []) if x.get("@id") != "dft-parity"]
    # insert before the literature entries, which tools/export_literature_data.py always rewrites at the end
    at = next((i for i, x in enumerate(dist) if str(x.get("@id", "")).startswith("literature")), len(dist))
    dist.insert(at, {"@type": "cr:FileObject", "@id": "dft-parity", "name": "dft_parity.json",
                 "description": "Literature formation enthalpies (kJ/mol-atom) paired with the Materials Project "
                                "and OQMD entry of the named structure, else the ground state, at the same "
                                "composition (never averaged), "
                                "with per-source statistics. Viewer: https://orchid-project.org/dft-parity/",
                 "contentUrl": "https://orchid-project.org/data/literature/dft_parity.json",
                 "encodingFormat": "application/json", "sha256": hashlib.sha256(data).hexdigest()})
    d["distribution"] = dist
    return json.dumps(d, indent=1, ensure_ascii=False).encode()


def main(argv=None):
    global OUT, LIT
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--open-tdb", default=OPEN_TDB, help="open-tdb checkout holding hulls/snapshot*/ (read only)")
    ap.add_argument("--sign-check", default=os.path.join(OPEN_TDB, "build", "sgte-corpus", "dft-check", "sign-check",
                                                          "sign_check.json"))
    ap.add_argument("--hulls", default=HULLS, help="released hulls root holding dft-mp/ and dft-oqmd/ (read only)")
    ap.add_argument("--lit", default=LIT, help="literature release directory read (default: %(default)s)")
    ap.add_argument("--no-croissant", action="store_true", help="do not touch data/croissant.json")
    ap.add_argument("--out", default=OUT, help="output dft_parity.json (default: %(default)s)")
    ap.add_argument("--check", action="store_true", help="exit 1 if the files on disk are stale")
    a = ap.parse_args(argv)
    OUT, LIT = a.out, a.lit
    data = build(a.open_tdb, a.sign_check, a.hulls)
    cro = croissant_update(CROISSANT, data) if os.path.exists(CROISSANT) and not a.no_croissant else None
    if a.check:
        bad = []
        if not os.path.exists(OUT) or open(OUT, "rb").read() != data:
            bad.append(os.path.relpath(OUT, SITE))
        if cro is not None and open(CROISSANT, "rb").read() != cro:
            bad.append(os.path.relpath(CROISSANT, SITE))
        if bad:
            print("stale: " + ", ".join(bad), file=sys.stderr)
            return 1
        print("dft_parity.json is current")
        return 0
    with open(OUT, "wb") as fh:
        fh.write(data)
    if cro is not None:
        with open(CROISSANT, "wb") as fh:
            fh.write(cro)
    d = json.loads(data)
    print(json.dumps(d["counts"], indent=1))
    for s in SOURCES:
        x = d["stats"][s]["experimental"]
        print(s, x)
    return 0


if __name__ == "__main__":
    sys.exit(main())
