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

Matching: the record's phase label is parsed into a formula (Greek / polymorph prefixes, sublattice brackets
and formula groups handled; tokens outside the record's system dropped); if no label parses, a phase
composition derived from stoichiometry is used. The reduced composition is matched to the lowest-energy
non-deprecated entry at exactly that composition (4-decimal mole fractions) in each DFT snapshot. MP and
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

SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OPEN_TDB = os.path.expanduser("~/projects/open-tdb")
LIT = os.path.join(SITE, "data", "literature")
OUT = os.path.join(LIT, "dft_parity.json")
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
    notes = []
    if re.search(r"\((g|gas|l|liq)\)", s, re.I):
        return None, "gas or liquid state marker in phase label"
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


def formula_str(counts):
    """Reduced formula: divide by the GCD of integer counts (non-integer counts kept as printed ratios)."""
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

    def lookup(self, key):
        cands = self.by_comp.get(key)
        if not cands:
            return None
        best = min(cands, key=lambda e: (float(e["formation_energy_per_atom"]), e["material_id"]))
        eh = best.get("energy_above_hull")
        return {"id": best["material_id"], "ef_eV": round(float(best["formation_energy_per_atom"]), 6),
                "dft": round(float(best["formation_energy_per_atom"]) * KJ_PER_EV, 3),
                "sg": best.get("spacegroup_symbol"), "ehull": None if eh is None else round(float(eh) * KJ_PER_EV, 3),
                "npoly": len(cands), "theoretical": best.get("theoretical")}


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
           "within_band": sum(1 for r in res if abs(r) <= BAND) / n, "r": None}
    if n >= 3:
        mx = sum(l for l, _ in pairs) / n
        my = sum(d for _, d in pairs) / n
        sxx = sum((l - mx) ** 2 for l, _ in pairs)
        syy = sum((d - my) ** 2 for _, d in pairs)
        sxy = sum((l - mx) * (d - my) for l, d in pairs)
        out["r"] = sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else None
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


def build(open_tdb, sign_check_path):
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
        unsigned = bool(flags.get("unsigned_in_print")) or sv == "no_minus_printed"
        sign_lost = sv == "minus_printed" and v > 0 and not flags.get("sign_pdf_page_verified")
        cited = [{"ref": c.get("reference") or "", "doi": c.get("doi") or ""} for c in (r.get("cited_sources") or [])][:3]
        loc = r.get("locator") or {}
        base = {
            "rid": r["id"], "sys": r.get("system"), "formula": formula_str(counts),
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
        }
        matched = {}
        for s, src in dft.items():
            m = src.lookup(key)
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
    DFT_KEYS = ("src", "dft", "dft_eV", "mid", "sg", "ehull", "npoly", "theoretical", "res")
    for base, matched in rows:
        used_papers.add(base["paper"])
        records_out[base["rid"]] = {k: v for k, v in base.items() if k != "rid"}
        for s in SOURCES:
            if s in matched:
                m = matched[s]
                p = dict(base)
                p.update({"src": s, "dft": m["dft"], "dft_eV": m["ef_eV"], "mid": m["id"], "sg": m["sg"],
                          "ehull": m["ehull"], "npoly": m["npoly"], "theoretical": m["theoretical"],
                          "res": round(m["dft"] - base["lit"], 3)})
                points.append(p)
    points.sort(key=lambda p: (p["src"], p["rid"]))

    # ---- aggregate statistics
    def default(p, classes=("experimental",)):
        return (p["ev"] in classes and p["acc"] and not p["unsigned"] and not p["sign_lost"]
                and "magnitude_suspect" not in p["suspect"])

    subsets = {
        "experimental": ("Experimental, accepted, signed, not magnitude-suspect (default)", lambda p: default(p)),
        "experimental_model": ("+ model values", lambda p: default(p, ("experimental", "model"))),
        "experimental_model_unspecified": ("+ model + unspecified", lambda p: default(p, CLASSES)),
        "experimental_with_unsigned": ("Experimental incl. unsigned-in-print",
                                       lambda p: p["ev"] == "experimental" and p["acc"]
                                       and "magnitude_suspect" not in p["suspect"]),
        "experimental_with_magnitude_suspect": ("Experimental incl. magnitude-suspect extractions",
                                                lambda p: p["ev"] == "experimental" and p["acc"]
                                                and not p["unsigned"] and not p["sign_lost"]),
        "experimental_incl_not_accepted": ("Experimental incl. not-accepted extractions",
                                           lambda p: p["ev"] == "experimental" and not p["unsigned"]
                                           and not p["sign_lost"] and "magnitude_suspect" not in p["suspect"]),
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
                g[p["formula"]].append(p)
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
        "description": "Literature formation enthalpies (HM_FORM, kJ/mol-atom) against the lowest-energy DFT "
                       "entry at the same composition. Materials Project and OQMD are separate sources with "
                       "different energy references and are never averaged. Residual = DFT - literature.",
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
                     "DFT entries; the lowest formation energy at that composition is used; npoly = number of "
                     "entries at that composition. Binary compounds only (snapshots hold unary/binary entries).",
            "units": "DFT eV/atom x 96.485 = kJ/mol-atom; ehull = upstream energy above hull in kJ/mol-atom "
                     "(OQMD: derived by the snapshot, not an upstream flag).",
            "default_subset": "evidence class experimental, extraction accepted, not unsigned-in-print, no "
                              "uncorrected sign loss found by the page-image sign check, not flagged "
                              "magnitude_suspect by the extraction checks (a flag set before and independently "
                              "of this DFT comparison).",
            "layout": "records{rid: literature fields} + points[] (one per record x DFT source, joined by rid) "
                      "+ papers{paper id: citation, DOI}.",
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
            "exclusions": dict(sorted(excl.items(), key=lambda kv: (-kv[1], kv[0]))),
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
                 "description": "Literature formation enthalpies (kJ/mol-atom) paired with the lowest-energy "
                                "Materials Project and OQMD entry at the same composition (never averaged), "
                                "with per-source statistics. Viewer: https://orchid-project.org/dft-parity/",
                 "contentUrl": "https://orchid-project.org/data/literature/dft_parity.json",
                 "encodingFormat": "application/json", "sha256": hashlib.sha256(data).hexdigest()})
    d["distribution"] = dist
    return json.dumps(d, indent=1, ensure_ascii=False).encode()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--open-tdb", default=OPEN_TDB, help="open-tdb checkout holding hulls/snapshot*/ (read only)")
    ap.add_argument("--sign-check", default=os.path.join(OPEN_TDB, "build", "sgte-corpus", "dft-check", "sign-check",
                                                          "sign_check.json"))
    ap.add_argument("--check", action="store_true", help="exit 1 if the files on disk are stale")
    a = ap.parse_args(argv)
    data = build(a.open_tdb, a.sign_check)
    cro = croissant_update(CROISSANT, data) if os.path.exists(CROISSANT) else None
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
