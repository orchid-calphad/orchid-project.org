#!/usr/bin/env python3
"""Sanity check of the model phase that the SGTE-corpus builder attached to each printed literature phase label
(ORCHID audit round 2, A2 class 1, 2026-10-24; called by tools/export_literature_data.py after the literature gates).

The corpus mapping sometimes sends an intermediate phase named only by a Greek letter (eps, beta, delta, phi ...) or a
Greek-prefixed compound (betaNi3Ge, eps-Cu3Sn) through an element-polymorph table to a terminal-solution model
(BCC_A2, TETRAGONAL_A6, GRAPHITE_A9 ...) and writes that model's Pearson symbol / prototype beside the label, which the
page does not print. Rules (a rejected mapping is removed; the model phase, its structure and the reason are kept in
`model_phase_rejected`; the printed label is never changed):

  R1 element-structure model (FCC_A1, BCC_A2, HCP_A3, DIAMOND_A4 ... GRAPHITE_A9, DHCP, A1_FCC ...): kept only when the
     printed label names a terminal solution or pure element ("(Fe)", "(alphaTi)", "gamma-Fe", "(Cu,Ni)"), names the
     structure itself ("beta (bcc)", "A2", "fcc"), or is an interstitial compound M X_n (X = C, N, H, O, B) that the
     target model holds on an interstitial (X, VA) sublattice (SGTE convention, e.g. TiC -> FCC_A1, V2C -> HCP_A3).
     A bare Greek letter is kept only when it is the conventional name of a polymorph of one of the system's elements
     with that structure (no other Greek -> structure fallback; Hume-Rothery names such as beta brass are left unmapped) (alpha-Ti hcp, beta-Zr bcc, gamma-Fe fcc, delta-Pu fcc ...; POLYMORPH); a Greek-prefixed
     compound, a primed / numbered Greek label or any other compound is never mapped to an element structure; an
     element label whose Greek polymorph contradicts the model (beta-B -> BCC_A2) is rejected.
  R2 a structure word in the label that contradicts the model (label "bcc", model FCC_A1) rejects the mapping.
  R3 a label that names a liquid ("L(phi2)", "L + ...") is never mapped to a non-liquid model.
  R3b a Greek / compound label mapped to the same model as a terminal-solution label of the same row ("(deltaMn)"
     and "eps" both BCC_A2) is rejected: two printed phases of one equilibrium are not one model phase.
  R4 a printed mole-fraction composition of the phase must lie inside the model phase's calculated stability range in
     the system (data/phases stability rows, any assessment) +- 0.1; otherwise the mapping is rejected.

Structure fields: every literature `pearson_symbol` / `prototype` comes from the mapped model phase (none is printed in
the table row), so each phase that keeps one gets `structure_basis: "mapped_model_phase"`.
Standard library only; deterministic.
"""
from __future__ import annotations

import re
from collections import Counter

RELEASE = "phase-label-check 2026-10-24 (audit round 2, A2 class 1)"
X_TOL = 0.10
GREEK = "αβγδεζηθικλμνξοπρστυφχψωϕϵ"
# disordered element-structure (terminal solution) models only; ordered variants (BCC_B2, FCC_L12, BCC_4SL ...) are not
ELEMENT_MODEL = re.compile(r"^(?:(?:FCC|BCC|HCP|DIAMOND|BCT|TETRAGONAL|RHOMBOHEDRAL|HEXAGONAL|GRAPHITE|CBCC|CUB|ORTHORHOMBIC)"
                           r"(?:_A\d{1,2})?|DHCP|GRAPHITE|DIAMOND|A\d{1,2}_(?:FCC|BCC|HCP|DIA|DIAMOND|BCT|GRAPHITE))$")
# Greek polymorph names of the elements (Massalski / Okamoto convention) -> structure families of that polymorph
POLYMORPH = {}
for _els, _tab in (
        ("FE", {"α": "BCC", "γ": "FCC", "δ": "BCC", "ε": "HCP"}),
        ("CO", {"α": "FCC HCP", "ε": "HCP"}),
        ("TI ZR HF SC Y GD TB DY HO LU TL BE", {"α": "HCP", "β": "BCC"}),
        ("ER TM", {"α": "HCP"}),
        ("LA", {"α": "DHCP", "β": "FCC", "γ": "BCC"}),
        ("CE", {"α": "FCC", "β": "DHCP", "γ": "FCC", "δ": "BCC"}),
        ("PR ND PM", {"α": "DHCP", "β": "BCC"}),
        ("SM", {"α": "RHOMBOHEDRAL", "β": "HCP", "γ": "BCC"}),
        ("YB", {"α": "HCP", "β": "FCC", "γ": "BCC"}),
        ("CA", {"α": "FCC", "β": "BCC"}), ("SR", {"α": "FCC", "β": "BCC HCP", "γ": "BCC"}),
        ("TH", {"α": "FCC", "β": "BCC"}), ("AM", {"α": "DHCP", "β": "FCC", "γ": "BCC"}),
        ("U NP", {"α": "ORTHORHOMBIC", "β": "TETRAGONAL", "γ": "BCC"}),
        ("PU", {"δ": "FCC", "δ'": "TETRAGONAL BCT", "ε": "BCC"}),
        ("MN", {"α": "CBCC", "β": "CUB", "γ": "FCC", "δ": "BCC"}),
        ("SN", {"α": "DIAMOND", "β": "BCT"}), ("B", {"α": "RHOMBO", "β": "RHOMBO"}), ("C", {"α": "GRAPHITE"})):
    for _e in _els.split():
        POLYMORPH[_e] = {g: set(f.split()) for g, f in _tab.items()}
LIQUID_MODEL = re.compile(r"^(?:LIQUID|LIQ|IONIC_LIQ|LIQUID_IONIC|ION|IONIC_LIQUID|KAP_LIQ|ASSOCI|OXIDE_L|SLAG)(?:_.*)?$")
STRUCT_WORDS = {"cph": "HCP", "dia": "DIAMOND", "rhom": "RHOMBOHEDRAL", "rhombo": "RHOMBOHEDRAL", "cub": "CUBIC",
                "cbcc": "CBCC", "rho": "RHOMBOHEDRAL", "rhomb": "RHOMBOHEDRAL", "rhombohedral": "RHOMBOHEDRAL", "hex": "HEX",
                "hexagonal": "HEX", "tetragonal": "TETRAGONAL", "tet": "TETRAGONAL", "orthorhombic": "ORTHORHOMBIC",
                "ortho": "ORTHORHOMBIC", "orth": "ORTHORHOMBIC", "ort": "ORTHORHOMBIC", "cubic": "CUBIC", "graph": "GRAPHITE", "gr": "GRAPHITE", "fcc": "FCC", "bcc": "BCC", "hcp": "HCP", "dhcp": "DHCP", "bct": "BCT", "diamond": "DIAMOND",
                "graphite": "GRAPHITE", "a1": "FCC", "a2": "BCC", "a3": "HCP", "a3'": "DHCP", "a4": "DIAMOND",
                "a5": "BCT", "a6": "TETRAGONAL", "a7": "RHOMBOHEDRAL", "a8": "HEXAGONAL", "a9": "GRAPHITE",
                "a12": "CBCC", "a13": "CUB", "a20": "ORTHORHOMBIC", "cf4": "FCC", "ci2": "BCC", "hp2": "HCP",
                "hp4": "DHCP", "ti4": "BCT", "cf8": "DIAMOND"}
MODEL_FAMILY = [("CBCC", r"CBCC|^A12(_|$)"), ("CUB", r"^CUB(_|$)|^A13(_|$)"), ("DHCP", r"DHCP"),
                ("ORTHORHOMBIC", r"ORTHORHOMBIC|^A20(_|$)|^A11(_|$)"), ("RHOMBOHEDRAL", r"^RHOMB|^A7(_|$)"),
                ("FCC", r"FCC|^A1(_|$)"), ("BCC", r"BCC|^A2(_|$)"), ("HCP", r"HCP|^A3(_|$)"),
                ("DIAMOND", r"DIAMOND|^A4(_|$)|_DIA$"), ("BCT", r"BCT|^A5(_|$)"), ("TETRAGONAL", r"TETRAGONAL|^A6(_|$)"),
                ("HEXAGONAL", r"HEXAGONAL|^A8(_|$)"), ("GRAPHITE", r"GRAPHITE|^A9(_|$)")]
INTERSTITIAL = {"C", "N", "H", "O", "B"}
_STRIP = re.compile(r"[\s()\[\]'’′″*\-–‐,.:;/]+")
_ELEM = re.compile(r"[A-Z][a-z]?")


def family(model):
    m = model.upper()
    for fam, pat in MODEL_FAMILY:
        if re.search(pat, m):
            return fam
    return None


GREEK_NAMES = {"alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε", "zeta": "ζ", "eta": "η",
               "theta": "θ", "kappa": "κ", "lambda": "λ", "mu": "μ", "sigma": "σ", "tau": "τ", "phi": "φ", "chi": "χ",
               "psi": "ψ", "omega": "ω"}
_GREEK_WORD = re.compile(r"(?i)(?<![a-z])(" + "|".join(sorted(GREEK_NAMES, key=len, reverse=True)) + r")(?![a-z])|"
                         r"^(" + "|".join(sorted(GREEK_NAMES, key=len, reverse=True)) + r")(?=[A-Z])")
MAGNETIC_SUFFIX = re.compile(r"(?<=\))\s*(?:ferro|para|pm|fm|af|afm|mag|nm)\s*$|\s*(?:FM|PM|AFM)\s*$", re.I)


def normalize_label(label):
    """Spelled Greek letters -> symbols (gammaFe -> γFe), micro sign -> mu, magnetic suffixes and underscores dropped."""
    s = str(label or "").replace("µ", "μ").replace("ϕ", "φ").replace("ϵ", "ε")
    s = _GREEK_WORD.sub(lambda m: GREEK_NAMES[(m.group(1) or m.group(2)).lower()], s)
    s = MAGNETIC_SUFFIX.sub("", s)
    s = re.sub(r"(?<=[%s])_" % GREEK, "", s)
    return s.strip()


def _greek_free(label):
    return re.sub("[%s]" % GREEK, "", label)


def struct_words(label):
    """Structure families named by words / Strukturbericht / A-type Pearson tokens in the label."""
    s = re.sub(r"(?i)\b([bfh])\.\s*c\.\s*([cpt])\.", r"\1c\2", str(label))      # b.c.c., h.c.p. -> bcc, hcp
    toks = re.findall(r"[A-Za-z]+\d*'?", s)
    out = set()
    low = s.lower()
    for w, fam_ in (("dhcp", "DHCP"), ("cbcc", "CBCC"), ("fcc", "FCC"), ("bcc", "BCC"), ("hcp", "HCP"), ("cph", "HCP"),
                    ("bct", "BCT")):
        if w in low:
            out.add(fam_)
            low = low.replace(w, " ")
    for t in toks:
        k = t.lower()
        if k.endswith("'") and k != "a3'":
            k = k[:-1]
        if k in STRUCT_WORDS and (not re.fullmatch(r"a\d+'?", k) or t[0] == "A"):
            out.add(STRUCT_WORDS[k])
    if "HEX" in out:          # "hex" names any hexagonal stacking
        out |= {"HCP", "DHCP", "HEXAGONAL"}
    if "CUBIC" in out:
        out |= {"FCC", "BCC", "DIAMOND", "CBCC", "CUB"}
    if "TETRAGONAL" in out:
        out.add("BCT")
    return out


def formula_counts(label, elements):
    """{El: n} when the Greek-free label is a formula of system elements only, else None."""
    s = _greek_free(label)
    s = re.sub(r"\b(?:phase|ht|lt|rt|hp|ss|s\.s\.)\b", "", s, flags=re.I)
    s = _STRIP.sub("", s)
    s = re.sub(r"_\{?(\d+(?:\.\d+)?)\}?", r"\1", s)
    s = re.sub(r"(?<=[A-Za-z0-9])[xδ](?![a-z])|(?<=\d)y(?![a-z])", "", s)     # VCx, TiC1-x: stoichiometry variables
    s = re.sub(r"(?<=\d)1(?=$)", "", s) if False else s
    if not s or not re.fullmatch(r"(?:[A-Z][a-z]?\d*(?:\.\d+)?)+", s):
        return None
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)", s):
        if el == "M":         # generic metal (M2C, MC1-x): any metal of the system
            el = "M"
        elif el.upper() not in elements:
            return None
        out[el.upper()] = out.get(el.upper(), 0.0) + float(n or 1)
    return out


def is_liquid_label(label):
    s = str(label or "").strip()
    return bool(re.match(r"(?i)^(?:l|liq|liquid|melt)\s*(?:$|[(\[+#\d'′,]|\s)", s)) and not re.match(r"^L\d{2}", s) \
        and not re.match(r"^La|^Li|^Lu", s)


def label_kind(label, elements):
    """element | solution | structure | interstitial:<X> | compound | greek | other"""
    s = normalize_label(label)
    if not s:
        return "other"
    if is_liquid_label(s):
        return "liquid"
    if struct_words(s) and not formula_counts(re.sub(r"\(.*?\)", "", s), elements):
        return "structure"
    s0 = re.sub(r"^[%s]['’′\d]*\s*[-–]?\s*(?=\()" % GREEK, "", s)       # beta(Fe,Zr): Greek name of a solution
    inner = s0[1:-1] if s0.startswith("(") and s0.endswith(")") else None
    if inner is not None:
        parts = [p for p in re.split(r"[,\s]+", _greek_free(inner).strip(" -–")) if p]
        if parts and all(re.fullmatch(r"[A-Z][a-z]?", p) and p.upper() in elements for p in parts):
            return "solution"
    fc = formula_counts(s, elements)
    if fc:
        if list(fc) == ["M"]:
            return "other"
        if len(fc) == 1:
            return "element"
        nonm = [e for e in fc if e in INTERSTITIAL]
        if len(nonm) == 1 and len(fc) == 2:
            return "interstitial:" + nonm[0]
        return "compound"
    if re.fullmatch(r"[%s]" % GREEK, s):
        return "greek"
    if re.fullmatch(r"[\s%s'’′″\d()_\-]+" % GREEK, s):
        return "greek_variant"
    return "other"


def interstitial_model(entry, x):
    """True when some model of the registry entry has a last sublattice holding X and VA."""
    for m in (entry or {}).get("models") or []:
        subs = m.get("sublattices") or []
        if subs and x in (subs[-1].get("constituents") or []) and "VA" in (subs[-1].get("constituents") or []):
            return True
    return False


def stable_range(entry, model_name):
    """(x_element, lo, hi) over the stable rows of the registry entry (all assessments), or None."""
    rows = [s for s in (entry or {}).get("stability") or [] if s.get("stable") and s.get("x_min") is not None]
    if not rows:
        return None
    xe = rows[0]["x_element"]
    rows = [s for s in rows if s["x_element"] == xe]
    return xe, min(s["x_min"] for s in rows), max(s["x_max"] for s in rows)


GREEK_X_TOL = 0.01


def majority_inside(entry, comps, model, T, rule):
    """True when, at the record's temperature, the model phase is present in the calculated diagram within
    GREEK_X_TOL of every printed composition in at least half of the system's assessments that have a diagram;
    False without a printed composition, a temperature or a diagram."""
    if not comps or not isinstance(T, (int, float)) or rule is None or not hasattr(rule, "in_field"):
        return False
    stems = sorted({s["tdb"][:-4] if s["tdb"].lower().endswith(".tdb") else s["tdb"]
                    for s in (entry or {}).get("stability") or [] if s.get("tdb")})
    for c in comps:
        votes = []
        for st in stems:
            G = rule.grids.get(st)
            if not G:
                continue
            xe = G["x_component"]
            x = c.get(xe)
            if x is None and len(c) == 1:
                (e1, v1), = c.items()
                x = 1.0 - v1
            if not isinstance(x, (int, float)):
                continue
            res = [rule.in_field(st, model, min(1.0, max(0.0, x + d)), T) for d in (0.0, -GREEK_X_TOL, GREEK_X_TOL)]
            if all(v is None for v in res):
                continue
            votes.append(any(res))
        if not votes or 2 * sum(votes) < len(votes):
            return False
    return True


def check(label, model, system, elements, entry, comps, T=None, rule=None):
    """None when the mapping is kept, else the rejection code. comps: printed mole-fraction compositions {El: x} of
    this phase in the record."""
    label = normalize_label(label)
    mp = str(model).upper().split("#")[0]
    kind = label_kind(label, elements)
    if kind == "liquid":
        return None if LIQUID_MODEL.match(mp) else "liquid_label_mapped_to_non_liquid_model"
    fam = family(mp) if ELEMENT_MODEL.match(mp) else None
    sw = struct_words(str(label))
    if fam and sw and fam not in sw:
        return "structure_word_contradicts_model"
    if fam:
        if kind == "element" or kind == "solution":
            m_ = re.search(r"([%s])\s*(['’′]?)" % GREEK, str(label))
            g = (m_.group(1) + ("'" if m_.group(2) else "")) if m_ else None
            els = list(formula_counts(re.sub(r"[()]", "", str(label)), elements) or {})
            tab = POLYMORPH.get(els[0]) if len(els) == 1 else None
            key = g if tab and g in tab else (g[0] if tab and g and g[0] in tab and len(g) == 1 else None)
            if key and fam not in tab[key]:
                return "element_polymorph_contradicts_model"
        elif kind == "structure":
            pass
        elif kind == "greek":
            fams = set().union(*[POLYMORPH[e].get(label, set()) for e in elements if e in POLYMORPH] or [set()])
            if fam not in fams:
                # no Greek -> structure fallback (A2 audit): Hume-Rothery names (beta brass bcc, zeta hcp) are lost with
                # the wrong ones (eps-Cu3Sn, eps-Fe2Ta); a calculated-field vote at the printed composition was tried
                # and rejected because invariant compositions sit on boundaries where the neighbouring phase votes
                return "greek_label_not_an_element_polymorph_of_this_structure"
        elif kind.startswith("interstitial:") and interstitial_model(entry, kind.split(":")[1]):
            pass
        elif kind in ("greek_variant", "compound") or kind.startswith("interstitial:"):
            return "greek_or_compound_label_mapped_to_element_structure"
        elif sw:
            pass
        else:
            return "label_not_a_solution_or_structure_name"
    if entry is not None and comps and len(elements) == 2:
        rng = stable_range(entry, mp)
        if rng:
            xe, lo, hi = rng
            for c in comps:
                x = c.get(xe)
                if x is None and len(c) == 1:
                    (e1, v1), = c.items()
                    x = 1.0 - v1 if e1 != xe else v1
                if isinstance(x, (int, float)) and 0 <= x <= 1 and not (lo - X_TOL <= x <= hi + X_TOL):
                    return "composition_outside_model_phase_range"
    return None


PRINTED_SOURCES = {"printed", "printed_range"}


def apply(records, rule):
    """rule: phase_rule.PhaseRule (registry by system / TDB name). Mutates the records; returns a report dict."""
    counts, items = Counter(), []
    for r in records:
        system = r.get("system") or ""
        elements = [e.upper() for e in r.get("elements") or system.split("-")]
        names = rule.by_tdbname.get(system, {}) if rule is not None else {}
        for p in r.get("phases") or []:
            mp = p.get("model_phase")
            if not mp:
                continue
            lab = p.get("label")
            comps = [pc.get("composition") or {} for pc in r.get("phase_compositions") or []
                     if pc.get("phase") == lab and pc.get("composition_source") in PRINTED_SOURCES
                     and str(pc.get("basis") or "").startswith("mole_fraction")]
            entry = names.get(str(mp).upper().split("#")[0])
            twin = [q.get("label") for q in r.get("phases") or [] if q is not p and q.get("model_phase") == mp
                    and q.get("label") != lab and label_kind(q.get("label"), elements) in ("element", "solution", "structure")]
            if twin and ELEMENT_MODEL.match(str(mp).upper().split("#")[0]) and \
                    label_kind(lab, elements) in ("greek", "greek_variant", "compound", "other"):
                code = "same_model_as_another_phase_of_the_row"
            else:
                code = None
            code = code or check(lab, mp, system, elements, entry, comps, r.get("temperature_K"), rule)
            if code is None:
                if p.get("pearson_symbol") or p.get("prototype"):
                    p["structure_basis"] = "mapped_model_phase"
                counts["kept"] += 1
                continue
            rej = {"model_phase": mp, "reason": code}
            for k in ("pearson_symbol", "prototype", "mapping_confidence"):
                if p.get(k) is not None:
                    rej[k] = p.pop(k)
            p.pop("model_phase")
            p["model_phase_rejected"] = rej
            r.setdefault("corrections", []).append({
                "gate": "PL", "field": "phases.model_phase", "before": mp, "after": None, "rule": code,
                "evidence": {"label": lab, "release": RELEASE}})
            counts["rejected:" + code] += 1
            counts["rejected"] += 1
            if len(items) < 20000:
                items.append({"id": r["id"], "system": system, "label": lab, "model_phase": mp, "reason": code})
    return {"release": RELEASE, "rules": __doc__.split("Rules")[1].split("Structure fields")[0].strip(),
            "counts": dict(sorted(counts.items())), "items": items}


def report_md(rep):
    """QUALITY_FIXES.md section."""
    c = rep.get("counts") or {}
    L = ["", "## Model-phase mapping check (audit round 2, A2 class 1)", "",
         "`tools/phase_label_check.py` (%s). A model phase attached by the corpus mapping is removed when the printed "
         "label cannot name it; the printed label is unchanged, the removed model phase, its Pearson symbol / prototype "
         "and the reason are kept in `phases[].model_phase_rejected`, and the record gets a `corrections` entry "
         "(gate `PL`). Pearson symbols and prototypes that remain come from the mapped model phase, not from the page: "
         "`phases[].structure_basis = mapped_model_phase`." % rep.get("release"), "",
         "Rules: " + " ".join(rep.get("rules", "").split()), "",
         "| result | phase mentions |", "|---|---:|"]
    for k, n in c.items():
        L.append("| %s | %d |" % (k, n))
    by = Counter((i["system"], i["label"], i["model_phase"], i["reason"]) for i in rep.get("items") or [])
    L += ["", "Most frequent rejected mappings:", "", "| system | label | model phase | reason | n |", "|---|---|---|---|---:|"]
    for (s_, l_, m_, r_), n in by.most_common(40):
        L.append("| %s | %s | %s | %s | %d |" % (s_, str(l_).replace("|", "\\|"), m_, r_, n))
    return "\n".join(L) + "\n"
