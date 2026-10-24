#!/usr/bin/env python3
"""Literature-level consistency gates for the ORCHID literature release (2026-10-23).

Runs inside tools/export_literature_data.py after literature_quality.apply() and before the ESPEI fields, on the
exported record dicts. Twelve gates (G1-G12) check each value against the rest of its table, its paper, the physics of
the quantity and the calculated SGTE diagrams of its system. A gate never deletes a record and never alters a printed
value or quote:

* a *correction* is made only when it is deterministic and evidenced (a printed header / caption / heading, or
  internal consistency of a whole table confirmed on the rendered PDF page). The stored value is kept
  (`value_stored`, `temperature_K_stored`, ...) and the record gets a `corrections` entry
  {gate, field, before, after, rule, evidence[, visual_check]}.
* otherwise the record is *demoted* or *flagged*: `gate_failures` [{gate, code, effect, detail[, visual_check]}],
  where effect `not_accepted` sets flags.extraction_accepted false and `not_espei_ready` sets flags.espei_ready
  (and espei_certified) false; `flag_only` changes nothing else.

Visual checks (rule of 2026-10-23: "if there are errors you need to check the source visually"): every table-level
correction, every system reassignment and the non-mechanical demotions were checked on the rendered PDF page; the
verdicts live in tools/literature_visual_checks.json ({key -> verdict, page, image sha256, what was seen}). No image
or page text is stored there. A correction that needs a visual verdict is applied only when the verdict confirms it;
without one the record is flagged instead.

Standard library only; deterministic.
"""
from __future__ import annotations

import json
import os
import re
import unicodedata
from collections import Counter, defaultdict

import build_dft_parity as parity  # formula parser shared with the DFT benchmark
import phase_rule  # liquid labels, SGTE liquid bounds, element melting points
from espei_constraints import ELEMENT_TM_TB

HERE = os.path.dirname(os.path.abspath(__file__))
VISUAL_FILE = os.path.join(HERE, "literature_visual_checks.json")
GATE_RELEASE = "literature-gates 2026-10-24 (cleanup)"

TEMPERATURE_PROPS = {"INVARIANT_TEMPERATURE", "TRANSITION_TEMPERATURE"}
EQUILIBRIUM_PROPS = TEMPERATURE_PROPS | {"PHASE_COMPOSITION", "PHASE_ASSEMBLAGE"}
NUM = re.compile(r"^\s*[-+]?\d+(?:\.\d+)?\s*$")
C_TO_K = 273.15

GATES = {
    "G1": "temperatures", "G2": "energies, entropies, heat capacities, activities", "G3": "signs",
    "G4": "system assignment", "G5": "property semantics", "G6": "invariant phase lists",
    "G7": "duplicates", "G8": "evidence class", "G9": "compositions", "G10": "reader tier / agreement",
    "G11": "CALPHAD plausibility fields", "G12": "pressure basis", "page": "page checks (visual)",
}


# ----------------------------------------------------------------------------------------------- helpers

def load_visual(path=VISUAL_FILE):
    if not os.path.exists(path):
        return {}
    d = json.load(open(path, encoding="utf-8"))
    return {c["key"]: c for c in d.get("checks", [])}


def vc_of(check):
    """The visual-check provenance stored on a record (no image, no page text)."""
    if not check:
        return None
    return {k: check[k] for k in ("page", "pdf_sha256_prefix", "image_sha256", "verdict", "seen", "checked")
            if check.get(k) is not None}


def fail(r, gate, code, effect, detail=None, check=None):
    """Record a gate failure. effect: not_accepted | not_espei_ready | not_certified | flag_only."""
    gf = r.setdefault("gate_failures", [])
    if any(x["code"] == code for x in gf):
        return False
    e = {"gate": gate, "code": code, "effect": effect}
    if detail:
        e["detail"] = detail
    if check:
        e["visual_check"] = vc_of(check)
    gf.append(e)
    fl = r["flags"]
    if effect == "not_accepted":
        fl["extraction_accepted"] = False
    if effect == "not_certified":
        fl["espei_certified"] = False
    if effect in ("not_accepted", "not_espei_ready"):
        if fl.get("espei_ready"):
            fl["espei_ready_before_gates"] = True
        fl["espei_ready"] = False
        fl["espei_certified"] = False
    return True


def correct(r, gate, field, before, after, rule, evidence, check=None):
    c = {"gate": gate, "field": field, "before": before, "after": after, "rule": rule, "evidence": evidence}
    if check:
        c["visual_check"] = vc_of(check)
    r.setdefault("corrections", []).append(c)


def labels(r):
    return [p.get("label") for p in r.get("phases") or []]


def is_liq(p):
    return phase_rule.is_liquid_label(p.get("label"), p.get("model_phase"))


def asset_of(r):
    loc = r.get("locator") or {}
    return loc.get("asset") or loc.get("ref")


def printed_number(s):
    if s is None:
        return None
    s = str(s)
    if NUM.match(s):
        return float(s)
    return None


def tm(el):
    v = ELEMENT_TM_TB.get(el)
    return v[0] if v else None


def norm_label(l):
    l = unicodedata.normalize("NFKC", l or "")
    l = l.replace("′", "'").replace("″", '"').replace("’", "'").replace("−", "-").replace("–", "-")
    return re.sub(r"\s+", "", l).lower()


GREEK_NAMES = {"alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε", "eps": "ε", "zeta": "ζ", "eta": "η",
               "theta": "θ", "kappa": "κ", "lambda": "λ", "mu": "μ", "sigma": "σ", "tau": "τ", "phi": "φ", "chi": "χ",
               "psi": "ψ", "omega": "ω", "xi": "ξ", "nu": "ν", "pi": "π", "rho": "ρ", "iota": "ι", "upsilon": "υ"}
_GREEK_WORD = re.compile(r"\\?\b(" + "|".join(sorted(GREEK_NAMES, key=len, reverse=True)) + r")(?!(?-i:[a-z]))", re.I)
_KEY_GLYPH = str.maketrans({"−": "-", "–": "-", "—": "-", "‐": "-", "′": "'", "’": "'", "″": '"', "ϕ": "φ", "ϵ": "ε"})
_KEY_WORDS = re.compile(r"(?i)\b(?:phase|sol\.?|solution|solid solution|ss)(?![a-z])")


def phase_key(label):
    """Spelling-independent key of a printed phase label: Unicode compatibility forms (subscripts), LaTeX and spelled
    Greek letters (\\sigma, gamma -> σ, γ), dash / prime glyphs (″ = two primes), case, white space, '_', '-', '.',
    outer parentheses and the words 'phase' / 'sol.' / 'solution' are ignored. Used to find a phase printed (or read)
    twice under two spellings; never to rename a stored label."""
    s = unicodedata.normalize("NFKC", label or "")
    s = _GREEK_WORD.sub(lambda m: GREEK_NAMES[m.group(1).lower()], s)
    s = s.replace("\\", "").replace("{", "").replace("}", "").replace("$", "").translate(_KEY_GLYPH).replace("''", '"')
    s = _KEY_WORDS.sub("", s)
    s = re.sub(r"[\s_\-·.]+", "", s).lower()
    while len(s) > 2 and s.startswith("(") and s.endswith(")") and s.count("(") == 1:
        s = s[1:-1]
    return s


_G6_LIQ = re.compile(r"^\(?\s*(?:L|Liq|Liquid|liq|liquid|LIQUID)\s*[#_]?\s*(?:\d|'+|\"|′+|″|I{1,3})?\s*#?\s*\)?$")


def is_liq_g6(p):
    """Liquid test of the invariant gate: phase_rule.is_liquid_label plus numbered liquids (L1, L2, L#1, L', LI)."""
    return is_liq(p) or bool(_G6_LIQ.match((p.get("label") or "").strip()))


def formula_elements(label):
    """Element set of a chemical formula in a phase label (any element), or None."""
    c, _ = parity.parse_phase(label, None)
    return {e.upper() for e in c} if c else None


_TOKEN = re.compile(r"^(?:[A-Z][a-z]?\d*(?:\.\d+)?)+$")
_PART = re.compile(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)")
STRUKTUR = re.compile(r"^(?:[A-DL]\d{1,2}[a-z]?|D\d_?\d{1,2}|L\d_?\d|[A-C]\d+'?)$")


def foreign_formula_elements(label, system_elements):
    """Elements outside the system named by a chemical formula token of the label: the token must be a complete
    formula (>= 2 element parts), contain at least one system element, and not be a structure designation
    (Strukturbericht symbol, '<prototype>-type', a polymorph letter after the formula) or a capitalisation slip of
    a formula of the system's own elements (FeHF2 for FeHf2)."""
    sysel = {e.upper() for e in system_elements}
    sym = {e.capitalize() for e in sysel}
    out = set()
    lab = unicodedata.normalize("NFKC", label or "")
    for tok in re.split(r"[\s,;:/()\[\]{}+=~'′’\"`*]+", lab):
        if not tok or tok.lower().endswith("type") or "-type" in tok.lower():
            continue
        tok = re.sub(r"^[^A-Z]+", "", tok.split("-")[0] if re.search(r"-[A-Za-z]", tok) and not _TOKEN.match(tok) else tok)
        if not tok or STRUKTUR.match(tok) or not _TOKEN.match(tok) or re.fullmatch(r"[A-Z]{3,}", tok):
            continue
        parts = _PART.findall(tok)
        if len(parts) < 2:
            continue
        els = {p[0].upper() for p in parts}
        if not els & sysel or not all(p[0].capitalize() in ELEMENT_SYMBOLS for p in parts):
            continue
        # capitalisation slip: the same letters spell a formula of system elements only
        if re.fullmatch("(?:%s)+" % "|".join(sorted((re.escape(x) + r"\d*(?:\.\d+)?" for x in sym), key=len, reverse=True)),
                        tok, re.I):
            continue
        out |= els - sysel
    return out


def text_of(r):
    loc = r.get("locator") or {}
    return " ".join(str(x) for x in (r.get("property_printed"), r.get("quote"), loc.get("caption")) if x)


class Report:
    def __init__(self):
        self.counts = Counter()
        self.items = defaultdict(list)     # gate code -> [entries]

    def add(self, gate, code, r, action, **kw):
        self.counts[(gate, code, action)] += 1
        e = {"id": r["id"], "system": r["system"], "property": r.get("property"), "action": action}
        e.update({k: v for k, v in kw.items() if v is not None})
        self.items[(gate, code)].append(e)


# ----------------------------------------------------------------------------------------------- G7 duplicates

_GLYPH = str.maketrans({"−": "-", "–": "-", "—": "-", "‐": "-", "′": "'", "’": "'", "″": '"'})


def _label_tokens(label):
    s = unicodedata.normalize("NFKC", label or "")
    s = _GREEK_WORD.sub(lambda m: GREEK_NAMES[m.group(1).lower()], s)
    s = s.replace("\\", "").replace("{", "").replace("}", "").replace("$", "").translate(_KEY_GLYPH).lower()
    toks = {t.strip("-_.'") for t in re.split(r"[\s,;()\[\]/]+", s)}
    return {t for t in toks if t and t not in ("phase", "sol", "sol.", "solution", "ss", "to", "and")}


def _cosmetic_phases(r):
    """[(key, label)] of the named phases of a record: the liquid under any spelling is one phase."""
    out = {}
    for p in r.get("phases") or []:
        lab = (p.get("label") or "").translate(_GLYPH).strip()
        if lab.upper() in phase_rule.NOT_PRINTED:
            continue
        k = "liquid" if is_liq_g6(p) and not re.search(r"\d", lab) else phase_key(lab)
        out.setdefault(k, lab)
    return sorted(out.items())


def _formula_fractions(label):
    c, _ = parity.parse_phase(label, None)
    if not c:
        return None
    tot = float(sum(c.values()))
    return {e.upper(): round(float(n) / tot, 2) for e, n in c.items()} if tot else None


def _same_phase(a, b):
    """Two spellings / readings of one phase label: the same key, one label's words contained in the other's
    ('Fe3O4' / 'Magnetite (Fe3O4)', 'C14 Laves' / 'C14 Laves TiMn2'), or the same formula composition
    ('UC1.90' / 'C0.655U0.345')."""
    (ka, la), (kb, lb) = a, b
    if ka == kb:
        return True
    ta, tb = _label_tokens(la), _label_tokens(lb)
    if ta and tb and (ta <= tb or tb <= ta):
        return True
    fa, fb = _formula_fractions(la), _formula_fractions(lb)
    return bool(fa and fb and len(fa) > 1 and set(fa) == set(fb) and all(abs(fa[e] - fb[e]) <= 0.01 for e in fa))


def _phases_compatible(pa, pb):
    if not pa or not pb:
        return True
    return all(any(_same_phase(x, y) for y in pb) for x in pa) and all(any(_same_phase(y, x) for x in pa) for y in pb)


def _specificity(r):
    named = sum(1 for p in r.get("phases") or [] if (p.get("label") or "").strip().upper() not in phase_rule.NOT_PRINTED)
    mapped = sum(1 for p in r.get("phases") or [] if p.get("model_phase"))
    # last tie-break: a label that names the formula ('NiAl (B2)') over one that does not ('B2')
    formula = any(_formula_fractions(p.get("label")) for p in r.get("phases") or [])
    toks = sum(len(_label_tokens(p.get("label"))) for p in r.get("phases") or [])
    return (bool(r["flags"].get("extraction_accepted")), bool(r["flags"].get("espei_ready")), named, mapped,
            len(r.get("readers") or []),
            len(r.get("phase_compositions") or []), bool(r.get("overall_composition")), bool(r.get("reference_state")),
            formula, toks)


def _compatible(a, b):
    """Two readings of one printed value differ only cosmetically: the same named phases (or none on one side) and the
    same compositions (or none on one side)."""
    if not _phases_compatible(_cosmetic_phases(a), _cosmetic_phases(b)):
        return False
    for f in ("overall_composition", "phase_compositions"):
        if a.get(f) and b.get(f) and json.dumps(a[f], sort_keys=True) != json.dumps(b[f], sort_keys=True):
            return False
    return True


def _compatible_but_phases(a, b):
    for f in ("overall_composition", "phase_compositions"):
        if a.get(f) and b.get(f) and json.dumps(a[f], sort_keys=True) != json.dumps(b[f], sort_keys=True):
            return False
    return bool(_cosmetic_phases(a)) and bool(_cosmetic_phases(b))


def g7_dedupe(records, rep):
    """Records identical in paper / property / value / unit / T / component / table / row (or quote) / evidence class
    whose only differences are cosmetic (an UNKNOWN or missing phase beside a named one, the liquid printed as L or
    liquid, minus / prime glyph variants, a composition or reference state given by one reading only, reader lists) are
    merged into the most specified record; merged ids are kept."""
    groups = defaultdict(list)
    for r in records:
        loc = r.get("locator") or {}
        if r.get("value") is None and r.get("property") not in EQUILIBRIUM_PROPS:
            continue
        key = (r.get("paper"), r["system"], r.get("property"), json.dumps(r.get("value")), r.get("unit"),
               json.dumps(r.get("temperature_K")), asset_of(r), loc.get("row"), r.get("quote"),
               r.get("evidence_class"), r.get("component"))
        groups[key].append(r)
    drop = set()
    for key, g in groups.items():
        if len(g) < 2:
            continue
        kept = []    # [(best, [others])]
        for r in sorted(g, key=lambda x: tuple(-int(v) for v in _specificity(x)) + (x["id"],)):
            fits = [k for k in kept if _compatible(k[0], r) and all(_compatible(o, r) for o in k[1])]
            if len(fits) == 1:
                fits[0][1].append(r)
            else:
                kept.append((r, []))     # a new reading, or ambiguous between two (kept apart)
        if len(kept) > 1:
            # one printed value read with different phase labels: at most one reading is right
            for best, _ in kept:
                if all(_compatible_but_phases(best, k[0]) for k in kept if k[0] is not best) and _cosmetic_phases(best):
                    fail(best, "G7", "phase_label_readings_disagree", "not_espei_ready",
                         "the same printed value is also stored with phase(s) %s" % "; ".join(
                             ", ".join(labels(k[0])) for k in kept if k[0] is not best)[:200])
                    rep.add("G7", "phase_label_readings_disagree", best, "flagged")
        for best, others in kept:
            if not others:
                continue
            ids = set(best.get("merged_ids") or [])
            pids = set(best.get("merged_pids") or [])
            for x in others:
                ids.add(x["id"])
                ids.update(x.get("merged_ids") or [])
                if x.get("pid") and x.get("pid") != best.get("pid"):
                    pids.add(x["pid"])
                pids.update(x.get("merged_pids") or [])
                drop.add(id(x))
            best["merged_ids"] = sorted(ids)
            # the other readings' spellings of a single phase stay with the kept record (`label_variants`)
            bp = best.get("phases") or []
            if len(bp) == 1:
                var = sorted({(x.get("phases") or [{}])[0].get("label") for x in others if len(x.get("phases") or []) == 1}
                             - {bp[0].get("label"), None, "UNKNOWN"})
                if var:
                    bp[0]["label_variants"] = sorted(set(bp[0].get("label_variants") or []) | set(var))
            if pids - {best.get("pid")}:
                best["merged_pids"] = sorted(pids - {best.get("pid")})
            rep.add("G7", "cosmetic_duplicate_merged", best, "merged", merged=len(others),
                    kept_phases=labels(best), dropped_phases=[labels(x) for x in others][:3])
    out = [r for r in records if id(r) not in drop]
    rep.counts[("G7", "records_removed", "merged")] = len(records) - len(out)
    return out


# ----------------------------------------------------------------------------------------------- G12 / G11

def g12_pressure(r, rep):
    b = r.get("pressure_basis")
    if b in ("stated_default_1_atm", "unstated_assumed") and r.get("pressure_Pa") == 101325:
        r["pressure_basis_stored"] = b
        r["pressure_basis"] = "assumed_default_1_atm_not_printed"
        rep.counts[("G12", "pressure_basis_relabelled", "relabelled")] += 1


def g11_plausibility(r, plaus_store, rep):
    """dx (mole-fraction differences per phase) and the calculation behind the CALPHAD check, from the record
    store (master_store.py), so the tolerance can be judged per value; fails with |dT| < 5 K are explained."""
    pl = r["flags"].get("plausibility")
    if not pl:
        return
    st = plaus_store.get(r.get("pid")) if plaus_store else None
    if st:
        if st.get("dx") is not None and "dx" not in pl:
            pl["dx"] = {k: round(v, 4) for k, v in sorted(st["dx"].items()) if isinstance(v, (int, float))}
        if st.get("calculation") and "calculation" not in pl:
            pl["calculation"] = st["calculation"]
    pl.setdefault("tolerance", {"dT_K": 25.0, "dx": 0.05})
    if pl.get("result") == "fail" and pl.get("dT_K") is not None and abs(pl["dT_K"]) < 5:
        dx = pl.get("dx") or {}
        big = sorted(k for k, v in dx.items() if abs(v) > pl["tolerance"]["dx"])
        if big:
            why = "composition"
            # a pure-element phase label compared on the wrong axis: |dx| equals the other element's fraction
            if all(re.fullmatch(r"\(?[A-Z][a-z]?\)?|[αβγδε]?[A-Z][a-z]?", k) for k in big):
                why = "composition_pure_element_label"
            pl["fail_reason"] = why
            rep.counts[("G11", "fail_small_dT_" + why, "explained")] += 1
        else:
            pl["fail_reason"] = "unexplained"
            rep.counts[("G11", "fail_small_dT_unexplained", "explained")] += 1


# ----------------------------------------------------------------------------------------------- G4 system

FOREIGN_SYSTEM_CAPTION = re.compile(r"\b([A-Z][a-z]?(?:\s*[–\-—]\s*[A-Z][a-z]?){1,4})\b(?:\s*[–\-—]\s*\w+)?\s+(?:ternary\s+|quaternary\s+|binary\s+)?(?:system|alloys?|melts?|phase diagram)")
CEMENT = re.compile(r"\bC\s*[-–]\s*S\b.*(?:system|cement)|cement", re.I)
ELEMENT_SYMBOLS = {e.capitalize() for e in ELEMENT_TM_TB} | {"O", "N", "C", "B", "H", "S", "F", "P", "Cl", "Br", "I", "Se",
                                                             "Te", "As", "Si", "Ge", "Ar", "He", "Ne", "Kr", "Xe"}


def caption_systems(text):
    """Element sets named as '<A>-<B>[-<C>] system' in a caption / quote."""
    out = []
    text = re.sub(r"\bAI(?=\s*[–\-—])|(?<=[–\-—])\s*AI\b", "Al", text or "")   # OCR: 'AI' for 'Al'
    for m in FOREIGN_SYSTEM_CAPTION.finditer(text):
        els = [e.strip() for e in re.split(r"[–\-—]", m.group(1))]
        if all(e in ELEMENT_SYMBOLS for e in els) and len(set(els)) >= 2:
            out.append(frozenset(e.upper() for e in els))
    return out


def g4_system(records, visual, rep):
    """Values whose phase formulas name an element outside the system: reassigned to the system named by all the
    phase formulas and the caption when that is unambiguous and confirmed on the page; otherwise not accepted."""
    by_asset = defaultdict(list)
    for r in records:
        by_asset[asset_of(r)].append(r)
    for r in records:
        chk = visual.get("record:" + r["id"])
        seen_sys = sorted(set((chk or {}).get("system", "").upper().split("-")) - {""})
        if chk and seen_sys and set(seen_sys) != set(r["elements"]) and chk.get("verdict") in ("misread", "higher_order", "confirmed", "as_stored"):
            new = "-".join(seen_sys)
            correct(r, "G4", "system", r["system"], new, "system_read_from_page", {"seen": chk.get("seen")}, chk)
            r["system_stored"], r["system"], r["elements"] = r["system"], new, seen_sys
            r["flags"]["system_reassigned"] = True
            r["_system_from_page"] = True
            rep.add("G4", "system_reassigned_from_page", r, "corrected", to=new)
    for r in records:
        if r.pop("_system_from_page", None):
            continue
        sysel = set(r["elements"])
        foreign, fels = set(), set()
        for lab in labels(r):
            fe = foreign_formula_elements(lab, sysel)
            if fe:
                foreign |= fe
                fels |= (formula_elements(lab) or set()) | fe
        cement = r["system"] == "C-S" and bool(CEMENT.search(text_of(r)) or any(
            re.search(r"Ca\d*Si\w*O|CaO|SiO2|Ca2SiO4|Ca3SiO5", l or "") for l in labels(r)))
        if not foreign and not cement:
            continue
        loc = r.get("locator") or {}
        caps = caption_systems((loc.get("caption") or "") + " " + (r.get("quote") or ""))
        target = None
        if cement:
            target = frozenset({"CA", "O", "SI"}) | frozenset(fels & {"AL", "FE", "MG"})
        else:
            cand = frozenset(sysel | foreign)
            if any(c == cand for c in caps):
                target = cand
            elif not caps and len(foreign) <= 2 and fels and fels >= sysel:
                target = frozenset(fels | sysel)
        key = "system:" + r["id"]
        akey = "system:" + str(asset_of(r))
        chk = visual.get(key) or visual.get(akey)
        seen_sys = set((chk or {}).get("system", "").upper().split("-")) - {""}
        if chk and seen_sys and seen_sys == sysel and chk.get("verdict") in ("confirmed", "misread", "as_stored"):
            # the page confirms the stored system: the foreign element is a misprint / misread in a phase label
            r["flags"]["phase_label_names_foreign_element"] = True
            fail(r, "G4", "phase_label_names_foreign_element", "flag_only",
                 "label names %s; the page confirms the system %s (label misprint or misread)" % (", ".join(sorted(foreign)) or "-", r["system"]), chk)
            rep.add("G4", "system_confirmed_label_typo", r, "flagged")
            continue
        if chk and seen_sys and chk.get("verdict") in ("confirmed", "misread", "higher_order"):
            target = frozenset(seen_sys)
        if target and chk and chk.get("verdict") in ("confirmed", "misread", "higher_order") and seen_sys == set(target):
            new = "-".join(sorted(target))
            correct(r, "G4", "system", r["system"], new, "system_from_phase_formulas_and_caption",
                    {"foreign_elements": sorted(foreign), "caption_systems": sorted("-".join(sorted(c)) for c in caps),
                     "cement_notation": cement or None}, chk)
            r["system_stored"] = r["system"]
            r["system"] = new
            r["elements"] = sorted(target)
            r["flags"]["system_reassigned"] = True
            rep.add("G4", "system_reassigned", r, "corrected", to=new)
        else:
            detail = "phase formula elements outside the system: %s" % ", ".join(sorted(foreign)) if foreign else \
                "cement notation C-S (CaO-SiO2) stored as the C-S element system"
            if target:
                detail += "; candidate system %s (not visually confirmed)" % "-".join(sorted(target))
            r["flags"]["system_mismatch"] = True
            fail(r, "G4", "system_mismatch", "not_accepted", detail, chk)
            rep.add("G4", "system_mismatch", r, "demoted", foreign=sorted(foreign), candidate="-".join(sorted(target)) if target else None)


NOBLE_GASES = {"HE", "NE", "AR", "KR", "XE", "RN"}
PRESSURE_UNIT = re.compile(r"\d\s*(?:[kMG]?Pa|hPa)\b")
PROTACTINIUM = re.compile(r"\bPa(?:\d|[A-Z]|\s*[–\-—]\s*[A-Z])|protactinium", re.I)
GAS_ATMOSPHERE = re.compile(r"\b(?:Ar|He|argon|helium|N2|nitrogen)\b[^.;]{0,30}\b(?:flow|atmosphere|gas|cover|purge|stream)|"
                            r"\b(?:under|in|flowing|purified|high[- ]purity)\s+(?:an?\s+)?(?:Ar|He|argon|helium)\b", re.I)


def g4_unit_gas_elements(records, rep):
    """James Saal review 2026-10-24: a system element taken from a unit (Pa in kPa / MPa), a cover or flow gas (Ar,
    He) or a garbled phase prefix (H-CuTe read as HeCuTe) is not a component. A noble gas is never a component of the
    condensed systems ORCHID covers; Pa counts only when the text names protactinium (Pa followed by a formula digit or
    element, or the word). Records whose system the page check already reassigned (G4 above) pass; the rest are not
    accepted (`system_element_from_unit_or_gas`) until a page check reassigns them."""
    for r in records:
        if r["flags"].get("system_reassigned"):
            continue
        els = set(r.get("elements") or [])
        txt = text_of(r) + " " + " ".join(l or "" for l in labels(r))
        bad = []
        for e in sorted(els & NOBLE_GASES):
            why = "cover / flow gas" if GAS_ATMOSPHERE.search(txt) else "noble gas (not a condensed-system component)"
            bad.append("%s: %s" % (e.capitalize(), why))
        if "PA" in els and PRESSURE_UNIT.search(txt) and not PROTACTINIUM.search(PRESSURE_UNIT.sub(" ", txt)):
            bad.append("Pa: pressure unit (Pa/kPa/MPa), not protactinium")
        if bad:
            r["flags"]["system_mismatch"] = True
            fail(r, "G4", "system_element_from_unit_or_gas", "not_accepted", "; ".join(bad))
            rep.add("G4", "system_element_from_unit_or_gas", r, "demoted", detail="; ".join(bad))


# ----------------------------------------------------------------------------------------------- G1 temperatures

def _anchors(rule, elements):
    a = [tm(e) for e in elements if tm(e)]
    els = sorted(elements)
    for i in range(len(els)):
        for j in range(i + 1, len(els)):
            a += list(rule.inv_T.get(els[i] + "-" + els[j], []))
    return a


def _median(xs):
    xs = sorted(xs)
    n = len(xs)
    return None if not n else (xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2]))


def _pure_element_row(r):
    """The element whose pure state the row describes (all printed phase compositions at x=0 or 1 of one element)."""
    pcs = [pc.get("composition") or {} for pc in r.get("phase_compositions") or []]
    if not pcs:
        return None
    pure = set()
    for c in pcs:
        if len(c) == 1:
            (el, v), = c.items()
            if v >= 0.999:
                pure.add(el)
            elif v <= 0.001:
                others = [e for e in r["elements"] if e != el]
                if len(others) == 1:
                    pure.add(others[0])
                else:
                    return None
            else:
                return None
        else:
            big = [e for e, v in c.items() if v >= 0.999]
            if len(big) != 1:
                return None
            pure.add(big[0])
    return pure.pop() if len(pure) == 1 else None


def table_temperature_evidence(rule, rows, tunits=None):
    """Evidence about the printed temperature unit of one table: what the stored temperatures assume and what the
    table's own numbers say (pure-element melting points, SGTE invariants of the system, reader unit votes)."""
    els = sorted({e for r in rows for e in r["elements"]})
    anc = _anchors(rule, els)
    modes = Counter()
    for r in rows:
        v = printed_number(r.get("temperature_printed"))
        T = r.get("temperature_K")
        if v is None or T is None:
            continue
        modes["K" if abs(T - v) < 0.01 else "C" if abs(T - v - C_TO_K) < 0.6 else "other"] += 1
    ev = {"stored_as": dict(modes)}
    pure_c, pure_k = [], []
    for r in rows:
        v = printed_number(r.get("temperature_printed"))
        el = _pure_element_row(r)
        if v is None or not el or not tm(el) or tm(el) < 350:
            continue
        if abs(v - (tm(el) - C_TO_K)) <= 2.0:
            pure_c.append("%s %g" % (el, v))
        elif abs(v - tm(el)) <= 2.0:
            pure_k.append("%s %g" % (el, v))
    if pure_c:
        ev["pure_element_rows_match_celsius"] = pure_c
    if pure_k:
        ev["pure_element_rows_match_kelvin"] = pure_k
    inv = [printed_number(r.get("temperature_printed")) for r in rows if r.get("property") in TEMPERATURE_PROPS]
    inv = [v for v in inv if v is not None]
    if anc and len(inv) >= 2:
        dK = _median([min(abs(a - v) for a in anc) for v in inv])
        dC = _median([min(abs(a - (v + C_TO_K)) for a in anc) for v in inv])
        ev["invariant_anchor_median_dK"] = {"as_K": round(dK, 1), "as_C": round(dC, 1), "n": len(inv)}
    if tunits:
        ev["reader_unit_votes"] = dict(tunits)
    # decision
    vote = None
    a = ev.get("invariant_anchor_median_dK")
    if pure_c and not pure_k:
        vote = "C"
    elif pure_k and not pure_c:
        vote = "K"
    if a:
        if a["as_K"] > 40 and a["as_C"] <= 12:
            vote = vote or "C"
            if vote == "K":
                vote = "conflict"
        elif a["as_C"] > 40 and a["as_K"] <= 12:
            vote = vote or "K"
            if vote == "C":
                vote = "conflict"
    if tunits and vote is None:
        tot = sum(tunits.values())
        if tunits.get("C", 0) >= 0.75 * tot and modes.get("K"):
            vote = "C?"
    ev["inferred_unit"] = vote
    return ev


def g1_temperatures(records, rule, visual, store_units, rep):
    # ---- (a) table-level unit inference
    tables = defaultdict(list)
    for r in records:
        if r["record_kind"] != "running_text" and r.get("temperature_K") is not None:
            tables[asset_of(r)].append(r)
    for asset, rows in sorted(tables.items(), key=lambda x: str(x[0])):
        tu = Counter()
        for r in rows:
            for k, n in ((store_units.get(r.get("pid")) or {}).get("votes") or {}).items():
                tu[k] += n
        ev = table_temperature_evidence(rule, rows, tu)
        modes = ev["stored_as"]
        inferred = ev["inferred_unit"]
        mixed = modes.get("K", 0) > 0 and modes.get("C", 0) > 0
        suspect = mixed or (inferred in ("C", "C?") and modes.get("K")) or (inferred == "K" and modes.get("C")) \
            or inferred == "conflict"
        chk = visual.get("table_T_unit:" + str(asset))
        if not suspect and not chk:
            continue
        unit = (chk or {}).get("unit") if (chk or {}).get("verdict") == "confirmed" else None
        if unit in ("C", "K"):
            n_changed = 0
            for r in rows:
                v = printed_number(r.get("temperature_printed"))
                if v is None:
                    # a footnote mark after the number ('287*', '271†') does not change the cell's unit
                    v = printed_number(re.sub(r"\s*[*†‡§]+\s*$", "", str(r.get("temperature_printed") or "")) or None)
                if v is None:
                    continue
                newT = round(v + C_TO_K, 4) if unit == "C" else v
                if abs((r.get("temperature_K") or 0) - newT) < 0.01:
                    continue
                before = r.get("temperature_K")
                r["temperature_K_stored"] = before
                r["temperature_K"] = newT
                if r.get("property") in TEMPERATURE_PROPS and r.get("value") == before:
                    r["value_stored"] = r.get("value")
                    r["value"] = newT
                correct(r, "G1", "temperature_K", before, newT, "table_temperature_unit_%s" % unit,
                        {k: v for k, v in ev.items() if v}, chk)
                pl = r["flags"].pop("plausibility", None)
                if pl:
                    r["flags"]["plausibility_before_correction"] = pl
                    r["flags"]["calphad_check"] = "not_available"
                    if r["flags"].get("espei_certified"):
                        r["flags"]["espei_certified"] = False
                n_changed += 1
                rep.add("G1", "table_temperature_unit_corrected", r, "corrected", asset=asset, unit=unit,
                        before=before, after=newT)
            if n_changed == 0 and suspect:
                rep.counts[("G1", "table_temperature_unit_confirmed_as_stored", "checked")] += 1
        elif chk and chk.get("verdict") == "as_stored":
            rep.counts[("G1", "table_temperature_unit_confirmed_as_stored", "checked")] += 1
        else:
            for r in rows:
                fail(r, "G1", "temperature_unit_table_suspect", "not_espei_ready",
                     "table temperatures stored as %s; table evidence: %s" % (
                         "/".join(sorted(modes)), json.dumps({k: v for k, v in ev.items() if v and k != "stored_as"},
                                                             ensure_ascii=False)), chk)
                rep.add("G1", "temperature_unit_table_suspect", r, "flagged", asset=asset)
    # ---- (a2) sibling rows: the same reaction / phase set of the same table, evidence class and paper stored once at
    # T and once at T + 273.15, the first without a printed temperature of its own: a °C number stored as K
    for asset, rows in sorted(tables.items(), key=lambda x: str(x[0])):
        sig = defaultdict(list)
        for r in rows:
            k = (tuple(sorted({"liquid" if is_liq(p) else norm_label(p.get("label")) for p in r.get("phases") or []})),
                 r.get("evidence_class"), r.get("property"))
            if k[0]:
                sig[k].append(r)
        for k, grp in sig.items():
            if len(grp) < 2:
                continue
            for r in grp:
                T = r.get("temperature_K")
                if T is None or any(c.get("field") == "temperature_K" for c in r.get("corrections") or []):
                    continue
                pv = printed_number(r.get("temperature_printed"))
                if pv is not None and abs(pv - T) > 0.01:
                    continue     # its own printed number was converted: not this case
                sib = [x for x in grp if x is not r and x.get("temperature_K") is not None
                       and abs(x["temperature_K"] - T - C_TO_K) < 0.6 and re.search(r"°\s*C|\bC\b", str(x.get("temperature_printed") or ""))]
                if not sib:
                    continue
                chk = visual.get("table_T_unit:" + str(asset))
                if chk and chk.get("verdict") == "confirmed" and chk.get("unit") == "C":
                    newT = round(T + C_TO_K, 4)
                    correct(r, "G1", "temperature_K", T, newT, "sibling_row_celsius_stored_as_K",
                            {"sibling": sib[0]["id"], "sibling_printed": sib[0].get("temperature_printed")}, chk)
                    r["temperature_K_stored"], r["temperature_K"] = T, newT
                    if r.get("property") in TEMPERATURE_PROPS and r.get("value") == T:
                        r["value_stored"], r["value"] = T, newT
                    pl = r["flags"].pop("plausibility", None)
                    if pl:
                        r["flags"]["plausibility_before_correction"] = pl
                        r["flags"]["calphad_check"] = "not_available"
                    rep.add("G1", "sibling_row_celsius_corrected", r, "corrected", asset=asset, before=T, after=newT)
                else:
                    fail(r, "G1", "temperature_unit_sibling_row", "not_accepted",
                         "a sibling row (%s, printed '%s') holds this temperature + 273.15 K: the number is °C" % (
                             sib[0]["id"], sib[0].get("temperature_printed")), chk)
                    rep.add("G1", "temperature_unit_sibling_row", r, "demoted", asset=asset)
    # ---- (b) non-temperatures, (c) impossible liquids, (d) ranges stored as points
    for r in records:
        T = r.get("temperature_K")
        if T is None:
            continue
        tp = r.get("temperature_printed")
        tps = str(tp).strip() if tp is not None else ""
        prop = r.get("property")
        key = "record:" + r["id"]
        chk = visual.get(key)
        # (b) zero / dash / composition read as a temperature
        if prop in EQUILIBRIUM_PROPS or r["record_kind"] != "running_text" and r.get("evidence_class") != "dft":
            why = None
            if tps in ("–", "-", "—", "−", "…", "...") and T is not None:
                why = "temperature printed as '%s' (no value)" % tps
            elif T <= phase_rule.MIN_T_K and prop in EQUILIBRIUM_PROPS and not re.search(r"\d\s*m?K\b", tps):
                why = "temperature %g K (zero or blank read as a temperature)" % T
            else:
                comps = [v for pc in r.get("phase_compositions") or [] if pc.get("composition_source") == "printed"
                         for v in (pc.get("composition") or {}).values()] + list((r.get("overall_composition") or {}).values())
                if prop in EQUILIBRIUM_PROPS and tp is None:
                    for v in comps:
                        if isinstance(v, (int, float)) and 0 < v < 1 and (abs(T - 100 * v) < 0.051 or abs(T - C_TO_K - 100 * v) < 0.051):
                            why = "temperature equals a printed composition (%g at.%%)" % round(100 * v, 3)
                            break
            if why:
                if chk and chk.get("verdict") == "corrected" and chk.get("temperature_printed"):
                    newT = float(chk["temperature_K"])
                    correct(r, "G1", "temperature_K", T, newT, "temperature_read_from_page", {"printed": chk["temperature_printed"]}, chk)
                    r["temperature_K_stored"] = T
                    r["temperature_K"] = newT
                    r["temperature_printed_page"] = chk["temperature_printed"]
                    if prop in TEMPERATURE_PROPS:
                        r["value_stored"], r["value"] = r.get("value"), newT
                    r["flags"].pop("plausibility", None)
                    rep.add("G1", "non_temperature_corrected_from_page", r, "corrected", before=T, after=newT)
                else:
                    fail(r, "G1", "temperature_not_a_temperature", "not_accepted", why, chk)
                    rep.add("G1", "temperature_not_a_temperature", r, "demoted", T=T, printed=tp)
                continue
        # (d) a printed range stored as a point
        m = re.match(r"^\s*([-−–]?\d+(?:\.\d+)?)\s*(?:°?\s*[CK])?\s*(?:to|–|—|-|…|\.\.\.|~)\s*([-−–]?\d+(?:\.\d+)?)\s*(°?\s*[CK])?\s*$", tps)
        if m and not r.get("temperature_range_K") and not r.get("temperature_max_K"):
            lo, hi = (float(x.replace("−", "-").replace("–", "-")) for x in m.groups()[:2])
            if hi > lo:
                fail(r, "G1", "temperature_range_stored_as_point", "not_accepted" if prop in TEMPERATURE_PROPS else "not_espei_ready",
                     "printed temperature '%s' is a range; the stored point %g K is one end" % (tps, T))
                rep.add("G1", "temperature_range_stored_as_point", r, "demoted" if prop in TEMPERATURE_PROPS else "flagged", printed=tps)
                continue
        # (c) a liquid below any physically possible temperature
        liq = [p for p in r.get("phases") or [] if is_liq(p)]
        if liq:
            bound, basis = liquid_bound_any(rule, r["elements"])
            if bound is not None and T < bound:
                tms = [tm(e) for e in r["elements"] if tm(e)]
                lowest = (rule.liq_min.get(r["system"]) or (None,))[0]
                gross = prop in EQUILIBRIUM_PROPS and (
                    (tms and len(tms) == len(r["elements"]) and T < phase_rule.EQUILIBRIUM_T_FRACTION * min(tms))
                    or (lowest is not None and not set(r["elements"]) & VOLATILE and T < 0.5 * lowest))
                if chk and chk.get("verdict") == "corrected" and chk.get("temperature_K") is not None:
                    newT = float(chk["temperature_K"])
                    correct(r, "G1", "temperature_K", T, newT, "temperature_read_from_page", {"printed": chk.get("temperature_printed")}, chk)
                    r["temperature_K_stored"], r["temperature_K"] = T, newT
                    if prop in TEMPERATURE_PROPS:
                        r["value_stored"], r["value"] = r.get("value"), newT
                    r["flags"].pop("plausibility", None)
                    rep.add("G1", "liquid_temperature_corrected_from_page", r, "corrected", before=T, after=newT)
                    continue
                fail(r, "G1", "liquid_below_physical_minimum", "not_accepted" if gross else "not_espei_ready",
                     "liquid-bearing value at %g K, below %g K (%s)" % (T, round(bound, 1), basis), chk)
                rep.add("G1", "liquid_below_physical_minimum", r, "demoted" if gross else "flagged", T=T, bound=round(bound, 1))


def liquid_bound_any(rule, elements):
    """Lowest temperature at which a liquid can exist in this system: the SGTE bound of the system (twice the
    phase-rule margin), else the lowest bound over its binary subsystems, else 0.40 x the lowest melting point."""
    els = sorted(elements)
    sysn = "-".join(els)
    if set(els) & VOLATILE:
        # oxides, nitrides, carbides, phosphides ...: the calculated SGTE diagram rarely spans the whole system and
        # molecular / volatile components melt far below the metals; only the element melting points apply
        tms = [tm(e) for e in els if tm(e)]
        if not tms or len(tms) != len(els):
            return None, None
        return phase_rule.LIQ_TM_FRACTION * min(tms), "0.40 x lowest element melting point"
    if sysn in rule.liq_min:
        b, basis = rule.liquid_bound(sysn, phase_rule.EQUILIBRIUM_MARGIN_FACTOR)
        return b, "calculated SGTE diagram %s, lowest liquid %s K" % (basis.get("tdb"), basis.get("lowest_liquid_T_K"))
    subs = []
    for i in range(len(els)):
        for j in range(i + 1, len(els)):
            s = els[i] + "-" + els[j]
            if s in rule.liq_min:
                b, basis = rule.liquid_bound(s, phase_rule.EQUILIBRIUM_MARGIN_FACTOR)
                subs.append((b, s))
    tms = [tm(e) for e in els if tm(e)]
    if not tms or len(tms) != len(els):
        return None, None
    if subs and len(els) > 2:
        b, s = min(subs)
        # higher-order systems melt lower than their binaries: never above 0.40 x min Tm, and 25 % below the binary bound
        return min(0.75 * b, phase_rule.LIQ_TM_FRACTION * min(tms)), "binary subsystems (lowest %s) and element melting points" % s
    return phase_rule.LIQ_TM_FRACTION * min(tms), "0.40 x lowest element melting point"


# ----------------------------------------------------------------------------------------------- G2 magnitudes

VOLATILE = {"O", "N", "F", "CL", "BR", "I", "H", "S", "P", "SE", "AS", "C"}
NONMETAL = {"O", "N", "F", "CL", "BR", "I", "H", "S", "SE", "TE", "P", "C", "B", "AS", "SI", "GE"}
PER_ATOM_UNITS = re.compile(r"atom|at\.|g\.?\s*at|per mol-atom", re.I)
J_NOT_KJ = re.compile(r"^\s*(?:J(?!\s*/\s*\(?\s*g\b)|joule)", re.I)
KJ = re.compile(r"\bkJ\b", re.I)
PER_FORMULA_PRINTED = re.compile(r"formula|f\.\s*u\.|\bf\.?u\b|/\s*mol\s*f\b|mol\s+of\s+(?:compound|formula)|per\s+mole?\s+of\s+compound", re.I)
MASS_BASIS = re.compile(r"J\s*/\s*\(?\s*g\s*[·.\s]?\s*K\b|J\s*g\s*[-−⁻]\s*1\s*K|J\s*/\s*\(?\s*K\s*[·.\s]?\s*g\b(?!\s*[-.·]?\s*at)", re.I)
SCALE_1E3 = re.compile(r"(?:[×x·*]\s*)?10\s*\^?\s*\{?\s*3\s*\}?(?![\d.])|10³")


def physical_range(r):
    """(lo, hi, per_atom) for the value's property; None when no range applies."""
    p = r.get("property") or ""
    u = r.get("unit") or ""
    metallic = not (set(r["elements"]) & NONMETAL)
    per_atom = bool(PER_ATOM_UNITS.search(u))
    if p in ("HM_FORM", "GM_FORM"):
        return (-200.0, 60.0, per_atom) if metallic else (-650.0, 120.0, per_atom)
    if p in ("HM_MIX", "GM_MIX"):
        return (-150.0, 60.0, per_atom)
    if p == "PARTIAL_ENTHALPY":
        return (-450.0, 250.0, per_atom)
    if p == "SM_FORM":
        return (-60.0, 60.0, per_atom)
    if p == "SM":
        return (0.0, 150.0, per_atom)
    if p == "CPM" and "mJ" not in u:
        return (0.0, 120.0, per_atom)
    if p.startswith("ACR_") or p == "ACTIVITY":
        return (0.0, 1.0001, True)
    return None


def single_formula_atoms(r):
    cs = [c for c in (parity.parse_phase(l, None)[0] for l in labels(r) if l) if c]
    if len({tuple(sorted(c.items())) for c in cs}) != 1:
        return None, None
    c = cs[0]
    return float(sum(c.values())), parity.formula_str(c) if hasattr(parity, "formula_str") else None


def in_range(v, lo, hi):
    return lo <= v <= hi


def _scale_groups(records):
    """(asset, property) groups whose out-of-range values all fall in range after one factor 1e-3 or 1e3."""
    groups = defaultdict(list)
    for r in records:
        v = r.get("value")
        if isinstance(v, (int, float)) and not isinstance(v, bool) and physical_range(r) and r["record_kind"] != "running_text" \
                and not _j_printed_as_kj(r):
            groups[(asset_of(r), r.get("property"))].append(r)
    out = {}
    for k, g in groups.items():
        bad = [r for r in g if not _value_ok(r, r["value"])]
        if not bad:
            continue
        if all(_value_ok(r, -r["value"]) for r in bad):
            continue     # a missing minus sign (G3 checks it), not a J / kJ scale
        for f in (1e-3, 1e3):
            if len(bad) >= 2 and all(_value_ok(r, r["value"] * f) for r in bad) and len(bad) >= 0.8 * len(g):
                out[k] = (f, len(bad), len(g))
                break
    return out


def _j_printed_as_kj(r):
    v = r.get("value")
    up = str(r.get("unit_printed") or "")
    vp = r.get("value_printed")
    return isinstance(v, (int, float)) and KJ.search(r.get("unit") or "") and J_NOT_KJ.search(up) and not KJ.search(up) \
        and isinstance(vp, (int, float)) and abs(vp - v) < 1e-6 * max(1.0, abs(v)) and not _value_ok(r, v) and _value_ok(r, v / 1000.0)


def _value_ok(r, v):
    rg = physical_range(r)
    if not rg:
        return True
    lo, hi, per_atom = rg
    n_at, _ = single_formula_atoms(r)
    return in_range(v, lo, hi) or (not per_atom and n_at and in_range(v / n_at, lo, hi)) or \
        (not per_atom and not n_at and in_range(v, 12 * lo, 12 * hi))


def g2_magnitudes(records, visual, rep):
    scale = _scale_groups(records)
    for r in records:
        sk = (asset_of(r), r.get("property"))
        v = r.get("value")
        if sk in scale and isinstance(v, (int, float)) and not _value_ok(r, v):
            f, nb, ng = scale[sk]
            chk = visual.get("table_scale:%s:%s" % sk)
            if chk and chk.get("verdict") == "as_stored":
                fail(r, "G2", "value_outside_typical_range_as_printed", "not_espei_ready",
                     "%s = %g %s is outside the typical range; the page prints it so" % (r.get("property"), v, r.get("unit")), chk)
                rep.add("G2", "value_outside_typical_range_as_printed", r, "flagged")
                continue
            if chk and chk.get("verdict") == "unit_not_printed":
                r["unit_stored"] = r.get("unit")
                r["unit"] = "not printed"
                correct(r, "G2", "unit", r["unit_stored"], "not printed", "unit_not_printed_on_page",
                        {"magnitudes": "out of the %s range by x%g" % (r.get("property"), 1 / f)}, chk)
                fail(r, "G2", "unit_not_printed", "not_espei_ready", "no unit printed in the table; the stored unit was assumed", chk)
                rep.add("G2", "unit_not_printed", r, "corrected")
                continue
            if chk and chk.get("verdict") == "confirmed" and float(chk.get("factor", 0)) == f:
                new = round(v * f, 10)
                correct(r, "G2", "value", v, new, "table_unit_scale_%g" % f,
                        {"unit_printed": r.get("unit_printed"), "property_printed": r.get("property_printed"),
                         "out_of_range_in_table": "%d of %d" % (nb, ng), "header_seen": chk.get("seen")}, chk)
                r["value_stored"], r["value"] = v, new
                rep.add("G2", "table_unit_scale_applied", r, "corrected", before=v, after=new, factor=f)
            else:
                fail(r, "G2", "table_unit_scale_suspect", "not_accepted",
                     "%d of %d %s values of this table are out of range and all fall in range x%g (J/kJ scale)" % (
                         nb, ng, r.get("property"), f), chk)
                rep.add("G2", "table_unit_scale_suspect", r, "demoted", factor=f)
            continue
    for r in records:
        if any(g["code"] == "table_unit_scale_suspect" for g in r.get("gate_failures") or []) or \
                any(c.get("rule", "").startswith("table_unit_scale") for c in r.get("corrections") or []):
            continue
        v = r.get("value")
        if not isinstance(v, (int, float)) or isinstance(v, bool):
            continue
        p = r.get("property") or ""
        up = str(r.get("unit_printed") or "")
        pp = str(r.get("property_printed") or "")
        tkey = "table_scale:%s:%s" % (asset_of(r), p)
        chk = visual.get("record:" + r["id"]) or visual.get(tkey)
        # (a) printed x10^3 scale on an activity header: the printed numbers are a x 10^3
        if (p.startswith("ACR_") or p == "ACTIVITY") and SCALE_1E3.search(pp) and v > 1.0001:
            if chk and chk.get("verdict") == "confirmed":
                new = round(v / 1000.0, 10)
                correct(r, "G2", "value", v, new, "printed_scale_factor_1e3", {"property_printed": pp}, chk)
                r["value_stored"], r["value"] = v, new
                rep.add("G2", "activity_scale_1e3_applied", r, "corrected", before=v, after=new)
                continue
            fail(r, "G2", "activity_scale_factor_not_applied", "not_espei_ready",
                 "header '%s' states a x10^3 scale; value %g stored unscaled" % (pp, v), chk)
            rep.add("G2", "activity_scale_factor_not_applied", r, "flagged")
            continue
        # (b) a mass-basis heat capacity / entropy printed per gram stored per mole
        if p in ("CPM", "SM", "HEAT_CONTENT") and MASS_BASIS.search(up) and "g" not in (r.get("unit") or ""):
            r["unit_stored"] = r.get("unit")
            r["unit"] = "J/g/K"
            correct(r, "G2", "unit", r["unit_stored"], "J/g/K", "unit_from_printed_header", {"unit_printed": up},
                    visual.get("table_unit:%s" % asset_of(r)))
            fail(r, "G2", "mass_basis_unit", "not_espei_ready", "printed per gram (%s); not a molar quantity" % up)
            rep.add("G2", "mass_basis_unit_relabelled", r, "corrected", unit_printed=up)
            continue
        rg = physical_range(r)
        if not rg:
            continue
        lo, hi, per_atom = rg
        if p.startswith("ACR_") or p == "ACTIVITY":
            if v == 0:
                fail(r, "G9", "activity_rounded_to_zero", "not_espei_ready", "activity stored as 0 (printed text kept)")
                rep.add("G9", "activity_rounded_to_zero", r, "flagged")
            elif not in_range(v, lo, hi):
                fail(r, "G2", "value_out_of_physical_range", "not_accepted", "activity %g outside (0, 1]" % v, chk)
                rep.add("G2", "value_out_of_physical_range", r, "demoted", value=v)
            continue
        n_at, _ = single_formula_atoms(r)
        ok = in_range(v, lo, hi) or (not per_atom and n_at and in_range(v / n_at, lo, hi)) or \
            (not per_atom and not n_at and in_range(v, 12 * lo, 12 * hi))
        if ok:
            continue
        unit = r.get("unit") or ""
        # (c) a J value under a printed J unit stored as kJ (no conversion made)
        if _j_printed_as_kj(r):
            new = round(v / 1000.0, 10)
            correct(r, "G2", "value", v, new, "printed_unit_J_stored_as_kJ", {"unit_printed": up})
            r["value_stored"], r["value"] = v, new
            rep.add("G2", "J_to_kJ_applied", r, "corrected", before=v, after=new, unit_printed=up)
            continue
        # (d) a per-formula-unit value labelled per atom (printed header says per formula)
        if per_atom and n_at and n_at > 1 and PER_FORMULA_PRINTED.search(up + " " + pp) and in_range(v / n_at, lo, hi):
            new = round(v / n_at, 10)
            correct(r, "G2", "value", v, new, "per_formula_printed_divided_by_atoms",
                    {"unit_printed": up, "atoms_per_formula_unit": n_at})
            r["value_stored"], r["value"] = v, new
            r["basis_correction"] = {"atoms_per_formula_unit": n_at, "operation": "divide", "gate": "G2"}
            rep.add("G2", "per_formula_divided", r, "corrected", before=v, after=new)
            continue
        if unit in ("unknown", "unspecified", "") or unit.startswith("unknown"):
            fail(r, "G2", "unit_not_read_value_out_of_range", "not_espei_ready",
                 "value %g with no unit read lies outside the %s range of %s" % (v, "per-atom" if per_atom else "", p), chk)
            rep.add("G2", "unit_not_read_value_out_of_range", r, "flagged", value=v)
            continue
        if chk and chk.get("verdict") == "corrected" and chk.get("value") is not None:
            new = chk["value"]
            correct(r, "G2", "value", v, new, "value_read_from_page", {"printed": chk.get("value_printed_page")}, chk)
            r["value_stored"], r["value"] = v, new
            rep.add("G2", "value_corrected_from_page", r, "corrected", before=v, after=new)
            continue
        if chk and chk.get("verdict") == "as_stored" and chk.get("kind") in ("record_value", "values"):
            # read correctly: the page prints this value; it stays accepted but is not used for fitting
            fail(r, "G2", "value_outside_typical_range_as_printed", "not_espei_ready",
                 "%s = %g %s is outside [%g, %g]; the page prints it so" % (p, v, unit, lo, hi), chk)
            r["flags"]["magnitude_suspect"] = True
            rep.add("G2", "value_outside_typical_range_as_printed", r, "flagged", value=v)
            continue
        effect = "not_accepted"
        if r.get("evidence_class") == "dft" and re.search(r"relative|ΔE|\bE\b", pp + " " + str((r.get("locator") or {}).get("caption"))):
            effect = "not_espei_ready"   # relative energies of hypothetical structures can be large
        fail(r, "G2", "value_out_of_physical_range", effect,
             "%s = %g %s outside [%g, %g]%s" % (p, v, unit, lo, hi, "" if per_atom else " (per atom of the formula)"), chk)
        r["flags"]["magnitude_suspect"] = True
        rep.add("G2", "value_out_of_physical_range", r, "demoted" if effect == "not_accepted" else "flagged", value=v, unit=unit)


# ----------------------------------------------------------------------------------------------- G3 signs

NEG_HEADING = re.compile(r"^\s*[(\[]?\s*[−\-–‐]\s*\(?\s*Δ")


def g3_signs(records, visual, rep):
    by_table = defaultdict(list)
    for r in records:
        pp = str(r.get("property_printed") or "")
        v = r.get("value")
        if NEG_HEADING.match(pp) and isinstance(v, (int, float)) and v > 0 and not r["flags"].get("sign_from_heading") \
                and (r.get("property") or "").startswith(("HM", "GM", "SM")):
            by_table[(asset_of(r), pp)].append(r)
    for (asset, pp), rows in sorted(by_table.items(), key=lambda x: str(x[0])):
        chk = visual.get("table_sign:%s" % asset)
        if chk and chk.get("verdict") == "confirmed":
            for r in rows:
                v = r["value"]
                correct(r, "G3", "value", v, -v, "sign_from_printed_heading", {"heading": pp}, chk)
                r["value_stored"], r["value"] = v, -v
                r["flags"]["sign_from_heading"] = True
                r["sign_heading_printed"] = pp
                rep.add("G3", "sign_from_heading_applied", r, "corrected", before=v, after=-v, heading=pp)
        else:
            for r in rows:
                fail(r, "G3", "sign_heading_not_applied", "not_espei_ready",
                     "heading '%s' prints minus the quantity; value stored %g" % (pp, r["value"]), chk)
                rep.add("G3", "sign_heading_not_applied", r, "flagged")
    for r in records:
        if r["flags"].get("unsigned_in_print"):
            if fail(r, "G3", "unsigned_in_print", "not_espei_ready", "minus sign absent in print"):
                rep.add("G3", "unsigned_in_print_blocks_espei", r, "flagged")
        p = r.get("property")
        v = r.get("value")
        if p in ("HM_FORM", "GM_FORM") and isinstance(v, (int, float)) and r.get("evidence_class") != "dft":
            n_at, _ = single_formula_atoms(r)
            per_atom = bool(PER_ATOM_UNITS.search(r.get("unit") or ""))
            va = v if per_atom else (v / n_at if n_at else None)
            if n_at and n_at > 1 and va is not None and va > 10.0 and not r["flags"].get("sign_from_heading"):
                chk = visual.get("record:" + r["id"])
                if chk and chk.get("verdict") == "corrected" and chk.get("value") is not None:
                    new = chk["value"]
                    correct(r, "G3", "value", v, new, "sign_read_from_page", {"printed": chk.get("value_printed_page")}, chk)
                    r["value_stored"], r["value"] = v, new
                    r["flags"]["sign_pdf_page_verified"] = True   # the sign was read on the rendered page
                    rep.add("G3", "sign_corrected_from_page", r, "corrected", before=v, after=new)
                elif chk and chk.get("verdict") == "as_stored":
                    rep.counts[("G3", "positive_formation_confirmed_on_page", "checked")] += 1
                else:
                    fail(r, "G3", "formation_positive_check_sign", "not_espei_ready",
                         "formation %s of compound %s is +%g kJ/mol-atom" % (p, labels(r), round(va, 2)), chk)
                    rep.add("G3", "formation_positive_check_sign", r, "flagged", value=v)


# ----------------------------------------------------------------------------------------------- G5 semantics

FORMATION_WORD = re.compile(r"Δ\s*[_ ]?f\b|Δ_?f[GHS]|Δ[GHS]\s*[_ ]?f\b|[GHS]\s*f\b|formation", re.I)
MIX_WORD = re.compile(r"mix|solution|Δ\s*[_ ]?mix|excess|partial|integral", re.I)
FUSION_WORD = re.compile(r"fusion|melting|\bfus\b|Δ\s*_?fus|ΔH_?m\b|ΔS_?m\b", re.I)
LATTICE_WORD = re.compile(r"lattice stabilit|transformation|→|->|relative to|[A-Z]{3}_[A-Z]\d|fcc|bcc|hcp", re.I)
INCREMENT = re.compile(r"S\s*\(?\s*T\s*\)?\s*[-−–]\s*S\s*\(?\s*298|S\s*°?\s*_?T\s*[-−–]\s*S\s*°?\s*_?298|S_T\s*[-−–]\s*S_298", re.I)
RK_PARAM = re.compile(r"^\s*(?:\^?\s*\d\s*L|L\s*[_^]?\s*\{?\s*\d|L\s*\(|[ab]\s*$|[ab]_?\d\s*$|interaction\s+param|"
                      r"Redlich|RK\b|\d\s*L\s*$|.*\bcoefficient\b|.*\bcoeff\.)", re.I)


ROOM_T_LIQUID_PROPS = {"HM_FORM", "GM_FORM", "SM_FORM", "SM", "CPM", "ENTROPY", "HEAT_CONTENT"}
GAS_ELEMENTS = {"O", "N", "H", "F", "CL", "AR", "HE", "NE", "KR", "XE"}


def g5_semantics(records, rep):
    for r in records:
        p = r.get("property") or ""
        pp = str(r.get("property_printed") or "")
        txt = text_of(r)
        # absolute-entropy code on an excess / magnetic contribution or an entropy of formation
        if p in ("SM", "ENTROPY") and not INCREMENT.search(pp):
            cap = str((r.get("locator") or {}).get("caption") or "")
            new = None
            if re.search(r"\bmag", pp, re.I):
                new, rule = "SM_MAGNETIC", "printed_magnetic_contribution"
            elif re.search(r"\^?\s*ex\b|excess|\bE\b|S\s*\^?\s*E\b", pp):
                new, rule = "SM_EXCESS", "printed_excess_entropy"
            elif FORMATION_WORD.search(pp) or (re.fullmatch(r"\s*[sS]\s*(?:\(.*\))?\s*", pp) and re.search(r"entrop\w* \(?s?\)? ?of formation|formation", cap, re.I)):
                new, rule = "SM_FORM", "printed_entropy_of_formation"
            if new:
                correct(r, "G5", "property", p, new, rule, {"property_printed": pp, "caption": cap[:120] or None})
                r["property_code_before_gates"], r["property"] = p, new
                r["property_group"] = "formation" if new == "SM_FORM" else "heat capacity / entropy / enthalpy"
                if new != "SM_FORM":
                    fail(r, "G5", "entropy_contribution_not_absolute", "not_espei_ready", "printed %s" % pp)
                rep.add("G5", rule, r, "corrected", printed=pp)
                p = new
        # MIX code but printed as formation
        if p in ("GM_MIX", "HM_MIX", "SM_MIX"):
            head = pp or ""
            if FORMATION_WORD.search(head) and not MIX_WORD.search(head):
                new = p.replace("_MIX", "_FORM")
                correct(r, "G5", "property", p, new, "printed_quantity_is_formation", {"property_printed": pp})
                r["property_code_before_gates"] = p
                r["property"] = new
                r["property_group"] = "formation"
                if any(is_liq(x) for x in r.get("phases") or []):
                    fail(r, "G5", "formation_on_liquid", "not_espei_ready", "formation quantity of a liquid")
                rep.add("G5", "mix_to_form", r, "corrected", printed=pp)
                continue
        if p in ("HM_FORM", "GM_FORM", "SM_FORM", "HM_MIX", "GM_MIX", "SM_MIX") and len(r["elements"]) == 1:
            if FUSION_WORD.search(pp + " " + txt):
                new = {"H": "HM_FUSION", "G": "GM_FUSION", "S": "SM_FUSION"}[p[0]]
                rule = "unary_fusion_quantity"
            elif LATTICE_WORD.search(pp + " " + " ".join(str(x) for x in labels(r))):
                new = {"H": "HM_LATTICE_STABILITY", "G": "GM_LATTICE_STABILITY", "S": "SM_LATTICE_STABILITY"}[p[0]]
                rule = "unary_lattice_stability"
            else:
                new = None
            if new:
                correct(r, "G5", "property", p, new, rule, {"property_printed": pp, "system": r["system"]})
                r["property_code_before_gates"], r["property"] = p, new
                r["property_group"] = "heat capacity / entropy / enthalpy"
                fail(r, "G5", rule, "not_espei_ready", "unary %s (not a formation / mixing quantity)" % new)
                rep.add("G5", rule, r, "corrected", printed=pp)
            else:
                fail(r, "G5", "formation_code_in_unary_system", "not_accepted",
                     "%s in the one-element system %s (printed '%s')" % (p, r["system"], pp))
                rep.add("G5", "formation_code_in_unary_system", r, "demoted", printed=pp)
            continue
        if p in ("HM_FORM", "GM_FORM", "SM_FORM") and r.get("phases") and all(is_liq(x) for x in r["phases"]):
            fail(r, "G5", "formation_on_liquid", "not_espei_ready", "formation quantity assigned to a liquid")
            rep.add("G5", "formation_on_liquid", r, "flagged")
        T = r.get("temperature_K")
        if p in ROOM_T_LIQUID_PROPS and r.get("phases") and all(is_liq(x) for x in r["phases"]) and \
                isinstance(T, (int, float)) and T <= 310:
            tms = [tm(e) for e in r["elements"] if e not in GAS_ELEMENTS and tm(e)]
            if tms and min(tms) > 350:
                fail(r, "G5", "liquid_property_at_room_temperature", "not_espei_ready",
                     "%s of a liquid / amorphous phase at %g K, far below the melting points of the system's "
                     "metals: a hypothetical (model) state, not an equilibrium datum" % (p, T))
                rep.add("G5", "liquid_property_at_room_temperature", r, "flagged")
        if p in ("SM", "ENTROPY") and INCREMENT.search(pp):
            correct(r, "G5", "property", p, "SM_INCREMENT", "printed_entropy_increment", {"property_printed": pp})
            r["property_code_before_gates"], r["property"] = p, "SM_INCREMENT"
            fail(r, "G5", "entropy_increment_not_absolute", "not_espei_ready", "printed S(T)-S(298)")
            rep.add("G5", "entropy_increment", r, "corrected", printed=pp)
            continue
        if p.startswith(("SM", "HM", "GM", "ENTROPY", "CPM", "HEAT_CONTENT")) and RK_PARAM.match(pp):
            fail(r, "G5", "model_parameter_not_property", "not_accepted", "printed '%s' is a model (Redlich-Kister) parameter" % pp)
            rep.add("G5", "model_parameter_not_property", r, "demoted", printed=pp)


# ----------------------------------------------------------------------------------------------- G6 invariants

NON_PHASE = re.compile(r"^(?:unknown|other|liquidus|solidus|polymorphic|eutectic|peritectic|eutectoid|peritectoid|"
                       r"monotectic|congruent|max(?:imum)?|min(?:imum)?|un|c|e|p|u|\?|n/?a|none|-|–)$", re.I)
COMPOSITION_POINT = re.compile(r"^\s*[\d.]+\s*(?:at\.?\s*%|%|wt)|^\s*x\s*[=(]|^\s*[\d.]+\s*$", re.I)
THREE_PHASE = {"eutectic": 1, "peritectic": 1, "metatectic/catatectic": 1, "monotectic": 2, "syntectic": 2,
               "eutectoid": 0, "peritectoid": 0, "monotectoid": 0}
TWO_PHASE = {"congruent melting": 1, "congruent transformation": 0, "polymorphic": 0}


def _liq_key(p):
    n = norm_label(p.get("label"))
    if is_liq_g6(p):
        m = re.match(r"^\(?(?:liquid|liq\.?|melt|l)[#_]?(\d|'+|\"+|i{1,3}(?=\)?$))?", n)
        return "liquid" + ((m.group(1) or "") if m else n)
    return None


def normalize_phases(r):
    """Drop non-phase labels and repeated phases (the same phase under two spellings: phase_key, or the same
    unnumbered liquid printed twice); the dropped labels are kept in `phases_dropped`."""
    ph = r.get("phases") or []
    if len(ph) < 2:
        return []
    named = [p for p in ph if not NON_PHASE.match((p.get("label") or "").strip()) and not COMPOSITION_POINT.match(p.get("label") or "")]
    dropped = []
    if named and len(named) < len(ph):
        for p in ph:
            if p not in named:
                dropped.append({"label": p.get("label"), "reason": "not a phase (%s)" % (
                    "composition point" if COMPOSITION_POINT.match(p.get("label") or "") else "reaction / line / placeholder label")})
    else:
        named = list(ph)
    seen, keep = {}, []
    for p in named:
        k = _liq_key(p) or phase_key(p.get("label"))
        if k in seen:
            prev = seen[k]
            if not prev.get("model_phase") and p.get("model_phase"):
                keep[keep.index(prev)] = p
                seen[k] = p
                p, prev = prev, p
            dropped.append({"label": p.get("label"), "reason": "repeats phase '%s'" % prev.get("label")})
            continue
        seen[k] = p
        keep.append(p)
    if dropped:
        r.setdefault("phases_stored", [x.get("label") for x in ph])
        r["phases"] = keep
        r.setdefault("phases_dropped", []).extend(dropped)
    return dropped


REACTION_TYPES = set(THREE_PHASE) | set(TWO_PHASE) | {"other", "quasi-peritectic", "critical point", "magnetic"}


def _table_check(r, visual):
    return visual.get("invariant_table:" + str(asset_of(r)))


def _record_action(chk, r):
    return ((chk or {}).get("records") or {}).get(r["id"]) or {}


def g6_page_phase_list(r, chk, rep):
    """Phase-list corrections read on the page (tools/literature_visual_checks.json, key invariant_table:<asset>):
    stored labels that are not phases of the row (column headers, composition labels, a phase repeated under another
    spelling) are dropped, phases the printed reaction names but the reading missed are added, and the reaction type is
    set to the printed one. The stored list stays in `phases_stored`; compositions printed under a dropped label move
    to `phase_compositions_dropped`."""
    if not chk:
        return False
    act = _record_action(chk, r)
    # the page check judged the listed rows (`rows`); its dropped labels apply to those rows only
    drop = chk.get("drop_labels") or {} if (not chk.get("rows") or r["id"] in chk["rows"]) else {}
    ph = r.get("phases") or []
    changed = False
    gone = [p for p in ph if p.get("label") in drop]
    if gone and len(ph) - len(gone) >= 1:
        before = [p.get("label") for p in ph]
        r.setdefault("phases_stored", list(before))
        r["phases"] = [p for p in ph if p.get("label") not in drop]
        for p in gone:
            r.setdefault("phases_dropped", []).append({"label": p.get("label"), "reason": "not a phase of the row on the page: %s" % drop[p.get("label")]})
        pcs = r.get("phase_compositions") or []
        moved = [pc for pc in pcs if pc.get("phase") in drop]
        if moved:
            r["phase_compositions_dropped"] = moved
            r["phase_compositions"] = [pc for pc in pcs if pc.get("phase") not in drop]
        correct(r, "G6", "phases", before, [p.get("label") for p in r["phases"]], "non_phase_labels_dropped_on_page",
                {"dropped": {p.get("label"): drop[p.get("label")] for p in gone}}, chk)
        rep.add("G6", "phase_list_corrected_on_page", r, "corrected", dropped=[p.get("label") for p in gone])
        changed = True
    if act.get("action") == "add_labels" and act.get("labels"):
        have = {phase_key(p.get("label")) for p in r.get("phases") or []}
        add = [l for l in act["labels"] if phase_key(l) not in have]
        if add:
            before = [p.get("label") for p in r.get("phases") or []]
            r.setdefault("phases_stored", list(before))
            r["phases"] = list(r.get("phases") or []) + [{"label": l, "label_source": "page_check"} for l in add]
            correct(r, "G6", "phases", before, [p.get("label") for p in r["phases"]], "printed_reaction_phase_added_on_page",
                    {"added": add, "reaction_printed": (r.get("reaction") or {}).get("printed")}, chk)
            rep.add("G6", "phase_added_on_page", r, "corrected", added=add)
            changed = True
    if act.get("action") == "reaction_type" and act.get("type") in REACTION_TYPES and act["type"] != r.get("transition_kind"):
        before = r.get("transition_kind")
        r["transition_kind_stored"] = before
        r["transition_kind"] = act["type"]
        if isinstance(r.get("reaction"), dict) and r["reaction"].get("type") == before:
            r["reaction"]["type"] = act["type"]
        correct(r, "G6", "transition_kind", before, act["type"], "reaction_type_read_from_page",
                {"printed": act.get("printed")}, chk)
        rep.add("G6", "reaction_type_corrected_on_page", r, "corrected", before=before, after=act["type"])
        changed = True
    return changed


def g6_invariants(records, rep, visual=None):
    visual = visual or {}
    for r in records:
        if r.get("property") not in TEMPERATURE_PROPS and r.get("observation_type") not in ("invariant_reaction", "phase_transition"):
            continue
        d = normalize_phases(r)
        if d:
            rep.add("G6", "phase_list_normalized", r, "corrected", dropped=[x["label"] for x in d])
        chk = _table_check(r, visual)
        page_fixed = g6_page_phase_list(r, chk, rep)
        act = _record_action(chk, r)
        ph = r.get("phases") or []
        n = len(ph)
        if n == 0:
            continue
        nl = sum(1 for p in ph if is_liq_g6(p))
        tk = r.get("transition_kind")
        nel = len(r["elements"])
        why = None
        if nel == 2:
            if n > 3:
                why = ("binary_more_than_three_phases", "%d distinct phases in a binary invariant (%s): likely a "
                       "higher-order system or several reactions in one row" % (n, ", ".join(str(x) for x in labels(r))))
            elif tk in THREE_PHASE and n == 3 and nl != THREE_PHASE[tk]:
                why = ("reaction_liquid_count_inconsistent", "%s with %d liquid phase(s)" % (tk, nl))
            elif tk in THREE_PHASE and n == 2 and r.get("property") == "INVARIANT_TEMPERATURE":
                why = ("reaction_phase_count_inconsistent", "%s with 2 phases" % tk)
            elif tk in TWO_PHASE and n == 2 and nl != TWO_PHASE[tk]:
                why = ("reaction_liquid_count_inconsistent", "%s with %d liquid phase(s)" % (tk, nl))
        elif nel == 3 and n > 4:
            why = ("ternary_more_than_four_phases", "%d distinct phases in a ternary invariant" % n)
        if act.get("action") == "misread":
            fail(r, "G6", "phase_list_misread_on_page", "not_accepted", "page check: %s" % (act.get("note") or chk.get("seen") or "")[:200], chk)
            rep.add("G6", "phase_list_misread_on_page", r, "demoted")
            continue
        if act.get("action") == "higher_order" and not r["flags"].get("system_reassigned"):
            fail(r, "page", "higher_order_system_row", "not_accepted",
                 "page check: a reaction of a higher-order system; the page names no unambiguous system", chk)
            rep.add("page", "higher_order_system_row", r, "demoted")
            continue
        if why and act.get("action") == "as_printed":
            code = "binary_degenerate_reaction_as_printed" if why[0] == "binary_more_than_three_phases" else why[0] + "_as_printed"
            why = (code, why[1] + "; the page prints it so")
        if why:
            fail(r, "G6", why[0], "not_espei_ready", why[1], chk)
            rep.add("G6", why[0], r, "flagged", kind=tk, phases=labels(r))
            if why[0].startswith("binary_more"):
                for lab in labels(r):
                    if foreign_formula_elements(lab, r["elements"]):
                        r["flags"]["system_mismatch"] = True
                        break
        elif page_fixed:
            r["_g6_page_fixed"] = True


# ----------------------------------------------------------------------------------------------- G9 compositions

def formula_fraction(label, elements):
    c, _ = parity.parse_phase(label, [e.capitalize() for e in elements])
    if not c:
        return None
    tot = float(sum(c.values()))
    out = {e.upper(): float(n) / tot for e, n in c.items()}
    return out if set(out) <= set(elements) else None


def g9_compositions(records, visual, rep):
    for r in records:
        els = r["elements"]
        oc = r.get("overall_composition") or {}
        if len(oc) == len(els) and len(els) > 1 and all(isinstance(v, (int, float)) for v in oc.values()) and abs(sum(oc.values()) - 1) > 0.02:
            fail(r, "G9", "overall_composition_sum_not_one", "not_espei_ready", "overall composition sums to %.3f" % sum(oc.values()))
            rep.add("G9", "overall_composition_sum_not_one", r, "flagged")
        for pc in r.get("phase_compositions") or []:
            comp = pc.get("composition") or {}
            if len(comp) == len(els) and len(els) > 1 and str(pc.get("basis", "")).startswith("mole_fraction") \
                    and abs(sum(comp.values()) - 1) > 0.02:
                fail(r, "G9", "phase_composition_sum_not_one", "not_espei_ready",
                     "composition of %s sums to %.3f" % (pc.get("phase"), sum(comp.values())))
                rep.add("G9", "phase_composition_sum_not_one", r, "flagged")
                break
    # compositions copied from the first row of a table: the same printed composition on rows whose formulas differ
    tables = defaultdict(list)
    for r in records:
        if r["record_kind"] != "running_text" and r.get("phase_compositions"):
            tables[asset_of(r)].append(r)
    for asset, rows in tables.items():
        sig = Counter()
        for r in rows:
            for pc in r["phase_compositions"]:
                if pc.get("composition_source") == "printed":
                    sig[json.dumps(pc.get("composition"), sort_keys=True)] += 1
        for r in rows:
            for pc in r["phase_compositions"]:
                if pc.get("composition_source") != "printed" or not str(pc.get("basis", "")).startswith("mole_fraction"):
                    continue
                ff = formula_fraction(pc.get("phase"), r["elements"])
                if not ff or len(ff) < 2:
                    continue
                comp = pc.get("composition") or {}
                if not all(e in comp for e in ff):
                    continue
                dev = max(abs(comp.get(e, 0) - ff.get(e, 0)) for e in set(comp) | set(ff))
                if dev <= 0.15:
                    continue
                copied = sig[json.dumps(comp, sort_keys=True)] >= 3
                chk = visual.get("composition:%s" % r["id"]) or visual.get("composition_table:%s" % asset)
                if chk and chk.get("verdict") == "confirmed" and copied:
                    new = {e: round(v, 6) for e, v in sorted(ff.items())}
                    correct(r, "G9", "phase_compositions[%s]" % pc.get("phase"), comp, new,
                            "composition_from_printed_formula", {"formula": pc.get("phase"), "same_composition_rows": sig[json.dumps(comp, sort_keys=True)]}, chk)
                    pc["composition_stored"] = comp
                    pc["composition"] = new
                    pc["composition_source"] = "stoichiometry_of_printed_formula"
                    rep.add("G9", "composition_from_formula", r, "corrected", phase=pc.get("phase"))
                else:
                    code = "composition_copied_across_rows" if copied else "composition_label_mismatch"
                    # a printed composition that deviates from the label's formula is, on the page sample, a real
                    # solubility / homogeneity range as printed: flagged only; copied compositions block ESPEI
                    fail(r, "G9", code, "not_espei_ready" if copied else "flag_only",
                         "phase %s printed composition %s deviates %.2f from its formula" % (pc.get("phase"), comp, dev), chk)
                    rep.add("G9", code, r, "flagged", phase=pc.get("phase"))
                break


# ----------------------------------------------------------------------------------------------- G8 evidence class

CALC_WORDS = re.compile(r"\b(calculat\w*|assess\w*|predict\w*|optimi[sz]\w*|thermodynamic (?:model|description|dataset|database|parameters?)|"
                        r"CALPHAD|comput\w*|modell?ed|model|simulat\w*|Thermo-?Calc|Pandat|FactSage|fitted|estimated|evaluat\w*|"
                        r"extrapolat\w*|fictitious|agree\w*|compar\w*|close to|consistent with|this work|present (?:work|study|"
                        r"description|parameters|model)|correction|reproduc\w*)\b", re.I)
EXP_WORDS = re.compile(r"\b(measured|observed|determined|detected|DTA|DSC|calorimetr\w*|e\.?m\.?f\.?|EPMA|metallograph\w*|"
                       r"dilatometr\w*|thermal analysis|reported by)\b", re.I)
DFT_WORDS = re.compile(r"\b(first[- ]principles|ab[- ]initio|DFT|density functional|VASP|GGA|PBE|LDA)\b", re.I)
EVALUATION_TITLE = re.compile(r"^\s*(?:The\s+)?[A-Z][a-z]?(?:\s*[-–]\s*[A-Z][a-z]?){1,3}\s*\((?:[A-Za-z]+\s*[-–]\s*)+[A-Za-z]+\)")
DFT_EXCL = re.compile(r"CALPHAD|assess|optimi|thermodynamic|parameters|compar|agree|close to|consistent|this work|present|"
                      r"correction|Miedema|experiment", re.I)


def g8_evidence(records, visual, describe, rep):
    for r in records:
        if r["record_kind"] != "running_text":
            continue
        q = r.get("quote") or ""
        ev = r.get("evidence_class")
        chk = visual.get("evidence:" + r["id"])
        new, rule = None, None
        if ev == "model" and DFT_WORDS.search(q) and not DFT_EXCL.search(q):
            new, rule = "dft", "quote_wording_first_principles"
        elif ev == "model" and EXP_WORDS.search(q) and not CALC_WORDS.search(q) and not DFT_WORDS.search(q) \
                and (r.get("cited_sources") or r.get("citations_printed")):
            new, rule = "experimental", "quote_wording_experimental_cited"
        if new:
            if chk and chk.get("verdict") == "rejected":
                rep.counts[("G8", "evidence_change_rejected_on_check", "checked")] += 1
                continue
            if not (chk and chk.get("verdict") == "confirmed"):
                # wording alone is not enough (a page sample found assessed / handbook values and comparisons): only a
                # change confirmed on the page is made; the others are flagged for review
                fail(r, "G8", "evidence_class_review", "flag_only", "quote wording suggests %s (not checked on the page)" % new)
                rep.add("G8", "evidence_class_review", r, "flagged")
                continue
            correct(r, "G8", "evidence_class", ev, new, rule, {"quote_keywords": sorted({m.group(0).lower() for m in (DFT_WORDS if new == 'dft' else EXP_WORDS).finditer(q)})}, chk)
            r["evidence_class_before_gates"] = ev
            r["evidence_class"] = new
            r.setdefault("origin", {})["basis_gate"] = rule
            if new == "experimental" and not (r.get("origin") or {}).get("own_work_of_paper"):
                r["origin"]["class_gate"] = "experimental (cited)"
            fail(r, "G8", "evidence_class_reclassified", "not_espei_ready", "%s -> %s by quote wording" % (ev, new), chk)
            rep.add("G8", "%s_to_%s" % (ev, new), r, "corrected", quote=q[:160])
    # assessed / evaluated sources (ASM / JPE phase-diagram evaluations, titled "<A>-<B> (<Name>-<Name>)"): every value
    # read from them is flagged; a value given as the evaluation's own experimental work is not a measurement
    for r in records:
        if not describe:
            break
        meta = describe(r.get("paper")) or {}
        title = meta.get("title") or ""
        if not (EVALUATION_TITLE.search(title) and re.search(r"phase|JPE|Phs|alloy", meta.get("container") or "", re.I)):
            continue
        r["flags"]["from_evaluation_paper"] = True
        rep.counts[("G8", "from_evaluation_paper_flagged", "flagged")] += 1
        if r.get("evidence_class") == "experimental" and (r.get("origin") or {}).get("own_work_of_paper"):
            r["evidence_class_before_gates"] = "experimental"
            r["evidence_class"] = "unspecified"
            r["flags"]["assessed_source"] = True
            correct(r, "G8", "evidence_class", "experimental", "unspecified", "assessed_or_evaluated_source",
                    {"paper_title": title[:160]}, visual.get("evidence:" + r["id"]))
            rep.add("G8", "assessed_source_not_experimental", r, "corrected", title=title[:120])


# ----------------------------------------------------------------------------------------------- G4b flags

def g4_flags(records, rep):
    """system_mismatch / magnitude_suspect / source_disputed block ESPEI (phase_rule.BLOCKING_FLAGS also does);
    a system mismatch confirmed by a foreign element in a phase formula is not accepted (G4 above)."""
    for r in records:
        fl = r["flags"]
        for f in ("system_mismatch", "magnitude_suspect", "source_disputed"):
            if fl.get(f) and fl.get("espei_ready"):
                fail(r, "G4", f, "not_espei_ready", "record-builder flag %s" % f)
                rep.add("G4", f + "_blocks_espei", r, "flagged")


# ----------------------------------------------------------------------------------------------- G10 readers

READER_AGREEMENT_LABELS = {
    "multi_reader": "confirmed (two or more readers)",
    "single_reader": "unconfirmed (single reader)",
    "single_reader_text_layer": "single reader + PDF text layer",
}


def g10_readers(r, rep):
    """`reader_tier` is a pipeline status (which reader passes ran on the table / passage), not an agreement claim.
    The claim is `reader_agreement`; where the record builder did not classify a value (invariant tables, running
    text) it is derived from the readers whose readings carry the value (settled rows only)."""
    if r.get("reader_agreement"):
        return
    vis = sorted({x for x in r.get("readers") or [] if x != "text_layer"})
    if r["record_kind"] == "invariant_table" and r["flags"].get("row_state") not in (None, "settled"):
        return
    if r["record_kind"] == "running_text" and r["flags"].get("second_reader") not in (None, "done"):
        code = "single_reader" if len(vis) == 1 else None
    elif len(vis) >= 2:
        code = "multi_reader"
    elif len(vis) == 1:
        code = "single_reader_text_layer" if "text_layer" in (r.get("readers") or []) else "single_reader"
    else:
        code = None
    if code:
        r["reader_agreement"], r["reader_agreement_label"] = code, READER_AGREEMENT_LABELS[code]
        r["reader_agreement_basis"] = "derived_from_readers"
        rep.counts[("G10", "reader_agreement_derived_" + code, "derived")] += 1


# ----------------------------------------------------------------------------------------------- page verdicts

def g_page_verdicts(records, visual, rep):
    """Record-level page checks that found the stored value not to be what is printed (and no single correct value
    could be read): the value is not accepted. A row of a higher-order system stored under a lower-order one, without
    a system the page names unambiguously, is not accepted either."""
    for r in records:
        chk = visual.get("record:" + r["id"])
        if not chk or r["flags"].get("system_reassigned"):
            continue
        v = chk.get("verdict")
        if r.get("_g6_page_fixed") and chk.get("kind") == "record_phases":
            # the misreading was the phase list (header labels read as phases), corrected on the page check of the table
            rep.counts[("page", "phase_list_misread_corrected", "corrected")] += 1
            continue
        if v == "misread":
            fail(r, "page", "misread_on_page", "not_accepted", "page check: %s" % (chk.get("seen") or "")[:200], chk)
            rep.add("page", "misread_on_page", r, "demoted")
        elif v == "higher_order":
            fail(r, "page", "higher_order_system_row", "not_accepted", "page check: %s" % (chk.get("seen") or "")[:200], chk)
            rep.add("page", "higher_order_system_row", r, "demoted")


# ----------------------------------------------------------------------------------------------- page fixes / T_ref

STD_T_PRINTED = re.compile(r"^\s*~?\s*298(?:\.[012]0?)?\s*(?:K)?\s*(?:\((?:caption|subscript[^)]*)\)?.*)?$")
PAGE_FIX_FIELDS = ("value", "unit", "temperature_K", "evidence_class", "property", "uncertainty", "reference_state",
                   "temperature_basis", "cited_sources", "quote_page")


def g0_page_fixes(records, visual, rep):
    """Record-level corrections read on the rendered page (key page_fix:<id>, verdict corrected): the stored field is
    kept (`value_stored`, `unit_stored`, `temperature_K_stored`, `evidence_class_before_gates`,
    `property_code_before_gates`) and the record gets a `corrections` entry with the page check. Applied before the
    other gates, which then judge the corrected value."""
    for r in records:
        chk = visual.get("page_fix:" + r["id"])
        if not chk or chk.get("verdict") != "corrected":
            continue
        for f in PAGE_FIX_FIELDS:
            if f not in (chk.get("fields") or {}):
                continue
            new = chk["fields"][f]
            before = r.get(f)
            if before == new:
                continue
            store = {"value": "value_stored", "unit": "unit_stored", "temperature_K": "temperature_K_stored",
                     "evidence_class": "evidence_class_before_gates", "property": "property_code_before_gates",
                     "uncertainty": "uncertainty_stored", "reference_state": "reference_state_stored",
                     "temperature_basis": "temperature_basis_stored", "cited_sources": None,
                     "quote_page": None}[f]   # quote_page: the unreadable stored page is in quality_fixes.json changes
            if store:   # cited sources: the unresolved entries stay in the correction's `before`
                r.setdefault(store, before)
            r[f] = new
            if f == "temperature_K":
                if r.get("property") in TEMPERATURE_PROPS and r.get("value") == before:
                    r.setdefault("value_stored", r.get("value"))
                    r["value"] = new
                pl = r["flags"].pop("plausibility", None)
                if pl:
                    r["flags"]["plausibility_before_correction"] = pl
                    r["flags"]["calphad_check"] = "not_available"
                    r["flags"]["espei_certified"] = False
            if f == "evidence_class":
                r.setdefault("origin", {})["class_gate"] = chk.get("origin_class") or new
            if f == "quote_page" and new is not None:
                r["flags"]["page_unreadable"] = False   # the page was read on the rendered PDF
            if f in ("value", "unit") and r["flags"].get("magnitude_suspect"):
                # the magnitude was the suspicion; the page settles the value
                r["flags"]["magnitude_suspect_before_page_fix"] = True
                r["flags"]["magnitude_suspect"] = False
            ev = {"printed": chk.get("printed")}
            if f in (chk.get("conversions") or {}):
                # a value converted on the page's own evidence (e.g. to the standard reference state): the
                # conversion (source of the auxiliary data, numbers) is part of the correction
                ev["conversion"] = chk["conversions"][f]
            correct(r, "page", f, before, new, (chk.get("rules") or {}).get(f) or chk.get("rule") or "read_from_page", ev, chk)
            rep.add("page", "page_fix_" + f, r, "corrected", before=before, after=new)
    for r in records:
        chk = visual.get("flag:" + r["id"])
        if not chk:
            continue
        for fl in chk.get("flags") or ([{"code": chk["code"], "effect": chk.get("effect"), "detail": chk.get("detail")}] if chk.get("code") else []):
            fail(r, "page", fl["code"], fl.get("effect") or "not_espei_ready", (fl.get("detail") or chk.get("seen") or "")[:240], chk)
            rep.add("page", fl["code"], r, "flagged")
            if fl["code"] == "pressure_not_1_atm_on_page" and r.get("pressure_Pa") == 101325:
                # the page prints a pressure other than 1 atm: the assumed 1 atm is withdrawn (audit A6 H10)
                r["pressure_Pa_stored"], r["pressure_Pa"] = r["pressure_Pa"], None
                r.setdefault("pressure_basis_stored", r.get("pressure_basis"))
                r["pressure_basis"] = "printed_not_1_atm_see_gate_failures"
                rep.counts[("page", "pressure_1_atm_withdrawn", "corrected")] += 1
        if chk.get("related_records"):
            r["related_records"] = chk["related_records"]


def g0_standard_temperature(r, rep):
    """PI rule (2026-10-24) '298 K is 298.15 K': a temperature printed as 298 K / 298.2 K (or read from a '298'
    subscript) denotes the standard reference temperature and is stored as 298.15 K; the printed string is kept."""
    T = r.get("temperature_K")
    if not isinstance(T, (int, float)) or T not in (298.0, 298.1, 298.2) or r.get("temperature_range_K"):
        return
    tp = r.get("temperature_printed")
    if tp is not None and not STD_T_PRINTED.match(str(tp)):
        return
    r["temperature_K_stored"] = T
    r["temperature_K"] = 298.15
    if r.get("property") in TEMPERATURE_PROPS and r.get("value") == T:
        r["value_stored"], r["value"] = r.get("value"), 298.15
    correct(r, "G1", "temperature_K", T, 298.15, "standard_reference_temperature",
            {"printed": tp if tp is not None else (r.get("property_printed") or (r.get("locator") or {}).get("caption") or "")[:120]})
    rep.counts[("G1", "standard_reference_temperature_298.15", "corrected")] += 1


# ----------------------------------------------------------------------------------------------- driver

def expand_visual(visual):
    """Table-level page checks name single rows of a higher-order system (with the system the page names): they are
    turned into record-level checks (record:<id>, verdict higher_order, system) so that G4 reassigns the system."""
    out = dict(visual)
    for k, c in visual.items():
        if not k.startswith("invariant_table:"):
            continue
        for rid, act in (c.get("records") or {}).items():
            if act.get("action") != "higher_order" or not act.get("system"):
                continue
            prev = dict(out.get("record:" + rid) or {})
            e = {x: c[x] for x in ("page", "pdf_sha256_prefix", "image_sha256", "checked") if c.get(x) is not None}
            e.update({"key": "record:" + rid, "kind": "record_phases", "verdict": "higher_order", "system": act["system"],
                      "seen": (act.get("note") or c.get("seen") or "")[:300], "from_table_check": k})
            if prev and not prev.get("system"):
                prev.update({"system": act["system"], "from_table_check": k})
                e = prev
            elif prev:
                continue
            out["record:" + rid] = e
    return out


def apply(records, rule, visual=None, store_units=None, plaus_store=None, describe=None):
    """Mutates the records; returns (records, report dict)."""
    _visual_given = visual is not None
    visual = expand_visual(load_visual() if visual is None else visual)
    rep = Report()
    acc0 = sum(1 for r in records if r["flags"].get("extraction_accepted"))
    ready0 = sum(1 for r in records if r["flags"].get("espei_ready"))
    n0 = len(records)
    g0_page_fixes(records, visual, rep)
    g4_system(records, visual, rep)
    g4_unit_gas_elements(records, rep)
    for r in records:
        g0_standard_temperature(r, rep)
        g12_pressure(r, rep)
        g11_plausibility(r, plaus_store or {}, rep)
        g10_readers(r, rep)
    g1_temperatures(records, rule, visual, store_units or {}, rep)
    g5_semantics(records, rep)
    g2_magnitudes(records, visual, rep)
    g3_signs(records, visual, rep)
    g6_invariants(records, rep, visual)
    g9_compositions(records, visual, rep)
    g8_evidence(records, visual, describe, rep)
    g4_flags(records, rep)
    g_page_verdicts(records, visual, rep)
    records = g7_dedupe(records, rep)
    for r in records:
        r.pop("_g6_page_fixed", None)
        if r.get("gate_failures"):
            r["gate_failures"].sort(key=lambda x: (x["gate"], x["code"]))
            r["flags"]["gate_failures"] = [x["code"] for x in r["gate_failures"]]
    report = {
        "release": GATE_RELEASE,
        "gates": GATES,
        "records_in": n0, "records_out": len(records),
        "accepted_before": acc0, "accepted_after": sum(1 for r in records if r["flags"].get("extraction_accepted")),
        "espei_ready_before": ready0, "espei_ready_after": sum(1 for r in records if r["flags"].get("espei_ready")),
        "espei_ready_after_note": ("espei_ready after the gates and before the ESPEI phase rule and blocking flags "
                                   "(export_literature_data.espei_fields); the released count is quality_fixes.json "
                                   "espei.espei_ready"),
        "records_corrected": sum(1 for r in records if r.get("corrections")),
        "records_with_gate_failures": sum(1 for r in records if r.get("gate_failures")),
        "visual_checks": len(visual),
        "visual_check_entries": (len(json.load(open(VISUAL_FILE, encoding="utf-8")).get("checks", []))
                                 if not _visual_given and os.path.exists(VISUAL_FILE) else None),
        "visual_checks_note": ("visual_checks counts the keys the gates look up (a table-level ledger entry that names "
                               "single rows of a higher-order system expands to one key per row); visual_check_entries "
                               "counts the entries of tools/literature_visual_checks.json (one rendered page each)"),
        "visual_check_summary": [{"kind": k, "verdict": v, "n": n} for (k, v), n in sorted(Counter(
            (c.get("kind") or k.split(":")[0], c.get("verdict")) for k, c in visual.items()).items())],
        "counts": [{"gate": g, "code": c, "action": a, "n": n} for (g, c, a), n in sorted(rep.counts.items())],
        "items": {"%s:%s" % k: v for k, v in sorted(rep.items.items())},
    }
    return records, report


def gates_md(report):
    """QUALITY_FIXES.md section for the gates."""
    L = ["## 7. Literature-level consistency gates (%s)" % report["release"], "",
         "Every value was checked against its table, its paper, the physics of its quantity and the calculated SGTE "
         "diagrams of its system (tools/literature_gates.py; tests tools/tests/test_literature_gates.py). No record is "
         "deleted and no printed value or quote is altered. A correction is made only when it is deterministic and "
         "evidenced (a printed header, caption or heading, or the internal consistency of a whole table) and, for "
         "table-level unit / sign / scale changes and system reassignments, confirmed on the rendered PDF page "
         "(`visual_check` {page, image sha256, verdict, what was seen}; tools/literature_visual_checks.json). The "
         "stored value is kept (`value_stored`, `temperature_K_stored`, `system_stored`, `property_code_before_gates`, "
         "`evidence_class_before_gates`, `phases_stored`) and the record carries `corrections`. Otherwise the record "
         "is demoted or flagged: `gate_failures` [{gate, code, effect, detail}] (codes also in `flags.gate_failures`); "
         "effect `not_accepted` sets `extraction_accepted` false, `not_espei_ready` sets `espei_ready` and "
         "`espei_certified` false, `flag_only` changes nothing else.", "",
         "| | before | after |", "|---|---:|---:|",
         "| records | %d | %d |" % (report["records_in"], report["records_out"]),
         "| extraction_accepted | %d | %d |" % (report["accepted_before"], report["accepted_after"]),
         "| espei_ready (after the gates, before the ESPEI phase rule and blocking flags; released count in `espei.espei_ready`) | %d | %d |" % (report["espei_ready_before"], report["espei_ready_after"]),
         "| records with corrections | | %d |" % report["records_corrected"],
         "| records with gate failures | | %d |" % report["records_with_gate_failures"],
         "| page-check ledger entries (tools/literature_visual_checks.json) | | %s |" % (report.get("visual_check_entries") if report.get("visual_check_entries") is not None else "-"),
         "| page-check keys looked up by the gates (table entries expanded to rows) | | %d |" % report["visual_checks"], "",
         "| gate | code | action | n |", "|---|---|---|---:|"]
    for c in report["counts"]:
        L.append("| %s %s | `%s` | %s | %d |" % (c["gate"], GATES.get(c["gate"], ""), c["code"], c["action"], c["n"]))
    L.append("")
    if report.get("visual_check_summary"):
        L += ["Visual checks (rendered PDF pages inspected one by one: every table-level group, every system reassignment, "
              "and samples of at least 50 per demotion category where the category is larger; `verdict_reader` keeps the "
              "reader's verdict where a decision was curated):", "", "| kind | verdict | n |", "|---|---|---:|"]
        for c in report["visual_check_summary"]:
            L.append("| %s | %s | %d |" % (c["kind"], c["verdict"], c["n"]))
        L.append("")
    L += ["Gate rules:", "",
          "* G1 temperatures. (a) Per table, the stored temperatures (as K, or converted from °C) are compared with "
          "the table's own numbers: rows describing a pure element against its melting point in K and °C, invariant "
          "rows against the system's SGTE invariants and melting points (median distance, K vs °C reading), and the "
          "readers' unit votes. A table whose stored convention disagrees with that evidence is converted as a whole "
          "from the printed numbers once the page confirms the header unit; unconfirmed: `temperature_unit_table_suspect`. "
          "The CALPHAD cross-check of a corrected value was computed on the old temperature and is moved to "
          "`flags.plausibility_before_correction`. (b) A temperature printed as a dash, T <= 1 K for a phase "
          "equilibrium, or a temperature equal to a printed composition: `temperature_not_a_temperature` (not accepted). "
          "(c) A liquid-bearing value below the lowest temperature a liquid can have (SGTE diagram of the system minus "
          "twice the phase-rule margin; binary subsystems for higher orders; else 0.40 x the lowest melting point): "
          "`liquid_below_physical_minimum` (not ESPEI-ready; not accepted below 0.25 x the lowest melting point). "
          "(d) A printed range stored as one temperature: `temperature_range_stored_as_point`.",
          "* G2 magnitudes. Physical ranges (kJ/mol-atom unless stated): formation enthalpy / Gibbs energy [-200, 60] "
          "(metallic) or [-650, 120] (with O, N, halogens, chalcogens, pnictogens, C, B, Si, Ge); mixing [-150, 60]; "
          "partial [-450, 250]; entropy of formation [-60, 60] J/(mol-atom K); absolute entropy [0, 150]; heat capacity "
          "[0, 120]; activity (0, 1]. Values of unstated basis are tested per atom of their single printed formula. "
          "Corrections: a printed x10^3 activity header (values divided by 1000 when the page confirms the header); a "
          "J value under a printed J unit stored as kJ; a per-formula value under a printed per-formula unit labelled "
          "per atom; a per-gram unit relabelled. Unexplained values out of range: `value_out_of_physical_range` (not "
          "accepted; relative DFT energies only not ESPEI-ready); no unit read: `unit_not_read_value_out_of_range`.",
          "* G3 signs. A positive value under a printed (-Δ...) heading is negated when the page confirms the heading "
          "(`sign_from_heading`); `unsigned_in_print` blocks ESPEI; a positive formation enthalpy / Gibbs energy above "
          "+10 kJ/mol-atom of a compound is flagged `formation_positive_check_sign` unless the page confirms it.",
          "* G4 system. A phase formula naming an element outside the system (or cement notation C-S = CaO-SiO2) moves "
          "the value to the system named by its phase formulas and caption when the page confirms it; otherwise "
          "`system_mismatch` (not accepted). system_mismatch, magnitude_suspect and source_disputed block ESPEI. "
          "A system element that is a pressure unit (Pa in kPa / MPa), a cover or flow gas (Ar, He) or a garbled phase "
          "prefix is never a component: such records are reassigned from the rendered page (record checks of "
          "2026-10-24) or, without a page check, not accepted (`system_element_from_unit_or_gas`).",
          "* G5 semantics. A mixing code whose printed quantity is a formation quantity becomes the formation code; "
          "formation / mixing codes in one-element systems become fusion or lattice-stability codes when printed as "
          "such, otherwise are not accepted; S(T)-S(298) becomes `SM_INCREMENT`; formation quantities of a liquid are "
          "flagged; Redlich-Kister parameters stored as properties are not accepted.",
          "* G6 invariants. Phase lists drop non-phase labels (composition points, 'liquidus', 'eutectic', 'UNKNOWN' "
          "beside a named phase) and repeated phases (`phases_dropped`, `phases_stored`); binary three-phase reactions "
          "must have three phases with the right number of liquids, congruent melting two phases with one liquid, "
          "solid-state reactions no liquid; a binary row with four or more phases is flagged (likely a higher-order "
          "system). Repeats are found with a spelling-independent key (`phase_key`: subscripts, LaTeX and spelled "
          "Greek letters, prime glyphs, separators, outer parentheses); numbered liquids (L1, L2, Liquid#2) count as "
          "liquids. Page checks of whole tables (cleanup 2026-10-24, key `invariant_table:<asset>`): labels that are "
          "not phases of the row (column headers such as 'Phase 2', 'X^M', 'T/K', composition labels, a phase from a "
          "neighbouring column) are dropped (compositions printed under them move to `phase_compositions_dropped`); "
          "a phase the printed reaction names but the reading missed is added (`label_source: page_check`); the "
          "reaction type is set to the printed one (`transition_kind_stored`); a row of a higher-order system moves "
          "to the system the page names (else not accepted); a row whose odd phase count or type is printed so stays "
          "accepted and not ESPEI-ready (`*_as_printed`, `binary_degenerate_reaction_as_printed`).",
          "* G7 duplicates. Values identical in paper, property, value, unit, temperature, table and row (or quote) "
          "that differ only cosmetically (UNKNOWN beside a named phase, glyph variants, reader lists) are merged into "
          "the most specified record (`merged_ids`). Labels are compared with `phase_key`, by contained words "
          "('Fe3O4' / 'Magnetite (Fe3O4)') and by formula composition ('UC1.90' / 'C0.655U0.345'); two readings of one "
          "printed value with different phases are kept apart and flagged `phase_label_readings_disagree` (not "
          "ESPEI-ready).",
          "* Display fix (2026-10-24, not a data change): the compound pages' DFT table labelled the energy above hull "
          "'eV/atom', but `dft_parity.json` stores it in kJ/mol-atom (`ehull`; Parquet `ehull_kJ_per_mol_atom`). The "
          "column now reads 'E above hull (kJ/mol-atom)' with an eV/atom column beside it; the MCP server's "
          "`e_above_hull_eV` is now converted (and `e_above_hull_kJ_per_mol_atom` added).",
          "* Standard reference temperature (PI rule 2026-10-24: '298 K is 298.15 K'). A temperature stored as 298 / "
          "298.1 / 298.2 K and printed as 298 K (or read from a '298' subscript) is stored as 298.15 K, printed string "
          "kept (`standard_reference_temperature`).",
          "* Page fixes (cleanup 2026-10-24, keys `page_fix:<id>`, `flag:<id>`): a value, unit, temperature, property "
          "code or evidence class read on the rendered page replaces the stored one (kept in `*_stored` / "
          "`*_before_gates`, `corrections` with the page check) before the other gates run; values the page prints as "
          "stored but outside the physical range stay accepted and are not ESPEI-ready "
          "(`value_outside_typical_range_as_printed`).",
          "* G8 evidence class. Running-text model values whose sentence names only first-principles methods become "
          "`dft`; those whose sentence names only experimental methods (and cite a source) become `experimental` "
          "(cited); both only when the page check confirms the sentence context (otherwise `evidence_class_review`). Values "
          "read from ASM / JPE phase-diagram evaluations (title '<A>-<B> (<Name>-<Name>)') carry "
          "`flags.from_evaluation_paper`; an evaluation's own value labelled experimental becomes `unspecified` with "
          "`flags.assessed_source` (cited measurements stay experimental). CALPHAD 'thermodynamic assessment' papers with "
          "their own experiments are not evaluations in this sense.",
          "* G9 compositions. Compositions not summing to 1, printed compositions far from the printed formula (and "
          "compositions repeated from another row) are flagged; recomputed from the printed formula when the page "
          "confirms the formula; activities stored as 0 are not ESPEI-ready. A printed composition far from its "
          "label's formula is, on the page sample, a real solubility / homogeneity range as printed: flag only "
          "(`composition_label_mismatch`); compositions copied across rows block ESPEI.",
          "* G10 readers. `reader_tier` is the pipeline status of the table or passage (which reader passes ran), "
          "not a claim about the value; a row read only by a later reader can carry the table's earlier tier. The "
          "agreement claim is `reader_agreement`; where the record builder did not classify a value (invariant tables, "
          "running text) it is derived from `readers` (`reader_agreement_basis: derived_from_readers`; settled rows "
          "only).",
          "* G11 plausibility. `flags.plausibility` carries dx (mole-fraction differences per phase) and the "
          "calculation; fails with |dT| < 5 K have `fail_reason` (composition; pure-element label compared on the "
          "other axis).",
          "* G12 pressure. 1 atm defaults (stated_default_1_atm, unstated_assumed) are relabelled "
          "`assumed_default_1_atm_not_printed` (`pressure_basis_stored` kept).", ""]
    return "\n".join(L) + "\n"


def load_store_evidence(corpus):
    """(pid -> reader temperature-unit votes, pid -> record-builder plausibility) from the record store (read only)."""
    import glob
    db = os.path.join(corpus, "db-v2.1")
    units, plaus = {}, {}
    for sub in ("invariants/master", "general/records"):
        for f in sorted(glob.glob(os.path.join(db, sub, "*.jsonl"))):
            reads, obs = {}, []
            with open(f, encoding="utf-8") as fh:
                for line in fh:
                    if '"record_type": "reading"' in line or '"record_type":"reading"' in line:
                        rd = json.loads(line)
                        pr = rd.get("printed") or {}
                        reads[rd["id"]] = (pr.get("temperature_unit"), pr.get("temperature_printed"))
                    elif '"observation"' in line:
                        o = json.loads(line)
                        if o.get("record_type") == "observation":
                            obs.append(o)
            for o in obs:
                c = Counter(reads[x][0] for x in o.get("readings") or [] if x in reads and reads[x][0])
                tp = sorted({str(reads[x][1]) for x in o.get("readings") or [] if x in reads and reads[x][1] is not None})
                if c or tp:
                    units[o.get("pid")] = {"votes": dict(c), "printed": tp}
                pl = ((o.get("quality") or {}).get("status") or {}).get("plausibility") or {}
                if pl.get("dx") is not None or pl.get("calculation"):
                    plaus[o.get("pid")] = {"dx": pl.get("dx"), "calculation": pl.get("calculation")}
    return units, plaus
