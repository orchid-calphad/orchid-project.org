#!/usr/bin/env python3
"""Deterministic phase mapping: a printed phase (handbook or paper table) -> a phase of a specific SGTE TDB file
(2026-10-21). Works from the file's own PHASE declarations and the 'Prototype X (pearson, sg)' comment before each.
Rules (first that applies; every mapping records its basis):
  liquid / gas        L, L#k, L1.. -> LIQUID (or IONIC_LIQ); G -> GAS
  solution by Pearson (X), (aX), (A,B): the TDB phase whose comment Pearson symbol equals the printed one (unique)
  solution by name    Strukturbericht A1/A2/A3... -> FCC_A1/BCC_A2/HCP_A3... if the TDB declares it
  compound by formula normalised formula == TDB phase name without suffix; Pearson must agree when both are known;
                      several variants (e.g. ALAU2, ALAU2_BETA) -> the Greek prefix chooses, else open
  rock salt (B1)       own NaCl-type TDB phase (prototype NaCl, name with B1/HALITE/ROCKSALT), else FCC_A1 (M)(C|N,VA)
                      for carbides/nitrides (SGTE convention); 2026-10-24
Anything else stays open with its candidates. Usable for handbook phases and for table phases (invariant/general).

Audit round 2 (2026-10-24, A2 classes 2-3; under version control in orchid-project.org/tools since then):
  * a printed Pearson symbol or Strukturbericht that contradicts the formula-matched TDB phase blocks the match
    ("polymorph not modelled"); "+pearson" is written only when the printed and the TDB Pearson symbols agree
    (rhombohedral hR counts compared on both the rhombohedral and the hexagonal cell, x3);
  * liquid, glass and amorphous TDB phases (LIQUID, IONIC_LIQ, ION, *_L, comment "liquid"/"glass") never match a
    formula; ION is a liquid for L labels;
  * the Greek table covers alpha..omega; a Greek prefix selects a variant by whole name token (eta never matches BETA
    or ZETA);
  * a solution "(El)" matched by Pearson symbol must have El among the TDB phase's constituents;
  * Pearson-only standard names (hP4 -> DHCP, hP12 -> C14 ...) apply only when no phase of the printed formula
    exists; a printed Strukturbericht key still applies first;
  * a label containing '?' is never mapped.
"""
import glob, json, os, re, sys
REPO = os.environ.get("OPEN_TDB", os.path.expanduser("~/projects/open-tdb"))
STD = {"A1": "FCC_A1", "A2": "BCC_A2", "A3": "HCP_A3", "A4": "DIAMOND_A4", "A5": "BCT_A5", "A6": "TETRAGONAL_A6",
       "A7": "RHOMBOHEDRAL_A7", "A8": "HEXAGONAL_A8", "A9": "GRAPHITE_A9", "A12": "CBCC_A12", "A13": "CUB_A13", "A20": "ORTHORHOMBIC_A20"}
# standard CALPHAD names for non-A structures (Strukturbericht or Pearson symbol -> name prefix)
STD_COMPOUND = {"C14": "C14_LAVES", "hP12": "C14_LAVES", "C15": "C15_LAVES", "cF24": "C15_LAVES", "C36": "C36_LAVES", "hP24": "C36_LAVES",
                "D8B": "SIGMA", "tP30": "SIGMA", "D85": "MU_D85", "hR13": "MU_D85", "hR39": "MU_D85", "A12": "CBCC_A12", "cI58": "CBCC_A12",
                "A13": "CUB_A13", "cP20": "CUB_A13", "B3": "ZINCBLENDE_B3", "A3'": "DHCP", "hP4": "DHCP",
                "hR105": "BETA_RHOMBO_B", "hR12": "ALPHA_RHOMBO_B"}
# 2026-10-24 (James Saal review): "cF8": "DIAMOND_A4" removed. cF8 is shared by diamond (A4), rock salt (B1, NaCl) and
# sphalerite (B3, ZnS), and the handbook tables also print cF8 for FeSi-type (B20, really cP8) compounds, so every NaCl-type
# carbide/telluride, sphalerite SiC and the FeSi-type silicides were named DIAMOND_A4. A Pearson-only key now applies only
# when no contradicting Strukturbericht is printed; rock-salt phases go through map_rocksalt().
INTERSTITIAL_B1 = ("C", "N")   # SGTE convention: MC / MN rock-salt carbides and nitrides are FCC_A1 (M)(C,VA) / (M)(N,VA)
GREEK = {"α": "ALPHA", "β": "BETA", "γ": "GAMMA", "δ": "DELTA", "ε": "EPSILON", "ζ": "ZETA", "η": "ETA", "θ": "THETA",
         "ι": "IOTA", "κ": "KAPPA", "λ": "LAMBDA", "μ": "MU", "ν": "NU", "ξ": "XI", "ο": "OMICRON", "π": "PI", "ρ": "RHO",
         "σ": "SIGMA", "τ": "TAU", "υ": "UPSILON", "φ": "PHI", "ϕ": "PHI", "χ": "CHI", "ψ": "PSI", "ω": "OMEGA"}
GREEK_CHARS = "".join(GREEK)
LIQUID_NAMES = ("LIQUID", "IONIC_LIQ", "ION", "IONIC_LIQUID", "LIQUID_IONIC")
LIQUID_LIKE = re.compile(r"^(?:LIQUID|LIQ|IONIC_LIQ|IONIC_LIQUID|LIQUID_IONIC|ION|KAP_LIQ|ASSOCI|OXIDE_L|SLAG|GLASS)$|"
                         r"(?:^|_)(?:LIQ|LIQUID|GLASS|AMORPH\w*)$")   # not '_L' alone: CU6SN5_L is the low-T phase


def liquid_like(name, info=None):
    """A TDB phase that is a liquid, glass or amorphous phase (never matched to a printed solid compound)."""
    return bool(LIQUID_LIKE.search(name.upper())) or bool(re.search(r"(?i)liquid|glass|amorph|\bmelt\b",
                                                                    (info or {}).get("comment") or ""))


def pearson_equal(a, b):
    """Printed and TDB Pearson symbols denote the same cell: equal after '*' is stripped, or rhombohedral hR counts that
    differ by the factor 3 between the rhombohedral and the hexagonal setting (the handbook prints hexagonal counts)."""
    a, b = (a or "").rstrip("*").strip(), (b or "").rstrip("*").strip()
    if not a or not b:
        return None
    if a == b:
        return True
    ma, mb = re.fullmatch(r"(hR)(\d+)", a), re.fullmatch(r"(hR)(\d+)", b)
    if ma and mb:
        na, nb = int(ma.group(2)), int(mb.group(2))
        return na == 3 * nb or nb == 3 * na
    return False


def name_tokens(p):
    return [t for t in re.split(r"[_#]+", p.upper()) if t]

def tdb_phases(path):
    """{phase: {'pearson':..,'prototype':..}} from PHASE lines and the comment just before each."""
    out, last, comment = {}, None, None
    for line in open(path, encoding="latin-1"):
        m = re.match(r"\s*\$\s*(?:Prototype\s+)?(\S.*?)\s*\((\w+\*?),\s*([^)]*)\)", line)
        if m:
            last = {"prototype": m.group(1), "pearson": m.group(2), "space_group": m.group(3)}
            continue
        mc = re.match(r"\s*\$+\s*(\S.*)$", line)
        if mc:
            comment = mc.group(1).strip()
            continue
        m = re.match(r"\s*PHASE\s+([A-Z0-9_#]+)(?::\S+)?\s", line, re.I)
        if m:
            out[m.group(1).upper()] = dict(last or {}, **({"comment": comment} if comment else {})); last = None
            comment = None
    txt = open(path, encoding="latin-1").read()
    for m in re.finditer(r"CONST\w*\s+([A-Z0-9_#]+)(?::\S+)?\s*:(.*?)!", txt, re.I | re.S):
        ph = m.group(1).upper()
        if ph in out:
            out[ph]["constituents"] = [sorted({x.strip().upper().rstrip("%") for x in re.split(r"[,\s]+", sl) if x.strip()})
                                       for sl in m.group(2).split(":") if sl.strip()]
    return out


def is_rocksalt(info):
    """Printed NaCl type: Strukturbericht B1 or prototype NaCl (the Pearson symbol cF8 alone is ambiguous)."""
    sb = norm_sb((info or {}).get("strukturbericht"))
    proto = re.sub(r"\s", "", (info or {}).get("prototype") or "")
    return sb == "B1" or proto.upper() == "NACL"


def map_rocksalt(n, info, tdb, comps, fd, elements):
    """A printed rock-salt (B1) compound with no phase of its own formula: (1) the TDB phase whose comment prototype is NaCl
    or whose name carries B1 / HALITE / ROCKSALT; else (2) for a carbide or nitride, FCC_A1 when the TDB models it as the
    interstitial solution (M)(X,VA) of that non-metal X (SGTE convention, e.g. HfC, NbC, TiC, ZrC, VC, TaC, beta-WC).
    Composition ranges, when printed and computable, must overlap. Returns (phase, basis, candidates) or None."""
    rng = (info or {}).get("composition_at_pct")

    def overlaps(p):
        if not rng or comps is None or p not in comps:
            return True
        return comps[p][0] <= rng[1] / 100 + 0.03 and comps[p][1] >= rng[0] / 100 - 0.03
    own = [p for p, v in tdb.items() if re.sub(r"\s", "", v.get("prototype") or "").upper() == "NACL"
           or "B1" in p.split("_") or any(t in p for t in ("HALITE", "ROCKSALT"))]
    own = [p for p in own if overlaps(p)]
    if len(own) == 1:
        return own[0], "rocksalt_b1_by_tdb_phase", own
    if own:
        return None, "rock salt: several candidates", own
    nonmetal = [e for e in (fd or {}) if e in INTERSTITIAL_B1] if fd and len(fd) >= 2 else \
               [e for e in (elements or []) if e in INTERSTITIAL_B1]
    fcc = tdb.get("FCC_A1", {}).get("constituents") or []
    if len(nonmetal) == 1 and len(fcc) == 2 and nonmetal[0] in fcc[1] and "VA" in fcc[1] and overlaps("FCC_A1"):
        return "FCC_A1", "rocksalt_b1_as_interstitial_fcc_a1 (SGTE convention: (M)(%s,VA))" % nonmetal[0], ["FCC_A1"]
    return None

def tdb_compositions(path, elements):
    """{phase: (x_min, x_max)} of the system's second element, from PHASE site ratios and CONST(ITUENT) lines; the range
    spans the phase's end members (a single-species phase gives one point). VA is ignored."""
    import itertools
    txt = open(path, encoding="latin-1").read()
    ratios = {m.group(1).upper(): [float(x) for x in m.group(3).split()] for m in
              re.finditer(r"PHASE\s+([A-Z0-9_#]+)(:\S+)?\s+\S+\s+\d+\s+([\d.\s]+)!", txt, re.I)}
    out = {}
    for m in re.finditer(r"CONST\w*\s+([A-Z0-9_#]+)(:\S+)?\s*:(.*?)!", txt, re.I | re.S):
        ph = m.group(1).upper(); subl = [x.strip() for x in m.group(3).split(":") if x.strip()]
        r = ratios.get(ph)
        if not r or len(r) != len(subl):
            continue
        sets = []
        for sp in subl:
            sp_ = [x.strip().upper().rstrip("%") for x in re.split(r"[,\s]+", sp) if x.strip()]
            sets.append([x for x in sp_ if x in elements] + (["VA"] if "VA" in sp_ else []))
        if any(not st for st in sets) or len(list(itertools.islice(itertools.product(*sets), 5000))) >= 5000:
            continue
        xs = []
        for combo in itertools.product(*sets):
            tot = {e: 0.0 for e in elements}
            for ratio, sp in zip(r, combo):
                if sp != "VA":
                    tot[sp] += ratio
            if sum(tot.values()) > 0:
                xs.append(tot[elements[1]] / sum(tot.values()))
        if xs:
            out[ph] = (min(xs), max(xs))
    return out

def formula_dict(s):
    s = re.sub(r"[()\s'’*]", "", s); s = re.sub(r"^[%s]+" % GREEK_CHARS, "", s)
    d = {}
    for e, n in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", s):
        d[e.upper()] = d.get(e.upper(), 0) + float(n or 1)
    t = sum(d.values()) or 1
    return {k: round(v / t, 4) for k, v in d.items()}

def norm_formula(s):
    s = re.sub(r"[()\s'’*]", "", s); s = re.sub(r"^[%s]+" % GREEK_CHARS, "", s)
    parts = re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", s)
    return "".join(f"{e.upper()}{n if n not in ('', '1') else ''}" for e, n in parts), [e for e, _ in parts]

def norm_sb(sb):
    return re.sub(r"[_\s]", "", sb or "").upper()

def map_phase(name, info, tdb, comps=None, elements=None):
    """name: printed phase; info: {'pearson','strukturbericht','prototype'} (any may be None); tdb: from tdb_phases."""
    n = (name or "").strip(); pear = (info or {}).get("pearson") or (info or {}).get("pearson_symbol"); sb = norm_sb((info or {}).get("strukturbericht"))
    if "?" in n:
        return None, "label uncertain ('?' printed): never mapped", []
    if re.fullmatch(r"L\d?|L#\d", n):
        for c in LIQUID_NAMES:
            if c in tdb: return c, "liquid", [c]
        return None, "no liquid phase in TDB", []
    if n == "G":
        return ("GAS", "gas", ["GAS"]) if "GAS" in tdb else (None, "no gas phase in TDB", [])
    solution = n.startswith("(") and n.endswith(")")
    if solution:
        named = {e.upper() for e in re.findall(r"[A-Z][a-z]?", re.sub("[%s]" % GREEK_CHARS, "", n))}
        if pear:
            c = [p for p, v in tdb.items() if pearson_equal(v.get("pearson"), pear) and not liquid_like(p, v)
                 and (not v.get("constituents") or named <= {x for sl in v["constituents"] for x in sl})]
            std = [p for p in c if p in STD.values()]
            if len(std) == 1: return std[0], "solution_by_pearson", c
            if len(c) == 1: return c[0], "solution_by_pearson", c
        if sb in STD and STD[sb] in tdb and (not tdb[STD[sb]].get("constituents") or
                                             named <= {x for sl in tdb[STD[sb]]["constituents"] for x in sl}):
            return STD[sb], "solution_by_strukturbericht", [STD[sb]]
        pref = STD_COMPOUND.get(pear or "") or STD_COMPOUND.get(sb or "")
        if pref == "DHCP" and "DHCP" in tdb and (not tdb["DHCP"].get("constituents") or
                                                  named <= {x for sl in tdb["DHCP"]["constituents"] for x in sl}):
            return "DHCP", "solution_by_standard_structure_name (hP4 / A3' -> DHCP)", ["DHCP"]
        return None, "solution unresolved", [p for p in STD.values() if p in tdb]
    def standard_names(keys):
        for key in keys:   # standard compound names (Laves, sigma, mu ...)
            pref = STD_COMPOUND.get(key or "")
            if pref and key == pear and sb and STD_COMPOUND.get(sb) != pref and STD.get(sb) != pref:
                continue          # a Pearson-only key never overrides a printed Strukturbericht of another structure
            if pref:
                hit = [p for p in tdb if (p == pref or p.startswith(pref)) and not liquid_like(p, tdb[p])]
                if len(hit) == 1:
                    return hit[0], "standard_name_by_structure", hit
        return None
    if sb:   # a printed Strukturbericht names the structure directly
        hit = standard_names([sb])
        if hit:
            return hit
    if re.fullmatch(r"[PFICRA][a-z]?\d*", n):
        return None, "not a phase (space-group fragment from the table)", []
    f, els = norm_formula(n)
    fd = formula_dict(n)
    c = [p for p in tdb if fd and len(fd) >= 2 and formula_dict(p.split("_")[0]) == fd] or \
        [p for p in tdb if re.sub(r"_[A-Z0-9]+$", "", p) == f or p == f or p.split("_")[0] == f]
    c = [p for p in c if not liquid_like(p, tdb[p])]      # a printed solid is never the liquid / glass phase
    if c and (pear or sb):
        # polymorphs: a printed Pearson symbol / Strukturbericht that contradicts the TDB phase blocks the match
        def contradicts(p):
            toks = [t for t in name_tokens(p)[1:] if re.fullmatch(r"[A-L]\d{1,2}[A-Z]?|[A-L]\d\d?P?", t)]
            if sb and toks:   # the TDB name carries a Strukturbericht: it decides (a printed Pearson may be misread)
                return sb not in toks and sb.rstrip("'") not in toks
            tp = tdb[p].get("pearson")
            return bool(pear and tp and pearson_equal(tp, pear) is False)
        ok = [p for p in c if not contradicts(p)]
        if not ok:
            return None, ("compound_by_formula blocked: printed structure %s contradicts %s (polymorph not modelled)"
                          % ("/".join(x for x in (pear, sb) if x), ", ".join("%s %s" % (p, tdb[p].get("pearson") or "")
                                                                         for p in c))), c
        c = ok
    if len(c) == 1:
        agree = bool(pear and pearson_equal(tdb[c[0]].get("pearson"), pear))
        return c[0], "compound_by_formula" + ("+pearson" if agree else ""), c
    if len(c) > 1:
        g = next((GREEK[ch] for ch in n[:2] if ch in GREEK), None)
        pick = [p for p in c if g and g in name_tokens(p)] if g else \
            [p for p in c if not any(x in name_tokens(p) for x in GREEK.values())]
        if len(pick) == 1:
            return pick[0], "compound_variant_by_greek_prefix" if g else "compound_plain_variant", c
        return None, "compound: several variants", c
    if is_rocksalt(info):
        hit = map_rocksalt(n, info, tdb, comps, fd, elements)
        if hit:
            return hit
    if not (fd and len(fd) >= 2):   # Pearson-only standard names only for a phase not named by a formula
        hit = standard_names([pear])
        if hit:
            return hit
    # a phase named by letter (eta, sigma, gamma0 ...): Pearson symbol + composition range against the TDB phases
    rng = (info or {}).get("composition_at_pct")
    if comps is not None and rng and (len(fd) < 2):
        lo, hi = rng[0] / 100 - 0.03, rng[1] / 100 + 0.03
        cand = [p for p, (x0, x1) in comps.items() if x0 <= hi and x1 >= lo and (x1 - x0) < 0.5   # overlapping, not a full solution
                and p not in STD.values() and p not in ("LIQUID", "GAS") and not liquid_like(p, tdb.get(p))]
        if pear:
            cand = [p for p in cand if pearson_equal(tdb.get(p, {}).get("pearson"), pear)] or \
                   ([] if any(tdb.get(p, {}).get("pearson") for p in cand) else cand)
        if len(cand) == 1:
            agree = bool(pear and pearson_equal(tdb.get(cand[0], {}).get("pearson"), pear))
            return cand[0], "letter_name_by_composition" + ("+pearson" if agree else ""), cand
        if cand:
            return None, "letter name: several candidates", cand
    return None, "compound not in TDB", []

def sgte_files(system):
    a, b = system.split("-")[:2]; cap = lambda s: s[0] + s[1:].lower()
    return sorted(glob.glob(f"{REPO}/databases/sgte-binary/{cap(a)}{cap(b)}-*.tdb") + glob.glob(f"{REPO}/databases/sgte-binary/{cap(b)}{cap(a)}-*.tdb"))

if __name__ == "__main__":
    import collections
    hb = json.load(open(f"{REPO}/build/sgte-corpus/db-v2.1/okamoto/handbook_labels.json"))["systems"]
    stats = collections.Counter(); per_basis = collections.Counter(); rows = []; nsys = 0
    axis = {}
    RUN = os.environ.get("PM_RUN", "run-20261021c-full")
    # 2026-10-21: phase tables from the tracer's own reading of each diagram (fields.json), not the older
    # handbook_labels tables (junk rows / missing phases in 813 of 2,273 diagrams)
    own = {}
    for cf_ in glob.glob(f"{REPO}/build/okamoto-annotation/{RUN}/systems/*/calibrated.json"):
        d_ = os.path.dirname(cf_); s_, pg_ = os.path.basename(d_).rsplit("-p", 1)
        axis[s_] = (json.load(open(cf_)).get("component") or "").upper()
        if os.path.exists(f"{d_}/fields.json"):
            fj = json.load(open(f"{d_}/fields.json")); t_ = {}
            for ph, rng in (fj.get("phase_table") or {}).items():
                info = dict((fj.get("phase_structures") or {}).get(ph) or {}); info["composition_at_pct"] = list(rng); t_[ph] = info
            own.setdefault(s_, {})[pg_] = t_
    for system, rec in sorted(hb.items()):
        files = sgte_files(system)
        if not files: stats["system without SGTE TDB"] += 1; continue
        nsys += 1
        tabs = own.get(system) or rec.get("phase_tables") or {}
        phases = {}
        for pg, t in tabs.items():
            phases.update(t)
        A, B = system.split("-")[:2]
        if axis.get(system) == A.upper():   # composition ranges are in at.% of the axis element
            A, B = B, A
        for f in files:
            tdb = tdb_phases(f); comps = tdb_compositions(f, [A.upper(), B.upper()])
            for ph, info in phases.items():
                m, basis, cands = map_phase(ph, info, tdb, comps, [A.upper(), B.upper()])
                rows.append({"system": system, "tdb": os.path.basename(f), "handbook_phase": ph, "tdb_phase": m, "basis": basis, "candidates": cands})
                stats["mapped" if m else "open"] += 1; per_basis[basis] += 1
    os.makedirs(f"{REPO}/build/okamoto-annotation/phase-map", exist_ok=True)
    json.dump(rows, open(f"{REPO}/build/okamoto-annotation/phase-map/handbook_to_sgte.json", "w"), indent=0, ensure_ascii=False)
    print("systems with SGTE TDBs:", nsys, "| pairs (handbook phase x TDB file):", len(rows), dict(stats))
    for k, v in per_basis.most_common(): print(f"  {v:6d}  {k}")
