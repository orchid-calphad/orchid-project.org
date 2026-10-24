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
"""
import glob, json, os, re, sys
REPO = "/Users/raymundo.arroyave/projects/open-tdb"
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
GREEK = {"α": "ALPHA", "β": "BETA", "γ": "GAMMA", "δ": "DELTA", "ε": "EPSILON", "ζ": "ZETA", "η": "ETA"}

def tdb_phases(path):
    """{phase: {'pearson':..,'prototype':..}} from PHASE lines and the comment just before each."""
    out, last = {}, None
    for line in open(path, encoding="latin-1"):
        m = re.match(r"\s*\$\s*(?:Prototype\s+)?(\S.*?)\s*\((\w+\*?),\s*([^)]*)\)", line)
        if m:
            last = {"prototype": m.group(1), "pearson": m.group(2), "space_group": m.group(3)}
            continue
        m = re.match(r"\s*PHASE\s+([A-Z0-9_#]+)(?::\S+)?\s", line, re.I)
        if m:
            out[m.group(1).upper()] = dict(last or {}); last = None
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
    s = re.sub(r"[()\s'’*]", "", s); s = re.sub(r"^[αβγδεζηθ]+", "", s)
    d = {}
    for e, n in re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", s):
        d[e.upper()] = d.get(e.upper(), 0) + float(n or 1)
    t = sum(d.values()) or 1
    return {k: round(v / t, 4) for k, v in d.items()}

def norm_formula(s):
    s = re.sub(r"[()\s'’*]", "", s); s = re.sub(r"^[αβγδεζηθ]+", "", s)
    parts = re.findall(r"([A-Z][a-z]?)(\d*\.?\d*)", s)
    return "".join(f"{e.upper()}{n if n not in ('', '1') else ''}" for e, n in parts), [e for e, _ in parts]

def norm_sb(sb):
    return re.sub(r"[_\s]", "", sb or "").upper()

def map_phase(name, info, tdb, comps=None, elements=None):
    """name: printed phase; info: {'pearson','strukturbericht','prototype'} (any may be None); tdb: from tdb_phases."""
    n = (name or "").strip(); pear = (info or {}).get("pearson") or (info or {}).get("pearson_symbol"); sb = norm_sb((info or {}).get("strukturbericht"))
    if re.fullmatch(r"L\d?|L#\d", n):
        for c in ("LIQUID", "IONIC_LIQ"):
            if c in tdb: return c, "liquid", [c]
        return None, "no liquid phase in TDB", []
    if n == "G":
        return ("GAS", "gas", ["GAS"]) if "GAS" in tdb else (None, "no gas phase in TDB", [])
    solution = n.startswith("(") and n.endswith(")")
    if solution:
        if pear:
            c = [p for p, v in tdb.items() if (v.get("pearson") or "").rstrip("*") == pear.rstrip("*")]
            std = [p for p in c if p in STD.values()]
            if len(std) == 1: return std[0], "solution_by_pearson", c
            if len(c) == 1: return c[0], "solution_by_pearson", c
        if sb in STD and STD[sb] in tdb:
            return STD[sb], "solution_by_strukturbericht", [STD[sb]]
        return None, "solution unresolved", [p for p in STD.values() if p in tdb]
    for key in (norm_sb((info or {}).get("strukturbericht")), pear):   # standard compound names (Laves, sigma, mu ...)
        pref = STD_COMPOUND.get(key or "")
        if pref and key == pear and sb and STD_COMPOUND.get(sb) != pref and STD.get(sb) != pref:
            continue          # a Pearson-only key never overrides a printed Strukturbericht of another structure
        if pref:
            hit = [p for p in tdb if p == pref or p.startswith(pref)]
            if len(hit) == 1:
                return hit[0], "standard_name_by_structure", hit
    if re.fullmatch(r"[PFICRA][a-z]?\d*", n):
        return None, "not a phase (space-group fragment from the table)", []
    f, els = norm_formula(n)
    fd = formula_dict(n)
    c = [p for p in tdb if fd and len(fd) >= 2 and formula_dict(p.split("_")[0]) == fd] or \
        [p for p in tdb if re.sub(r"_[A-Z0-9]+$", "", p) == f or p == f or p.split("_")[0] == f]
    if pear:
        agree = [p for p in c if not tdb[p].get("pearson") or tdb[p]["pearson"].rstrip("*") == pear.rstrip("*")]
        if agree: c = agree
    if len(c) == 1:
        return c[0], "compound_by_formula" + ("+pearson" if pear and tdb[c[0]].get("pearson") else ""), c
    if len(c) > 1:
        g = next((GREEK[ch] for ch in n[:2] if ch in GREEK), None)
        pick = [p for p in c if g and g in p] if g else [p for p in c if not any(x in p for x in GREEK.values())]
        if len(pick) == 1:
            return pick[0], "compound_variant_by_greek_prefix" if g else "compound_plain_variant", c
        return None, "compound: several variants", c
    if is_rocksalt(info):
        hit = map_rocksalt(n, info, tdb, comps, fd, elements)
        if hit:
            return hit
    # a phase named by letter (eta, sigma, gamma0 ...): Pearson symbol + composition range against the TDB phases
    rng = (info or {}).get("composition_at_pct")
    if comps is not None and rng and (len(fd) < 2):
        lo, hi = rng[0] / 100 - 0.03, rng[1] / 100 + 0.03
        cand = [p for p, (x0, x1) in comps.items() if x0 <= hi and x1 >= lo and (x1 - x0) < 0.5   # overlapping, not a full solution
                and p not in STD.values() and p not in ("LIQUID", "GAS", "IONIC_LIQ")]
        if pear:
            cand = [p for p in cand if (tdb.get(p, {}).get("pearson") or "").rstrip("*") == pear.rstrip("*")] or \
                   ([] if any(tdb.get(p, {}).get("pearson") for p in cand) else cand)
        if len(cand) == 1:
            return cand[0], "letter_name_by_composition" + ("+pearson" if pear and tdb.get(cand[0], {}).get("pearson") else ""), cand
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
