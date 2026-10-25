"""ESPEI bundle: mechanical encoding fixes (E1-E3) and optional, tagged phase-constraint datasets.

Used by tools/build_espei_bundle.py; standard library only; deterministic.

E1  `VA` in `components` when the dataset's target TDB (the assessment the folder is named after) has a phase model
    with a vacancy sublattice; where the folder is not a registry TDB: when any registry TDB of the system has one,
    or the system is all-metallic. ESPEI 0.9.1 drops phases whose constituents are not all in the components of a
    dataset (e.g. FCC_A1/BCC_A2/HCP_A3 with their VA interstitial sublattice, and with them activity reference states).
E2  a ZPF vertex of a phase that is stoichiometric (line compound) in the target TDB is written with composition
    null. ESPEI 0.9.1 + pycalphad 0.11.2 return NaN for the isolated-phase energy of a line compound at a fixed
    composition (any NaN driving force makes the whole ZPF likelihood -inf); with null, ESPEI evaluates the compound
    at its own fixed stoichiometry, which is the formula composition the vertex carried. Where this leaves a region
    with no known composition and no hyperplane, a __HYPERPLANE__ is added at the mean of the (distinct) compound
    compositions, the bundle's existing rule for compound-only tie-lines.
E3  identical ZPF regions (same system, same TDB phase names, |dT| <= 0.05 K, every composition known in both within
    0.001 and the same set of known compositions) that appear in several datasets (several assessment folders or
    several papers quoting the same invariant) are written once; the kept region lists every source (observation
    ids, papers, DOIs) in `orchid.observations[i].also_reported_as`.

Phase constraints (optional; never mixed with measured data; folder and tag say what they are):
  derived/invariant-brackets/   tag derived_bracket            (from experimental invariant temperatures)
  derived/liquid-stability/     tag derived_liquid_stability   (from complete liquid miscibility in the sources)
  assessed/                     tag assessed_phase_diagram     (assessed-phase-assemblages in each TDB's phase names)
"""
from __future__ import annotations

import glob
import gzip
import json
import os
import re
from collections import Counter, defaultdict

T_TOL_DUP = 0.05
X_TOL_DUP = 1e-3
DEFAULT_T_SIGMA = 5.0         # K, invariant temperature uncertainty when none is printed
MIN_OFFSET = 10.0             # K
HYPER_DELTA = 0.01            # mole fraction: hyperplane estimated at x0 -/+ delta
LIQ_MARGIN = 50.0             # K above the highest liquidus evidence
LIQ_X = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
SKIP_TYPES = {"critical point", "magnetic", "order-disorder"}
ARROWS = ["<=>", "<->", "⇌", "↔", "⇄", "⇆", "->", "→", "="]
NONMETALLIC = {"H", "B", "C", "N", "O", "F", "SI", "P", "S", "SE", "TE", "AS", "SB", "CL", "BR", "I", "AT", "GE",
               "HE", "NE", "AR", "KR", "XE"}

# Melting and normal boiling points (K) at 1 atm of the elements (SGTE unary data for melting points, A.T. Dinsdale,
# CALPHAD 15 (1991) 317; boiling points CRC Handbook of Chemistry and Physics). Only elements that are condensed liquids
# at 1 atm; elements that sublime or are gases at 1 atm are absent (no liquid-stability points for their systems).
ELEMENT_TM_TB = {
    "AG": (1234.93, 2435), "AL": (933.47, 2792), "AU": (1337.33, 3129), "B": (2348, 4273), "BA": (1000, 2118),
    "BE": (1560, 2742), "BI": (544.55, 1837), "CA": (1115, 1757), "CD": (594.22, 1040), "CE": (1071, 3716),
    "CO": (1768, 3200), "CR": (2180, 2944), "CS": (301.6, 944), "CU": (1357.77, 2835), "DY": (1685, 2840),
    "ER": (1802, 3141), "EU": (1095, 1802), "FE": (1811, 3134), "GA": (302.91, 2477), "GD": (1586, 3546),
    "GE": (1211.4, 3106), "HF": (2506, 4876), "HO": (1747, 2973), "IN": (429.75, 2345), "IR": (2719, 4701),
    "K": (336.53, 1032), "LA": (1193, 3737), "LI": (453.65, 1603), "LU": (1925, 3675), "MG": (923, 1363),
    "MN": (1519, 2334), "MO": (2896, 4912), "NA": (370.94, 1156), "NB": (2750, 5017), "ND": (1297, 3347),
    "NI": (1728, 3186), "OS": (3306, 5285), "PB": (600.61, 2022), "PD": (1828.05, 3236), "PR": (1208, 3403),
    "PT": (2041.4, 4098), "RB": (312.46, 961), "RE": (3459, 5869), "RH": (2237, 3968), "RU": (2607, 4423),
    "SB": (903.78, 1860), "SC": (1814, 3109), "SI": (1687, 3538), "SM": (1345, 2067), "SN": (505.08, 2875),
    "SR": (1050, 1655), "TA": (3290, 5731), "TB": (1629, 3396), "TE": (722.66, 1261), "TH": (2023, 5061),
    "TI": (1941, 3560), "TL": (577, 1746), "TM": (1818, 2223), "U": (1408, 4404), "V": (2183, 3680),
    "W": (3695, 5828), "Y": (1799, 3609), "YB": (1097, 1469), "ZN": (692.68, 1180), "ZR": (2128, 4682),
    "SE": (494, 958), "PU": (912.5, 3501), "NP": (917, 4273),
}

BRACKET_RULE = (
    "For each binary invariant reaction (eutectic, peritectic, monotectic, eutectoid, peritectoid, metatectic, "
    "syntectic, congruent melting or transformation, polymorphic) with printed temperature T and printed temperature "
    "uncertainty u (u = sigma_effective of the value's uncertainty budget: the printed sigma or the method floor, "
    "whichever is larger; else the printed sigma; else %g K, documented default), two ZPF points are written at the invariant composition x0: "
    "at T + max(2u, %g K) the phase set stable above the reaction and at T - max(2u, %g K) the set stable below. The "
    "reaction is read as printed with the cooling convention (left side = high temperature; a reaction written on "
    "heating, i.e. a single liquid formed from solids, is turned round). x0 is the composition of the phase alone on "
    "its side (eutectic: the liquid; peritectic: the product; congruent: the common composition); the other side's "
    "phases are written with their printed compositions (null where not printed or for line compounds, E2). The "
    "target hyperplane is estimated at x0 - d and x0 + d (d = min(%g, half the distance to the nearest printed "
    "composition or to a pure element)) so that a line compound at x0 cannot make the hyperplane undefined (ESPEI "
    "0.9.1 returns zero driving force when the equilibrium at the hyperplane composition is a single phase without "
    "internal degrees of freedom). A point is not written when another invariant of the system that shares a "
    "non-liquid phase (or a liquid within 0.05 in composition of x0) lies between T and the bracket temperature. "
    "These points encode the semantics of the printed reaction (which phases are stable on which side); they are "
    "DERIVED, not measured." % (DEFAULT_T_SIGMA, MIN_OFFSET, MIN_OFFSET, HYPER_DELTA))
LIQUID_RULE = (
    "Single-phase LIQUID points at x = 0.1 ... 0.9 at T = T_liq_max + %g K, where T_liq_max is the highest of the "
    "melting points of the two elements (SGTE unary data), every accepted printed temperature of an experimental "
    "literature value of the system that involves the liquid (invariants, liquidus and ZPF points) and every assessed "
    "(handbook) two-phase field point that contains the liquid. Written only where the evidence shows complete liquid "
    "miscibility: the system has an assessed (handbook) diagram with a single-phase liquid field and no second liquid "
    "(no L#2 field), no literature value of the system has two liquids (monotectic, syntectic or two liquid labels), "
    "both elements are condensed liquids at 1 atm (metals and metalloids) and T stays below both normal boiling points. "
    "Purpose: forbid spurious liquid miscibility gaps (and solids stable above the liquidus) in a fit. DERIVED, not "
    "measured." % LIQ_MARGIN)


# ------------------------------------------------------------------------------------------- registry --
class Registry:
    def __init__(self, data):
        self.models = defaultdict(dict)        # tdb stem -> phase_name -> model (with 'phase_id')
        self.tdb_va = {}                       # tdb stem -> any phase model has VA
        self.sys_tdbs = defaultdict(set)       # system -> {tdb stem}
        self.by_pref = defaultdict(dict)       # system -> preferred_name -> phase entry
        self.phase_tdb_names = defaultdict(lambda: defaultdict(list))   # phase id -> stem -> [phase_name]
        for f in sorted(glob.glob(os.path.join(data, "phases", "systems", "*.json"))):
            d = json.load(open(f, encoding="utf-8"))
            sysn = d["system"]
            for p in d["phases"]:
                self.by_pref[sysn][p["preferred_name"]] = p
                for m in p.get("models") or []:
                    stem = m["tdb"][:-4] if m["tdb"].lower().endswith(".tdb") else m["tdb"]
                    mm = dict(m, phase_id=p["id"])
                    self.models[stem][m["phase_name"]] = mm
                    va = any("VA" in [str(c).upper().rstrip("%") for c in s.get("constituents") or []]
                             for s in m.get("sublattices") or [])
                    self.tdb_va[stem] = self.tdb_va.get(stem, False) or va
                    self.sys_tdbs[sysn].add(stem)
                    self.phase_tdb_names[p["id"]][stem].append(m["phase_name"])

    def stoichiometric(self, stem, phase):
        m = self.models.get(stem, {}).get(phase)
        return bool(m and m.get("stoichiometric"))

    def formula_x(self, stem, phase, el):
        m = self.models.get(stem, {}).get(phase)
        if m and m.get("stoichiometric") and m.get("composition"):
            c = {k.upper(): v for k, v in m["composition"].items()}
            return c.get(el.upper())
        return None

    def phase_id(self, stem, phase):
        m = self.models.get(stem, {}).get(phase)
        return m["phase_id"] if m else None

    def va_needed(self, system, stem):
        """(bool, basis)"""
        if stem in self.tdb_va:
            return self.tdb_va[stem], "target_tdb"
        if any(self.tdb_va.get(t) for t in self.sys_tdbs.get(system, ())):
            return True, "system_tdbs"
        els = system.split("-")
        if not (set(els) & NONMETALLIC):
            return True, "metallic_system"
        return False, "no_va_model"


def stem_of(rel):
    """dataset path <folder>/<SYS>/<stem>.<file>.json -> stem"""
    base = os.path.basename(rel)
    return base.rsplit(".", 2)[0]


# ------------------------------------------------------------------------------------------- E1 / E2 --
def apply_e1(d, system, stem, R, counts):
    need, basis = R.va_needed(system, stem)
    if need and "VA" not in [c.upper() for c in d["components"]]:
        d["components"] = sorted(set(d["components"]) | {"VA"})
        counts["E1_va_added"] += 1
        counts["E1_basis_" + basis] += 1
        return True
    return False


def apply_e2(d, stem, R, counts):
    """null the composition of line-compound vertices; returns number of vertices changed."""
    if d.get("output") != "ZPF":
        return 0
    n = 0
    hyper_added = 0
    for reg in d["values"]:
        formula = []
        for vx in reg:
            if vx[0] == "__HYPERPLANE__":
                continue
            if R.stoichiometric(stem, vx[0]) and any(x is not None for x in vx[2]):
                formula.append(tuple(vx[2]))
                vx[2] = [None] * len(vx[2])
                n += 1
        has_known = any(vx[0] != "__HYPERPLANE__" and any(x is not None for x in vx[2]) for vx in reg)
        has_hyper = any(vx[0] == "__HYPERPLANE__" for vx in reg)
        distinct = sorted(set(formula))
        if formula and not has_known and not has_hyper and len(distinct) >= 2 and len(distinct[0]) == 1:
            comps = [vx[1] for vx in reg if vx[0] != "__HYPERPLANE__"][0]
            xm = sum(x[0] for x in distinct) / len(distinct)
            reg.append(["__HYPERPLANE__", list(comps), [round(xm, 10)]])
            hyper_added += 1
    if n:
        d["orchid_zpf_encoding"] = {
            "line_compound_vertices_null": n, "hyperplanes_added": hyper_added,
            "rule": "E2: vertices of phases that are stoichiometric in the target TDB are written with composition "
                    "null (ESPEI evaluates a line compound at its own formula composition; a fixed composition "
                    "gives NaN in ESPEI 0.9.1 / pycalphad 0.11.2); the printed composition stays in the literature "
                    "release (orchid.observations ids)"}
        counts["E2_vertices_nulled"] += n
        counts["E2_datasets"] += 1
        counts["E2_hyperplanes_added"] += hyper_added
    return n


def zpf_inert(region, stem, R, e2=True):
    """True when a region has (after E2) no composition and no hyperplane: ESPEI 0.9.1 driving force 0. With e2,
    stoichiometric vertices count as null and the E2 hyperplane (two distinct binary line-compound formulas) counts."""
    formula, known = [], False
    for vx in region:
        if vx[0] == "__HYPERPLANE__":
            return False
        if any(x is not None for x in vx[2]):
            if e2 and R.stoichiometric(stem, vx[0]):
                formula.append(tuple(vx[2]))
            else:
                known = True
    if known:
        return False
    distinct = sorted(set(formula))
    return not (e2 and len(distinct) >= 2 and len(distinct[0]) == 1)


# ------------------------------------------------------------------------------------------- E3 --
def _region_sig(reg):
    out = {}
    for vx in reg:
        if vx[0] == "__HYPERPLANE__":
            continue
        out.setdefault(vx[0], []).append(tuple(vx[2]))
    return out


def _same_region(a, b):
    sa, sb = _region_sig(a), _region_sig(b)
    if set(sa) != set(sb) or any(len(sa[k]) != len(sb[k]) for k in sa):
        return False
    for k in sa:
        for xa, xb in zip(sorted(sa[k], key=str), sorted(sb[k], key=str)):
            if len(xa) != len(xb):
                return False
            for p, q in zip(xa, xb):
                if (p is None) != (q is None):
                    return False
                if p is not None and abs(p - q) > X_TOL_DUP:
                    return False
    return True


def dedup_zpf(C, folders=("experimental", "invariants")):
    """E3 across datasets of the same evidence class (experimental ZPF + invariants). Mutates C."""
    groups = defaultdict(list)    # (system, frozenset phases) -> [(rel, i)]
    for rel in sorted(C["datasets"]):
        if rel.split("/")[0] not in folders:
            continue
        d = C["datasets"][rel]
        if d["output"] != "ZPF":
            continue
        Ts = d["conditions"]["T"] if isinstance(d["conditions"]["T"], list) else [d["conditions"]["T"]] * len(d["values"])
        sysn = rel.split("/")[1]
        for i, reg in enumerate(d["values"]):
            ph = frozenset(vx[0] for vx in reg if vx[0] != "__HYPERPLANE__")
            groups[(sysn, ph)].append((rel, i, Ts[i]))
    order = {rel: k for k, rel in enumerate(sorted(C["datasets"], key=lambda r: (0 if r.startswith("experimental/") else 1, r)))}
    remove = defaultdict(set)
    merged = []
    for key, items in groups.items():
        items.sort(key=lambda t: (order[t[0]], t[1]))
        kept = []
        for rel, i, T in items:
            reg = C["datasets"][rel]["values"][i]
            hit = None
            for (r2, i2, T2) in kept:
                if abs(T - T2) <= T_TOL_DUP + 1e-9 and _same_region(reg, C["datasets"][r2]["values"][i2]):
                    hit = (r2, i2)
                    break
            if hit is None:
                kept.append((rel, i, T))
                continue
            remove[rel].add(i)
            merged.append({"kept": {"path": hit[0], "index": hit[1]}, "removed": {"path": rel, "index": i}})
    # attach sources to the kept observation, then cut
    for m in merged:
        kd = C["datasets"][m["kept"]["path"]]
        rd = C["datasets"][m["removed"]["path"]]
        ko = kd["orchid"]["observations"][m["kept"]["index"]]
        ro = rd["orchid"]["observations"][m["removed"]["index"]]
        src = {k: v for k, v in (("id", ro.get("id")), ("paper", ro.get("paper")), ("measured_in", ro.get("measured_in")),
                                  ("dataset", m["removed"]["path"])) if v}
        ko.setdefault("also_reported_as", []).append(src)
        m["kept_id"], m["removed_id"] = ko.get("id"), ro.get("id")
        for doi in rd["orchid"].get("dois") or []:
            if doi not in kd["orchid"]["dois"]:
                kd["orchid"]["dois"].append(doi)
        kd["orchid"]["dois"].sort()
        for pk in rd["orchid"].get("papers") or []:
            if pk not in kd["orchid"]["papers"]:
                kd["orchid"]["papers"].append(pk)
        kd["orchid"]["e3_merged"] = kd["orchid"].get("e3_merged", 0) + 1
    changed = set()
    for rel, idx in remove.items():
        d = C["datasets"][rel]
        keep = [i for i in range(len(d["values"])) if i not in idx]
        Ts = d["conditions"]["T"] if isinstance(d["conditions"]["T"], list) else [d["conditions"]["T"]] * len(d["values"])
        for i in sorted(idx):
            C["drops"].append({"source": C["rows"][rel]["source_dataset"].replace("db-v2.1/", ""), "index": i,
                               "reason": "duplicate_identical_zpf",
                               "orchid_id": d["orchid"]["observations"][i].get("id"),
                               "kept_in": next(m["kept"] for m in merged if m["removed"] == {"path": rel, "index": i})})
        if not keep:
            C["dropped_datasets"].append({"source": C["rows"][rel]["source_dataset"], "reason": ["duplicate_identical_zpf"],
                                          "values": len(d["values"])})
            C["datasets"].pop(rel)
            C["rows"].pop(rel)
            C["_recs"].pop(rel, None)
            continue
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
        changed.add(rel)
    for m in merged:
        changed.add(m["kept"]["path"])
    return merged, changed


# ------------------------------------------------------------------------------------------- brackets --
def _norm_label(s):
    return re.sub(r"\s+", "", (s or "")).lower()


def parse_reaction(printed):
    p = (printed or "").strip()
    for a in ARROWS:
        if a in p:
            parts = p.split(a)
            if len(parts) != 2:
                return None
            lhs = [t.strip() for t in parts[0].split("+") if t.strip()]
            rhs = [t.strip() for t in parts[1].split("+") if t.strip()]
            if lhs and rhs:
                return lhs, rhs
            return None
    return None


def _t_sigma(rec):
    """u of the bracket rule: sigma_effective of the invariant temperature's published uncertainty budget (method
    floor, never below the printed sigma; PI decision 2026-10-24), else the printed sigma, else DEFAULT_T_SIGMA."""
    q = ((rec.get("uncertainty_budget") or {}).get("quantities") or {}).get("value") or {}
    if isinstance(q.get("sigma_effective"), (int, float)) and q["sigma_effective"] > 0 and q.get("unit") == "K":
        m = q.get("sigma_method") or {}
        return float(q["sigma_effective"]), "sigma_effective of the uncertainty budget (method class %s, %s)" % (
            m.get("class"), m.get("basis"))
    tu = rec.get("temperature_uncertainty") or {}
    if tu.get("basis") in ("printed", "printed_range") and isinstance(tu.get("value"), (int, float)) and tu["value"] > 0:
        return float(tu["value"]), tu["basis"]
    return DEFAULT_T_SIGMA, "default (no printed temperature uncertainty)"


def build_brackets(C, R, folders=("experimental", "invariants")):
    """Returns (datasets {rel: d}, skipped Counter, per-system invariant list)."""
    inv = defaultdict(list)    # system -> [(T, ids frozenset, {id: x or None}, rel, i)]
    cand = []
    for rel in sorted(C["datasets"]):
        if rel.split("/")[0] not in folders:
            continue
        d = C["datasets"][rel]
        if d["output"] != "ZPF":
            continue
        sysn = rel.split("/")[1]
        if len(sysn.split("-")) != 2:
            continue
        stem = stem_of(rel)
        Ts = d["conditions"]["T"] if isinstance(d["conditions"]["T"], list) else [d["conditions"]["T"]] * len(d["values"])
        recs = C["_recs"][rel]
        for i, (reg, rec) in enumerate(zip(d["values"], recs)):
            if rec.get("observation_type") != "invariant_reaction":
                continue
            comp = {}
            for vx in reg:
                if vx[0] == "__HYPERPLANE__":
                    continue
                x = vx[2][0] if vx[2] else None
                if x is None:
                    x = R.formula_x(stem, vx[0], vx[1][0])
                comp[vx[0]] = x
            ids = {}
            for ph, x in comp.items():
                ids[R.phase_id(stem, ph) or (stem + ":" + ph)] = x
            inv[sysn].append((Ts[i], frozenset(ids), ids))
            cand.append((rel, i, Ts[i], reg, rec, comp))
    out, skipped, made = {}, Counter(), Counter()
    for rel, i, T, reg, rec, comp in cand:
        d = C["datasets"][rel]
        sysn = rel.split("/")[1]
        stem = stem_of(rel)
        rtype = (rec.get("reaction") or {}).get("type")
        if rtype in SKIP_TYPES:
            skipped["reaction_type_" + rtype.replace(" ", "_")] += 1
            continue
        pr = parse_reaction((rec.get("reaction") or {}).get("printed"))
        if not pr:
            skipped["reaction_not_parsed"] += 1
            continue
        lab2ph = {}
        for p in rec.get("phases") or []:
            lab2ph[_norm_label(p.get("label"))] = p.get("model_phase")
        lhs, rhs = pr
        try:
            L = [lab2ph[_norm_label(t)] for t in lhs]
            Rh = [lab2ph[_norm_label(t)] for t in rhs]
        except KeyError:
            skipped["reaction_label_not_matched"] += 1
            continue
        if any(p is None for p in L + Rh) or len(set(L + Rh)) != len(L + Rh):
            skipped["phase_unmapped_or_repeated"] += 1
            continue
        if set(L + Rh) != set(comp):
            skipped["reaction_phases_differ_from_dataset"] += 1
            continue
        # cooling convention; a reaction written on heating (solids -> one liquid) is turned round
        if len(Rh) == 1 and Rh[0].upper().startswith("LIQ") and not any(p.upper().startswith("LIQ") for p in L):
            L, Rh = Rh, L
        if not ((len(L) == 1 and len(Rh) in (1, 2)) or (len(Rh) == 1 and len(L) == 2)):
            skipped["not_a_binary_invariant (%d+%d phases)" % (len(L), len(Rh))] += 1
            continue
        single = L if len(L) == 1 else Rh
        x0 = comp[single[0]]
        if len(L) == 1 and len(Rh) == 1 and x0 is None:
            x0 = comp[Rh[0]]
        if x0 is None:
            skipped["invariant_composition_unknown"] += 1
            continue
        if x0 <= 0.002 or x0 >= 0.998:
            skipped["unary_reaction"] += 1
            continue
        other = Rh if single is L else L
        known = [comp[p] for p in other if comp[p] is not None]
        lo, hi = 0.0, 1.0
        if len(other) == 2:
            if len(known) == 2:
                lo, hi = min(known), max(known)
                if not (lo + 1e-4 < x0 < hi - 1e-4):
                    skipped["composition_not_between_the_two_phases"] += 1
                    continue
            elif len(known) == 1:
                if abs(known[0] - x0) <= 1e-4:
                    skipped["composition_not_between_the_two_phases"] += 1
                    continue
                if known[0] < x0:
                    lo = known[0]
                else:
                    hi = known[0]
        delta = min(HYPER_DELTA, (x0 - lo) / 2, (hi - x0) / 2)
        u, ubasis = _t_sigma(rec)
        off = max(2 * u, MIN_OFFSET)
        myids = {R.phase_id(stem, p) or (stem + ":" + p) for p in comp}
        liq_ids = {k for k in myids if k.endswith(":liquid") or k.split(":")[-1].upper().startswith("LIQ")}

        def interferes(T1, T2):
            for (To, ido, xo) in inv[sysn]:
                if ido == frozenset(myids) or not (min(T1, T2) < To < max(T1, T2)):
                    continue
                shared = ido & myids
                if shared - liq_ids:
                    return True
                if shared & liq_ids:
                    for k in shared & liq_ids:
                        if xo.get(k) is not None and abs(xo[k] - x0) <= 0.05:
                            return True
            return False
        comps = [vx[1] for vx in reg if vx[0] != "__HYPERPLANE__"][0]

        def region(phases):
            r = []
            for p in phases:
                x = comp[p] if len(phases) > 1 else x0
                if R.stoichiometric(stem, p):
                    x = None
                r.append([p, list(comps), [x]])
            r.append(["__HYPERPLANE__", list(comps), [round(x0 - delta, 10)]])
            r.append(["__HYPERPLANE__", list(comps), [round(x0 + delta, 10)]])
            return r
        oid = d["orchid"]["observations"][i].get("id")
        for side, phases, Tb in (("above", L, T + off), ("below", Rh, T - off)):
            if interferes(T, Tb):
                skipped["another_invariant_within_offset (%s)" % side] += 1
                continue
            o = out.setdefault(rel, {"T": [], "values": [], "deriv": []})
            o["T"].append(round(Tb, 4))
            o["values"].append(region(phases))
            o["deriv"].append({"source_observation": oid, "side": side, "T_invariant_K": T, "offset_K": off,
                               "T_sigma_K": u, "T_sigma_basis": ubasis, "x0": x0, "hyperplane_delta": delta,
                               "reaction_printed": (rec.get("reaction") or {}).get("printed"),
                               "reaction_type": rtype, "phases_above": L, "phases_below": Rh})
            made[side] += 1
    datasets = {}
    for rel, o in sorted(out.items()):
        src = C["datasets"][rel]
        sysn = rel.split("/")[1]
        name = "derived/invariant-brackets/%s/%s" % (sysn, os.path.basename(rel))
        ids = sorted({x["source_observation"] for x in o["deriv"] if x["source_observation"]})
        datasets[name] = {
            "components": list(src["components"]),
            "phases": sorted({vx[0] for reg in o["values"] for vx in reg if vx[0] != "__HYPERPLANE__"}),
            "conditions": {"P": 101325, "T": o["T"]}, "output": "ZPF", "values": o["values"],
            "reference": src.get("reference"), "bibtex": src.get("bibtex"),
            "tags": ["derived", "derived_bracket", "orchid"] + [t for t in ("via_assessment",) if t in (src.get("tags") or [])],
            "comment": ("DERIVED phase constraints (not measured): invariant-temperature brackets from the printed "
                        "invariant reactions of %s (ORCHID literature values, Version 1, machine-extracted). Rule: "
                        "README.md, 'Phase constraints'. Source ORCHID observation ids: %s" % (rel, " ".join(ids))),
            "orchid": {"evidence_class": "derived", "derived_from": rel, "derivation": "invariant_bracket",
                       "rule": BRACKET_RULE, "points": o["deriv"], "dois": src["orchid"].get("dois"),
                       "papers": src["orchid"].get("papers")},
        }
    return datasets, skipped, made


# ------------------------------------------------------------------------------------------- assessed per TDB --
def assessed_per_tdb(data, R, stems_by_sys):
    root = os.path.join(data, "assessed-phase-assemblages")
    conf, budget = {}, {}
    for rp in (os.path.join(root, "records.jsonl.gz"), os.path.join(root, "dense-field-grid", "records.jsonl.gz")):
        if os.path.exists(rp):
            with gzip.open(rp, "rt", encoding="utf-8") as fh:
                for line in fh:
                    r = json.loads(line)
                    conf[r["id"]] = r.get("confidence")
                    budget[r["id"]] = r.get("uncertainty_budget")
    out, skipped = {}, Counter()
    # primary points (one per field), then the dense field grid of the validated diagrams (separate datasets, suffix .dense)
    srcs = [(f, "espei-zpf", "") for f in sorted(glob.glob(os.path.join(root, "espei-zpf", "*.json")))] + \
           [(f, "espei-zpf-dense", ".dense") for f in sorted(glob.glob(os.path.join(root, "espei-zpf-dense", "*.json")))]
    for f, srcdir, suffix in srcs:
        d = json.load(open(f, encoding="utf-8"))
        diag = os.path.basename(f)[:-5]
        sysn = diag.rsplit("-p", 1)[0]
        ids = (d.get("orchid") or {}).get("record_ids") or []
        Ts = d["conditions"]["T"]
        for stem in sorted(stems_by_sys.get(sysn, ())):
            vals, T2, ids2 = [], [], []
            for reg, T, rid in zip(d["values"], Ts, ids):
                new, ok = [], True
                for vx in reg:
                    if vx[0] == "__HYPERPLANE__":
                        new.append(list(vx))
                        continue
                    p = R.by_pref.get(sysn, {}).get(vx[0])
                    names = (R.phase_tdb_names.get(p["id"], {}).get(stem) if p else None) or []
                    if len(names) != 1:
                        ok = False
                        skipped["phase_not_in_tdb" if not names else "phase_ambiguous_in_tdb"] += 1
                        break
                    x = vx[2]
                    if R.stoichiometric(stem, names[0]):
                        x = [None] * len(x)
                    new.append([names[0], list(vx[1]), list(x)])
                if ok and zpf_inert(new, stem, R, e2=False):
                    ok = False
                    skipped["region_inert_no_composition"] += 1
                if ok:
                    vals.append(new)
                    T2.append(T)
                    ids2.append(rid)
            if not vals:
                continue
            rel = "assessed/%s/%s.%s%s.json" % (sysn, stem, diag, suffix)
            out[rel] = {
                "components": list(d["components"]),
                "phases": sorted({vx[0] for reg in vals for vx in reg if vx[0] != "__HYPERPLANE__"}),
                "conditions": {"P": d["conditions"].get("P", 101325.0), "T": T2}, "output": "ZPF", "values": vals,
                "reference": d.get("reference"),
                "tags": list(d.get("tags") or ["assessed_phase_diagram", "orchid"]),
                "comment": d.get("comment", "") + " Phase names mapped to %s.tdb through the ORCHID phase registry "
                                                  "(points whose phases are not all in that TDB are left out)." % stem,
                "orchid": {"evidence_class": "assessed_phase_diagram", "diagram": diag, "target_tdb": stem + ".tdb",
                           "record_ids": ids2, "confidence": [conf.get(r) for r in ids2],
                           "points_left_out": len(d["values"]) - len(vals),
                           "source": "assessed-phase-assemblages/%s/%s.json" % (srcdir, diag)},
            }
            if suffix:
                o_ = d["orchid"]
                keep = [k for k, (rid, v) in enumerate(zip(ids, d["values"])) if rid in set(ids2)]
                out[rel]["orchid"].update(point_method=o_.get("point_method"), dense_method_version=o_.get("dense_method_version"),
                                          field=[o_["field"][k] for k in keep], field_weight_share=[o_["field_weight_share"][k] for k in keep])
            ubud = assessed_budget_points(ids2, budget)
            if ubud:
                out[rel]["orchid_uncertainty_budget"] = ubud
    return out, skipped


def assessed_budget_points(ids, budget):
    """Per point (value order) the component values of the assessed record's uncertainty budget (temperature and
    composition), or None when the records carry no budget."""
    if not any(budget.get(r) for r in ids):
        return None
    keys = ("sigma_printed", "sigma_reading", "sigma_method", "sigma_boundary_margin", "sigma_effective")

    def comp(q):
        q = q or {}
        return {k: (q.get(k) or {}).get("value") if isinstance(q.get(k), dict) else q.get(k) for k in keys}
    first = next(budget[r] for r in ids if budget.get(r))
    pts = []
    for r in ids:
        q = (budget.get(r) or {}).get("quantities") or {}
        pts.append({"record_id": r, "temperature_K": comp(q.get("temperature_K")),
                    "composition_mole_fraction": comp(q.get("composition"))})
    return {"sigma_rule": first.get("sigma_rule"), "calibration": first.get("calibration"), "points": pts}


# ------------------------------------------------------------------------------------------- liquid stability --
def liquid_stability(data, R, pub, stems_by_sys):
    root = os.path.join(data, "assessed-phase-assemblages")
    per = defaultdict(lambda: {"liq_single": False, "gap": False, "Tmax_l2": None, "Tmax_solid": None})
    rp = os.path.join(root, "records.jsonl.gz")
    if os.path.exists(rp):
        with gzip.open(rp, "rt", encoding="utf-8") as fh:
            for line in fh:
                r = json.loads(line)
                s = per[r["system"]]
                labs = [a.get("label") or "" for a in r["assemblage"]]
                names = [(a.get("canonical_name") or "").upper() for a in r["assemblage"]]
                isliq = [n == "LIQUID" or lb.startswith("L") and (lb == "L" or lb.startswith("L#")) for lb, n in zip(labs, names)]
                if any(lb.startswith("L#") for lb in labs):
                    s["gap"] = True
                if r["n_phases"] == 1 and isliq[0] and not labs[0].startswith("L#"):
                    s["liq_single"] = True
                if any(isliq) and r["n_phases"] >= 2 and r.get("T_K"):
                    s["Tmax_l2"] = max(s["Tmax_l2"] or 0, r["T_K"])
                if not all(isliq) and r.get("T_K"):
                    s["Tmax_solid"] = max(s["Tmax_solid"] or 0, r["T_K"])
    lit_T = defaultdict(float)
    lit_gap = set()
    for r in pub.values():
        els = r.get("elements") or []
        if len(els) != 2 or not (r.get("flags") or {}).get("extraction_accepted"):
            continue
        sysn = "-".join(sorted(e.upper() for e in els))
        mp = [(p.get("model_phase") or "").upper() for p in r.get("phases") or []]
        labs = [(p.get("label") or "") for p in r.get("phases") or []]
        rtype = (r.get("reaction") or {}).get("type")
        if rtype in ("monotectic", "syntectic") or sum(1 for m in mp if m.startswith("LIQ")) >= 2 or \
                any(re.fullmatch(r"L\s*[12'\"]|L_?[12]|Liquid\s*[12]|liq\s*[12]", lb or "") for lb in labs):
            lit_gap.add(sysn)
        if r.get("evidence_class") == "experimental" and any(m.startswith("LIQ") for m in mp) and \
                isinstance(r.get("temperature_K"), (int, float)) and \
                r.get("observation_type") in ("invariant_reaction", "phase_region_vertex_set", "phase_transition"):
            lit_T[sysn] = max(lit_T[sysn], float(r["temperature_K"]))
    out, why = {}, Counter()
    for sysn in sorted(stems_by_sys):
        els = sysn.split("-")
        if len(els) != 2:
            continue
        s = per.get(sysn)
        if s is None:
            why["no_assessed_diagram"] += 1
            continue
        if s["gap"] or sysn in lit_gap:
            why["liquid_miscibility_gap_in_sources"] += 1
            continue
        if not s["liq_single"]:
            why["no_single_phase_liquid_field"] += 1
            continue
        if not all(e in ELEMENT_TM_TB for e in els):
            why["element_not_a_condensed_liquid_at_1_atm"] += 1
            continue
        tm = max(ELEMENT_TM_TB[e][0] for e in els)
        Tl = max(tm, lit_T.get(sysn, 0.0), s["Tmax_l2"] or 0.0)
        Tp = round(Tl + LIQ_MARGIN, 2)
        tb = min(ELEMENT_TM_TB[e][1] for e in els)
        if Tp >= tb:
            why["above_a_boiling_point"] += 1
            continue
        if s["Tmax_solid"] and s["Tmax_solid"] >= Tp:
            why["assessed_solid_field_above_T"] += 1
            continue
        xel = els[1]
        for stem in sorted(stems_by_sys[sysn]):
            if "LIQUID" not in R.models.get(stem, {}):
                continue
            comps = sorted(els + (["VA"] if R.va_needed(sysn, stem)[0] else []))
            rel = "derived/liquid-stability/%s/%s.liquid-stability.json" % (sysn, stem)
            out[rel] = {
                "components": comps, "phases": ["LIQUID"], "conditions": {"P": 101325, "T": [Tp] * len(LIQ_X)},
                "output": "ZPF", "values": [[["LIQUID", [xel], [x]]] for x in LIQ_X],
                "reference": "derived by ORCHID from the system's sources (see orchid.basis)",
                "tags": ["derived", "derived_liquid_stability", "orchid"],
                "comment": "DERIVED phase constraints (not measured): single-phase liquid above the highest liquidus "
                           "evidence, where the sources show complete liquid miscibility. Rule: README.md, 'Phase "
                           "constraints'.",
                "orchid": {"evidence_class": "derived", "derivation": "liquid_stability", "rule": LIQUID_RULE,
                           "target_tdb": stem + ".tdb",
                           "basis": {"T_liq_max_K": round(Tl, 2), "element_melting_K": {e: ELEMENT_TM_TB[e][0] for e in els},
                                     "highest_literature_liquid_T_K": lit_T.get(sysn) or None,
                                     "highest_assessed_liquid_two_phase_T_K": s["Tmax_l2"],
                                     "lowest_boiling_point_K": tb, "margin_K": LIQ_MARGIN}},
            }
        why["written"] += 1
    return out, why
