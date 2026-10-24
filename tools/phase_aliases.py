#!/usr/bin/env python3
"""Per-system aliases of conventional (Hume-Rothery / intermetallic) phase names (ORCHID, 2026-10-24).

The printed-label check (tools/phase_label_check.py, audit round 2) no longer maps a bare Greek letter to an element
structure (no global Greek -> structure rule): "beta" is bcc in Cu-Zn but "eps" is Cu3Sn in Cu-Sn, and a global rule sent
eps-Cu3Sn to BCC_A2. The classical names are restored here by a PER-SYSTEM table, tools/phase_aliases.json, keyed by
(chemical system, printed label: Greek letter with optional prime / index and optional formula). Each entry names the
registry identity (data/phases) and, per SGTE TDB of the system, the TDB phase that models it, with the evidence:

  * the TDB phase declared in each assessment (PHASE / CONSTITUENT lines, sublattice model, calculated stable range);
  * the handbook / registry structure (Pearson symbol, Strukturbericht, prototype) and the handbook composition and
    temperature range of the phase;
  * literature records of the system that print the label at a consistent composition / temperature.

An alias is added only when (a) the label is used for that phase in the system's own literature / handbook, (b) each
listed TDB has exactly one matching independent phase (composition range and structure consistent; disordered parts
of an order/disorder pair are never targets), and (c) the mapping is unambiguous at the record's composition: a printed
mole-fraction composition of the phase (or the overall composition of a one-phase value) must lie in the alias range
+- X_TOL, and a printed temperature in the alias temperature range +- T_TOL; a formula printed in the label
(eps-Cu3Sn, gamma' (Ni3Al)) must be the alias formula or lie in its composition range; otherwise the alias is not applied.
Uncertain cases stay unmapped and are listed in the file (`not_aliased`).

Applied by tools/phase_label_check.py BEFORE its Greek rule; a mapped phase carries
`model_phase_basis: system_alias:<id>` and the record a `corrections` entry (gate PL).

  python3 tools/phase_aliases.py --refresh   # recompute the generated evidence blocks in phase_aliases.json
  python3 tools/phase_aliases.py --check     # exit 1 if the evidence is stale or an alias fails validation

Standard library only; deterministic.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ALIAS_FILE = os.path.join(HERE, "phase_aliases.json")
DEFAULT_DATA = os.environ.get("ORCHID_DATA", os.path.expanduser("~/projects/orchid-data"))
DEFAULT_TDB_DIR = os.path.expanduser("~/projects/open-tdb/databases/sgte-binary")
RELEASE = "phase-aliases 2026-10-24"
X_TOL = 0.02          # mole fraction beyond the alias (handbook) range
T_TOL = 50.0          # K beyond the alias temperature range
GREEK = "αβγδεζηθικλμνξοπρστυφχψωΓΔΘΛΞΠΣΦΨΩ"
PRIMES = "'’′‘`´ʹ"
_GREEK_NAMES = {"alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε", "zeta": "ζ", "eta": "η",
                "theta": "θ", "kappa": "κ", "lambda": "λ", "mu": "μ", "sigma": "σ", "tau": "τ", "phi": "φ",
                "chi": "χ", "psi": "ψ", "omega": "ω"}
_GW = re.compile(r"(?i)^(%s)(?![a-z])" % "|".join(sorted(_GREEK_NAMES, key=len, reverse=True)))
_CORE = r"([%s])((?:[%s]{1,2})|\d{1,2})?" % (GREEK, re.escape(PRIMES))
_FORMULA = r"((?:[A-Z][a-z]?\d*(?:\.\d+)?)+)"


def _canon_core(g, mod):
    mod = mod or ""
    if mod and mod[0] in PRIMES:
        mod = "′" * len(mod)
    return g + mod


def label_key(label):
    """(core, formula) of a printed label: core = Greek letter + optional prime(s) / index ('ε', 'γ′', 'β1');
    formula = printed formula or None. None when the label is not of that form."""
    s = str(label or "").replace("µ", "μ").replace("ϵ", "ε").replace("ϕ", "φ").strip()
    s = _GW.sub(lambda m: _GREEK_NAMES[m.group(1).lower()], s)
    s = re.sub(r"(?i)\s*\bphase\b\s*", " ", s).strip()
    s = re.sub(r"\s+", " ", s)
    m = re.fullmatch(_CORE, s)
    if m:
        return _canon_core(m.group(1), m.group(2)), None
    m = re.fullmatch(_CORE + r"\s*[-–,]?\s*\(?\s*" + _FORMULA + r"\s*\)?", s)
    if m:
        return _canon_core(m.group(1), m.group(2)), m.group(3)
    m = re.fullmatch(_FORMULA + r"\s*[-–,]?\s*\(\s*" + _CORE + r"\s*\)", s)
    if m:
        return _canon_core(m.group(2), m.group(3)), m.group(1)
    return None


def reduced(formula):
    """{EL: mole fraction} of a formula string, or None."""
    if not formula:
        return None
    d = Counter()
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*(?:\.\d+)?)", formula):
        d[el.upper()] += float(n or 1)
    t = sum(d.values())
    return {k: round(v / t, 4) for k, v in sorted(d.items())} if t else None


PRINTED_SOURCES = ("printed", "printed_range")


def record_compositions(r, label):
    """Printed mole-fraction compositions of the phase `label` in record r; for a one-phase value without a printed
    phase composition, the record's overall composition (the composition of that phase)."""
    comps = [pc.get("composition") or {} for pc in r.get("phase_compositions") or []
             if pc.get("phase") == label and pc.get("composition_source") in PRINTED_SOURCES
             and str(pc.get("basis") or "").startswith("mole_fraction")]
    if not comps and len(r.get("phases") or []) == 1 and isinstance(r.get("overall_composition"), dict):
        comps = [r["overall_composition"]]
    return comps


def x_of(comp, xe):
    """Mole fraction of element xe in a printed composition {El: x} of a binary (one-key form = the other element)."""
    if not isinstance(comp, dict) or not comp:
        return None
    v = comp.get(xe)
    if isinstance(v, (int, float)):
        return float(v)
    if len(comp) == 1:
        (_, w), = comp.items()
        if isinstance(w, (int, float)):
            return 1.0 - float(w)
    return None


class Aliases:
    def __init__(self, path=ALIAS_FILE):
        self.path = path
        self.doc = json.load(open(path, encoding="utf-8"))
        self.aliases = self.doc["aliases"]
        self.by_key = defaultdict(list)
        self.by_id = {}
        for a in self.aliases:
            self.by_id[a["id"]] = a
            for lab in [a["label"]] + list(a.get("label_variants") or []):
                self.by_key[(a["system"], lab)].append(a)

    def lookup(self, system, label, comps=None, T=None, tdb=None):
        """{'alias', 'target', 'check'} when the label is an alias of the system and the record's composition /
        temperature agree; {'miss': reason} otherwise; None when the label is no alias key of the system."""
        k = label_key(label)
        if not k or (system, k[0]) not in self.by_key:
            return None
        core, formula = k
        cands = self.by_key[(system, core)]
        els = system.split("-")
        if formula and re.fullmatch(r"(?:[A-Z][a-z]?){2,}", formula) and \
                sorted(e.upper() for e in re.findall(r"[A-Z][a-z]?", formula)) == sorted(els):
            formula = None   # "kappa(CuSi)": the system's elements without indices name the system, not a formula
        if formula:      # a printed formula must be the alias formula, or lie in the alias composition range
            rf = reduced(formula)
            cands = [a for a in cands if (reduced(a["formula"]) == rf if a.get("formula") else
                                          rf and set(rf) <= set(els) and
                                          a["x_range"][0] - X_TOL <= rf.get(a["x_element"], 0.0) <= a["x_range"][1] + X_TOL)]
            if not cands:
                return {"miss": "printed_formula_differs_from_alias"}
        if len(cands) != 1:
            return {"miss": "several_aliases"}
        a = cands[0]
        xe, (lo, hi) = a["x_element"], a["x_range"]
        xs = [x_of(c, xe) for c in comps or []]
        xs = [x for x in xs if x is not None]
        if any(not (lo - X_TOL <= x <= hi + X_TOL) for x in xs):
            return {"miss": "composition_outside_alias_range", "alias": a}
        tr = a.get("T_range_K")
        tok = isinstance(T, (int, float)) and T > 1
        if tr and tok and not ((tr[0] is None or T >= tr[0] - T_TOL) and (tr[1] is None or T <= tr[1] + T_TOL)):
            return {"miss": "temperature_outside_alias_range", "alias": a}
        tdbs = a.get("tdb_phases") or {}
        target = tdbs.get(tdb) if tdb and tdb in tdbs else a.get("model_phase")
        if tdb and tdb not in tdbs and tdb in (a.get("excluded_tdbs") or {}):
            return {"miss": "alias_excluded_in_record_tdb", "alias": a}
        check = "+".join(x for x, ok in (("composition", bool(xs)), ("temperature", bool(tr and tok))) if ok) or \
            "label_only"
        return {"alias": a, "target": target, "check": check}

    @staticmethod
    def names(a):
        return set((a.get("tdb_phases") or {}).values())


# ------------------------------------------------------------------------------------------- evidence / validation --
def pearson_equal(a, b):
    a, b = (a or "").rstrip("*").strip(), (b or "").rstrip("*").strip()
    if not a or not b:
        return None
    if a == b:
        return True
    ma, mb = re.fullmatch(r"hR(\d+)", a), re.fullmatch(r"hR(\d+)", b)
    if ma and mb:
        na, nb = int(ma.group(1)), int(mb.group(1))
        return na == 3 * nb or nb == 3 * na
    return False


def model_range(m, xe):
    """(x_min, x_max) of element xe spanned by a TDB model's end members (VA ignored), from site ratios."""
    import itertools
    subs = m.get("sublattices") or []
    sets = [[c for c in s.get("constituents") or []] for s in subs]
    if not sets or any(not s for s in sets):
        return None
    xs = []
    for combo in itertools.islice(itertools.product(*sets), 4000):
        tot = Counter()
        for s, c in zip(subs, combo):
            if c != "VA":
                tot[c] += float(s.get("site_ratio") or 0)
        t = sum(tot.values())
        if t:
            xs.append(tot.get(xe, 0.0) / t)
    return (round(min(xs), 4), round(max(xs), 4)) if xs else None


def field_overlap(s, a):
    """True when a calculated stable field {x_min, x_max, T_min_K, T_max_K} overlaps the alias x / T range."""
    if not s:
        return False
    lo, hi = a["x_range"]
    if s["x_max"] < lo - X_TOL or s["x_min"] > hi + X_TOL:
        return False
    tr = a.get("T_range_K")
    if tr and ((tr[0] is not None and s["T_max_K"] < tr[0] - T_TOL) or (tr[1] is not None and s["T_min_K"] > tr[1] + T_TOL)):
        return False
    return True


def page_checks(path=os.path.join(HERE, "literature_visual_checks.json")):
    """Page checks of alias mappings (ledger kind phase_alias; key phase_alias:<alias id>:<record id>)."""
    if not os.path.exists(path):
        return []
    return sorted((c for c in json.load(open(path, encoding="utf-8")).get("checks", []) if c.get("kind") == "phase_alias"),
                  key=lambda c: c["key"])


def load_registry(data, systems):
    out = {}
    for s in systems:
        p = os.path.join(data, "phases", "systems", s + ".json")
        out[s] = json.load(open(p, encoding="utf-8")) if os.path.exists(p) else None
    return out


def tdb_files_of(tdb_dir, system):
    a, b = system.split("-")
    cap = lambda x: x[0] + x[1:].lower()
    return sorted(f for f in os.listdir(tdb_dir) if f.lower().endswith(".tdb") and
                  (f.startswith(cap(a) + cap(b) + "-") or f.startswith(cap(b) + cap(a) + "-")))


def build_evidence(a, reg, tdb_dir, lit):
    """Generated part of an alias: per-TDB phases, declarations, stability; registry structure; literature records.
    Returns (fields, problems)."""
    problems = []
    phases = {p["id"]: p for p in (reg or {}).get("phases") or []}
    xe = a["x_element"]
    lo, hi = a["x_range"]
    per_tdb = defaultdict(list)
    fields = {}
    decl, struct = [], []
    for rid in a["registry_ids"]:
        p = phases.get(rid)
        if p is None:
            problems.append("registry id %s not in the registry" % rid)
            continue
        st = p.get("structure") or {}
        struct.append({"registry_id": rid, "preferred_name": p.get("preferred_name"), "pearson": st.get("pearson"),
                       "strukturbericht": st.get("strukturbericht"), "prototype": st.get("prototype"),
                       "structure_basis": st.get("basis")})
        stab = {(s["tdb"]): s for s in p.get("stability") or [] if s.get("stable") and s.get("x_min") is not None
                and s.get("x_element") == xe}
        for m in p.get("models") or []:
            if m.get("disordered_part_of"):
                continue          # a DIS_PART of an ordered phase is never an independent target
            if m["tdb"] in (a.get("excluded_tdbs") or {}):
                continue
            per_tdb[m["tdb"]].append(m["phase_name"])
            fields[(m["tdb"], m["phase_name"])] = stab.get(m["tdb"])
            v = (m.get("tdb_definition") or {}).get("verbatim") or {}
            s = stab.get(m["tdb"])
            mr = model_range(m, xe)
            decl.append({"tdb": m["tdb"], "phase": m["phase_name"], "registry_id": rid,
                         "phase_line": (v.get("phase") or {}).get("text"),
                         "constituent_line": (v.get("constituent") or {}).get("text"),
                         "line": (v.get("phase") or {}).get("line_start"),
                         "sublattice_model": " : ".join(",".join(sl.get("constituents") or []) + "/%g" % sl["site_ratio"]
                                                         for sl in m.get("sublattices") or []),
                         "model_x_range": list(mr) if mr else None,
                         "calculated_stable": {"x_min": round(s["x_min"], 4), "x_max": round(s["x_max"], 4),
                                               "T_min_K": s["T_min_K"], "T_max_K": s["T_max_K"]} if s else None})
            # composition consistency: the model must reach the alias range; its calculated field must overlap it
            if mr and (mr[1] < lo - X_TOL or mr[0] > hi + X_TOL):
                problems.append("%s %s model range %s misses the alias range %s" % (m["tdb"], m["phase_name"], mr, a["x_range"]))
        # structure consistency with the handbook structure of the alias
        hs = a.get("structure") or {}
        if pearson_equal(hs.get("pearson"), st.get("pearson")) is False:
            problems.append("registry %s Pearson %s differs from the handbook %s" % (rid, st.get("pearson"), hs.get("pearson")))
        sb1, sb2 = [re.sub(r"[_\s]", "", x or "").upper() for x in (hs.get("strukturbericht"), st.get("strukturbericht"))]
        if sb1 and sb2 and sb1 != sb2:
            problems.append("registry %s Strukturbericht %s differs from the handbook %s" % (rid, sb2, sb1))
    tdb_phases, not_in_field = {}, {}
    for t, ns in sorted(per_tdb.items()):
        ns = sorted(set(ns))
        # the listed phase whose calculated field (data/phase-diagrams) covers the alias range; a TDB whose listed
        # phases all have calculated fields outside it does not model the phase there (excluded, not a problem)
        with_field = [n for n in ns if fields.get((t, n))]
        in_field = [n for n in with_field if field_overlap(fields[(t, n)], a)]
        if in_field:
            ns = in_field
        elif with_field and len(with_field) == len(ns):
            not_in_field[t] = "calculated field of %s outside the alias range" % "/".join(ns)
            continue
        if len(ns) != 1:
            problems.append("TDB %s has several matching phases %s" % (t, ns))
        else:
            tdb_phases[t] = ns[0]
    if not tdb_phases:
        problems.append("no TDB phase")
    cnt = Counter(tdb_phases.values())
    model_phase = sorted(cnt, key=lambda n: (-cnt[n], n))[0] if cnt else None
    unmodelled = sorted(set(tdb_files_of(tdb_dir, a["system"])) - set(tdb_phases) - set(a.get("excluded_tdbs") or {}))
    for t, why in not_in_field.items():
        a.setdefault("excluded_tdbs", {})[t] = why
    recs = lit.get((a["system"], a["id"])) or []
    pcs = [{"record": c["key"].split(":")[-1], "page": c.get("page"), "pdf_sha256_prefix": c.get("pdf_sha256_prefix"),
            "verdict": c.get("verdict"), "checked": c.get("checked")}
           for c in page_checks() if c["key"].startswith("phase_alias:%s:" % a["id"])]
    if not recs:
        problems.append("no literature record prints the label at a consistent composition / temperature")
    return {"model_phase": model_phase, "tdb_phases": tdb_phases, "tdbs_without_the_phase": unmodelled,
            "evidence": {"tdb_declarations": decl, "registry_structure": struct,
                         "literature_records": recs[:5], "literature_records_total": len(recs),
                         "page_checks": pcs}}, problems


def literature_examples(path, table):
    """(system, alias id) -> [{id, paper, T_K, x, label}] of records whose label is that alias and whose printed
    composition / temperature agree (deterministic: sorted by id)."""
    out = defaultdict(list)
    if not os.path.exists(path):
        return out
    for line in gzip.open(path, "rt", encoding="utf-8"):
        r = json.loads(line)
        if r.get("system") not in {a["system"] for a in table.aliases}:
            continue
        for p in r.get("phases") or []:
            lab = p.get("label")
            comps = record_compositions(r, lab)
            hit = table.lookup(r["system"], lab, comps, r.get("temperature_K"))
            if not hit or "alias" not in hit or "miss" in hit or hit["check"] == "label_only":
                continue
            a = hit["alias"]
            out[(r["system"], a["id"])].append({
                "id": r["id"], "paper": r.get("paper"), "label": lab, "T_K": r.get("temperature_K"),
                "x": [round(x_of(c, a["x_element"]), 4) for c in comps if x_of(c, a["x_element"]) is not None],
                "observation_type": r.get("observation_type")})
    for k in out:
        out[k].sort(key=lambda e: e["id"])
    return out


def validate(table):
    """Uniqueness and internal consistency of the curated table. Returns a list of problems."""
    probs = []
    seen = Counter((a["system"], lab, (reduced(a.get("formula")) and json.dumps(reduced(a["formula"]))))
                   for a in table.aliases for lab in [a["label"]] + list(a.get("label_variants") or []))
    probs += ["duplicate alias key %s" % (k,) for k, n in seen.items() if n > 1]
    ids = Counter(a["id"] for a in table.aliases)
    probs += ["duplicate alias id %s" % k for k, n in ids.items() if n > 1]
    for a in table.aliases:
        for lab in [a["label"]] + list(a.get("label_variants") or []):
            if label_key(lab) != (lab, None):
                probs.append("%s: label %r is not a canonical core" % (a["id"], lab))
        lo, hi = a["x_range"]
        if not (0 <= lo <= hi <= 1):
            probs.append("%s: x_range %s" % (a["id"], a["x_range"]))
        if a.get("formula"):
            rf = reduced(a["formula"])
            x = rf.get(a["x_element"], 0.0) if rf else None
            if x is None or not (lo - X_TOL <= x <= hi + X_TOL):
                probs.append("%s: formula %s outside its own x_range" % (a["id"], a["formula"]))
        if not a.get("basis"):
            probs.append("%s: no handbook / literature basis" % a["id"])
    by_sys = defaultdict(list)
    for a in table.aliases:
        by_sys[a["system"]].append(a)
    for s, al in by_sys.items():      # two aliases of one system on one TDB phase only for disjoint composition ranges
        tgt = defaultdict(list)
        for a in al:
            for t, n in (a.get("tdb_phases") or {}).items():
                tgt[(t, n)].append(a)
        for k, als in tgt.items():
            for i in range(len(als)):
                for j in range(i + 1, len(als)):
                    (l1, h1), (l2, h2) = als[i]["x_range"], als[j]["x_range"]
                    if not (h1 + X_TOL < l2 - X_TOL or h2 + X_TOL < l1 - X_TOL):
                        probs.append("%s: TDB phase %s targeted by %s and %s with overlapping ranges"
                                     % (s, k, als[i]["id"], als[j]["id"]))
    return probs


def refresh(path=ALIAS_FILE, data=DEFAULT_DATA, tdb_dir=DEFAULT_TDB_DIR, write=True):
    table = Aliases(path)
    reg = load_registry(data, sorted({a["system"] for a in table.aliases}))
    lit = literature_examples(os.path.join(data, "literature", "orchid-literature.jsonl.gz"), table)
    problems = []
    for a in table.aliases:
        gen, pr = build_evidence(a, reg.get(a["system"]), tdb_dir, lit)
        problems += ["%s: %s" % (a["id"], p) for p in pr]
        a.update(gen)
    table.aliases.sort(key=lambda a: (a["system"], a["label"]))
    problems += validate(table)
    doc = table.doc
    doc["counts"] = {"aliases": len(table.aliases), "systems": len({a["system"] for a in table.aliases}),
                     "not_aliased": len(doc.get("not_aliased") or [])}
    text = json.dumps(doc, indent=1, ensure_ascii=False, sort_keys=False) + "\n"
    if write:
        open(path, "w", encoding="utf-8").write(text)
    return text, problems


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--refresh", action="store_true")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--data", default=DEFAULT_DATA)
    ap.add_argument("--tdb-dir", default=DEFAULT_TDB_DIR)
    A = ap.parse_args()
    sys.path.insert(0, HERE)
    text, problems = refresh(data=A.data, tdb_dir=A.tdb_dir, write=A.refresh)
    for p in problems:
        print("PROBLEM", p)
    if A.check:
        cur = open(ALIAS_FILE, encoding="utf-8").read()
        if cur != text:
            print("stale: phase_aliases.json evidence differs from a fresh refresh")
            sys.exit(1)
    t = Aliases(ALIAS_FILE)
    print("%d aliases in %d systems; %d problems" % (len(t.aliases), len({a["system"] for a in t.aliases}), len(problems)))
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()
