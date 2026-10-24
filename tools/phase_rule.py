#!/usr/bin/env python3
"""Physical phase rule for ESPEI use of literature values (shared by the literature export and the ESPEI bundle).

A value may be used in ESPEI with a phase only when

  1. the phase is what the paper says: a printed phase label of the value (or, when no label is printed, the table
     caption) names it; a label is mapped to a TDB phase through the phase registry (paper label -> phase, or the
     record's model phase); no default phase is ever assigned (in particular not LIQUID);
  2. it is physically possible at the value's temperature: a liquid is impossible below the lowest temperature at
     which any calculated SGTE diagram of the system has a liquid (data/phase-diagrams), less a margin; without a
     calculated diagram, below LIQ_TM_FRACTION x the lowest melting point of the system's elements;
  3. a property of a stoichiometric compound is usable only as a formation quantity that ESPEI 0.9.1 can fit as a
     non-equilibrium endmember dataset (HM_FORM, SM_FORM); other properties of a compound (CPM, SM, GM_FORM, a
     mixing quantity of a compound) are not usable in ESPEI (they stay in the literature data).

A solution-phase assignment at a composition / temperature outside that phase's field in the calculated diagram of
the target TDB is flagged only (`outside_calculated_field`), never excluded.

Standard library only; deterministic; reads data/phases/systems/*.json and data/phase-diagrams/sgte/*.json.
"""
from __future__ import annotations

import glob
import json
import os
import re
from collections import defaultdict

LIQ_MARGIN_K = 100.0          # a liquid is impossible below (lowest calculated liquid T) - max(margin, 10 %)
LIQ_MARGIN_FRACTION = 0.10
EQUILIBRIUM_MARGIN_FACTOR = 2.0   # printed invariant / phase-composition records: twice the margin (the printed
#                                 reaction is itself liquid evidence; only gross misreadings are rejected)
LIQ_TM_FRACTION = 0.40        # no calculated diagram: impossible below 0.40 x the lowest element melting point
INVARIANT_T_TOL = 50.0        # K: a binary three-phase ZPF region must lie this close to a calculated invariant
MIN_T_K = 1.0                 # T <= 1 K: zero or a composition read as a temperature
EQUILIBRIUM_T_FRACTION = 0.25  # a printed phase equilibrium below 0.25 x the lowest element melting point is a misreading
BLOCKING_FLAGS = ("unsigned_in_print", "magnitude_suspect", "system_mismatch", "source_disputed")
FORMATION_FITTABLE = ("HM_FORM", "SM_FORM")    # ESPEI 0.9.1 non-equilibrium properties for compound endmembers
NOT_PRINTED = {"", "UNKNOWN", "NONE", "N/A", "?"}

RULE_TEXT = (
    "Physical phase rule (tools/phase_rule.py): a literature value is assigned to a phase only when the paper names "
    "that phase (printed phase label of the value; when no label is printed, a table caption that names the liquid) "
    "and the phase registry maps the label to a model phase; no default phase is assigned. A liquid assignment must be "
    "physically possible: not below the lowest liquid temperature of any calculated SGTE diagram of the system minus "
    "max(%g K, %d %%) (twice that for printed invariant and phase-composition records; no calculated diagram: not below %.2f x the lowest element melting point). Formation enthalpies "
    "and entropies of a stoichiometric compound go to that compound's TDB endmember (ESPEI non-equilibrium dataset, "
    "`solver`); other properties of a compound are not used in ESPEI. A solution-phase assignment outside the phase's "
    "field in the calculated diagram of the target TDB is flagged (`outside_calculated_field`), not excluded. "
    "Temperatures: T <= %g K is not a measured state; a printed phase equilibrium (invariant, phase composition, "
    "transition) below %.2f x the lowest element melting point is a misreading. A value flagged %s is not ESPEI-ready."
    % (LIQ_MARGIN_K, int(LIQ_MARGIN_FRACTION * 100), LIQ_TM_FRACTION, MIN_T_K, EQUILIBRIUM_T_FRACTION, " / ".join("`%s`" % f for f in BLOCKING_FLAGS)))

# Melting points (K), SGTE unary data (A.T. Dinsdale, CALPHAD 15 (1991) 317); elements without a liquid at 1 atm absent.
from espei_constraints import ELEMENT_TM_TB  # noqa: E402
import build_dft_parity as parity  # noqa: E402  (formula parser shared with the DFT benchmark)

STANDARD_NAMES = {"fcc": "FCC_A1", "(fcc)": "FCC_A1", "bcc": "BCC_A2", "(bcc)": "BCC_A2", "hcp": "HCP_A3",
                  "(hcp)": "HCP_A3", "sigma": "SIGMA", "σ": "SIGMA", "a1": "FCC_A1", "a2": "BCC_A2", "a3": "HCP_A3"}
_LIQ_WORD = re.compile(r"(?i)liq(?!uidus)|\bmelt\b|molten")
_LIQ_EXACT = re.compile(r"\s*\(?\s*(?:L|L'|L#\d|L\s*\(\s*liquid\s*\))\s*\)?\s*")
_CAP_REF = re.compile(r"(?i)(?:pure|the)\s+liquid\s+(?:components?|elements?|metals?)|reference\s+states?[^.;)]*"
                      r"|relative\s+to\s+(?:the\s+)?(?:pure\s+)?liquid\s+\w+|with\s+respect\s+to[^.;)]*")
_CAP_LIQ = re.compile(r"(?i)\bliquid\b(?!us)|\bmelts?\b|\bmolten\b")
_CAP_SOLID = re.compile(r"(?i)\bsolid\b|intermetallic|compound|\boxides?\b|silicide|boride|carbide|nitride|"
                        r"crystal|\bphases\b")


def is_liquid_label(label, model_phase=None):
    if model_phase and str(model_phase).upper().split("#")[0] in ("LIQUID", "IONIC_LIQ", "LIQUID_IONIC"):
        return True
    lab = str(label or "")
    if not lab.strip() or lab.strip().upper() in NOT_PRINTED:
        return False
    if re.search(r"(?i)liquidus|solidus|solid", lab):
        return False
    return bool(_LIQ_WORD.search(lab) or _LIQ_EXACT.fullmatch(lab))


def caption_liquid(caption):
    """True when the caption names a liquid as the state measured (reference-state phrases removed; any solid /
    compound wording makes it not liquid)."""
    c = _CAP_REF.sub(" ", str(caption or ""))
    return bool(_CAP_LIQ.search(c)) and not _CAP_SOLID.search(c)


def x_of(rec, xel):
    comp = rec.get("overall_composition") or {}
    v = comp.get(xel)
    if isinstance(v, (int, float)):
        return float(v)
    if len(comp) == 1:
        (el, w), = comp.items()
        if isinstance(w, (int, float)):
            return 1.0 - float(w)
    return None


class PhaseRule:
    def __init__(self, data, open_tdb=None):
        self.data = data
        self.open_tdb = open_tdb or os.path.expanduser("~/projects/open-tdb")
        self.by_paper = defaultdict(list)   # "paper:<id>:<label>" -> [phase entries]
        self.alias_paper_names = set()     # (paper name, phase id) linked through a per-system alias
        self.by_tdbname = defaultdict(dict)  # system -> TDB phase name -> phase entry
        self.by_formula = defaultdict(lambda: defaultdict(list))   # system -> reduced composition -> compounds
        self.phases = {}
        for f in sorted(glob.glob(os.path.join(data, "phases", "systems", "*.json"))):
            d = json.load(open(f, encoding="utf-8"))
            for p in d["phases"]:
                self.phases[p["id"]] = p
                for n in p.get("names") or []:
                    if n.get("namespace") == "paper":
                        self.by_paper[n["name"]].append(p)
                        if str(n.get("basis") or "").startswith("system_alias:"):
                            self.alias_paper_names.add((n["name"], p["id"]))
                for m in p.get("models") or []:
                    self.by_tdbname[d["system"]].setdefault(m["phase_name"].upper(), p)
                if self.compound(p):
                    keys = {tuple(sorted((k.upper(), round(v, 4)) for k, v in (m.get("composition") or {}).items()))
                            for m in p["models"]}
                    if len(keys) == 1:
                        self.by_formula[d["system"]][keys.pop()].append(p)
        self.liq_min, self.grids, self.inv_T = {}, {}, defaultdict(list)
        idx = os.path.join(data, "phase-diagrams", "index.json")
        systems = json.load(open(idx, encoding="utf-8"))["systems"] if os.path.exists(idx) else {}
        for sysn, files in systems.items():
            for e in files:
                p = os.path.join(data, "phase-diagrams", "sgte", e["file"])
                if not os.path.exists(p):
                    continue
                d = json.load(open(p, encoding="utf-8"))
                g = d.get("grid")
                if d.get("status") != "ok" or not g:
                    continue
                stem = e["tdb"][:-4] if e["tdb"].lower().endswith(".tdb") else e["tdb"]
                liq = [i for i, lb in enumerate(g["labels"]) if any(is_liquid_label(None, ph) for ph in lb.split("+"))]
                tmin = None
                for k, row in enumerate(g["rows"]):
                    if any(i in liq for i, _ in row):
                        tmin = g["T0"] + k * g["dT"]
                        break
                for inv in d.get("invariants") or []:
                    if any(is_liquid_label(None, pt.get("phase")) for pt in inv.get("points") or []):
                        tmin = inv["T"] if tmin is None else min(tmin, inv["T"])
                self.inv_T[sysn].extend(inv["T"] for inv in d.get("invariants") or [])
                if tmin is not None:
                    cur = self.liq_min.get(sysn)
                    self.liq_min[sysn] = (tmin, e["tdb"]) if cur is None or tmin < cur[0] else cur
                self.grids[stem] = {"x_component": d.get("x_component", "").upper(), "grid": g}
        self._ser = {}

    # ------------------------------------------------------------------ liquid temperature bound --
    def liquid_bound(self, system, factor=1.0):
        """(T_K below which a liquid is impossible, basis dict) or (None, basis)."""
        if system in self.liq_min:
            t, tdb = self.liq_min[system]
            return t - factor * max(LIQ_MARGIN_K, LIQ_MARGIN_FRACTION * t), {
                "basis": "calculated_sgte_diagram", "lowest_liquid_T_K": round(t, 2), "tdb": tdb}
        tms = [ELEMENT_TM_TB[e][0] for e in system.split("-") if e in ELEMENT_TM_TB]
        if tms and len(tms) == len(system.split("-")):
            return LIQ_TM_FRACTION * min(tms), {"basis": "element_melting_points", "lowest_melting_K": min(tms)}
        return None, {"basis": "no_bound (element without a liquid at 1 atm, no calculated diagram)"}

    # ------------------------------------------------------------------ label -> registry phase --
    def map_label(self, rec, label, model_phase):
        pk = rec.get("paper") or ""
        ns = ("paper:" + pk[4:]) if pk.startswith("doi:") else ("paper:sha256-" + pk[7:] if pk.startswith("sha256:")
                                                                  else None)
        key = "%s:%s" % (ns, label)
        # a paper label linked through a per-system alias (tools/phase_aliases.json) names a phase of that system
        # only: it never maps the same printed label in another system of the paper
        cands = {p["id"]: p for p in (self.by_paper.get(key, []) if ns else [])
                 if p.get("system") == rec.get("system") or (key, p["id"]) not in self.alias_paper_names}
        if not cands and model_phase:
            p = self.by_tdbname.get(rec.get("system"), {}).get(str(model_phase).upper().split("#")[0])
            if p:
                cands[p["id"]] = p
        if not cands:   # the printed label is the TDB phase name or a standard structure name of one
            names = self.by_tdbname.get(rec.get("system"), {})
            key = re.sub(r"[\s-]+", "_", str(label).strip()).upper()
            for n in (key, STANDARD_NAMES.get(str(label).strip().lower())):
                if n and n in names:
                    cands[names[n]["id"]] = names[n]
                    break
        if not cands:   # the printed label is a formula: the one registry compound with that composition
            system = [e.capitalize() for e in (rec.get("system") or "").split("-")]
            c, _ = parity.parse_phase(label, system)
            if c:
                tot = sum(float(v) for v in c.values())
                key = tuple(sorted((k.upper(), round(float(v) / tot, 4)) for k, v in c.items()))
                for p in self.by_formula.get(rec.get("system"), {}).get(key, []):
                    cands[p["id"]] = p
        return list(cands.values())

    @staticmethod
    def compound(p):
        ms = p.get("models") or []
        return bool(ms) and all(m.get("stoichiometric") for m in ms)

    # ------------------------------------------------------------------ the rule --
    def classify(self, rec):
        """{'result': pass|fail|not_applicable, 'reason', 'phases': [{label, kind, registry_id}], 'basis'}"""
        ot = rec.get("observation_type")
        T = rec.get("temperature_K")
        system = rec.get("system") or ""
        labels = [(p.get("label"), p.get("model_phase")) for p in rec.get("phases") or []]
        printed = [(lb, mp) for lb, mp in labels if str(lb or "").strip().upper() not in NOT_PRINTED or mp]
        out = {"result": "pass", "phases": []}
        liquid = False
        if ot != "equilibrium_property":
            liquid = any(is_liquid_label(lb, mp) for lb, mp in printed)
            out["basis"] = "printed_labels"
        else:
            if not printed:
                cap = (rec.get("locator") or {}).get("caption")
                if caption_liquid(cap):
                    out["phases"] = [{"label": None, "kind": "liquid", "basis": "caption"}]
                    out["basis"] = "caption_names_liquid"
                    liquid = True
                else:
                    return {"result": "fail", "reason": "phase_not_printed",
                            "detail": "no phase label printed and the caption does not name the liquid"}
            else:
                out["basis"] = "printed_labels"
                for lb, mp in printed:
                    if is_liquid_label(lb, mp):
                        liquid = True
                        out["phases"].append({"label": lb, "kind": "liquid"})
                        continue
                    c = self.map_label(rec, lb, mp)
                    if len(c) != 1:
                        return {"result": "fail", "reason": "phase_label_not_mapped" if not c
                                else "phase_label_ambiguous", "detail": "label %r" % lb}
                    p = c[0]
                    if not p.get("models"):
                        return {"result": "fail", "reason": "phase_without_tdb_model", "detail": p["id"]}
                    out["phases"].append({"label": lb, "kind": "compound" if self.compound(p) else "solution",
                                          "registry_id": p["id"]})
                kinds = [p["kind"] for p in out["phases"]]
                if len(kinds) == 1 and kinds[0] == "compound" and rec.get("property") not in FORMATION_FITTABLE:
                    return {"result": "fail", "reason": "compound_property_not_fittable",
                            "detail": "%s of a stoichiometric compound (ESPEI 0.9.1 fits only %s of compound "
                                      "endmembers)" % (rec.get("property"), "/".join(FORMATION_FITTABLE))}
                if "compound" in kinds and len(kinds) == 1:
                    out["encoding"] = "compound_endmember"
        if isinstance(T, (int, float)) and rec.get("evidence_class") != "dft":
            tms = [ELEMENT_TM_TB[e][0] for e in system.split("-") if e in ELEMENT_TM_TB]
            if T <= MIN_T_K:
                return {"result": "fail", "reason": "temperature_implausible",
                        "detail": "T = %g K is not a temperature of a measured state (zero or a composition)" % T}
            if ot != "equilibrium_property" and tms and T < EQUILIBRIUM_T_FRACTION * min(tms):
                return {"result": "fail", "reason": "temperature_implausible",
                        "detail": "phase equilibrium at %g K, below %.2f x the lowest element melting point (%g K)"
                                  % (T, EQUILIBRIUM_T_FRACTION, min(tms))}
        if liquid and isinstance(T, (int, float)):
            tb, basis = self.liquid_bound(system, 1.0 if ot == "equilibrium_property" else EQUILIBRIUM_MARGIN_FACTOR)
            out["liquid_bound"] = dict(basis, impossible_below_K=round(tb, 2) if tb is not None else None)
            if tb is not None and T < tb:
                return {"result": "fail", "reason": "liquid_impossible_at_T",
                        "detail": "liquid at %g K, below %.0f K (%s)" % (T, tb, basis["basis"]),
                        "liquid_bound": out["liquid_bound"]}
        return out

    def near_invariant(self, system, T, tol=INVARIANT_T_TOL):
        """True when a calculated invariant of the system lies within tol of T, or the system has no diagram."""
        if system not in self.inv_T:
            return True
        return any(abs(T - t) <= tol for t in self.inv_T[system])

    # ------------------------------------------------------------------ diagram field (flag only) --
    def in_field(self, stem, phase, x, T):
        """True / False (phase present in the calculated diagram cell at x, T of the TDB), or None (no diagram)."""
        G = self.grids.get(stem)
        if G is None or x is None or T is None:
            return None
        g = G["grid"]
        k = int(round((T - g["T0"]) / g["dT"]))
        if k < 0 or k >= len(g["rows"]):
            return None
        j = int(round(x / g["x_step"]))
        pos, lab = 0, None
        for i, n in g["rows"][k]:
            if pos <= j < pos + n:
                lab = g["labels"][i]
                break
            pos += n
        if lab is None:
            return None
        return phase.upper() in {s.upper() for s in lab.split("+")}

    # ------------------------------------------------------------------ SER phases of a TDB --
    def ser_phases(self, model_path):
        if model_path not in self._ser:
            p = os.path.join(self.open_tdb, model_path)
            out = {}
            if os.path.exists(p):
                for line in open(p, encoding="utf-8", errors="replace"):
                    m = re.match(r"\s*ELEM(?:ENT)?\s+(\S+)\s+(\S+)", line, re.I)
                    if m:
                        out[m.group(1).upper()] = m.group(2).upper()
            self._ser[model_path] = out
        return self._ser[model_path]
