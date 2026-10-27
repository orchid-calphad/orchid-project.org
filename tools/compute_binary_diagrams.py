#!/usr/bin/env python3
"""Compute binary phase diagrams (ZPF boundaries + invariants) from TDB files with pycalphad's BinaryStrategy and
write one compact JSON per TDB for the ORCHID phase-diagram viewer (2026-10-22).
  compute_binary_diagrams.py OUT_DIR TDB [TDB ...] [--workers N] [--tmax 2500]"""
import argparse, hashlib, json, math, os, re, sys, time, multiprocessing as mp
import warnings; warnings.filterwarnings("ignore")

def elements_from_name(fn):
    stem = os.path.basename(fn).split("-")[0]
    return re.findall(r"[A-Z][a-z]?", stem)

def r(a, n):  # round arrays for compact JSON
    return [round(float(x), n) for x in (a if getattr(a, "shape", None) != () else [a])]

def prepare_source(raw):
    """Decode legacy text and apply documented, bounded compatibility repairs.

    Original bytes remain the source of sha256. Every removed statement is
    recorded, and active order/disorder directives are never changed.
    """
    try:
        text = raw.decode("utf-8")
        encoding = "utf-8"
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
        encoding = "latin-1"
    clean = re.sub(r"(?m)\$[^\n]*", "", text)
    phases = re.findall(r"\bPHASE\s+(\S+)\s+(\S+)", clean, re.I)
    names = {name.split(":")[0].upper() for name, _ in phases}
    used = set("".join(code for _, code in phases))
    removed = []
    # COST 507 pp. 224-226 and the frozen COST507 database both support two
    # substitutional sublattices for MgY. This exact Hallstedt file has an
    # extra third VA constituent group, incompatible with its own PHASE and
    # PARAMETER declarations. Do not generalize this repair to other files.
    if hashlib.sha256(raw).hexdigest() == "68740131a132fee3651ec3eb1fd044d3a56428ed2828054a4818d2079e1115cc":
        original = " CONST MGY_B2 : MG : MG Y% : VA : !"
        replacement = " CONST MGY_B2 : MG : MG Y% : !"
        if text.count(original) != 1:
            raise ValueError("Verified MgY repair no longer matches source statement")
        text = text.replace(original, replacement, 1)
        removed.append({"rule": "verified_mgy_extra_vacancy_sublattice",
                        "statement": original.strip(), "replacement": replacement.strip(),
                        "evidence": "review/mgy-1998-20261026/provenance.json"})
    pattern = re.compile(r"(?im)^[ \t]*TYPE[-_]DEF(?:INITION)?\s+(\S)\s+GES\s+(?:AMEND_PHASE_DESCRIPTION|A_P_D)\s+(\S+)\s+DIS_PART\s+(\S+)\s*!")
    def replace(match):
        char, ordered, disordered = match.groups()
        if char not in used and ordered.upper() not in names:
            removed.append({"rule": "unused_order_directive_for_absent_phase", "statement": match.group(0).strip()})
            return ""
        return match.group(0)
    text = pattern.sub(replace, text)
    return text, encoding, removed

def validate_geometry(rec):
    """Reject nonfinite and nonphysical coordinates before publishing success."""
    if not rec["regions"] or not any(branch["x"] for region in rec["regions"] for branch in region):
        raise ValueError("No mapped phase boundaries; success cannot be established")
    lo, hi = rec["T_range"]
    def coordinate(x, t):
        if not (math.isfinite(x) and -1e-5 <= x <= 1 + 1e-5
                and math.isfinite(t) and lo - 0.02 <= t <= hi + 0.02):
            raise ValueError("Invalid phase-diagram coordinate")
    for region in rec["regions"]:
        for branch in region:
            if len(branch["x"]) != len(branch["T"]):
                raise ValueError("Mismatched boundary coordinates")
            for x, t in zip(branch["x"], branch["T"]):
                coordinate(x, t)
    for inv in rec["invariants"]:
        for point in inv["points"]:
            coordinate(point["x"], inv["T"])

def source_phase_selection(text, phase_names):
    """Apply only explicit source DEFAULT_COM REJECT_PHASE instructions.

    This is opt-in because it changes the phase domain of the calculation.
    Unknown phase names fail closed instead of silently ignoring a directive.
    """
    clean = re.sub(r"(?m)\$[^\n]*", "", text)
    statements = re.findall(r"(?im)^[ \t]*DEFAULT[-_]COM(?:MAND)?\s+REJECT_PHASE\s+([^!]+)!", clean)
    excluded = set()
    for statement in statements:
        names = {n.upper().split(":")[0] for n in statement.split()}
        if not names.issubset(set(phase_names)):
            raise ValueError("Unknown phase in source REJECT_PHASE: " + str(sorted(names-set(phase_names))))
        excluded.update(names)
    active = sorted(set(phase_names)-excluded)
    if not active:
        raise ValueError("Source defaults exclude every phase")
    return active, {"rule": "source_DEFAULT_COM_REJECT_PHASE", "excluded": sorted(excluded),
                    "statements": ["DEFAULT_COM REJECT_PHASE " + s.strip() + " !" for s in statements]}

REFRACTORY = {"W","Re","Os","Ta","Mo","Nb","Ir","Ru","Hf","Tc","Rh","C","B","Zr","Cr","V","Th","Pu","U"}
def _alarm(signum, frame): raise TimeoutError("diagram exceeded time limit")
def one(args):
    import signal; signal.signal(signal.SIGALRM, _alarm); signal.alarm(600)
    try:
        return _one(args)
    finally:
        signal.alarm(0)
def _one(args):
    tdb, out, tmax = args[:3]
    interior_seeds = len(args) > 3 and args[3]
    source_defaults = len(args) > 4 and args[4]
    if out and os.path.exists(os.path.join(out, os.path.basename(tdb).rsplit(".", 1)[0] + ".json")):
        return os.path.basename(tdb), "exists", 0
    if tmax <= 0: tmax = 4000 if set(elements_from_name(tdb)) & REFRACTORY else 2500
    from pycalphad import Database, variables as v
    from pycalphad.mapping import BinaryStrategy
    els = elements_from_name(tdb); rec = {"tdb": os.path.basename(tdb), "elements": els}
    try:
        raw = open(tdb, "rb").read()
        h = hashlib.sha256(raw).hexdigest()
        source_text, encoding, adjustments = prepare_source(raw)
        rec["source_adjustments"] = adjustments
        rec["parsed_text_sha256"] = hashlib.sha256(source_text.encode("utf-8")).hexdigest()
        rec.update(sha256=h, source_encoding=encoding,
                   engine="pycalphad " + __import__("pycalphad").__version__)
        if len(els) != 2:
            rec.update(status="not_applicable", error="Expected a binary assessment; filename does not identify two elements")
            json.dump(rec, open(os.path.join(out, os.path.basename(tdb).rsplit(".", 1)[0] + ".json"), "w"), separators=(",", ":"))
            return rec["tdb"], rec["status"], None
        db = Database(source_text); comps = [e.upper() for e in els] + (["VA"] if "VA" in db.elements else [])
        a, b = comps[0], comps[1]
        phases = sorted(db.phases.keys()); t0 = time.time()
        if source_defaults:
            phases, rec["phase_selection"] = source_phase_selection(source_text, phases)
        x_limits = (1e-5, 1-1e-5, 0.02) if interior_seeds else (0, 1, 0.02)
        s = BinaryStrategy(db, comps, phases, {v.X(b): x_limits, v.T: (300, tmax, 10), v.P: 101325, v.N: 1})
        if interior_seeds:
            # Use the supported custom-start API without changing the energy model.
            # A seed failure is recorded, not treated as successful equilibrium.
            attempts = []
            def seed_interior():
                for temperature in sorted(set(t for t in [400, 800, 1200, 1800, tmax-100] if 300 < t < tmax)):
                    for composition in [.05, .25, .5, .75, .95]:
                        trial = dict(T_K=temperature, x=composition)
                        try:
                            trial["added"] = s.add_nodes_from_conditions(
                                {v.X(b): composition, v.T: temperature, v.P: 101325, v.N: 1}, force_add=False)
                        except (ZeroDivisionError, ValueError, AssertionError) as exc:
                            trial["error"] = repr(exc)
                        attempts.append(trial)
            s.generate_automatic_starting_points = seed_interior
            rec["numerical_protocol"] = {"initialization": "interior_equilibria", "attempts": attempts,
                                         "x_range": list(x_limits[:2]), "pure_endpoints_excluded": True,
                                         "validation_status": "requires_independent_equilibrium_checks",
                                         "coverage": "not_established_by_seed_sampling"}
        s.do_map()
        regions = []
        for sd in s.get_tieline_data(v.X(b), v.T):
            regions.append([{"phase": p.phase, "x": r(p.x, 5), "T": r(p.y, 2)} for p in sd.data])
        inv = []
        for sd in s.get_invariant_data(v.X(b), v.T):
            pts = [{"phase": p.phase, "x": round(float(p.x), 5)} for p in sd.data]
            inv.append({"T": round(float(sd.data[0].y), 2), "points": sorted(pts, key=lambda q: q["x"])})
        rec.update(sha256=h, x_component=b, components=[a, b], T_range=[300, tmax], phases=phases,
                   regions=regions, invariants=sorted(inv, key=lambda q: -q["T"]), seconds=round(time.time() - t0, 1),
                   engine="pycalphad " + __import__("pycalphad").__version__, status="ok")
        validate_geometry(rec)
    except Exception as e:
        rec.pop("regions", None)
        rec.pop("invariants", None)
        rec.update(status="failed", error=f"{type(e).__name__}: {str(e)[:300]}")
    json.dump(rec, open(os.path.join(out, os.path.basename(tdb).rsplit(".", 1)[0] + ".json"), "w"), separators=(",", ":"))
    return rec["tdb"], rec["status"], rec.get("seconds")

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("out"); ap.add_argument("tdbs", nargs="+")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--tmax", type=float, default=0, help="0 = auto (4000 K with a refractory element, else 2500 K)")
    ap.add_argument("--interior-seeds", action="store_true", help="Experimental numerical recovery using interior equilibria; validate results independently")
    ap.add_argument("--respect-default-rejections", action="store_true", help="Use source DEFAULT_COM REJECT_PHASE domain, recorded in output")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    with mp.Pool(a.workers, maxtasksperchild=4) as pool:
        for res in pool.imap_unordered(one, [(t, a.out, a.tmax, a.interior_seeds, a.respect_default_rejections) for t in a.tdbs]):
            print(*res, flush=True)

if __name__ == "__main__":
    main()
