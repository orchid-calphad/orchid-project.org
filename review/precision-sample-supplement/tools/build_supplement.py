#!/usr/bin/env python3
"""Supplementary precision sample (PI approval 2026-10-24): ~150 items drawn from the FINAL ORCHID release, for the
student reviewers, alongside the frozen 300-value Dataset 1 sample (review/precision-sample, never changed here).

Three sub-samples, each stratified, drawn with one fixed seed (recorded):

  S-REG    phase-registry entries (Dataset 5, orchid-data/phases). Strata: structure tier (A standard name, B DFT
           unique match, C assessment-declared, liquid or gas, unassigned) x chemistry (metal-metal, or with a
           metalloid / non-metal). Reviewers check name and aliases, the structure (only where one is asserted), the TDB
           PHASE / CONSTITUENT statements against the file, and that the example DOIs refer to the phase.
  S-ESPEI  values of the ESPEI-ready datasets that come from the literature (orchid-data/espei: experimental/,
           model/, invariants/, dft/literature/; one item = one value of one dataset). Strata: data type (invariant
           ZPF, phase-boundary ZPF, formation enthalpy, mixing / partial / entropy, heat capacity, activity) x evidence
           class (experimental, model, dft). Reviewers check against the printed page the value, units, T / P /
           composition, phases, reference state, any conversion (e.g. bcc -> hcp Zr) and the uncertainty budget
           (printed sigma as printed, components taken correctly).
  S-ASM    assessed phase-assemblage points (Dataset 6, orchid-data/assessed-phase-assemblages; handbook diagrams).
           Strata: naming basis x boundary flag: the release's `boundary_basis` (printed_solid / mixed /
           dashed_estimated; frame_only when null) when the records carry it, else a heuristic (`dashed`: a dashed,
           i.e. estimated, boundary segment of the drawing lies within 8 raster cells of the field; `solid` otherwise). Same blind protocol and interface as the
           earlier expert sample (open-tdb build/okamoto-annotation/tier5/blind-sample): the reviewer names the
           tinted field from the drawing without seeing the machine's name, basis or neighbours; one field per
           diagram; the 48 diagrams of that sample and every diagram with an expert (tier-5) decision are excluded.

Allocation per sub-sample: the size is split over the primary strata proportionally (largest remainder) with a minimum
per primary stratum, then over the secondary levels proportionally with at least one item per non-empty level when
there is room. Candidates are sorted by id within a stratum cell and drawn with random.Random(seed) (S-ASM: shuffled,
then taken in order skipping diagrams already used). Items are then shuffled within each sub-sample (so that the order
does not reveal an S-ASM item's stratum) and numbered R001.., E001.., A001...

Packets: K reviewer packets (default 3). Packet k owns every item at position i with i % K == k (each packet mixes the
three sub-samples) and also gets, from the next packet's own items, round(10 %) of each sub-sample for inter-reviewer
agreement. The owner's answer is the one scored for precision; the overlap answers measure agreement.

Files in review/precision-sample-supplement/ (committed: sample.json, manifest.json, INSTRUCTIONS.md, tools/; local
only, .gitignore: index.html, packets/, img/, key.json, answers/):

  sample.json    what the page shows: S-REG and S-ESPEI items with their released content and stratum; S-ASM items
                 blind (no record id, machine name, basis or stratum). Public release data only; no page crops.
  key.json       full design: every item's stratum, source ids and (S-ASM) the hidden machine name; packet plan;
                 per-item view hashes for the drift check. LOCAL ONLY (blind key).
  img/           page crops: S-ESPEI the publisher PDF page of the value, S-ASM the handbook diagram (whole + zoom).
                 RESTRICTED: never committed, never public; sha256 recorded in manifest.json.
  index.html     self-contained review page (all items, crops embedded). LOCAL ONLY.
  packets/       one self-contained page per reviewer + packets.json. LOCAL ONLY.
  manifest.json  seed, design, allocation, population snapshot (sha256 of every release file read, orchid-data commit),
                 sha256 of sample.json, key.json, index.html, every packet and every image.

Commands (python with numpy and Pillow; poppler's pdftoppm on PATH):

  build_supplement.py draw [--out DIR] [--seed S] [--k 3] [--allow-dirty] [--no-images]
      draw, render and freeze in one step. Refuses a directory that already holds sample.json / key.json /
      manifest.json, and (without --allow-dirty) an orchid-data working tree with uncommitted changes.
  build_supplement.py check [--out DIR] [--drift-json PATH]
      verify every file against manifest.json and report, per item, whether the current release still carries the
      reviewed content (unchanged / changed / missing). Default command.

Serve: python3 -m http.server 8765 -d review/precision-sample-supplement  ->  http://localhost:8765/
(packets: http://localhost:8765/packets/ORCHID_supplement_packet_1.html). The pages are self-contained and also open
from disk. Results are scored by tools/score_supplement.py.
"""
import argparse
import base64
import gzip
import hashlib
import html
import io
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
KIT = HERE.parent
SITE = KIT.parent.parent
DATA = SITE.parent / "orchid-data"
ANNOT = SITE.parent / "open-tdb" / "build" / "okamoto-annotation"
BLIND_KEY = ANNOT / "tier5" / "blind-sample" / "key.json"
PDROOT = Path.home() / "data" / "PD-extraction"
HANDBOOK = PDROOT / "references" / "handbooks" / "okamoto-desk-handbook-2010.pdf"
PDF_ROOTS = ("references/pdfs", "references/pdfs-cited-elsevier", "references/pdfs-cited-publishers", "references/ill",
             "references/inbox", "references/books", "papers")
SEED = 20261024
SIZES = {"S-REG": 50, "S-ESPEI": 50, "S-ASM": 50}
MIN_PRIMARY = {"S-REG": 6, "S-ESPEI": 7, "S-ASM": 6}
PREFIX = {"S-REG": "R", "S-ESPEI": "E", "S-ASM": "A"}
K_DEFAULT = 3
OVERLAP_FRACTION = 0.10
DASH_RADIUS = 8
SCHEMA = "orchid-precision-supplement/1"
KEY_SCHEMA = "orchid-precision-supplement-key/1"
MANIFEST_SCHEMA = "orchid-precision-supplement-manifest/1"
RESULT_SCHEMA = "orchid-precision-supplement-review/1"
ESPEI_FOLDERS = ("experimental", "model", "invariants", "dft/literature")
TIER_ORDER = ("A_standard_name", "B_dft_unique_match", "C_assessment_declared", "liquid_or_gas", "unassigned")
CHEM_ORDER = ("metal-metal", "with_metalloid_or_nonmetal")
DTYPE_ORDER = ("invariant_zpf", "phase_boundary_zpf", "formation_enthalpy", "mixing_partial_entropy", "heat_capacity",
               "activity")
EV_ORDER = ("experimental", "model", "dft")
BASIS_ORDER = ("deterministic_printed_label", "deterministic_phase_rules", "reader_consensus", "deterministic_adjacency",
               "adjudicator", "dft_hull_corroborated")
DASH_ORDER = ("printed_solid", "mixed", "dashed_estimated", "frame_only", "solid", "dashed")
# metalloids and non-metals (the registry's chemistry split for S-REG)
NONMETALS = {"H", "B", "C", "N", "O", "F", "SI", "P", "S", "CL", "GE", "AS", "SE", "BR", "SB", "TE", "I"}


def canonical(obj):
    """Sorted keys, 2-space indent, UTF-8, trailing newline (same bytes as the page's canonicalJSON)."""
    return json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False) + "\n"


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def view_sha(view):
    return sha256_bytes(json.dumps(view, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8"))


# ---------------------------------------------------------------- allocation

def allocate(pop, n, minimum):
    """Proportional allocation (largest remainder) with a per-stratum minimum, capped at stratum size.
    Strata with population 0 get 0; ties are broken by stratum key."""
    pop = {k: v for k, v in pop.items() if v > 0}
    n = min(n, sum(pop.values()))
    fixed = {}
    while True:
        free = sorted(k for k in pop if k not in fixed)
        left = n - sum(fixed.values())
        if left < 0:
            raise ValueError("minimum per stratum exceeds the sample size")
        if not free:
            quota = {}
            break
        total = sum(pop[k] for k in free)
        quota = {k: left * pop[k] / total for k in free}
        floor = {k: min(minimum, pop[k]) for k in free}
        low = {k for k in free if quota[k] < floor[k] or quota[k] > pop[k]}
        if not low:
            break
        for k in low:
            fixed[k] = floor[k] if quota[k] < floor[k] else pop[k]
    base = {k: math.floor(q) for k, q in quota.items()}
    rest = n - sum(fixed.values()) - sum(base.values())
    for k in sorted(quota, key=lambda k: (-(quota[k] - base[k]), str(k)))[:rest]:
        base[k] += 1
    base.update(fixed)
    return base


def allocate2(cells, n, min_primary):
    """cells {(primary, secondary): population} -> {(primary, secondary): sampled}. Primary strata get a proportional
    share with a minimum; within a primary stratum the secondary levels are proportional with a minimum of one per
    non-empty level when the primary allocation has room for it."""
    prim = Counter()
    for (p, _), v in cells.items():
        prim[p] += v
    ap = allocate(dict(prim), n, min_primary)
    out = {c: 0 for c in cells}
    for p, k in ap.items():
        sub = {s: v for (pp, s), v in cells.items() if pp == p and v > 0}
        if not k:
            continue
        a = allocate(sub, k, 1 if k >= len(sub) else 0)
        for s, x in a.items():
            out[(p, s)] = x
    return out


def order_key(orders):
    def key(cell):
        return tuple((o.index(c) if c in o else len(o), str(c)) for c, o in zip(cell, orders))
    return key


# ---------------------------------------------------------------- populations

def reg_tier(p):
    st = p.get("structure") or {}
    if p.get("status") == "assigned":
        if st.get("basis") in ("liquid", "gas"):
            return "liquid_or_gas"
        return {"A": "A_standard_name", "B": "B_dft_unique_match", "C": "C_assessment_declared"}.get(st.get("tier"),
                                                                                                    "unassigned")
    return "unassigned"


def reg_chem(elements):
    return "with_metalloid_or_nonmetal" if any(e.upper() in NONMETALS for e in elements) else "metal-metal"


def reg_population(data):
    reg = json.loads((data / "phases" / "registry.json").read_text(encoding="utf-8"))
    cand = defaultdict(list)
    for p in reg["phases"]:
        cand[(reg_tier(p), reg_chem(p["elements"]))].append(p["id"])
    return cand, {p["id"]: p for p in reg["phases"]}


def espei_dtype(folder, output):
    if folder == "invariants":
        return "invariant_zpf"
    if output == "ZPF":
        return "phase_boundary_zpf"
    if output == "HM_FORM":
        return "formation_enthalpy"
    if output.startswith("ACR"):
        return "activity"
    if output == "CPM":
        return "heat_capacity"
    return "mixing_partial_entropy"


def espei_value_ids(d):
    o = d.get("orchid") or {}
    obs = o.get("observations")
    if obs:
        return [x.get("id") for x in obs]
    return [m.get("id") for m in o.get("mapping") or []]


def espei_population(data):
    """{(dtype, evidence): [(path, index, orchid id)]} over the literature-derived ESPEI datasets; and a digest of the
    dataset files read."""
    ix = json.loads((data / "espei" / "index.json").read_text(encoding="utf-8"))
    cand = defaultdict(list)
    dig = hashlib.sha256()
    for r in sorted(ix["datasets"], key=lambda r: r["path"]):
        if r["folder"] not in ESPEI_FOLDERS:
            continue
        b = (data / "espei" / r["path"]).read_bytes()
        dig.update(r["path"].encode() + b"\0" + hashlib.sha256(b).digest())
        d = json.loads(b)
        ev = (d.get("orchid") or {}).get("evidence_class") or r.get("evidence_class")
        dt = espei_dtype(r["folder"], d["output"])
        for i, oid in enumerate(espei_value_ids(d)):
            cand[(dt, ev)].append(("%s#%d" % (r["path"], i), r["path"], i, oid))
    return cand, dig.hexdigest()


def load_grid(sd):
    import numpy as np
    g = np.load(sd / "grid.npz")
    return g["lab"], g["x"], g["T"]


def dashed_fields(fields, grid):
    """Field ids within DASH_RADIUS raster cells of a dashed (estimated) boundary segment of the drawing."""
    import numpy as np
    lab, xs, Ts = grid
    out = set()
    for j in fields.get("endpoint_joins") or []:
        if not j.get("dashed"):
            continue
        for x, T in j.get("points") or []:
            i = int(np.searchsorted(xs, x))
            t = int(np.searchsorted(Ts, T))
            win = lab[max(0, t - DASH_RADIUS):t + DASH_RADIUS + 1, max(0, i - DASH_RADIUS):i + DASH_RADIUS + 1]
            out.update(int(v) for v in np.unique(win) if v > 0)
    return out


def asm_population(data, annot, blind_key):
    """{(basis, dashed|solid): [record]} after the exclusions; exclusion counts; per-diagram field geometry cache."""
    recs = [json.loads(l) for l in gzip.decompress((data / "assessed-phase-assemblages" / "records.jsonl.gz").read_bytes())
            .decode("utf-8").splitlines() if l.strip()]
    blind = json.loads(Path(blind_key).read_text(encoding="utf-8"))["items"]
    blind_diagrams = {v["diagram"] for v in blind.values()}
    expert_diagrams = {r["diagram"] for r in recs if r.get("review_state") != "not reviewed by the expert"}
    excl = Counter()
    keep = []
    for r in recs:
        if r["diagram"] in blind_diagrams:
            excl["diagram_in_earlier_blind_sample"] += 1
        elif r["diagram"] in expert_diagrams:
            excl["diagram_with_expert_tier5_decision"] += 1
        else:
            keep.append(r)
    cand = defaultdict(list)
    if any("boundary_basis" in r for r in keep):
        # the release's own flag (James Saal review 2026-10-24): boundary lines bounding the field along the point's
        # isotherm and isopleth are printed_solid, dashed_estimated or mixed; null = frame edges only
        for r in keep:
            cand[(r["basis"], r.get("boundary_basis") or "frame_only")].append(r)
        info = {"population_released": len(recs), "excluded": dict(excl), "blind_sample_diagrams": len(blind_diagrams),
                "expert_diagrams": len(expert_diagrams), "boundary_flag": "release field boundary_basis"}
        return cand, info
    dashed = {}
    for d in sorted({r["diagram"] for r in keep}):
        run = {r["method"]["annotation_run"] for r in keep if r["diagram"] == d}
        assert len(run) == 1, d
        sd = Path(annot) / run.pop() / "systems" / d
        fields = json.loads((sd / "fields.json").read_text(encoding="utf-8"))
        if any(j.get("dashed") for j in fields.get("endpoint_joins") or []):
            dashed[d] = dashed_fields(fields, load_grid(sd))
        else:
            dashed[d] = set()
    for r in keep:
        cand[(r["basis"], "dashed" if r["field"] in dashed[r["diagram"]] else "solid")].append(r)
    info = {"population_released": len(recs), "excluded": dict(excl), "blind_sample_diagrams": len(blind_diagrams),
            "expert_diagrams": len(expert_diagrams),
            "boundary_flag": "heuristic: dashed when a dashed boundary segment lies within %d raster cells of the field "
                             "(the release has no boundary_basis field)" % DASH_RADIUS}
    return cand, info


# ---------------------------------------------------------------- item views

def paper_info(papers, pid):
    p = papers.get(pid) or {}
    return {"id": pid, "citation": p.get("citation"), "doi": p.get("doi"), "doi_url": p.get("doi_url"),
            "title": p.get("title"), "year": p.get("year")}


def reg_view(full, seed):
    """Reviewed content of a registry entry: names, structure, one TDB model (drawn per item), first 3 references."""
    models = sorted(full.get("models") or [], key=lambda m: (m.get("tdb") or "", m.get("phase_name") or ""))
    model = None
    if models:
        m = random.Random("%s:%s" % (seed, full["id"])).choice(models)
        td = m.get("tdb_definition") or {}
        model = {"tdb": m.get("tdb"), "path": m.get("path"), "sha256": m.get("sha256"), "phase_name": m.get("phase_name"),
                 "sublattices": m.get("sublattices"), "stoichiometric": m.get("stoichiometric"),
                 "verbatim": td.get("verbatim"), "provenance": td.get("provenance"),
                 "tdb_comment_structure": m.get("tdb_comment_structure")}
    st = full.get("structure")
    if st:
        st = {k: st.get(k) for k in ("structure_type", "basis", "tier", "confidence", "label", "prototype",
                                     "strukturbericht", "pearson", "space_group", "caveats", "source")}
    names = [{k: n.get(k) for k in ("name", "namespace", "label", "basis") if n.get(k) is not None}
             for n in full.get("names") or []]
    refs = [{k: r.get(k) for k in ("doi", "citation", "roles", "labels", "values")} for r in (full.get("references") or [])[:3]]
    return {"phase_id": full["id"], "system": full["system"], "elements": full.get("elements"),
            "preferred_name": full.get("preferred_name"), "status": full.get("status"), "structure": st,
            "candidates": (full.get("candidates") or [])[:6], "review_reasons": full.get("review_reasons") or [],
            "names": names, "n_models": len(models), "model": model, "references": refs,
            "n_references": len(full.get("references") or []), "total_papers": full.get("total_papers"),
            "stoichiometry": full.get("stoichiometry")}


def flat(v):
    if isinstance(v, list):
        out = []
        for x in v:
            out.extend(flat(x))
        return out
    return [v]


def espei_value(d, i):
    """(value, conditions at the value, solver configuration) of value i of an ESPEI dataset."""
    c = d.get("conditions") or {}
    if d["output"] == "ZPF":
        T = c.get("T")
        cond = {"P": c.get("P"), "T": T[i] if isinstance(T, list) else T}
        return d["values"][i], cond, None
    Ps = c.get("P") if isinstance(c.get("P"), list) else [c.get("P")]
    Ts = c.get("T") if isinstance(c.get("T"), list) else [c.get("T")]
    vals = flat(d["values"])
    nX = len(vals) // (len(Ps) * len(Ts))
    p, t, x = i // (len(Ts) * nX), (i // nX) % len(Ts), i % nX
    cond = {"P": Ps[p], "T": Ts[t]}
    for k, v in c.items():
        if k.startswith("X_"):
            cond[k] = v[x] if isinstance(v, list) else v
    solver = None
    if d.get("solver"):
        s = d["solver"]
        cfg = s.get("sublattice_configurations") or []
        occ = s.get("sublattice_occupancies") or []
        solver = {"mode": s.get("mode"), "sublattice_site_ratios": s.get("sublattice_site_ratios"),
                  "configuration": cfg[x] if x < len(cfg) else None, "occupancies": occ[x] if x < len(occ) else None}
    return vals[i], cond, solver


LIT_FIELDS = ("id", "system", "property", "property_printed", "property_group", "value", "value_printed", "value_range",
              "value_expression", "unit", "unit_printed", "unit_basis", "temperature_K", "temperature_printed",
              "temperature_basis", "pressure_Pa", "pressure_basis", "overall_composition", "component", "phases",
              "phase_compositions", "reaction", "reference_state", "uncertainty", "uncertainty_budget", "locator",
              "quote", "quote_page", "corrections", "evidence_class", "origin", "cited_sources", "record_kind", "paper")
NOTE_KEYS = ("orchid_value_note", "orchid_standard_temperature", "orchid_temperature_note", "orchid_reference_state")


def espei_view(d, path, i, lit):
    value, cond, solver = espei_value(d, i)
    o = d.get("orchid") or {}
    ou = d.get("orchid_uncertainty") or {}
    uv = (ou.get("values") or [])
    obs = (o.get("observations") or o.get("mapping") or [])
    w = d.get("weight")
    notes = {k: d[k] for k in NOTE_KEYS if d.get(k) is not None}
    if o.get("temperature_note"):
        notes["dft_temperature_note"] = o["temperature_note"]
    m = re.search(r"eV/atom x [0-9.]+ = J/mol-atom", d.get("comment") or "")
    if m:
        notes["unit_conversion"] = m.group(0)
    rec = {k: lit.get(k) for k in LIT_FIELDS if lit.get(k) is not None} if lit else None
    return {"dataset": "espei/" + path, "value_index": i, "output": d["output"], "components": d.get("components"),
            "phases": d.get("phases"), "unit": ou.get("unit"), "value": value, "conditions": cond, "solver": solver,
            "reference_state": d.get("reference_state"), "reference_states": d.get("reference_states"),
            "tags": d.get("tags"), "reference": d.get("reference"), "notes": notes,
            "observation": obs[i] if i < len(obs) else None,
            "uncertainty": uv[i] if i < len(uv) else None,
            "uncertainty_rules": {k: ou.get(k) for k in ("default_sigma", "sigma_rule", "weight_rule")},
            "weight": w[i] if isinstance(w, list) else w, "record": rec}


def espei_reviewed(v):
    """The part of an S-ESPEI view compared with the current release by `check`."""
    return {k: v.get(k) for k in ("value", "conditions", "solver", "phases", "reference_state", "reference_states",
                                  "unit", "notes", "uncertainty", "observation")}


# ---------------------------------------------------------------- draw

def git_state(repo):
    try:
        head = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"], capture_output=True, text=True,
                              check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain"], capture_output=True, text=True,
                               check=True).stdout.strip()
        return head, bool(dirty)
    except (OSError, subprocess.CalledProcessError):
        return None, None


def draw_cells(cand, n, min_primary, orders, rng, sort_key, pick=None):
    cells = {c: len(v) for c, v in cand.items()}
    alloc = allocate2(cells, n, min_primary)
    chosen = []
    for c in sorted(cells, key=order_key(orders)):
        k = alloc.get(c, 0)
        if not k:
            continue
        pool = sorted(cand[c], key=sort_key)
        got = pick(pool, k) if pick else rng.sample(pool, k)
        if len(got) < k:
            alloc[c] = len(got)
        chosen += [(c, x) for x in got]
    strata = [{"stratum": "%s|%s" % c, "primary": c[0], "secondary": c[1], "population": cells[c],
               "sampled": alloc.get(c, 0)} for c in sorted(cells, key=order_key(orders))]
    return chosen, strata


def phase_buttons(fields):
    printed = {t.strip() for l in fields.get("diagram_labels") or [] for t in l.split("+")}
    phases = list(fields.get("phase_table") or {}) + ["L", "L#2", "G"]
    phases += sorted(p for p in printed if p and p not in phases and len(p) < 20 and any(c.isalpha() for c in p)
                     and p not in ("L1", "L2", "?"))
    return phases


def build(data=DATA, annot=ANNOT, blind_key=BLIND_KEY, seed=SEED, k=K_DEFAULT, sizes=None):
    """-> (sample, key, population snapshot). Deterministic for the same release files and seed."""
    sizes = sizes or SIZES
    rng = random.Random(seed)
    papers = json.loads((data / "literature" / "papers.json").read_text(encoding="utf-8"))
    snap = {"orchid_data": str(data), "files_sha256": {}}
    for rel in ("phases/registry.json", "espei/index.json", "literature/orchid-literature.jsonl.gz", "literature/papers.json",
                "assessed-phase-assemblages/records.jsonl.gz"):
        snap["files_sha256"][rel] = sha256_file(data / rel)
    snap["blind_key_sha256"] = sha256_file(blind_key)

    # S-REG
    rc, reg_compact = reg_population(data)
    reg_chosen, reg_strata = draw_cells(rc, sizes["S-REG"], MIN_PRIMARY["S-REG"], (TIER_ORDER, CHEM_ORDER), rng,
                                        sort_key=lambda x: x)
    sysfiles = {}
    reg_items = []
    for c, pid in reg_chosen:
        det = reg_compact[pid]["detail"]
        if det not in sysfiles:
            b = (data / "phases" / det).read_bytes()
            snap["files_sha256"]["phases/" + det] = sha256_bytes(b)
            sysfiles[det] = {p["id"]: p for p in json.loads(b)["phases"]}
        v = reg_view(sysfiles[det][pid], seed)
        reg_items.append({"kind": "S-REG", "stratum": "%s|%s" % c, "source_id": pid, "view": v})

    # S-ESPEI
    ec, digest = espei_population(data)
    snap["espei_literature_datasets_digest"] = digest
    esp_chosen, esp_strata = draw_cells(ec, sizes["S-ESPEI"], MIN_PRIMARY["S-ESPEI"], (DTYPE_ORDER, EV_ORDER), rng,
                                        sort_key=lambda x: x[0])
    want = {x[3] for _, x in esp_chosen}
    lit = {}
    for line in gzip.decompress((data / "literature" / "orchid-literature.jsonl.gz").read_bytes()).decode("utf-8").splitlines():
        if line.strip():
            r = json.loads(line)
            if r["id"] in want:
                lit[r["id"]] = r
            for mid in r.get("merged_ids") or []:
                if mid in want and mid not in lit:
                    lit[mid] = r
    esp_items = []
    for c, (uid, path, i, oid) in esp_chosen:
        d = json.loads((data / "espei" / path).read_text(encoding="utf-8"))
        v = espei_view(d, path, i, lit.get(oid))
        v["paper"] = paper_info(papers, (lit.get(oid) or {}).get("paper"))
        esp_items.append({"kind": "S-ESPEI", "stratum": "%s|%s" % c, "source_id": uid, "orchid_id": oid, "view": v})

    # S-ASM (one field per diagram across the sub-sample)
    ac, asm_info = asm_population(data, annot, blind_key)
    used = set()

    def pick(pool, kk):
        pool = pool[:]
        rng.shuffle(pool)
        got = []
        for r in pool:
            if len(got) == kk:
                break
            if r["diagram"] in used:
                continue
            used.add(r["diagram"])
            got.append(r)
        return got

    asm_chosen, asm_strata = draw_cells(ac, sizes["S-ASM"], MIN_PRIMARY["S-ASM"], (BASIS_ORDER, DASH_ORDER), rng,
                                        sort_key=lambda r: r["id"], pick=pick)
    asm_items = []
    for c, r in asm_chosen:
        sd = Path(annot) / r["method"]["annotation_run"] / "systems" / r["diagram"]
        fields = json.loads((sd / "fields.json").read_text(encoding="utf-8"))
        cal = json.loads((sd / "calibrated.json").read_text(encoding="utf-8"))
        import numpy as np
        lab, xs, Ts = load_grid(sd)
        jj, ii = np.nonzero(lab == r["field"])
        blind_view = {"system": r["system"], "diagram": r["diagram"], "handbook_page": cal.get("page"),
                      "axis_element": cal.get("component"), "composition_axis": r.get("composition_axis"),
                      "T_range_C": [round(float(Ts[jj].min())), round(float(Ts[jj].max()))] if len(jj) else None,
                      "x_range_at_pct": [round(float(xs[ii].min()) * 100, 1), round(float(xs[ii].max()) * 100, 1)] if len(ii) else None,
                      "point": {"T_C": round(r["T_C"], 1), "x_at_pct": round(r["x"] * 100, 2), "x_component": r["x_component"]},
                      "phase_buttons": phase_buttons(fields)}
        secret = {"record_id": r["id"], "field": r["field"], "machine_name": r["assemblage_label"],
                  "canonical_names": [p.get("canonical_name") for p in r.get("assemblage") or []],
                  "basis": r["basis"], "confidence": r.get("confidence"), "boundary": c[1], "flags": r.get("flags"),
                  "T_K": r["T_K"], "x": r["x"], "annotation_run": r["method"]["annotation_run"]}
        asm_items.append({"kind": "S-ASM", "stratum": "%s|%s" % c, "source_id": r["id"], "view": blind_view,
                          "secret": secret})

    # order, number, packets
    items = []
    for sub, lst in (("S-REG", reg_items), ("S-ESPEI", esp_items), ("S-ASM", asm_items)):
        lst = lst[:]
        rng.shuffle(lst)
        for j, it in enumerate(lst, 1):
            it["item"] = "%s%03d" % (PREFIX[sub], j)
        items += lst
    for n, it in enumerate(items, 1):
        it["n"] = n
    plan = packet_plan(items, k)

    design = {"seed": seed, "sizes": sizes, "min_per_primary_stratum": MIN_PRIMARY,
              "allocation": "per sub-sample: primary strata proportional (largest remainder) with a minimum; secondary "
                            "levels proportional within a primary stratum, at least one per non-empty level when there is room",
              "draw": "candidates sorted by id within a stratum cell; random.Random(seed).sample in fixed stratum order "
                      "(S-ASM: shuffle, then one field per diagram); items shuffled within each sub-sample, then numbered",
              "strata": {"S-REG": "structure tier x chemistry", "S-ESPEI": "data type x evidence class",
                         "S-ASM": "naming basis x boundary flag (%s)" % asm_info["boundary_flag"]},
              "packets": {"k": k, "own": "position i with i %% k == packet-1",
                          "overlap": "round(%.0f%%) of each sub-sample's own items of the next packet, evenly spaced"
                                     % (100 * OVERLAP_FRACTION)},
              "s_asm_exclusions": asm_info,
              "s_asm_protocol": "blind: the page shows the drawing with the field tinted and the released point marked; "
                                "no machine name, basis, stratum or record id"}
    strata = {"S-REG": reg_strata, "S-ESPEI": esp_strata, "S-ASM": asm_strata}

    def public(it):
        o = {"item": it["item"], "n": it["n"], "kind": it["kind"], "view": it["view"]}
        if it["kind"] != "S-ASM":
            o["stratum"] = it["stratum"]
        return o

    sample = {"schema": SCHEMA,
              "description": "Supplementary precision sample of the final ORCHID release (phase registry, ESPEI-ready "
                             "values, assessed phase assemblages) for a human precision check.",
              "design": {k_: design[k_] for k_ in ("seed", "sizes", "strata", "packets")},
              "strata": {s: v for s, v in strata.items() if s != "S-ASM"},
              "items": [public(it) for it in items]}
    key_items = []
    for it in items:
        o = {"item": it["item"], "n": it["n"], "kind": it["kind"], "stratum": it["stratum"], "source_id": it["source_id"],
             "owner_packet": plan["owner"][it["item"]], "overlap_packets": plan["overlap"].get(it["item"], [])}
        if it["kind"] == "S-ASM":
            o["secret"] = it["secret"]
            o["reviewed_sha256"] = view_sha({k_: it["secret"][k_] for k_ in ("machine_name", "T_K", "x")})
        elif it["kind"] == "S-ESPEI":
            o["orchid_id"] = it["orchid_id"]
            o["reviewed_sha256"] = view_sha(espei_reviewed(it["view"]))
        else:
            o["reviewed_sha256"] = view_sha(it["view"])
        key_items.append(o)
    key = {"schema": KEY_SCHEMA, "note": "LOCAL ONLY: holds the hidden machine names of the blind S-ASM items.",
           "design": design, "strata": strata, "packets": plan["packets"], "items": key_items}
    return sample, key, snap


def packet_plan(items, k):
    own = [[it for i, it in enumerate(items) if i % k == p] for p in range(k)]
    owner = {it["item"]: p + 1 for p in range(k) for it in own[p]}
    overlap = defaultdict(list)
    packets = []
    for p in range(k):
        nxt = own[(p + 1) % k] if k > 1 else []
        extra = []
        for sub in ("S-REG", "S-ESPEI", "S-ASM"):
            lst = [it for it in nxt if it["kind"] == sub]
            m = min(len(lst), max(1, int(round(OVERLAP_FRACTION * len(lst))))) if lst else 0
            extra += [lst[int((j + 0.5) * len(lst) / m)] for j in range(m)] if m else []
        for it in extra:
            overlap[it["item"]].append(p + 1)
        packets.append({"packet": p + 1, "file": "ORCHID_supplement_packet_%d.html" % (p + 1),
                        "own": [it["item"] for it in own[p]], "overlap_from_packet": (p + 1) % k + 1 if k > 1 else None,
                        "overlap": [it["item"] for it in extra]})
    return {"owner": owner, "overlap": dict(overlap), "packets": packets}


# ---------------------------------------------------------------- images

def pdf_index(root=PDROOT):
    idx = {}
    for sub in PDF_ROOTS:
        d = Path(root) / sub
        if not d.is_dir():
            continue
        for p in sorted(d.glob("*.pdf")):
            b = p.stem
            m = re.match(r"^\d{4}_(.*)$", b)
            idx.setdefault((m.group(1) if m else b).lower(), p)
    return idx


def to_jpeg(png, out, maxw, quality):
    from PIL import Image
    im = Image.open(png).convert("RGB")
    im.thumbnail((maxw, int(maxw * 1.6)))
    out.parent.mkdir(parents=True, exist_ok=True)
    im.save(out, "JPEG", quality=quality, optimize=True, progressive=True)


def render_page(pdf, page, out, cache):
    pref = Path(cache) / ("%s-p%d" % (sha256_bytes(str(pdf).encode())[:16], page))
    if not pref.with_suffix(".png").exists():
        subprocess.run(["pdftoppm", "-f", str(page), "-l", str(page), "-r", "130", "-png", "-singlefile", str(pdf),
                        str(pref)], check=True, capture_output=True)
    to_jpeg(pref.with_suffix(".png"), out, 1150, 62)


def render_asm(sd, field, point, out_whole, out_zoom, cache, handbook=HANDBOOK):
    """Handbook diagram with the field tinted orange and the released point marked (black ring); whole + zoom
    (same geometry as the earlier blind sample)."""
    import numpy as np
    from PIL import Image, ImageDraw
    cal = json.loads((sd / "calibrated.json").read_text(encoding="utf-8"))
    cx, cy = cal["calibration"]["x"], cal["calibration"]["y"]
    S, H = 300 / 72, 792.0
    pref = Path(cache) / ("hb-p%d" % cal["page"])
    if not pref.with_suffix(".png").exists():
        subprocess.run(["pdftoppm", "-f", str(cal["page"]), "-l", str(cal["page"]), "-r", "300", "-png", "-singlefile",
                        str(handbook), str(pref)], check=True, capture_output=True)
    base = Image.open(pref.with_suffix(".png")).convert("RGBA")
    lab, xs, Ts = load_grid(sd)

    def to(x, T):
        return (((x * 100 - cx["b"]) / cx["a"]) * S, (H - (T - cy["b"]) / cy["a"]) * S)

    jj, ii = np.nonzero(lab == field)
    ov = Image.new("RGBA", base.size, (0, 0, 0, 0))
    dr = ImageDraw.Draw(ov)
    for j, i_ in zip(jj, ii):
        X, Y = to(xs[i_], Ts[j])
        dr.rectangle((X - 2, Y - 2, X + 2, Y + 2), fill=(230, 60, 0, 110))
    PX, PY = to(point[0], point[1])
    for rad, wdt in ((13, 4),):
        dr.ellipse((PX - rad, PY - rad, PX + rad, PY + rad), outline=(0, 0, 0, 255), width=wdt)
    im = Image.alpha_composite(base, ov).convert("RGB")
    x0, y0, x1, y1 = cal["frame_pdf"]
    whole = im.crop((int((x0 - 45) * S), int((H - y1 - 25) * S), int((x1 + 15) * S), int((H - y0 + 40) * S)))
    xa, xb, Ta, Tb = xs[ii].min(), xs[ii].max(), Ts[jj].min(), Ts[jj].max()
    mx, mT = max(0.08, (xb - xa) * 0.6), max(0.15 * (Ts[-1] - Ts[0]), (Tb - Ta) * 0.6)
    (X0, Y0), (X1, Y1) = to(max(xs[0], xa - mx), min(Ts[-1], Tb + mT)), to(min(xs[-1], xb + mx), max(Ts[0], Ta - mT))
    z = im.crop((int(X0), int(Y0), int(X1), int(Y1)))
    sc = min(4.0, 900 / max(z.width, 1))
    z = z.resize((max(1, int(z.width * sc)), max(1, int(z.height * sc))))
    tmp = Path(cache)
    whole.save(tmp / "w.png")
    z.save(tmp / "z.png")
    to_jpeg(tmp / "w.png", out_whole, 760, 60)
    to_jpeg(tmp / "z.png", out_zoom, 620, 62)


def render_images(sample, key, out, annot=ANNOT, data=DATA):
    """Write img/<item>-*.jpg; return {item: {slot: relative path}} and the list of items without an image."""
    cache = Path(tempfile.mkdtemp(prefix="orchid-supp-render-"))
    idx = pdf_index()
    imgs, missing = {}, []
    secrets = {k["item"]: k for k in key["items"]}
    try:
        for it in sample["items"]:
            v = it["view"]
            if it["kind"] == "S-ESPEI":
                rec = v.get("record") or {}
                L = rec.get("locator") or {}
                page = rec.get("quote_page") if L.get("artifact") != "table" and rec.get("quote_page") else L.get("page")
                pid = rec.get("paper") or ""
                doi = pid[4:] if pid.startswith("doi:") else pid
                pdf = idx.get(doi.lower().replace("/", "_"))
                if not (pdf and page):
                    missing.append(it["item"])
                    continue
                rel = "img/%s-page.jpg" % it["item"]
                render_page(pdf, int(page), out / rel, cache)
                imgs[it["item"]] = {"page": rel, "page_number": int(page)}
            elif it["kind"] == "S-ASM":
                s = secrets[it["item"]]["secret"]
                sd = Path(annot) / s["annotation_run"] / "systems" / v["diagram"]
                pt = (s["x"], s["T_K"] - 273.15)
                rw, rz = "img/%s-whole.jpg" % it["item"], "img/%s-zoom.jpg" % it["item"]
                render_asm(sd, s["field"], pt, out / rw, out / rz, cache)
                imgs[it["item"]] = {"whole": rw, "zoom": rz}
    finally:
        shutil.rmtree(cache, ignore_errors=True)
    return imgs, missing


# ---------------------------------------------------------------- page

def instructions_html(md):
    """Tiny Markdown subset (headings, bullets, numbered lists, paragraphs, **bold**, `code`) -> HTML."""
    out, para, lst = [], [], None

    def inline(t):
        t = html.escape(t)
        t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
        return re.sub(r"`(.+?)`", r"<code>\1</code>", t)

    def flush():
        nonlocal para, lst
        if para:
            out.append("<p>" + inline(" ".join(para)) + "</p>")
            para = []
        if lst:
            out.append("</%s>" % lst)
            lst = None

    for line in md.splitlines():
        s = line.strip()
        m = re.match(r"^(#{1,3})\s+(.*)$", s)
        b = re.match(r"^[-*]\s+(.*)$", s)
        o = re.match(r"^\d+\.\s+(.*)$", s)
        if not s:
            flush()
        elif m:
            flush()
            lvl = min(4, len(m.group(1)) + 2)
            out.append("<h%d>%s</h%d>" % (lvl, inline(m.group(2)), lvl))
        elif b or o:
            if para:
                out.append("<p>" + inline(" ".join(para)) + "</p>")
                para = []
            tag = "ul" if b else "ol"
            if lst != tag:
                if lst:
                    out.append("</%s>" % lst)
                out.append("<%s>" % tag)
                lst = tag
            out.append("<li>" + inline((b or o).group(1)) + "</li>")
        elif lst and line.startswith("  "):
            out[-1] = out[-1][:-5] + " " + inline(s) + "</li>"
        else:
            if lst:
                out.append("</%s>" % lst)
                lst = None
            para.append(s)
    flush()
    return "\n".join(out)


PAGE = (HERE / "page_template.html")


def embed_json(obj):
    return json.dumps(obj, ensure_ascii=False, sort_keys=True).replace("</", "<\\/").replace("<!--", "<\\!--")


def data_uri(p):
    return "data:image/jpeg;base64," + base64.b64encode(Path(p).read_bytes()).decode("ascii")


def render_html(sample_text, sample, items, imgs, out, title, packet=None, instructions=""):
    sha = sha256_bytes(sample_text.encode("utf-8"))
    sub = dict(sample, items=items)
    if packet:
        sub["packet"] = packet
    images = {}
    for it in items:
        for slot, rel in (imgs.get(it["item"]) or {}).items():
            if slot == "page_number":
                continue
            images.setdefault(it["item"], {})[slot] = data_uri(out / rel)
    tpl = PAGE.read_text(encoding="utf-8")
    return (tpl.replace("__TITLE__", html.escape(title)).replace("__INSTRUCTIONS__", instructions)
            .replace("__SAMPLE__", embed_json(sub)).replace("__IMAGES__", embed_json(images)).replace("__SHA__", sha))


def write_pages(out, sample_text, sample, key, imgs, instructions_md):
    instr = instructions_html(instructions_md)
    by = {it["item"]: it for it in sample["items"]}
    (out / "index.html").write_text(render_html(sample_text, sample, sample["items"], imgs, out,
                                                "Precision sample supplement", instructions=instr), encoding="utf-8")
    pk = out / "packets"
    pk.mkdir(exist_ok=True)
    k = len(key["packets"])
    plan = []
    for p in key["packets"]:
        its = [by[i] for i in p["own"] + p["overlap"]]
        info = {"index": p["packet"], "of": k, "own": len(p["own"]), "overlap": len(p["overlap"]),
                "overlap_from_packet": p["overlap_from_packet"]}
        (pk / p["file"]).write_text(render_html(sample_text, sample, its, imgs, out,
                                                "ORCHID supplement packet %d of %d" % (p["packet"], k), packet=info,
                                                instructions=instr), encoding="utf-8")
        plan.append(dict(info, file=p["file"], items=len(its)))
    (pk / "packets.json").write_text(canonical({"k": k, "packets": plan}), encoding="utf-8")


def manifest(out, sample_bytes, key_bytes, snap, data_commit, data_dirty, site_commit, missing):
    files = {"sample.json": sha256_bytes(sample_bytes), "key.json": sha256_bytes(key_bytes)}
    for rel in ("index.html", "INSTRUCTIONS.md"):
        if (out / rel).exists():
            files[rel] = sha256_file(out / rel)
    for d in ("packets", "img"):
        if (out / d).is_dir():
            for p in sorted((out / d).iterdir()):
                if p.is_file():
                    files["%s/%s" % (d, p.name)] = sha256_file(p)
    key = json.loads(key_bytes)
    counts = {s: [{k: x[k] for k in ("stratum", "population", "sampled")} for x in v] for s, v in key["strata"].items()}
    return {"schema": MANIFEST_SCHEMA,
            "note": "Freeze record of the supplementary precision sample. Local-only files (index.html, packets/, img/, "
                    "key.json) are restricted or blind and never committed; their sha256 is recorded here. "
                    "tools/build_supplement.py check verifies them.",
            "sample_sha256": files["sample.json"], "key_sha256": files["key.json"], "files_sha256": files,
            "design": key["design"], "strata": counts,
            "population": dict(snap, orchid_data_commit=data_commit, orchid_data_dirty=data_dirty,
                               orchid_project_org_commit=site_commit),
            "items_without_page_image": missing,
            "n_items": len(key["items"]),
            "n_items_by_subsample": dict(Counter(it["kind"] for it in key["items"]))}


def cmd_draw(a):
    out = a.out
    if any((out / f).exists() for f in ("sample.json", "key.json", "manifest.json")):
        print("refusing to draw into %s: a sample exists there (samples are frozen; draw into a new directory)" % out,
              file=sys.stderr)
        return 1
    commit, dirty = git_state(a.data)
    if dirty and not a.allow_dirty:
        print("refusing to draw: %s has uncommitted changes (the sample must come from the committed final release; "
              "--allow-dirty for a dry run)" % a.data, file=sys.stderr)
        return 1
    site_commit, _ = git_state(SITE)
    sample, key, snap = build(a.data, a.annot, a.blind_key, a.seed, a.k)
    out.mkdir(parents=True, exist_ok=True)
    imgs, missing = ({}, []) if a.no_images else render_images(sample, key, out, a.annot, a.data)
    for it in sample["items"]:
        if it["item"] in imgs and "page_number" in imgs[it["item"]]:
            it["view"]["page_image_page"] = imgs[it["item"]]["page_number"]
    sample_text = canonical(sample)
    key_text = canonical(key)
    (out / "sample.json").write_text(sample_text, encoding="utf-8")
    (out / "key.json").write_text(key_text, encoding="utf-8")
    ip = out / "INSTRUCTIONS.md" if (out / "INSTRUCTIONS.md").exists() else KIT / "INSTRUCTIONS.md"
    instr = ip.read_text(encoding="utf-8") if ip.exists() else ""
    write_pages(out, sample_text, sample, key, imgs, instr)
    man = manifest(out, sample_text.encode("utf-8"), key_text.encode("utf-8"), snap, commit, dirty, site_commit, missing)
    (out / "manifest.json").write_text(canonical(man), encoding="utf-8")
    print("supplementary sample in %s: sample sha256 %s, %d items, seed %d, %d packets; orchid-data %s%s" % (
        out, man["sample_sha256"][:16], len(sample["items"]), a.seed, a.k, (commit or "?")[:9],
        " (DIRTY: dry run)" if dirty else ""))
    for s, v in key["strata"].items():
        print("  %s" % s)
        for x in v:
            print("    %-48s population %6d  sampled %3d" % (x["stratum"], x["population"], x["sampled"]))
    if missing:
        print("  items without a page image: %s" % ", ".join(missing))
    return 0


# ---------------------------------------------------------------- check

def current_views(key, data=DATA):
    """{item: (status, detail)} comparing each item's reviewed content with the current release."""
    out = {}
    reg = None
    asm = None
    sysfiles = {}
    for it in key["items"]:
        if it["kind"] == "S-REG":
            if reg is None:
                reg = {p["id"]: p for p in json.loads((data / "phases" / "registry.json").read_text(encoding="utf-8"))["phases"]}
            c = reg.get(it["source_id"])
            if c is None:
                out[it["item"]] = "missing"
                continue
            det = c["detail"]
            if det not in sysfiles:
                sysfiles[det] = {p["id"]: p for p in json.loads((data / "phases" / det).read_text(encoding="utf-8"))["phases"]}
            v = reg_view(sysfiles[det][it["source_id"]], key["design"]["seed"])
            out[it["item"]] = "unchanged" if view_sha(v) == it["reviewed_sha256"] else "changed"
        elif it["kind"] == "S-ESPEI":
            path, i = it["source_id"].rsplit("#", 1)
            p = data / "espei" / path
            if not p.exists():
                out[it["item"]] = "missing"
                continue
            d = json.loads(p.read_text(encoding="utf-8"))
            ids = espei_value_ids(d)
            if int(i) >= len(ids) or ids[int(i)] != it["orchid_id"]:
                out[it["item"]] = "missing"
                continue
            v = espei_view(d, path, int(i), None)
            out[it["item"]] = "unchanged" if view_sha(espei_reviewed(v)) == it["reviewed_sha256"] else "changed"
        else:
            if asm is None:
                asm = {}
                for l in gzip.decompress((data / "assessed-phase-assemblages" / "records.jsonl.gz").read_bytes()).decode("utf-8").splitlines():
                    if l.strip():
                        r = json.loads(l)
                        asm[r["id"]] = r
            r = asm.get(it["source_id"])
            if r is None:
                out[it["item"]] = "missing"
                continue
            v = {"machine_name": r["assemblage_label"], "T_K": r["T_K"], "x": r["x"]}
            out[it["item"]] = "unchanged" if view_sha(v) == it["reviewed_sha256"] else "changed"
    return out


def verify(out):
    """-> list of problems (files changed or missing, files not in the freeze record)."""
    mp = out / "manifest.json"
    if not mp.exists():
        return ["%s missing" % mp]
    man = json.loads(mp.read_text(encoding="utf-8"))
    bad = []
    for rel, sha in man["files_sha256"].items():
        p = out / rel
        if not p.exists():
            bad.append("frozen file missing: " + rel)
        elif sha256_file(p) != sha:
            bad.append("frozen file changed: " + rel)
    for d in ("packets", "img"):
        if (out / d).is_dir():
            for p in sorted((out / d).iterdir()):
                if p.is_file() and "%s/%s" % (d, p.name) not in man["files_sha256"]:
                    bad.append("file not in the freeze record: %s/%s" % (d, p.name))
    return bad


def cmd_check(a):
    bad = verify(a.out)
    if bad:
        print("supplement: " + "; ".join(bad[:10]) + (" ... (%d)" % len(bad) if len(bad) > 10 else ""), file=sys.stderr)
        return 1
    man = json.loads((a.out / "manifest.json").read_text(encoding="utf-8"))
    key = json.loads((a.out / "key.json").read_text(encoding="utf-8"))
    st = current_views(key, a.data)
    if a.drift_json:
        a.drift_json.write_text(canonical(st), encoding="utf-8")
    c = Counter(st.values())
    print("supplementary sample frozen and intact (sample %s, %d items); release check: %s" % (
        man["sample_sha256"][:16], len(st), ", ".join("%s %d" % kv for kv in sorted(c.items()))))
    ch = [i for i, s in st.items() if s != "unchanged"]
    if ch:
        print("  not unchanged: " + ", ".join("%s %s" % (i, st[i]) for i in ch[:20]) + (" ..." if len(ch) > 20 else ""))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("command", nargs="?", default="check", choices=["check", "draw"])
    ap.add_argument("--out", type=Path, default=KIT, help="kit directory (default: %(default)s)")
    ap.add_argument("--data", type=Path, default=DATA, help="orchid-data release checkout (default: %(default)s)")
    ap.add_argument("--annot", type=Path, default=ANNOT, help="open-tdb okamoto-annotation build dir (S-ASM geometry)")
    ap.add_argument("--blind-key", type=Path, default=BLIND_KEY, help="key.json of the earlier blind sample (excluded)")
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--k", type=int, default=K_DEFAULT, help="draw: number of reviewer packets")
    ap.add_argument("--allow-dirty", action="store_true", help="draw: allow uncommitted orchid-data changes (dry run)")
    ap.add_argument("--no-images", action="store_true", help="draw: skip page crops (tests)")
    ap.add_argument("--drift-json", type=Path, help="check: write the per-item release status here")
    a = ap.parse_args(argv)
    return {"draw": cmd_draw, "check": cmd_check}[a.command](a)


if __name__ == "__main__":
    sys.exit(main())
