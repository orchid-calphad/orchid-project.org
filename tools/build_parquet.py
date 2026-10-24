#!/usr/bin/env python3
"""Build analysis-ready Parquet tables of the ORCHID data repository (<orchid-data>/parquet/).

Reads (never writes) the published JSON of an orchid-data checkout:

  literature/orchid-literature.jsonl.gz   one record per literature value
  literature/papers.json                  bibliography
  literature/contexts.json                experimental-context items per paper
  literature/dft_parity.json              literature x DFT formation-enthalpy pairs
  literature/index.json                   release string
  phase-diagrams/index.json, sgte/*.json  calculated SGTE binary phase diagrams
  hulls/dft-mp/, hulls/dft-oqmd/          Materials Project and OQMD binary convex hulls

Writes <orchid-data>/parquet/:

  literature_values.parquet          one row per value, flat typed columns, split column (by paper)
  papers.parquet                     one row per paper / cited work
  experimental_context.parquet       one row per context item
  dft_parity.parquet                 one row per literature record x DFT source
  phase_diagram_index.parquet        one row per calculated diagram (incl. failed calculations)
  phase_diagram_invariants.parquet   one row per calculated invariant
  hull_vertices.parquet              one row per MP / OQMD hull vertex
  MANIFEST.json                      rows, schema, content SHA-256 and input hashes per table

Deterministic: rows are sorted by stable keys, zstd level and row-group size are fixed, no timestamps. The
content hash is the SHA-256 of the canonical JSON of every row (independent of Parquet encoder version), so
`--check` compares row counts, schema and content of the files on disk with a fresh in-memory build.

Train/validation/test split (literature_values.split): sha256("orchid-split-v1:" + paper_id) -> first 8 bytes
as an integer mod 100: < 80 train, < 90 validation, else test. No paper spans two splits.
Machine-learning folds (literature_values, dft_parity; tools/ml_folds.py): chemical_system, fold_lso =
int(sha256(chemical_system), 16) mod 5 (leave-chemical-system-out), fold_lpo = int(sha256(paper_id), 16) mod 5
(leave-paper-out). Subsystems are systems of their own.

Requires pyarrow. No network.

  python3 tools/build_parquet.py [--data ~/projects/orchid-data]           # write parquet/
  python3 tools/build_parquet.py [--data ...] --check                      # exit 1 if parquet/ is stale
"""
from __future__ import annotations

import argparse
import glob
import gzip
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orchid_cite  # noqa: E402
import ml_folds  # noqa: E402

# table -> citable Zenodo dataset it belongs to (tools/orchid_cite.py)
TABLE_DATASET = {"literature_values": "literature", "papers": "literature", "experimental_context": "literature",
                 "dft_parity": "dft-benchmark", "phase_diagram_index": "phase-diagrams",
                 "phase_diagram_invariants": "phase-diagrams", "hull_vertices": "hulls"}
KJ_PER_MOL_PER_EV = 96.485  # same factor as literature/dft_parity.json
SPLIT_SALT = "orchid-split-v1:"
ZSTD_LEVEL = 9
ROW_GROUP = {"literature_values": 8192}
DEFAULT_ROW_GROUP = 65536


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def split_of(paper_id):
    n = int.from_bytes(hashlib.sha256((SPLIT_SALT + paper_id).encode()).digest()[:8], "big") % 100
    return "train" if n < 80 else ("validation" if n < 90 else "test")


def fnum(v):
    if v is None or isinstance(v, bool):
        return None
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def page_int(v):
    """Page numbers outside 1..100000 (two corrupted quote_page values in the 2026-10-22 release) -> null."""
    return v if isinstance(v, int) and not isinstance(v, bool) and 0 < v <= 100000 else None


def L(t):
    """List type with the child name Parquet uses ("element"), so built and read-back schemas compare equal."""
    import pyarrow as pa
    return pa.list_(pa.field("element", t))


def jstr(v):
    return None if v is None else json.dumps(v, sort_keys=True, ensure_ascii=False)


# ---------------------------------------------------------------- literature values

BOUND = {"<": "upper", ">": "lower", "upper": "upper", "lower": "lower"}
BOOL_FLAGS = ["espei_ready", "judge_adjudicated", "magnitude_suspect", "source_disputed", "system_mismatch",
              "unsigned_in_print", "sign_pdf_page_verified", "has_dropped_fields",
              "composition_conflict_between_readers", "unit_conflict_between_readers", "page_unreadable",
              "basis_inferred_from_magnitude", "basis_corrected_from_formula", "sign_from_heading",
              "dft_outlier_unexplained", "reference_state_dropped_not_in_quote", "espei_certified",
              "espei_ready_upstream"]


def lit_schema(pa):
    s, f, b, i = pa.string(), pa.float64(), pa.bool_(), pa.int32()
    ls = L(s)
    comp = L(pa.struct([("phase", s), ("element", s), ("value", f), ("basis", s), ("source", s)]))
    phd = L(pa.struct([("label", s), ("model_phase", s), ("mapping_confidence", s),
                              ("pearson_symbol", s), ("prototype", s)]))
    cols = [
        ("id", s), ("pid", s), ("system", s), ("elements", ls), ("n_elements", i), ("record_kind", s),
        ("observation_type", s), ("property", s), ("property_group", s), ("property_printed", s),
        ("value", f), ("value_min", f), ("value_max", f), ("value_expression", s), ("unit", s),
        ("value_stored", f), ("unit_stored", s), ("value_printed", f), ("unit_printed", s),
        ("uncertainty_value", f), ("uncertainty_basis", s), ("uncertainty_note", s),
        ("sigma_printed", f), ("sigma_printed_kind", s), ("sigma_reading", f), ("sigma_method", f),
        ("sigma_method_class", s), ("sigma_method_basis", s), ("sigma_method_n_pairs", i), ("sigma_method_floor", f),
        ("sigma_method_floor_unit", s), ("sigma_composition", f), ("sigma_scatter", f), ("sigma_conversion", f),
        ("sigma_effective", f), ("sigma_unit", s), ("sigma_rule", s), ("uncertainty_calibration", s),
        ("sigma_x_effective", f), ("sigma_x_method", f), ("sigma_x_method_class", s),
        ("sigma_T_condition_effective", f), ("uncertainty_budget", s),
        ("bound", s), ("temperature_bound_K", f),
        ("temperature_K", f), ("temperature_min_K", f), ("temperature_max_K", f), ("temperature_printed", s),
        ("temperature_basis", s), ("temperature_approximate", b), ("temperature_uncertainty_value", f),
        ("temperature_uncertainty_basis", s), ("pressure_Pa", f), ("pressure_basis", s),
        ("phases", ls), ("model_phases", ls), ("phase_details", phd),
        ("reaction", s), ("reaction_type", s), ("transition_kind", s), ("transition_kind_printed", s),
        ("property_code_raw", s), ("unit_basis", s), ("per_formula_unit_divisor", f), ("value_before_quality_fix", f),
        ("merged_ids", ls),
        ("compositions", comp),
        ("overall_composition", L(pa.struct([("element", s), ("value", f)]))),
        ("overall_composition_basis", s), ("reference_state_json", s), ("component", s),
        ("structure_pearson_symbol", s), ("structure_prototype", s), ("structure_strukturbericht", s),
        ("structure_space_group", s), ("structure_space_group_number", i),
        ("evidence_class", s), ("origin", s), ("origin_class", s), ("origin_measured_in", s),
        ("origin_own_work_of_paper", b), ("accepted", b),
    ] + [("flag_" + k, b) for k in BOOL_FLAGS] + [
        ("flag_state_fully_specified", b), ("flag_value_in_quote", b), ("dropped_field_count", i),
        ("review_state", s), ("row_state", s), ("review_flags", ls), ("specification_gaps", ls),
        ("plausibility_result", s), ("plausibility_dT_K", f), ("plausibility_model", s),
        ("calphad_check", s), ("espei_certified_basis", s), ("espei_phase_rule_result", s),
        ("espei_phase_rule_reason", s), ("espei_blocked_by", ls),
        ("plausibility_dx_json", s), ("gate_failures", ls), ("corrections_json", s), ("temperature_K_stored", f),
        ("system_stored", s), ("evidence_class_before_gates", s), ("property_code_before_gates", s),
        ("reader_agreement_basis", s),
        ("method", s), ("method_basis", s), ("conditions_kind", s), ("conditions_paper", s),
        ("paper_id", s), ("paper_doi", s), ("paper_citation", s), ("paper_year", i),
        ("cited_dois", ls), ("cited_labels", ls),
        ("quote", s), ("quote_page", i), ("link_rule", s),
        ("artifact", s), ("artifact_label", s), ("artifact_caption", s), ("page", i), ("table_row", i),
        ("model_tdb", s), ("model_tdb_path", s), ("model_tdb_sha256", s), ("model_doi", s),
        ("model_reproduced_by_calculation", b), ("model_unlinked_reason", s),
        ("readers", ls), ("reader_tier", s), ("reader_agreement", s), ("reader_agreement_label", s),
        ("tool_release", s), ("release", s), ("split", s),
        ("chemical_system", s), ("fold_lso", i), ("fold_lpo", i),
    ]
    return pa.schema(cols)


def lit_row(r, papers, release):
    fl = r.get("flags") or {}
    unc = r.get("uncertainty") or {}
    vr = r.get("value_range") or {}
    tr = r.get("temperature_range_K") or {}
    tu = r.get("temperature_uncertainty") or {}
    tb = r.get("temperature_bound")
    if isinstance(tb, dict):
        bound, tb_k = BOUND.get(tb.get("op")), fnum(tb.get("value_K"))
    else:
        bound, tb_k = BOUND.get(tb), None  # "none" -> None
    ph = r.get("phases") or []
    rx = r.get("reaction") or {}
    rtype = rx.get("type")
    comps = []
    for pc in r.get("phase_compositions") or []:
        for el, x in sorted((pc.get("composition") or {}).items()):
            comps.append({"phase": pc.get("phase"), "element": el, "value": fnum(x), "basis": pc.get("basis"),
                          "source": pc.get("composition_source")})
    st = r.get("structure") or {}
    org = r.get("origin") or {}
    pl = fl.get("plausibility") or {}
    cond = r.get("conditions") or {}
    loc = r.get("locator") or {}
    mdl = r.get("model") or {}
    pid = r["paper"]
    pap = papers.get(pid) or {}
    cited = r.get("cited_sources") or []
    out = {
        "id": r["id"], "pid": r.get("pid"), "system": r["system"], "elements": list(r.get("elements") or []),
        "n_elements": len(r.get("elements") or []), "record_kind": r.get("record_kind"),
        "observation_type": r.get("observation_type"), "property": r.get("property"),
        "property_group": r.get("property_group"), "property_printed": r.get("property_printed"),
        "value": fnum(r.get("value")), "value_min": fnum(vr.get("min")), "value_max": fnum(vr.get("max")),
        "value_expression": r.get("value_expression"), "unit": r.get("unit"),
        "value_stored": fnum(r.get("value_stored")), "unit_stored": r.get("unit_stored"),
        "value_printed": fnum(r.get("value_printed")), "unit_printed": r.get("unit_printed"),
        "uncertainty_value": fnum(unc.get("value")), "uncertainty_basis": unc.get("basis"),
        "uncertainty_note": unc.get("note"), "bound": bound, "temperature_bound_K": tb_k,
        "temperature_K": fnum(r.get("temperature_K")), "temperature_min_K": fnum(tr.get("min")),
        "temperature_max_K": fnum(r.get("temperature_max_K", tr.get("max"))),
        "temperature_printed": r.get("temperature_printed"), "temperature_basis": r.get("temperature_basis"),
        "temperature_approximate": r.get("temperature_approximate"),
        "temperature_uncertainty_value": fnum(tu.get("value")), "temperature_uncertainty_basis": tu.get("basis"),
        "pressure_Pa": fnum(r.get("pressure_Pa")), "pressure_basis": r.get("pressure_basis"),
        "phases": [p.get("label") for p in ph], "model_phases": [p.get("model_phase") for p in ph],
        "phase_details": [{k: p.get(k) for k in ("label", "model_phase", "mapping_confidence", "pearson_symbol",
                                                 "prototype")} for p in ph],
        "reaction": rx.get("printed"), "reaction_type": rtype,
        "transition_kind": r.get("transition_kind"), "transition_kind_printed": r.get("transition_kind_printed"),
        "property_code_raw": r.get("property_code_raw"), "unit_basis": r.get("unit_basis"),
        "per_formula_unit_divisor": fnum(r.get("per_formula_unit_divisor")),
        "value_before_quality_fix": fnum(r.get("value_before_quality_fix")),
        "merged_ids": list(r.get("merged_ids") or []),
        "compositions": comps,
        "overall_composition": [{"element": k, "value": fnum(v)}
                                for k, v in sorted((r.get("overall_composition") or {}).items())],
        "overall_composition_basis": r.get("overall_composition_basis"),
        "reference_state_json": jstr(r.get("reference_state")), "component": r.get("component"),
        "structure_pearson_symbol": st.get("pearson_symbol"), "structure_prototype": st.get("prototype"),
        "structure_strukturbericht": st.get("strukturbericht"), "structure_space_group": st.get("space_group"),
        "structure_space_group_number": st.get("space_group_number"),
        "evidence_class": r.get("evidence_class"), "origin": org.get("basis"), "origin_class": org.get("class"),
        "origin_measured_in": org.get("measured_in"), "origin_own_work_of_paper": org.get("own_work_of_paper"),
        "accepted": bool(fl.get("extraction_accepted")),
        "flag_state_fully_specified": fl.get("state_fully_specified"),
        "flag_value_in_quote": fl.get("value_in_quote"), "dropped_field_count": fl.get("dropped_field_count"),
        "review_state": fl.get("review_state"), "row_state": fl.get("row_state"),
        "review_flags": [str(x) for x in fl.get("review_flags") or []],
        "specification_gaps": [str(x) for x in fl.get("specification_gaps") or []],
        "plausibility_result": pl.get("result"), "plausibility_dT_K": fnum(pl.get("dT_K")),
        "plausibility_model": pl.get("model"),
        "calphad_check": fl.get("calphad_check"), "espei_certified_basis": fl.get("espei_certified_basis"),
        "espei_phase_rule_result": (fl.get("espei_phase_rule") or {}).get("result"),
        "espei_phase_rule_reason": (fl.get("espei_phase_rule") or {}).get("reason"),
        "espei_blocked_by": [str(x) for x in fl.get("espei_blocked_by") or []],
        "plausibility_dx_json": jstr(pl.get("dx")),
        "gate_failures": [str(x.get("code")) for x in r.get("gate_failures") or []],
        "corrections_json": jstr(r.get("corrections")), "temperature_K_stored": fnum(r.get("temperature_K_stored")),
        "system_stored": r.get("system_stored"), "evidence_class_before_gates": r.get("evidence_class_before_gates"),
        "property_code_before_gates": r.get("property_code_before_gates"),
        "reader_agreement_basis": r.get("reader_agreement_basis"),
        "method": r.get("method"), "method_basis": r.get("method_basis"),
        "conditions_kind": cond.get("kind"), "conditions_paper": cond.get("paper"),
        "paper_id": pid, "paper_doi": pap.get("doi") or (pid[4:] if pid.startswith("doi:") else None),
        "paper_citation": pap.get("citation"), "paper_year": pap.get("year"),
        "cited_dois": sorted({c["doi"] for c in cited if c.get("doi")}),
        "cited_labels": [c.get("label") for c in cited if c.get("label") is not None],
        "quote": r.get("quote"), "quote_page": page_int(r.get("quote_page")), "link_rule": r.get("link_rule"),
        "artifact": loc.get("artifact"), "artifact_label": loc.get("label"), "artifact_caption": loc.get("caption"),
        "page": page_int(r.get("quote_page", loc.get("page"))), "table_row": loc.get("row"),
        "model_tdb": mdl.get("tdb"), "model_tdb_path": mdl.get("tdb_path"), "model_tdb_sha256": mdl.get("tdb_sha256"),
        "model_doi": mdl.get("model_doi"), "model_reproduced_by_calculation": mdl.get("reproduced_by_calculation"),
        "model_unlinked_reason": r.get("model_unlinked_reason"),
        "readers": list(r.get("readers") or []), "reader_tier": r.get("reader_tier"),
        "reader_agreement": r.get("reader_agreement"), "reader_agreement_label": r.get("reader_agreement_label"),
        "tool_release": r.get("tool_release"), "release": release, "split": split_of(pid),
        # machine-learning folds (tools/ml_folds.py), as in the JSON release
        "chemical_system": r["system"], "fold_lso": ml_folds.fold(r["system"]), "fold_lpo": ml_folds.fold(pid),
    }
    for k in BOOL_FLAGS:
        out["flag_" + k] = bool(fl.get(k))
    out.update(budget_columns(r.get("uncertainty_budget")))
    return out


def budget_columns(b):
    """Flat columns of the record's uncertainty budget (literature SCHEMA.md, `uncertainty_budget`): the value
    quantity's components, the phase-composition and equilibration-temperature effective sigmas, the full JSON."""
    b = b or {}
    qs = b.get("quantities") or {}
    q = qs.get("value") or {}
    qx = qs.get("composition") or {}
    qt = qs.get("temperature_K") or {}

    def v(c):
        return fnum((c or {}).get("value")) if isinstance(c, dict) else None
    m = q.get("sigma_method") or {}
    mx = qx.get("sigma_method") or {}
    n = m.get("n_pairs")
    return {"sigma_printed": v(q.get("sigma_printed")), "sigma_printed_kind": (q.get("sigma_printed") or {}).get("kind"),
            "sigma_reading": v(q.get("sigma_reading")), "sigma_method": v(m), "sigma_method_class": m.get("class"),
            "sigma_method_basis": m.get("basis"), "sigma_method_n_pairs": int(n) if isinstance(n, int) else None,
            "sigma_method_floor": fnum(m.get("floor")), "sigma_method_floor_unit": m.get("floor_unit"),
            "sigma_composition": v(q.get("sigma_composition")), "sigma_scatter": v(q.get("sigma_scatter")),
            "sigma_conversion": v(q.get("sigma_conversion")), "sigma_effective": fnum(q.get("sigma_effective")),
            "sigma_unit": q.get("unit") if q else None, "sigma_rule": b.get("sigma_rule"),
            "uncertainty_calibration": b.get("calibration"),
            "sigma_x_effective": fnum(qx.get("sigma_effective")), "sigma_x_method": v(mx),
            "sigma_x_method_class": mx.get("class"), "sigma_T_condition_effective": fnum(qt.get("sigma_effective")),
            "uncertainty_budget": jstr(b) if b else None}


# ---------------------------------------------------------------- other tables

def papers_schema(pa):
    s, ls = pa.string(), L(pa.string())
    return pa.schema([("id", s), ("doi", s), ("doi_url", s), ("citation", s), ("title", s), ("authors", ls),
                      ("year", pa.int32()), ("container", s), ("kind", s), ("roles", ls), ("source_kinds", ls),
                      ("systems", ls), ("metadata_source", s), ("pdf_sha256_prefix", s),
                      ("n_values", pa.int32())])


CTX_FIELDS = ["label", "purity", "starting_materials", "synthesis", "nominal_composition", "composition_basis",
              "measured_composition", "temperature_printed", "temperature_status", "duration",
              "atmosphere_or_container", "cooling", "purpose", "method", "instrument", "heating_or_cooling_rate",
              "atmosphere", "crucible", "calibration", "stated_uncertainty", "applies_to"]


def ctx_schema(pa):
    s, b = pa.string(), pa.bool_()
    return pa.schema([("paper_id", s), ("context_id", s), ("section_basis", s), ("context_reader_tier", s),
                      ("context_readers", L(s)), ("item_id", s), ("type", s)]
                     + [(k, s) for k in CTX_FIELDS]
                     + [("temperature_K", pa.float64()), ("duration_hours", pa.float64()), ("quote", s),
                        ("page", pa.int32()), ("extraction_accepted", b), ("second_reader_agrees", b),
                        ("has_dropped_fields", b), ("review_flags", L(s))])


def parity_schema(pa):
    s, f, b = pa.string(), pa.float64(), pa.bool_()
    return pa.schema([
        ("record_id", s), ("dft_source", s), ("system", s), ("elements", L(s)), ("formula", s),
        ("phase_label", s), ("composition", L(pa.struct([("element", s), ("x", f)]))),
        ("formula_match", s), ("match_note", s), ("lit_kJ_per_mol_atom", f), ("dft_kJ_per_mol_atom", f),
        ("dft_eV_per_atom", f), ("residual_kJ_per_mol_atom", f), ("ehull_kJ_per_mol_atom", f),
        ("dft_material_id", s), ("dft_spacegroup", s), ("dft_theoretical", b), ("dft_n_polymorphs", pa.int32()),
        ("evidence_class", s), ("accepted", b), ("own_work_of_paper", b), ("unsigned", b), ("sign_lost", b),
        ("sign_verified", b), ("sign_check", s), ("suspect", L(s)), ("chem_class", s),
        ("method_class", s), ("method", s), ("temperature_K", f), ("record_kind", s), ("paper_id", s),
        ("page", pa.int32()), ("table", s), ("quote", s), ("cited_dois", L(s)), ("default_subset", b),
        ("dft_release", s), ("quality_fix", s), ("unit_basis", s), ("basis_inferred_from_magnitude", b),
        ("chemical_system", s), ("fold_lso", pa.int32()), ("fold_lpo", pa.int32()),
        ("ehull_upstream_kJ_per_mol_atom", f), ("dft_structure_match", s), ("structure_designators", L(s)),
        ("dft_ground_state_id", s), ("x_offset_measured", f), ("sigma_composition_kJ_per_mol_atom", f),
        ("excluded_from_default_by", L(s))])


def pd_index_schema(pa):
    s = pa.string()
    return pa.schema([("system", s), ("elements", L(s)), ("file", s), ("tdb", s), ("tdb_sha256", s),
                      ("status", s), ("listed_in_index", pa.bool_()), ("has_grid", pa.bool_()),
                      ("x_component", s), ("T_min_K", pa.float64()), ("T_max_K", pa.float64()),
                      ("phases", L(s)), ("n_invariants", pa.int32()), ("engine", s),
                      ("citation", s), ("doi", s)])


def pd_inv_schema(pa):
    s, f = pa.string(), pa.float64()
    return pa.schema([("system", s), ("file", s), ("tdb", s), ("tdb_sha256", s), ("x_component", s),
                      ("invariant_index", pa.int32()), ("T_K", f), ("n_phases", pa.int32()),
                      ("phases", L(s)), ("compositions", L(pa.struct([("phase", s), ("x", f)]))),
                      ("x_min", f), ("x_max", f), ("engine", s), ("citation", s), ("doi", s)])


def hull_schema(pa):
    s, f, b = pa.string(), pa.float64(), pa.bool_()
    return pa.schema([("system", s), ("source", s), ("a_element", s), ("x_element", s), ("formula", s),
                      ("x", f), ("formation_energy_eV_per_atom", f), ("formation_energy_kJ_per_mol_atom", f),
                      ("energy_per_atom_eV", f), ("material_id", s), ("spacegroup_symbol", s),
                      ("spacegroup_number", pa.int32()), ("theoretical", b), ("collinear", b),
                      ("degenerate_ids", L(s)), ("is_endpoint", b), ("reference_kind", s),
                      ("release", s), ("thermo_type", s), ("entry_pool", s), ("snapshot_sha256", s)])


def build_tables(data, pa):
    lit_dir = os.path.join(data, "literature")
    inputs = {}

    def src(rel):
        p = os.path.join(data, rel)
        inputs[rel] = sha256_file(p)
        return p

    release = json.load(open(src("literature/index.json")))["release"]
    papers = json.load(open(src("literature/papers.json")))
    rows, n_by_paper = [], {}
    with gzip.open(src("literature/orchid-literature.jsonl.gz"), "rt", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                r = json.loads(line)
                rows.append(lit_row(r, papers, release))
                n_by_paper[r["paper"]] = n_by_paper.get(r["paper"], 0) + 1
    rows.sort(key=lambda x: (x["system"], x["id"]))
    tables = {"literature_values": (rows, lit_schema(pa))}

    prow = []
    for pid in sorted(papers):
        p = papers[pid]
        prow.append({"id": pid, "doi": p.get("doi"), "doi_url": p.get("doi_url"), "citation": p.get("citation"),
                     "title": p.get("title"), "authors": p.get("authors"), "year": p.get("year"),
                     "container": p.get("container"), "kind": p.get("kind"), "roles": p.get("roles"),
                     "source_kinds": p.get("source_kinds"), "systems": p.get("systems"),
                     "metadata_source": p.get("metadata_source"), "pdf_sha256_prefix": p.get("pdf_sha256_prefix"),
                     "n_values": n_by_paper.get(pid, 0)})
    tables["papers"] = (prow, papers_schema(pa))

    ctx = json.load(open(src("literature/contexts.json")))
    crow = []
    for pid in sorted(ctx):
        c = ctx[pid]
        for it in c.get("items") or []:
            t = it.get("temperature") or {}
            d = it.get("duration_parsed") or {}
            row = {"paper_id": c.get("paper", pid), "context_id": c.get("context_id"),
                   "section_basis": c.get("section_basis"), "context_reader_tier": c.get("reader_tier"),
                   "context_readers": c.get("readers"), "item_id": it.get("id"), "type": it.get("type"),
                   "temperature_K": fnum(t.get("T_K")), "duration_hours": fnum(d.get("hours")),
                   "quote": it.get("quote"), "page": it.get("page"),
                   "extraction_accepted": it.get("extraction_accepted"),
                   "second_reader_agrees": it.get("second_reader_agrees"),
                   "has_dropped_fields": bool(it.get("has_dropped_fields")),
                   "review_flags": [str(x) for x in it.get("review_flags") or []]}
            for k in CTX_FIELDS:
                row[k] = it.get(k)
            crow.append(row)
    tables["experimental_context"] = (crow, ctx_schema(pa))

    par = json.load(open(src("literature/dft_parity.json")))
    rel = {k: (v.get("release") or "") for k, v in par["dft_sources"].items()}
    drow = []
    for pt in par["points"]:
        r = par["records"][pt["rid"]]
        default = (pt["default"] if "default" in pt else
                   (r["ev"] == "experimental" and r["acc"] and not r["unsigned"] and not r["sign_lost"]
                    and "magnitude_suspect" not in (r.get("suspect") or [])))
        drow.append({
            "record_id": pt["rid"], "dft_source": pt["src"], "system": r["sys"], "elements": r["elements"],
            "formula": r["formula"], "phase_label": r["phase"],
            "composition": [{"element": k, "x": fnum(v)} for k, v in sorted(r["comp"].items())],
            "formula_match": r["match"], "match_note": r["match_note"], "lit_kJ_per_mol_atom": fnum(r["lit"]),
            "dft_kJ_per_mol_atom": fnum(pt["dft"]), "dft_eV_per_atom": fnum(pt["dft_eV"]),
            "residual_kJ_per_mol_atom": fnum(pt["res"]), "ehull_kJ_per_mol_atom": fnum(pt["ehull"]),
            "dft_material_id": pt["mid"], "dft_spacegroup": pt["sg"], "dft_theoretical": pt["theoretical"],
            "dft_n_polymorphs": pt["npoly"], "evidence_class": r["ev"], "accepted": r["acc"],
            "own_work_of_paper": r["own"], "unsigned": r["unsigned"], "sign_lost": r["sign_lost"],
            "sign_verified": r["sign_verified"], "sign_check": r["sign_check"],
            "suspect": [str(x) for x in r["suspect"]], "chem_class": r["ccls"], "method_class": r["mcls"],
            "method": r["method"], "temperature_K": fnum(r["T"]), "record_kind": r["kind"], "paper_id": r["paper"],
            "page": r["page"], "table": r["table"], "quote": r["quote"],
            "cited_dois": sorted({c["doi"] for c in r["cited"] if c.get("doi")}), "default_subset": default,
            "dft_release": rel.get(pt["src"]), "quality_fix": r.get("qfix"), "unit_basis": r.get("basis"),
            "basis_inferred_from_magnitude": bool(r.get("basis_inferred")),
            # the literature record's own ML keys (audit A1 C5: one fold per record_id in both tables)
            "chemical_system": r.get("chemical_system") or r["sys"],
            "fold_lso": r["fold_lso"] if r.get("fold_lso") is not None else ml_folds.fold(r.get("chemical_system") or r["sys"]),
            "fold_lpo": r["fold_lpo"] if r.get("fold_lpo") is not None else ml_folds.fold(r["paper"]),
            "ehull_upstream_kJ_per_mol_atom": fnum(pt.get("ehull_upstream")),
            "dft_structure_match": pt.get("smatch"), "structure_designators": sorted(r.get("designators") or []),
            "dft_ground_state_id": pt.get("ground_id"), "x_offset_measured": fnum(r.get("dx_measured")),
            "sigma_composition_kJ_per_mol_atom": fnum(pt.get("sigma_comp")),
            "excluded_from_default_by": [str(x) for x in pt.get("excluded_by") or []]})
    drow.sort(key=lambda x: (x["system"], x["record_id"], x["dft_source"]))
    tables["dft_parity"] = (drow, parity_schema(pa))

    pdi = json.load(open(src("phase-diagrams/index.json")))
    listed = {e["file"]: (sysn, e) for sysn, es in pdi["systems"].items() for e in es}
    irow, invrow = [], []
    sgte_digest = hashlib.sha256()
    for fp in sorted(glob.glob(os.path.join(data, "phase-diagrams", "sgte", "*.json"))):
        fn = os.path.basename(fp)
        sgte_digest.update(("%s  %s\n" % (sha256_file(fp), fn)).encode())
        j = json.load(open(fp))
        els = sorted(e.upper() for e in j.get("elements") or [])
        sysn, ent = listed.get(fn, ("-".join(els), {}))
        osrc = j.get("original_source") or {}
        cit, doi = ent.get("citation", osrc.get("citation")), ent.get("doi", osrc.get("doi"))
        tr = j.get("T_range") or [None, None]
        invs = j.get("invariants") or []
        irow.append({"system": sysn, "elements": els, "file": fn, "tdb": j.get("tdb"),
                     "tdb_sha256": j.get("sha256"), "status": j.get("status"), "listed_in_index": fn in listed,
                     "has_grid": bool(j.get("grid")), "x_component": j.get("x_component"),
                     "T_min_K": fnum(tr[0]), "T_max_K": fnum(tr[1]), "phases": j.get("phases") or [],
                     "n_invariants": len(invs), "engine": j.get("engine"), "citation": cit, "doi": doi})
        for k, inv in enumerate(invs):
            pts = inv.get("points") or []
            xs = [fnum(p.get("x")) for p in pts if p.get("x") is not None]
            invrow.append({"system": sysn, "file": fn, "tdb": j.get("tdb"), "tdb_sha256": j.get("sha256"),
                           "x_component": j.get("x_component"), "invariant_index": k, "T_K": fnum(inv.get("T")),
                           "n_phases": len(pts), "phases": [p.get("phase") for p in pts],
                           "compositions": [{"phase": p.get("phase"), "x": fnum(p.get("x"))} for p in pts],
                           "x_min": min(xs) if xs else None, "x_max": max(xs) if xs else None,
                           "engine": j.get("engine"), "citation": cit, "doi": doi})
    inputs["phase-diagrams/sgte/*.json (sha256 of 'sha256  name' lines)"] = sgte_digest.hexdigest()
    irow.sort(key=lambda x: (x["system"], x["file"]))
    invrow.sort(key=lambda x: (x["system"], x["file"], x["invariant_index"]))
    tables["phase_diagram_index"] = (irow, pd_index_schema(pa))
    tables["phase_diagram_invariants"] = (invrow, pd_inv_schema(pa))

    hrow = []
    for source, layer in (("mp", "dft-mp"), ("oqmd", "dft-oqmd")):
        cat_rel = "hulls/%s/catalog.json" % layer
        cat = json.load(open(src(cat_rel)))
        for s in cat["systems"]:
            hp = os.path.join(data, "hulls", s["hull_path"])
            h = json.load(open(hp))
            if sha256_file(hp) != s["hull_sha256"]:
                raise SystemExit("hull hash mismatch: " + s["hull_path"])
            pv, ref = h.get("provenance") or {}, h.get("reference") or {}
            els, xel = h["elements"], h["x_element"]
            ael = [e for e in els if e != xel][0]
            for iso in h["isotherms"]:
                for v in iso["vertices"]:
                    e = fnum(v.get("energy"))
                    hrow.append({
                        "system": h["system_id"], "source": source, "a_element": ael, "x_element": xel,
                        "formula": v.get("formula"), "x": fnum(v.get("x")), "formation_energy_eV_per_atom": e,
                        "formation_energy_kJ_per_mol_atom": None if e is None else round(e * KJ_PER_MOL_PER_EV, 6),
                        "energy_per_atom_eV": fnum(v.get("energy_per_atom")), "material_id": v.get("phase"),
                        "spacegroup_symbol": v.get("spacegroup_symbol"),
                        "spacegroup_number": v.get("spacegroup_number"), "theoretical": v.get("theoretical"),
                        "collinear": v.get("collinear"),
                        "degenerate_ids": [d.get("phase") for d in v.get("degenerate") or []],
                        "is_endpoint": fnum(v.get("x")) in (0.0, 1.0), "reference_kind": ref.get("kind"),
                        "release": pv.get("source_release", pv.get("mp_release")), "thermo_type": pv.get("thermo_type"),
                        "entry_pool": pv.get("entry_pool"), "snapshot_sha256": pv.get("snapshot_sha256")})
    hrow.sort(key=lambda x: (x["system"], x["source"], x["x"], x["material_id"] or ""))
    tables["hull_vertices"] = (hrow, hull_schema(pa))

    expect = {"literature_values": len(rows), "papers": len(papers),
              "experimental_context": sum(len(c.get("items") or []) for c in ctx.values()),
              "dft_parity": len(par["points"]),
              "phase_diagram_index": len(glob.glob(os.path.join(data, "phase-diagrams", "sgte", "*.json"))),
              "phase_diagram_invariants": len(invrow),
              "hull_vertices": None}
    return tables, inputs, release, expect


def content_hash(rows):
    h = hashlib.sha256()
    for r in rows:
        h.update(json.dumps(r, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode())
        h.update(b"\n")
    return h.hexdigest()


def table_metadata(name, release, data):
    """Parquet key-value metadata: release, licence and how to cite (orchid.cite_as = JSON of tools/orchid_cite.py)."""
    c = orchid_cite.cite_as(TABLE_DATASET[name], data)
    return {"orchid.release": release, "orchid.table": name, "orchid.license": "CC-BY-4.0",
            "orchid.generator": "orchid-project.org/tools/build_parquet.py",
            "orchid.dataset": TABLE_DATASET[name], "orchid.doi": c["doi"], "orchid.citation": c["text"],
            "orchid.cite_as": json.dumps(c, sort_keys=True, ensure_ascii=False)}


def write_table(pa, pq, name, rows, schema, path, release, data):
    md = table_metadata(name, release, data)
    t = pa.Table.from_pylist(rows, schema=schema.with_metadata(md))
    pq.write_table(t, path, compression="zstd", compression_level=ZSTD_LEVEL,
                   row_group_size=ROW_GROUP.get(name, DEFAULT_ROW_GROUP), write_statistics=True,
                   use_dictionary=True, data_page_version="2.0")
    return t


def manifest_cite(data):
    """MANIFEST.json cite_as: the umbrella record plus each table's dataset (tables name theirs in `dataset`)."""
    return {"umbrella": orchid_cite.cite_as(orchid_cite.UMBRELLA, data),
            "datasets": {n: orchid_cite.cite_as(n, data) for n in sorted(set(TABLE_DATASET.values()))}}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=os.environ.get("ORCHID_DATA", os.path.expanduser("~/projects/orchid-data")),
                    help="orchid-data checkout (reads its JSON, writes parquet/)")
    ap.add_argument("--check", action="store_true", help="exit 1 if parquet/ differs from a fresh build")
    a = ap.parse_args()
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        sys.exit("build_parquet.py needs pyarrow (pip install pyarrow); nothing written")

    tables, inputs, release, expect = build_tables(a.data, pa)
    out = os.path.join(a.data, "parquet")
    man_path = os.path.join(out, "MANIFEST.json")
    built = {}
    for name, (rows, schema) in tables.items():
        if expect.get(name) is not None and expect[name] != len(rows):
            sys.exit("%s: %d rows, source has %d" % (name, len(rows), expect[name]))
        # normalise through Arrow so the hash covers exactly what Parquet stores
        canon = pa.Table.from_pylist(rows, schema=schema).to_pylist()
        built[name] = {"rows": len(rows), "schema": str(schema), "content_sha256": content_hash(canon)}

    if a.check:
        bad = []
        man = json.load(open(man_path)) if os.path.exists(man_path) else {"tables": {}, "inputs": {}}
        if man.get("inputs") != inputs:
            bad.append("inputs changed since the last build")
        if man.get("cite_as") != manifest_cite(a.data):
            bad.append("MANIFEST.json cite_as differs")
        for name, b in built.items():
            p = os.path.join(out, name + ".parquet")
            if not os.path.exists(p):
                bad.append(name + ": missing")
                continue
            t = pq.read_table(p)
            disk = {"rows": t.num_rows, "schema": str(t.schema.remove_metadata()),
                    "content_sha256": content_hash(t.to_pylist())}
            for k in ("rows", "schema", "content_sha256"):
                if disk[k] != b[k] or (man["tables"].get(name) or {}).get(k) != b[k]:
                    bad.append("%s: %s differs" % (name, k))
            kv = {k.decode(): v.decode() for k, v in (t.schema.metadata or {}).items() if k.startswith(b"orchid.")}
            if kv != table_metadata(name, release, a.data):
                bad.append("%s: key-value metadata (release / cite_as) differs" % name)
            if (man["tables"].get(name) or {}).get("dataset") != TABLE_DATASET[name]:
                bad.append("%s: MANIFEST dataset differs" % name)
        print("\n".join(bad) if bad else "parquet/ is up to date (%d tables)" % len(built))
        sys.exit(1 if bad else 0)

    os.makedirs(out, exist_ok=True)
    man = {"description": "ORCHID analysis-ready Parquet tables; content_sha256 = SHA-256 of the canonical JSON "
                          "rows (sorted keys), independent of the Parquet encoder.",
           "release": release, "split_rule": "sha256('%s' + paper_id)[:8] mod 100: <80 train, <90 validation, "
                                              "else test" % SPLIT_SALT,
           "fold_rule": ml_folds.DOC,
           "kJ_per_mol_per_eV": KJ_PER_MOL_PER_EV, "compression": "zstd level %d" % ZSTD_LEVEL,
           "inputs": inputs, "cite_as": manifest_cite(a.data), "tables": {}}
    for name, (rows, schema) in tables.items():
        p = os.path.join(out, name + ".parquet")
        tmp = p + ".tmp"
        write_table(pa, pq, name, rows, schema, tmp, release, a.data)
        os.replace(tmp, p)
        man["tables"][name] = dict(built[name], file=name + ".parquet", bytes=os.path.getsize(p),
                                   dataset=TABLE_DATASET[name],
                                   file_sha256=sha256_file(p))
        print("%-28s %8d rows %10d bytes" % (name, built[name]["rows"], os.path.getsize(p)))
    with open(man_path + ".tmp", "w") as f:
        json.dump(man, f, indent=1, sort_keys=True, ensure_ascii=False)
        f.write("\n")
    os.replace(man_path + ".tmp", man_path)


if __name__ == "__main__":
    main()
