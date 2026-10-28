#!/usr/bin/env python3
"""Write the Zenodo deposition metadata of the ORCHID release: <orchid-data>/zenodo/<name>.metadata.json.

Ten records: the umbrella, the eight datasets (titles with the formal dataset names ORCHID-Values, ORCHID-Diagrams,
ORCHID-Hulls, ORCHID-Parity, ORCHID-Registry, ORCHID-Assemblages, ORCHID-Bibliography and the ESPEI datasets of
ORCHID-Values) and the software record orchid-examples. Every number in a description is READ from the release files
(JSON summaries and Parquet tables) at build time; zenodo/metadata_numbers.json lists each number with its source file,
so nothing is typed by hand and a rebuilt release regenerates the text.

  python3 tools/zenodo_release.py metadata [--check]      (or: python3 tools/build_zenodo_metadata.py [--check])

Single sources of truth (tools/orchid_cite.py): titles (DATASETS, SOFTWARE), the Zenodo community (ZENODO_COMMUNITY),
the funding sentence (FUNDING_NOTE), the licence-scope sentence (LICENCE_SCOPE), the authors. The counts of the
preliminary human review (not in the release files) are read from tools/human_review_summary_20261026.json.
Descriptions follow the X4 Zenodo-record review of 2026-10-27 (flags and checks are described as screening results,
not certification). The record of ORCHID-Assemblages is written by
tools/build_assessed_assemblages.py (its numbers come from that builder's run); this tool leaves it alone.

No `grants` entry: NSF award 2608825 is not in Zenodo's awards database and would be rejected; the funding is entered by
hand as a custom award (zenodo/ZENODO_FUNDING_MANUAL.md) and stated in the notes of every record.
Related identifier to the paper: set orchid_cite.PAPER when the paper has a DOI or arXiv id (see the runbook).

Needs pandas and pyarrow (--parquet-python of zenodo_release.py; run with that interpreter).
"""
from __future__ import annotations

import argparse
import html
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import orchid_cite as oc  # noqa: E402

WEBSITE, DATA_URL = oc.WEBSITE, oc.DATA_URL
TOOLS = os.path.dirname(os.path.abspath(__file__))
HUMAN_REVIEW = "human_review_summary_20261026.json"
DOI_SGTE, DOI_MP, DOI_OQMD1, DOI_OQMD2 = ("10.1016/j.calphad.2025.102833", "10.1063/1.4812323",
                                           "10.1007/s11837-013-0755-4", "10.1038/npjcompumats.2015.10")
DOI_KIM, DOI_KIM_DATA = "10.1038/sdata.2017.162", "10.6084/m9.figshare.7289723"
NAMES_ALL = list(oc.NAMES) + list(oc.SOFTWARE_NAMES)
BUILT_ELSEWHERE = "assessed-phase-assemblages"   # tools/build_assessed_assemblages.py

# ---------------------------------------------------------------------------------------------------- facts --
N, SRC = {}, {}


def put(key, value, source):
    N[key], SRC[key] = value, source


def f(key):
    v = N[key]
    if isinstance(v, bool) or isinstance(v, str):
        return str(v)
    if isinstance(v, int):
        return "{:,}".format(v)
    return "{:,.1f}".format(v)


def j(*p):
    return os.path.join(*p)


def load(data, *p):
    with open(j(data, *p), encoding="utf-8") as fh:
        return json.load(fh)


def collect(data, examples_manifest=None):
    import pandas as pd
    N.clear(), SRC.clear()
    # --- ORCHID-Values (literature) ---
    S = "parquet/literature_values.parquet"
    v = pd.read_parquet(j(data, "parquet", "literature_values.parquet"))
    put("values", len(v), S), put("accepted", int(v.accepted.sum()), S)
    put("lit_papers", int(v.paper_id.nunique()), S), put("lit_systems", int(v.system.nunique()), S)
    ev = v.evidence_class.value_counts()
    for k in ("experimental", "model", "dft", "unspecified"):
        put("ev_" + k, int(ev.get(k, 0)), S)
    put("tdb_linked", int(v.model_tdb.notna().sum()), S)
    pg = v.property_group.value_counts()
    for k, n in (("invariant reaction", "pg_invariant"), ("phase transition", "pg_transition"),
                 ("phase equilibrium", "pg_equilibrium"), ("formation", "pg_formation"), ("crystal structure", "pg_structure")):
        put(n, int(pg.get(k, 0)), S + " (property_group)")
    put("basis_unstated", int(((v.unit == "kJ/mol (basis unstated)") & (v.property_group == "formation")).sum()), S)
    put("basis_inferred", int(v.flag_basis_inferred_from_magnitude.sum()), S)
    put("lit_biblio", len(load(data, "literature", "papers.json")), "literature/papers.json")
    put("lit_context", len(load(data, "literature", "contexts.json")), "literature/contexts.json")
    put("lit_dups", load(data, "literature", "quality_fixes.json")["duplicates"]["records_removed"], "literature/quality_fixes.json")
    hr = load(TOOLS, HUMAN_REVIEW)
    S = "orchid-project.org tools/%s (review summary sha256 %s)" % (HUMAN_REVIEW, hr["source"]["sha256"][:16])
    put("hr_sample", hr["sample_size"], S), put("hr_reviewed", hr["unique_records_reviewed"], S)
    put("hr_adjudicated", hr["adjudication_queue_records"], S)
    # --- ORCHID-Diagrams ---
    S = "parquet/phase_diagram_index.parquet"
    pi = pd.read_parquet(j(data, "parquet", "phase_diagram_index.parquet"))
    put("pd_files", len(pi), S), put("pd_ok", int((pi.status == "ok").sum()), S), put("pd_failed", int((pi.status != "ok").sum()), S)
    iv = pd.read_parquet(j(data, "parquet", "phase_diagram_invariants.parquet"))
    put("pd_invariants", len(iv), "parquet/phase_diagram_invariants.parquet")
    put("pd_hits", int(iv.n_mapper_hits.sum()), "parquet/phase_diagram_invariants.parquet (sum of n_mapper_hits)")
    put("pd_engine", load(data, "phase-diagrams", "index.json")["engine"], "phase-diagrams/index.json")
    # --- ORCHID-Hulls ---
    for layer, n in (("dft-mp", "mp"), ("dft-oqmd", "oqmd")):
        cat = load(data, "hulls", layer, "catalog.json")
        s, pv, S = cat["summary"], cat["provenance"], "hulls/%s/catalog.json" % layer
        put("h_%s_systems" % n, s["systems"], S), put("h_%s_elements" % n, s["elements"], S)
        put("h_%s_vertices" % n, s["total_vertices"], S), put("h_%s_release" % n, pv["source_release"], S)
    s = load(data, "hulls", "dft-mp", "catalog.json")["summary"]
    put("h_mp_agree", s["systems"] - s["systems_disagreeing_with_mp"], "hulls/dft-mp/catalog.json")
    s = load(data, "hulls", "calphad-sgte", "catalog.json")["summary"]
    put("h_sgte_assess", s["assessments"], "hulls/calphad-sgte/catalog.json"), put("h_sgte_systems", s["systems"], "hulls/calphad-sgte/catalog.json")
    # --- ORCHID-Parity ---
    dp = load(data, "literature", "dft_parity.json")
    S, c = "literature/dft_parity.json", dp["counts"]
    put("par_points", c["points"], S), put("par_mp", c["points_by_source"]["mp"], S), put("par_oqmd", c["points_by_source"]["oqmd"], S)
    put("par_records", c["records_matched"], S), put("par_compounds", c["compounds"], S)
    for src in ("mp", "oqmd"):
        x = dp["stats"][src]["experimental"]
        put("par_%s_n" % src, x["n"], S), put("par_%s_mae" % src, round(x["mae"], 1), S)
        put("par_%s_spearman" % src, "%.2f" % x["spearman"], S)
    k = load(data, "literature", "kim2017_crosscheck.json")
    S = "literature/kim2017_crosscheck.json"
    put("kim_overlap", k["result"]["n_overlap_formulas"], S), put("kim_mae", "%.2f" % k["result"]["mae"], S)
    put("kim_sign", "%d" % round(100 * k["result"]["sign_agreement"]), S), put("kim_differ", k["result"]["n_differ_more_than_10"], S)
    # --- ORCHID-Registry ---
    r = load(data, "phases", "registry.json")["counts"]
    S = "phases/registry.json"
    for key in ("systems", "phases", "names", "structure_types", "liquid", "gas", "tdb_models", "review_queue"):
        put("reg_" + key, r[key], S)
    put("reg_tdb_names", r["names_by_namespace"]["tdb"], S)
    put("reg_tierA", r["assigned_tier_A_standard_name"] + r.get("assigned_tier_A_interstitial_compound", 0), S)
    put("reg_tierB", r["assigned_tier_B_dft_unique_match"], S)
    put("reg_tierC", r["assigned_tier_C_assessment_declared"] + r.get("assigned_tier_C_tdb_comment_structure", 0), S)
    # --- ORCHID-Assemblages (numbers also used by the umbrella list) ---
    a = load(data, "assessed-phase-assemblages", "summary.json")
    S = "assessed-phase-assemblages/summary.json"
    put("asm_points", a["points"], S), put("asm_diagrams", a["diagrams_with_points"], S)
    # --- ORCHID-Values ESPEI datasets ---
    idx = load(data, "espei", "index.json")
    val = load(data, "espei", "validation.json")
    S = "espei/index.json"
    rows = [x for x in idx.get("datasets", []) if x.get("evidence_class") not in ("derived", "assessed_phase_diagram")]
    put("es_datasets", len(rows), S), put("es_values", sum(int(x.get("n_values") or 0) for x in rows), S)
    put("es_systems", len({x.get("system") for x in rows}), S)
    from collections import Counter
    fc = Counter((x.get("folder") or "").split("/")[0] for x in rows)
    for kf in ("experimental", "model", "invariants", "dft"):
        put("es_" + kf, fc.get(kf, 0), S + " (datasets by folder)")
    allr = idx.get("datasets", [])
    put("es_derived", sum(1 for x in allr if x.get("evidence_class") == "derived"), S)
    put("es_assessed", sum(1 for x in allr if x.get("folder") == "assessed"), S)
    put("es_version", val["espei_version"], "espei/validation.json")
    # --- ORCHID-Bibliography ---
    b = load(data, "bibliography", "index.json")["counts"]
    S = "bibliography/index.json"
    put("bib_papers", b["papers"], S), put("bib_doi", b["papers_with_doi"], S), put("bib_nodoi", b["papers_without_doi"], S)
    for key, n in (("corpus_assessment_read", "bib_read"), ("cited_original_source", "bib_cited"),
                   ("tdb_source_assessment", "bib_tdb"), ("handbook_diagram_source", "bib_handbook"), ("discovery_lead", "bib_leads")):
        put(n, b["papers_by_role"][key], S)
    put("bib_edges", b["edges"], S)
    for key in ("cites", "models_system_in_tdb", "reassesses", "handbook_reproduces", "cited_by_crossref"):
        put("bib_e_" + key, b["edges_by_type"][key], S)
    dm = b["doi_match"]
    put("bib_verified", dm["verified"], S), put("bib_weak", dm["weak"], S), put("bib_unverifiable", dm["unverifiable"], S)
    put("bib_dropped", b["doi_mismatches_dropped"]["dois"], S)
    put("bib_crossref", ", ".join(sorted(set(b["crossref_retrieved"]))), S + " (crossref_retrieved)")
    # --- orchid-examples (software) ---
    if examples_manifest:
        S = "zenodo/orchid-examples.files.json"
        put("ex_files", examples_manifest["n_files"], S), put("ex_commit", examples_manifest["commit"][:10], S)
        put("ex_kb", int(round(examples_manifest["total_bytes"] / 1000)), S + " (total_bytes / 1000)")


# --------------------------------------------------------------------------------------------- html pieces --
STATUS_HTML = lambda name: "<p><b>Status.</b> %s</p>" % html.escape(oc.STATUS_BY.get(name, oc.STATUS), quote=False)
LICENCE = ("<h3>Licence</h3><p>ORCHID's data are released under CC BY 4.0, as are the SGTE (B. Hallstedt, Calphad 89 (2025) "
           "102833), Materials Project and OQMD data from which they are derived; attributions and statements of changes "
           "are in NOTICE.md, the licence of each folder in LICENSES.md (both in the bundle).%s</p><p>" +
           html.escape(oc.LICENCE_SCOPE, quote=False) + "</p>")
QUOTES = (" The quotations from papers that accompany values are short excerpts kept for verification; they remain under "
          "the publishers' copyright and are not relicensed. Source papers, page images and full texts are not redistributed.")


def files_html(name):
    return ("<p>Files: one zip (<code>%s</code>) holding the files at their repository paths, a <code>MANIFEST.sha256</code>, "
            "the README, CHANGELOG and licence files. Its file list with sizes and SHA-256 is <code>zenodo/%s.files.json</code> "
            "in the data repository.</p>" % (_bundle(name), name))


def _bundle(name):
    return "%s-%s.zip" % (name if name.startswith("orchid-") else "orchid-" + name, oc.RELEASE_TAG)  # as zenodo_release.bundle_name


def cite_html(name, data):
    e = oc.entry(name, data)
    h = "<h3>How to cite</h3><p>%s</p>" % html.escape(e["text"], quote=False)
    if name != oc.UMBRELLA:
        h += "<p>Part of: %s</p>" % html.escape(oc.entry(oc.UMBRELLA, data)["text"], quote=False)
    also = [oc.THIRD_PARTY[k] for k in oc.DATASETS[name][2]]
    h += "<p>Also cite %s.</p>" % html.escape("; ".join(also), quote=False)
    blocks = [e["bibtex"]] + ([oc.entry(oc.UMBRELLA, data)["bibtex"]] if name != oc.UMBRELLA else [])
    h += "<pre>%s</pre>" % html.escape("\n\n".join(blocks), quote=True)
    return h


import project_credits


def footer_html():
    return project_credits.html_block() + ("<p>Browse: <a href=\"%s\">orchid-project.org</a> · data site <a href=\"%s\">data.orchid-project.org</a> · source "
            "<a href=\"https://github.com/orchid-calphad/orchid-data\">github.com/orchid-calphad/orchid-data</a>.</p><p>%s</p>"
            % (WEBSITE, DATA_URL, oc.FUNDING_NOTE))


def lis(items):
    return "<ul>" + "".join("<li>%s</li>" % x for x in items) + "</ul>"


def I(text):   # fill {key} from the facts with thousands separators
    class D(dict):
        def __missing__(self, k):
            raise KeyError("no fact %r" % k)
    return text.format_map(D({k: f(k) for k in N}))


STATUS_V1 = "Version 1"


def description(name, data):
    t = oc.DATASETS[name][0] if name in oc.DATASETS else oc.SOFTWARE[name][0]
    part = "Part of the ORCHID data release %s." % oc.RELEASE_TAG
    if name == "literature":
        return I("<p><b>%s</b>: thermodynamic and phase-equilibrium values read from published papers. Records retain the "
                 "reported quantity, available conditions, phase labels, evidence class and source locator; some records "
                 "describe transitions, expressions or incompletely specified states. " % t + part +
                 " Short name in the paper: Values.</p>") + STATUS_HTML(name) + \
            "<h3>What</h3>" + I(lis([
                "{values} values from {lit_papers} papers covering {lit_systems} chemical systems; {accepted} carry the accepted "
                "flag of the implemented screening rules.",
                "Evidence classes, never mixed or averaged: experimental {ev_experimental}, model {ev_model} (CALPHAD or other "
                "models; {tdb_linked} linked to the exact TDB file by SHA-256), DFT {ev_dft}, unspecified {ev_unspecified}.",
                "Invariant reactions {pg_invariant}, phase transitions {pg_transition}, phase equilibria {pg_equilibrium}, "
                "formation quantities {pg_formation}, crystal structure {pg_structure}, and heat capacity, entropy, mixing, "
                "activity and EMF values.",
                "Bibliography of {lit_biblio} reporting papers and cited works; experimental context (samples, heat treatments, "
                "measurements) of {lit_context} papers, each item with its quote and page.",
                "Formats: JSON per chemical system (<code>literature/systems/</code>), full JSONL.gz and CSV.gz dumps, Parquet "
                "tables (<code>literature_values</code> with a by-paper train/validation/test split, <code>papers</code>, "
                "<code>experimental_context</code>); fields documented in <code>literature/SCHEMA.md</code>, datasheet in "
                "<code>DATASET_CARD.md</code>."])) + \
            "<h3>How</h3>" + I("<p>Tables and passages of lawfully held papers were read independently by language-model readers "
                               "(Gemini Flash, Claude Opus, Gemini Pro for confirmation) and the PDF text layer; the accepted flag "
                               "requires reader agreement plus deterministic checks: quote verification (the sentence must be verbatim "
                               "on one page and contain the printed value; released quotes are at most 500 characters, a longer one is clipped at export to a contiguous part that still contains the printed value, flag <code>quote_clipped</code>), the field-placement rule (fields not printed in the quote "
                               "or row are dropped), printed-page sign checks of formation enthalpies, and cross-checks against DFT "
                               "(Materials Project, OQMD) and Kim et al. 2017. {lit_dups} exact duplicates were merged; reaction types "
                               "and property codes mapped to controlled vocabularies; every correction is listed with its evidence in "
                               "<code>literature/QUALITY_FIXES.md</code>. The accepted flag records the result of the implemented screening "
                               "rules, and espei_ready records eligibility under the export rules; neither flag certifies that every "
                               "value, phase or measurement condition is correct.</p>") + \
            "<h3>Limitations</h3>" + I(lis([
                "Version 1 is machine-extracted. Preliminary human responses cover {hr_reviewed} distinct records from a frozen "
                "sample of {hr_sample}, followed by a separate AI-assisted adjudication of {hr_adjudicated} records selected for "
                "concerns or comments; these results do not establish population-wide accuracy or independently validate the "
                "repaired release.",
                "Values in figures are not included; coverage across systems is uneven.",
                "{basis_unstated} formation values have a per-atom vs per-formula basis the paper does not state (unit "
                "<code>kJ/mol (basis unstated)</code>, not converted); {basis_inferred} carry a per-atom basis inferred only by a "
                "magnitude heuristic.",
                "Uncertainties with basis <code>stated_default</code> or <code>reading_resolution</code> are placeholders, not "
                "measurements.",
                "<code>experimental</code> values cited from older work carry the citing paper's reading; <code>unspecified</code> "
                "values may be experimental or calculated; check <code>flag_source_disputed</code>, "
                "<code>flag_magnitude_suspect</code> and <code>flag_system_mismatch</code> before use.",
                "<code>flag_espei_ready</code> marks eligibility under the export rules, not a review verdict: do not feed "
                "values to an optimizer without review."])) + LICENCE % QUOTES + files_html(name) + cite_html(name, data) + footer_html()
    if name == "phase-diagrams":
        return I("<p><b>%s</b>: binary phase diagrams calculated with {pd_engine} from the SGTE binary collection (Hallstedt 2025, "
                 "CC BY 4.0), each tied to the exact TDB file by SHA-256 and to the original assessment it encodes (citation and "
                 "DOI from the TDB header). " % t + part + " Short name in the paper: Diagrams.</p>") + \
            "<h3>What</h3>" + I(lis([
                "{pd_files} SGTE binary databases: {pd_ok} calculated, {pd_failed} failed (recorded with status <code>failed</code>).",
                "Per diagram (<code>phase-diagrams/sgte/&lt;TDB&gt;.json</code>): phase-region boundaries, {pd_invariants} distinct "
                "invariant equilibria ({pd_hits} hits of the mapper, which finds an invariant from several starting points; the "
                "Parquet table has one row per distinct invariant) with temperatures, phases and phase compositions over all "
                "diagrams, an equilibrium phase grid, temperature range, TDB SHA-256, original-assessment citation and DOI.",
                "<code>phase-diagrams/index.json</code> lists diagrams per system; Parquet tables <code>phase_diagram_index</code> "
                "and <code>phase_diagram_invariants</code>."])) + \
            I("<h3>How</h3><p>{pd_engine} at 1 atm, run on every SGTE binary file of the ORCHID catalog (open-tdb).</p>") + \
            "<h3>Limitations</h3>" + I(lis([
                "Calculated values are model results, not measurements, and are as good as the assessment and the calculation grid.",
                "{pd_failed} files did not calculate; metastable and pressure-dependent features are not represented.",
                "The checks establish consistency between sampled mapped boundaries and equilibrium solutions from the same "
                "engine; they do not independently validate the assessments or all boundary points.",
                "Requested temperature ranges and the actual extents of stored boundaries are distinct, and a stored point "
                "outside the requested range carries no additional validation."])) + \
            LICENCE % "" + files_html(name) + cite_html(name, data) + footer_html()
    if name == "hulls":
        return I("<p><b>%s</b>: binary formation-energy convex hulls from two DFT databases and isotherms of assessed CALPHAD "
                 "descriptions, stored as ordered composition intervals (interval schema v1) so energy and the equilibrium phase "
                 "set at any composition are a lookup. " % t + part + " Short name in the paper: Hulls.</p>") + \
            "<h3>What</h3>" + I(lis([
                "<code>hulls/dft-mp/</code>: Materials Project release {h_mp_release} (GGA/GGA+U with the MP2020 compatibility "
                "corrections applied by the Materials Project; no further correction added), {h_mp_systems} binary systems of "
                "{h_mp_elements} elements, {h_mp_vertices} hull vertices; above-hull entries within 0.05 eV/atom.",
                "<code>hulls/dft-oqmd/</code>: OQMD {h_oqmd_release} (fit-partial elemental references applied by OQMD), "
                "{h_oqmd_systems} systems, {h_oqmd_vertices} vertices.",
                "<code>hulls/calphad-sgte/</code>: {h_sgte_assess} assessed SGTE binary descriptions ({h_sgte_systems} systems) at "
                "298.15 K, same interval schema, solution phases included.",
                "Per layer: <code>catalog.json</code> (provenance, summary, per-system SHA-256), <code>elements.json</code>, "
                "entries; Parquet table <code>hull_vertices</code>."])) + \
            I("<h3>How</h3><p>Deterministic lower-hull builder (standard library, open-tdb) over pinned snapshots; collinear and "
              "degenerate structures are retained, upstream stability disagreement is enumerated rather than reconciled; "
              "{h_mp_agree} of {h_mp_systems} MP systems reproduce the upstream stability flag exactly. Agreement in a hull "
              "consistency check tests the stated computational representation, not the experimental accuracy of the "
              "energies.</p>") + \
            "<h3>Limitations</h3>" + lis([
                "DFT hulls are at 0 K and carry formation enthalpy only.",
                "Materials Project and OQMD use different reference states and corrections; their energies are not comparable and "
                "must not be averaged.",
                "The MP, OQMD and assessed hulls retain their respective database versions, energy conventions and temperatures; "
                "comparisons across them require explicit reference-state treatment.",
                "OQMD stability is derived by the builder (no upstream cross-check)."]) + \
            LICENCE % "" + files_html(name) + cite_html(name, data) + footer_html()
    if name == "dft-benchmark":
        return I("<p><b>%s</b>: formation enthalpies of binary compounds read from published papers (ORCHID-Values) paired with "
                 "the lowest-energy Materials Project and OQMD entries at exactly the same composition, as separate sources that are "
                 "never averaged. " % t + part + " Short name in the paper: Parity.</p>") + STATUS_HTML(name) + \
            "<h3>What</h3>" + I(lis([
                "{par_points} literature-record x DFT-source points ({par_mp} Materials Project, {par_oqmd} OQMD) from {par_records} "
                "matched records covering {par_compounds} compounds; residual = DFT - literature, kJ/mol-atom (1 eV/atom = 96.485 "
                "kJ/mol-atom).",
                "Default subset (experimental, accepted, signed, not magnitude-suspect): MAE {par_mp_mae} kJ/mol-atom over {par_mp_n} "
                "Materials Project pairs (Spearman {par_mp_spearman}) and {par_oqmd_mae} over {par_oqmd_n} OQMD pairs (Spearman "
                "{par_oqmd_spearman}); statistics by chemistry, method and evidence class, exclusions counted by reason. The "
                "full pair table includes multiple evidence classes; the reported DFT–experiment statistics use the documented "
                "screened experimental subset.",
                "Cross-check against Kim et al. 2017 experimental formation enthalpies: {kim_overlap} formulas in common, MAE "
                "{kim_mae} kJ/mol-atom, sign agreement {kim_sign}%, {kim_differ} formulas differing by more than 10 kJ/mol-atom listed "
                "with record ids. The Kim comparison describes its recorded build and must be regenerated against the frozen "
                "release before its numerical results are used as current validation.",
                "Files: <code>literature/dft_parity.json</code> (rules, statistics, points, records, papers), "
                "<code>literature/kim2017_crosscheck.json|.md</code>, Parquet <code>dft_parity</code>."])) + \
            "<h3>How</h3><p>Formulas are parsed from printed phase labels; reduced compositions (4-decimal mole fractions) are matched " \
            "to non-deprecated DFT entries; large residuals were classified deterministically (per-formula values labelled per atom " \
            "converted, unsigned magnitudes flagged and excluded, the rest flagged and kept; <code>literature/QUALITY_FIXES.md</code>)." \
            "</p><h3>Limitations</h3>" + lis([
                "Literature values are machine-extracted (Version 1); each point keeps its quote and page so it can be checked.",
                "Binary compounds only; values whose per-atom basis is not stated are excluded.",
                "These discrepancies reflect DFT error together with residual uncertainty in extraction, composition, "
                "temperature, provenance and reference state, and are not independent measurements of DFT error alone.",
                "MP and OQMD energies use different references and are not comparable with each other."]) + \
            LICENCE % (QUOTES + " The Kim et al. 2017 values used in the cross-check come from the Figshare deposit %s under the MIT "
                       "licence (notice in NOTICE.md)." % DOI_KIM_DATA) + files_html(name) + cite_html(name, data) + footer_html()
    if name == "phase-registry":
        return I("<p><b>%s</b>, version 0: one entry per phase of each binary system, with every name the phase goes by, the CALPHAD "
                 "models that describe it, where calculated phase diagrams make it stable, example literature references, and a "
                 "crystal-structure type assigned when the stored evidence supports an assignment under the stated tier rules. " % t + part + " Short name in the paper: Registry.</p>") + \
            "<h3>What</h3>" + I(lis([
                "<code>phases/registry.json</code>: compact index (id, system, names, structure tier and prototype, counts) of "
                "{reg_phases} phases in {reg_systems} binary systems, with the full records (models with TDB definitions, stability, "
                "references) in <code>phases/systems/&lt;SYS&gt;.json</code>; id <code>phase:&lt;sorted elements&gt;:&lt;slug&gt;</code> "
                "(the SGTE-corpus canonical convention). Names in namespaces <code>tdb</code> (TDB file and phase name), "
                "<code>paper</code> (DOI and printed label), <code>corpus</code>, <code>mp</code>, <code>oqmd</code>: {reg_names} names.",
                "Models: {reg_tdb_models} TDB phase declarations (file SHA-256, sublattices, site ratios, constituents) from the bundled "
                "SGTE binary collection; stability ranges per calculated diagram; literature value counts; up to 10 example references "
                "per phase (assessments first) with DOI and role.",
                "Structure: {reg_tierA} phases assigned from the standard phase name (tier A), {reg_tierB} from a unique Materials Project "
                "+ OQMD on-hull match (tier B, \"DFT-supported\"), {reg_tierC} assessment-declared (tier C, assessor-stated: exactly one "
                "consistent structure in the Hallstedt 2025 crystallographic table for that phase name and system), {reg_liquid} liquid "
                "and {reg_gas} gas; {reg_review_queue} phases carry review candidates and are marked <code>needs_review</code>.",
                "<code>phases/structure_types.json</code>: {reg_structure_types} structure types (prototype, Strukturbericht, Pearson "
                "symbol, space group), each with its source; also <code>structure_types.csv</code> and "
                "<code>structure_types.parquet</code>; <code>registry.csv</code>, <code>registry.parquet</code>, "
                "<code>registry.jsonld</code> (schema.org DefinedTermSet), <code>review_queue.json</code>, <code>CONFLICTS.md</code>, "
                "<code>SCHEMA.md</code>."])) + \
            "<h3>How</h3><p>Deterministic builder (<code>tools/build_phase_registry.py</code>, standard library) over the bundled TDB " \
            "files, the Hallstedt 2025 crystallographic table, the ORCHID SGTE-corpus phase mappings and literature values, the " \
            "calculated binary phase diagrams and the Materials Project / OQMD hull layers. No model or network calls.</p>" \
            "<h3>Limitations</h3>" + lis([
                "A structure type is assigned when the stored evidence supports an assignment under the stated tier rules (for "
                "tier C, as stated by the assessment's crystallographic table); everything else is a review candidate, not a claim.",
                "Phase identities, aliases and structure assignments retain their evidence tier and review state; assignments "
                "supported by the implemented rules can still require source-specific adjudication.",
                "DFT entries are computed 0 K structures, not the CALPHAD phase; tier B requires agreement of both databases and a "
                "unique structure type for the space group and stoichiometry.",
                "Paper labels come from machine-extracted literature values (Version 1); mapping basis is recorded per name.",
                "Binary systems only; stability comes from calculated (not measured) diagrams."]) + \
            LICENCE % "" + files_html(name) + cite_html(name, data) + footer_html()
    if name == "espei":
        return I("<p><b>%s</b>: ESPEI {es_version} datasets built from the ORCHID literature values (ORCHID-Values) and from DFT, "
                 "one folder per evidence class. " % t + part + "</p>") + STATUS_HTML(name) + \
            "<h3>What</h3>" + I(lis([
                "{es_datasets} datasets with {es_values} values in {es_systems} systems, in folders <code>experimental/</code> "
                "(thermochemical and phase-composition ZPF tables; {es_experimental} datasets), <code>invariants/</code> (ZPF from "
                "invariant-reaction tables; {es_invariants}), <code>model/</code> (values the papers calculated with a CALPHAD model, "
                "TDB SHA-256; {es_model}) and <code>dft/</code> (<code>literature</code>, <code>mp</code>, <code>oqmd</code>: "
                "HM_FORM of TDB phase endmembers mapped through the ORCHID phase registry; {es_dft}).",
                "Literature-derived values retain their ORCHID identifiers and source provenance and are screened against the "
                "literature release; direct MP and OQMD entries retain their database provenance and follow the corresponding "
                "DFT export rules.",
                "Uncertainties are in physical units (printed; ESPEI's experimental default where none was printed; DFT sigma "
                "calibrated against experiment, leave-system-out). Weights follow the exported uncertainty policy; unequal "
                "thermochemical uncertainties are represented by separate fitting records, and ZPF likelihood tolerances are "
                "not obtained by converting temperature uncertainties into thermochemical weights.",
                "Optional phase constraints, tagged and kept in their own folders (never mixed with measured data): "
                "<code>derived/invariant-brackets/</code> and <code>derived/liquid-stability/</code> ({es_derived} datasets) and "
                "<code>assessed/</code> (assessed handbook phase-field points in each assessment's phase names; {es_assessed} datasets).",
                "Validated with ESPEI {es_version} <code>check_dataset</code> and <code>load_datasets</code> (<code>VALIDATION.md</code>); "
                "encoding fixes are counted there, literature records to correct are listed in <code>UPSTREAM_CORRECTIONS.md</code>. "
                "ESPEI validation checks the dataset format and loading behavior, while fitting suitability remains dependent on "
                "state identification, model-phase mapping and source quality."])) + \
            "<h3>Limitations</h3>" + lis([
                "Version 1, machine-extracted literature values.",
                "Model-calculated values carry a method floor calibrated against experimental values of the same state, a lower bound "
                "of the model error because many models were fitted to those experiments: choose your own weight.",
                "0 K DFT energies are used as 298.15 K enthalpies; DFT residuals are heavy-tailed.",
                "Use one DFT source per fit (different elemental references).",
                "Derived phase constraints encode the semantics of printed reactions; they are not measurements."]) + \
            LICENCE % QUOTES + files_html(name) + cite_html(name, data) + footer_html()
    if name == "bibliography":
        return I("<p><b>%s</b>, version 2026.10: papers known to the ORCHID CALPHAD data project, with Crossref bibliographic "
                 "metadata, the role each paper plays in ORCHID, counts of the evidence ORCHID released from it, and typed relations. "
                 "The bibliography combines works found through the documented corpus, citation, TDB and handbook sources; it is "
                 "not an exhaustive bibliography of the field. " % t + part + " Short name in the paper: Bibliography.</p>") + \
            STATUS_HTML(name) + \
            "<h3>What</h3>" + I(lis([
                "{bib_papers} works ({bib_doi} with DOI; {bib_nodoi} without, with stable ids from first author, year and title): "
                "{bib_read} read by ORCHID, {bib_cited} cited as the original source of extracted values, {bib_tdb} assessments behind "
                "TDB files, {bib_handbook} original evaluations behind assessed (handbook) diagrams, {bib_leads} discovery leads "
                "(Crossref higher-order assessment discovery, not read); a paper can have several roles.",
                "{bib_edges} typed edges with provenance: cites / used_in_assessment ({bib_e_cites} each; from the citation printed for "
                "each extracted value), models_system_in_tdb ({bib_e_models_system_in_tdb}), reassesses ({bib_e_reassesses}; citation + "
                "same exact system, never from dates), handbook_reproduces ({bib_e_handbook_reproduces}), cited_by_crossref "
                "({bib_e_cited_by_crossref}; publisher-deposited reference lists, kept separate).",
                "DOIs matched from printed references (handbook source lines, cited sources of values, TDB header citations) are "
                "checked against Crossref (first author, year, title, journal): {bib_verified} verified, {bib_weak} weak (kept, "
                "flagged), {bib_unverifiable} without a Crossref record; {bib_dropped} mismatched DOIs were dropped (works kept without "
                "DOI; listed in DOI_MISMATCHES.md). A verified citation link records the outcome of the stated bibliographic check "
                "and does not certify the measurements or claims of the linked work.",
                "Files: papers.jsonl/.csv/.parquet, edges.jsonl/.csv/.parquet, papers/&lt;key&gt;.json (one record per paper), "
                "systems/&lt;SYS&gt;.json and systems_index.json, BibTeX and CSL-JSON (Crossref metadata), papers.jsonld (schema.org "
                "ScholarlyArticle), SCHEMA.md, DOI_MISMATCHES.md, CITATION.cff."])) + \
            I("<h3>How</h3><p>Deterministic builder (<code>tools/build_bibliography.py</code>) over the ORCHID release files (literature "
              "values, ESPEI bundle, phase registry, assessed phase assemblages), the SGTE-corpus link tables and handbook source lines, "
              "and a cache of Crossref REST API responses (retrieved {bib_crossref}). No model calls; no abstracts, full text or "
              "PDFs.</p>") + \
            "<h3>Limitations</h3>" + lis([
                "Crossref metadata as deposited; doi_match grades DOIs matched from printed references (weak matches are kept and "
                "flagged; mismatches were dropped).",
                "Open-access status only from the Crossref license field.",
                "Discovery leads were not read; their inferred candidate systems are hypotheses.",
                "Version 1, machine-extracted: evidence counts derive from machine-extracted values."]) + \
            LICENCE % " Crossref bibliographic metadata: CC0. ORCHID-derived fields (roles, evidence counts, relations, checks): CC BY 4.0." + \
            files_html(name) + cite_html(name, data) + footer_html()
    if name == oc.UMBRELLA:
        rows = []
        heads = {"literature": "{values} values from {lit_papers} papers", "phase-diagrams": "{pd_ok} calculated SGTE binary diagrams",
                 "hulls": "{h_mp_systems} MP and {h_oqmd_systems} OQMD binary hulls, {h_sgte_assess} SGTE descriptions",
                 "dft-benchmark": "{par_points} DFT-vs-experiment points", "phase-registry": "{reg_phases} phases in {reg_systems} systems",
                 "assessed-phase-assemblages": "{asm_points} points in {asm_diagrams} assessed diagrams",
                 "espei": "{es_datasets} ESPEI datasets", "bibliography": "{bib_papers} works, {bib_edges} typed edges"}
        for n in ("literature", "espei", "phase-diagrams", "hulls", "dft-benchmark", "phase-registry",
                  "assessed-phase-assemblages", "bibliography"):
            e = oc.entry(n, data)
            rows.append("<a href=\"https://doi.org/%s\">%s</a> (doi:%s): %s" % (e["doi"], html.escape(e["title"], quote=False), e["doi"], I(heads[n])))
        ex = "<a href=\"https://doi.org/%s\">%s</a> (doi:%s)" % (_doi(data, "orchid-examples"), html.escape(oc.SOFTWARE["orchid-examples"][0], quote=False), _doi(data, "orchid-examples"))
        return ("<p><b>%s</b> (Open Repository for CALPHAD Hypotheses, Intelligence, and Discovery), release %s: umbrella record of "
                "eight citable datasets, each a separate Zenodo record with its own DOI (the ESPEI datasets are derived from "
                "ORCHID-Values and deposited apart because they are used apart):</p>" % ("ORCHID data", oc.RELEASE_TAG)) + lis(rows) + \
            "<p>The runnable workflows that use the release are archived as a software record: %s.</p>" % ex + \
            "<p>This umbrella links eight dataset records and the examples software archive within one versioned ORCHID release; " \
            "several records are derived from the same observations and should not be treated as independent evidence.</p>" + \
            STATUS_HTML("literature") + \
            "<p>This record holds the release documentation: README, CHANGELOG, DATASET_CARD (datasheet with loading examples for pandas, " \
            "DuckDB and Hugging Face <code>datasets</code>), licences and attributions, the Parquet manifest, and the Zenodo metadata of " \
            "all ten records with the list of every number in it (<code>zenodo/metadata_numbers.json</code>). The data files are in the " \
            "eight dataset records.</p><h3>Limitations</h3><p>See each dataset record and DATASET_CARD.md (known limitations). " \
            "Literature values are machine-extracted (Version 1); calculated and DFT values are model results, not measurements.</p>" + \
            LICENCE % (QUOTES + " Kim et al. 2017 cross-check values: MIT (NOTICE.md); Crossref metadata in bibliography/: CC0.") + \
            files_html(name) + cite_html(name, data) + footer_html()
    if name == "orchid-examples":
        return ("<p><b>%s</b>: runnable, documented examples of the workflows described in the ORCHID paper: tuning an MLIP-informed Cu-Zr CALPHAD description with gradient-based ESPEI against "
                "ORCHID data (PhaseForge+), the classic ESPEI MCMC route as a documented negative example, prevalence-weighted scoring of "
                "one-shot MLIP phase diagrams against ORCHID evidence, the saturation of MLIP formation-enthalpy accuracy at the DFT level, "
                "PhaseForge input structures, and a quickstart notebook. Cached examples can be inspected without rerunning the "
                "expensive calculations, while full model generation and fitting require the dependencies and resources documented "
                "for each workflow. Each example pins the ORCHID release commit and its software versions. Software record of release %s; the data it uses are in the ORCHID data records (umbrella: doi:%s).</p>"
                % ("ORCHID examples", oc.RELEASE_TAG, _doi(data, oc.UMBRELLA))) + \
            I("<p>Contents: {ex_files} files ({ex_kb} kB), the git archive of commit <code>{ex_commit}</code> of "
              "<a href=\"%s\">github.com/orchid-calphad/orchid-examples</a>. The archive is built from the commit recorded in "
              "its manifest; reproduction results apply only to the tested commit and the documented cached paths.</p>"
              % oc.REPO_EXAMPLES) + \
            "<h3>Licence</h3><p>Code (all <code>.py</code> and <code>.sh</code> files, notebook code cells, CI workflow): MIT (<code>LICENSE</code>). " \
            "Data in the repository keep the licences of their sources (<code>LICENSE-DATA.md</code>): values derived from ORCHID and the SGTE " \
            "TDB files are CC BY 4.0, with attribution to the ORCHID release and to Hallstedt (2025). The record's licence field states the " \
            "code licence (MIT).</p>" \
            "<h3>Limitations</h3>" + lis(["Cached outputs are small summaries; the bulky raw outputs of the full runs are not part of this record.",
                                          "Refined Cu-Zr TDBs are review-only demonstration models, not assessments."]) + \
            "<h3>How to cite</h3><p>%s. Please also cite the ORCHID data release you used, the ORCHID paper, the software each example credits and the original sources of the data.</p>" \
            % html.escape("R. Arróyave, J. Saal, D. Shin, %s (version %s), Zenodo (%s). https://doi.org/%s" % (t, oc.VERSION, oc.YEAR, _doi(data, "orchid-examples")), quote=False) + footer_html()
    raise KeyError(name)


def _doi(data, name):
    return oc.dois(data)[name]


# --------------------------------------------------------------------------------------------- record shell --
KEYWORDS = {
    "orchid-data": ["CALPHAD", "thermodynamics", "phase diagrams", "phase equilibria", "formation enthalpy", "DFT", "convex hull",
                    "materials science", "alloys", "AI-ready data"],
    "literature": ["CALPHAD", "thermodynamics", "phase equilibria", "invariant reactions", "formation enthalpy", "heat capacity",
                   "information extraction", "scientific literature", "alloys", "AI-ready data"],
    "phase-diagrams": ["CALPHAD", "phase diagrams", "pycalphad", "SGTE", "thermodynamic database", "TDB", "invariant reactions", "alloys"],
    "hulls": ["convex hull", "formation energy", "DFT", "Materials Project", "OQMD", "CALPHAD", "SGTE", "phase stability"],
    "dft-benchmark": ["formation enthalpy", "DFT", "Materials Project", "OQMD", "benchmark", "calorimetry", "experimental thermochemistry", "CALPHAD"],
    "phase-registry": ["phase registry", "crystal structure", "prototype", "Strukturbericht", "CALPHAD", "SGTE", "phase diagrams", "DFT",
                       "Materials Project", "OQMD", "AI-ready data"],
    "espei": ["ESPEI", "CALPHAD", "thermodynamics", "phase diagrams", "ZPF", "DFT", "formation enthalpy", "uncertainty", "AI-ready data"],
    "bibliography": ["bibliography", "citation graph", "CALPHAD", "thermodynamic assessment", "phase diagrams", "provenance", "Crossref",
                     "AI-ready data"],
    "orchid-examples": ["CALPHAD", "ESPEI", "pycalphad", "PhaseForge", "machine-learned interatomic potentials", "phase diagrams", "benchmark"],
}
URL = lambda u, rel, **kw: dict({"identifier": u, "relation": rel, "scheme": "url"}, **kw)
DOI = lambda d, rel, rt=None: dict({"identifier": d, "relation": rel, "scheme": "doi"}, **({"resource_type": rt} if rt else {}))


def related(name, data):
    d = oc.dois(data)
    U = d[oc.UMBRELLA]
    site = [URL(WEBSITE, "isDocumentedBy"), URL(DATA_URL, "isDocumentedBy")]
    if name == oc.UMBRELLA:
        r = [DOI(d[n], "hasPart", "dataset") for n in ("literature", "phase-diagrams", "hulls", "dft-benchmark", "phase-registry",
                                                       "assessed-phase-assemblages", "espei", "bibliography")]
        r += [DOI(d["orchid-examples"], "isReferencedBy", "software")] + site
        r += [DOI(x, "references") for x in (DOI_SGTE, DOI_MP, DOI_OQMD1, DOI_OQMD2)]
    elif name == "orchid-examples":
        r = [DOI(U, "references", "dataset"), URL(oc.REPO_EXAMPLES, "isSupplementTo"), URL(DATA_URL, "references")]
        r += [DOI(x, "references") for x in ("10.1016/j.actamat.2025.121453", "10.1016/j.actamat.2026.122412", "10.1557/mrc.2019.59",
                                             "10.5281/zenodo.15730912", "10.5281/zenodo.17519087")]
    else:
        r = [DOI(U, "isPartOf", "dataset")]
        r += {"literature": [], "phase-diagrams": [DOI(DOI_SGTE, "references")],
              "hulls": [DOI(x, "references") for x in (DOI_MP, DOI_OQMD1, DOI_OQMD2, DOI_SGTE)],
              "dft-benchmark": [DOI(d["literature"], "isDerivedFrom", "dataset"), DOI(d["hulls"], "references", "dataset")] +
                               [DOI(x, "references") for x in (DOI_MP, DOI_OQMD1, DOI_OQMD2, DOI_KIM, DOI_KIM_DATA)],
              "phase-registry": [DOI(x, "references") for x in (DOI_MP, DOI_OQMD1, DOI_OQMD2, DOI_SGTE)],
              "espei": [DOI(d["literature"], "isDerivedFrom", "dataset"), DOI(d["hulls"], "isDerivedFrom", "dataset"),
                        DOI(d["phase-registry"], "references", "dataset")] + [DOI(x, "references") for x in (DOI_MP, DOI_OQMD1, DOI_OQMD2)],
              "bibliography": [DOI(d["literature"], "isDerivedFrom", "dataset"), DOI(d["phase-registry"], "references", "dataset"),
                               DOI(DOI_SGTE, "references"), URL(DATA_URL + "/bibliography/", "isIdenticalTo", resource_type="dataset")],
              }[name]
        if name != "bibliography":
            r += site
    if oc.PAPER:   # orchid_cite.PAPER, set when the paper has an arXiv id or DOI
        r.append({"identifier": oc.PAPER["identifier"], "relation": "isSupplementTo", "scheme": oc.PAPER["scheme"]})
    existing = {(e["identifier"], e["relation"]) for e in r}
    r += [e for e in project_credits.related_identifiers() if (e["identifier"], e["relation"]) not in existing]
    return r


def record(name, data):
    software = name in oc.SOFTWARE
    e = None if software else oc.entry(name, data)
    title = oc.SOFTWARE[name][0] if software else e["title"]
    notes = oc.FUNDING_NOTE
    if not software and name in oc.EXTRACTED:
        notes += " " + oc.STATUS_BY.get(name, oc.STATUS)
    if software:
        notes += (" Code is MIT-licensed; data files keep their sources' licences (LICENSE-DATA.md). Software record of the "
                  "ORCHID release %s; zip built with git archive from a pinned commit of the repository." % oc.RELEASE_TAG)
    meta = {"title": title, "upload_type": "software" if software else "dataset", "version": oc.VERSION,
            "description": description(name, data),
            "creators": [{"name": "%s, %s" % (f_, g), "affiliation": af, "orcid": o} for f_, g, af, o in oc.AUTHORS],
            "keywords": KEYWORDS[name], "license": "mit" if software else "cc-by-4.0", "access_right": "open", "language": "eng",
            "communities": [{"identifier": oc.ZENODO_COMMUNITY}], "related_identifiers": related(name, data), "notes": notes}
    return {"_comment": "Zenodo deposition metadata (legacy REST format: PUT /api/deposit/depositions/<id> with this object as the body, "
                        "or enter the fields in the web form). Written by tools/build_zenodo_metadata.py (numbers: metadata_numbers.json); "
                        "do not edit by hand. No `grants` payload: enter both NSF awards in the Funding form; use a custom "
                        "award if no matching indexed award exists, see ZENODO_FUNDING_MANUAL.md. publication_date is set by Zenodo at publication. DOI placeholders "
                        "unassigned per-record DOI are replaced by tools/zenodo_release.py writeback.",
            "record": name, "doi": _doi(data, name), "metadata": meta}


def write_sheets(data, outdir):
    """One HTML sheet per record: every web-form field in the order of the form, the description as rich text to copy."""
    os.makedirs(outdir, exist_ok=True)
    E = html.escape
    for name in NAMES_ALL:
        p = j(data, "zenodo", name + ".metadata.json")
        m = json.load(open(p, encoding="utf-8"))
        md = m["metadata"]
        rel = "".join("<tr><td>%s</td><td><code>%s</code></td><td>%s</td></tr>" % (E(r["relation"]), E(r["identifier"]), E(r.get("resource_type", "")))
                      for r in md["related_identifiers"])
        body = ("<!doctype html><meta charset=utf-8><title>%s: Zenodo form sheet</title><style>body{font:15px/1.45 system-ui;max-width:56em;"
                "margin:2em auto;padding:0 1em}td,th{border:1px solid #bbb;padding:3px 8px;vertical-align:top;text-align:left}table{border-collapse:collapse}"
                "code,pre{background:#f3f3f3}.d{border:2px solid #36a;padding:0 1em;margin:1em 0}h2{margin-top:1.4em}</style>"
                "<h1>Zenodo form sheet: %s</h1><p>Generated from <code>zenodo/%s.metadata.json</code>. Copy field by field; select the "
                "blue box (the description) and paste it into the Description editor.</p>" % (E(name), E(name), E(name)))
        body += "<table><tr><th>Field</th><th>Value</th></tr>"
        body += "<tr><td>Reserved DOI</td><td><code>%s</code></td></tr>" % E(m["doi"])
        body += "<tr><td>Resource type</td><td>%s</td></tr>" % ("Software" if md["upload_type"] == "software" else "Dataset")
        body += "<tr><td>Title</td><td>%s</td></tr>" % E(md["title"])
        body += "<tr><td>Creators (this order)</td><td>%s</td></tr>" % "<br>".join(
            "%s; ORCID %s; %s" % (E(c["name"]), c["orcid"], E(c["affiliation"])) for c in md["creators"])
        body += "<tr><td>Publication date</td><td>TODO: the day of publication (set when you publish)</td></tr>"
        body += "<tr><td>Version</td><td>%s</td></tr>" % E(md["version"])
        body += "<tr><td>Licence</td><td>%s</td></tr>" % ("MIT License" if md["license"] == "mit" else "Creative Commons Attribution 4.0 International")
        body += "<tr><td>Keywords (one each)</td><td>%s</td></tr>" % "<br>".join(E(k) for k in md["keywords"])
        body += "<tr><td>Language</td><td>English</td></tr><tr><td>Access</td><td>Open (public, no embargo)</td></tr>"
        body += "<tr><td>Community</td><td><code>%s</code></td></tr>" % E(md["communities"][0]["identifier"])
        body += "<tr><td>Funding</td><td>NSF awards 2608825 and 2534344; see ZENODO_FUNDING_MANUAL.md</td></tr>"
        body += "<tr><td>Additional notes</td><td>%s</td></tr></table>" % E(md["notes"])
        body += "<h2>Related works</h2><table><tr><th>Relation</th><th>Identifier</th><th>Resource type</th></tr>%s</table>" % rel
        body += "<h2>Description</h2><div class=d>%s</div>" % md["description"]
        with open(j(outdir, name + ".form.html"), "w", encoding="utf-8") as fh:
            fh.write(body)
    print("wrote %d form sheets to %s" % (len(NAMES_ALL), outdir))


def dump(obj):
    return json.dumps(obj, indent=1, ensure_ascii=False) + "\n"


def run(data, check=False, export=None):
    mp = j(data, "zenodo", "orchid-examples.files.json")
    em = json.load(open(mp, encoding="utf-8")) if os.path.exists(mp) else None
    collect(data, em)
    if em is None:
        put("ex_files", 0, "n/a"), put("ex_commit", "pending", "n/a"), put("ex_kb", 0, "n/a")
    out = {}
    for name in NAMES_ALL:
        if name == BUILT_ELSEWHERE:
            continue
        out[j("zenodo", name + ".metadata.json")] = dump(record(name, data))
    nums = {"_comment": "Every number in the descriptions of zenodo/*.metadata.json (except ORCHID-Assemblages, written by "
                        "tools/build_assessed_assemblages.py) with the release file it was read from. Written by "
                        "tools/build_zenodo_metadata.py.",
            "numbers": {k: {"value": N[k], "source": SRC[k]} for k in sorted(N)}}
    out[j("zenodo", "metadata_numbers.json")] = dump(nums)
    bad = []
    for rel, text in out.items():
        p = j(data, rel)
        cur = open(p, encoding="utf-8").read() if os.path.exists(p) else None
        if cur != text:
            if check:
                bad.append(rel)
            else:
                with open(p + ".tmp", "w", encoding="utf-8") as fh:
                    fh.write(text)
                os.replace(p + ".tmp", p)
                print("wrote " + rel)
    if export:
        write_sheets(data, export)
    # consistency of the builder-written record
    # consistency of the builder-written record (title/community/grants)
    ap = j(data, "zenodo", BUILT_ELSEWHERE + ".metadata.json")
    if os.path.exists(ap):
        m = json.load(open(ap, encoding="utf-8"))["metadata"]
        if m["title"] != oc.DATASETS[BUILT_ELSEWHERE][0] or m.get("grants") or \
                m["communities"] != [{"identifier": oc.ZENODO_COMMUNITY}]:
            bad.append("zenodo/%s.metadata.json (title / community / grants differ from orchid_cite; rerun build_assessed_assemblages.py)"
                       % BUILT_ELSEWHERE)
    if check:
        if bad:
            print("stale Zenodo metadata: " + ", ".join(bad), file=sys.stderr)
            return 1
        print("zenodo metadata current (%d records, %d numbers traced)" % (len(NAMES_ALL), len(N)))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--data", default=oc.DEFAULT_DATA)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    return run(a.data, a.check)


if __name__ == "__main__":
    sys.exit(main())
