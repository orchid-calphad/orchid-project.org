#!/usr/bin/env python3
"""How to cite the ORCHID data release: one shared definition for every builder.

The release is archived on Zenodo as eight citable datasets plus an umbrella record ("ORCHID data") that links
them (DataCite isPartOf / hasPart). Their DOIs live in one file of the data repository,

  <orchid-data>/zenodo/dois.json      {"orchid-data": "10.5281/zenodo....", "literature": ..., ...}

which holds placeholders ("10.5281/zenodo.PENDING-<name>") until the PI reserves the DOIs on Zenodo and runs
`tools/zenodo_release.py writeback`. Builders call `cite_as(name)` and embed the result, so a rebuild after the
writeback reproduces exactly what the text replacement wrote and every `--check` keeps passing.

Standard library only; no network.
"""
from __future__ import annotations

import json
import os

SITE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_DATA = os.environ.get("ORCHID_DATA") or os.path.join(os.path.dirname(SITE), "orchid-data")

VERSION = "2026.10"
RELEASE_TAG = "v2026.10"
YEAR = "2026"
STATUS = ("Version 1, machine-extracted: values were read from papers by language-model readers constrained to quote "
          "their source and checked automatically (quote verification, field-placement rule, printed-page sign checks, "
          "DFT and Kim et al. 2017 cross-checks); a hand-checked precision sample is in progress.")
PLACEHOLDER_PREFIX = "10.5281/zenodo.PENDING-"
UMBRELLA = "orchid-data"
LICENSE = "CC-BY-4.0"
WEBSITE = "https://orchid-project.org"
DATA_URL = "https://data.orchid-project.org"

AUTHORS = [  # (family, given, affiliation, ORCID)
    ("Arróyave", "Raymundo", "Texas A&M University", "0000-0001-7548-8686"),
    ("Saal", "James", "Citrine Informatics", "0000-0003-0935-158X"),
    ("Shin", "Dongwon", "Oak Ridge National Laboratory", "0000-0002-5797-3423"),
]
AUTHORS_TEXT = "R. Arróyave, J. Saal, D. Shin"
AUTHORS_BIBTEX = "Arr{\\'o}yave, Raymundo and Saal, James and Shin, Dongwon"

THIRD_PARTY = {
    "sgte": "B. Hallstedt, The SGTE collection of binary datasets, Calphad 89 (2025) 102833, "
            "https://doi.org/10.1016/j.calphad.2025.102833 (CC BY 4.0), and the original assessment named in each file",
    "mp": "A. Jain et al., The Materials Project: A materials genome approach to accelerating materials innovation, "
          "APL Materials 1 (2013) 011002, https://doi.org/10.1063/1.4812323 (CC BY 4.0)",
    "oqmd": "J. E. Saal et al., JOM 65 (2013) 1501, https://doi.org/10.1007/s11837-013-0755-4 and S. Kirklin et al., "
            "npj Comput. Mater. 1 (2015) 15010, https://doi.org/10.1038/npjcompumats.2015.10 (CC BY 4.0)",
    "kim2017": "G. Kim, S. V. Meschel, P. Nash, W. Chen, Experimental formation enthalpies for intermetallic phases "
               "and other inorganic compounds, Scientific Data 4 (2017) 170162, https://doi.org/10.1038/sdata.2017.162",
    "assessed": "for every diagram used, its original evaluation (citations[role=original_evaluation], DOI where matched) "
                "and H. Okamoto, Desk Handbook: Phase Diagrams for Binary Alloys, 2nd ed., ASM International (2010) as the "
                "compilation",
    "papers": "the original source of every value you use (paper_citation / paper_doi, and cited_dois for values "
              "the paper cites from earlier work)",
    "papers_espei": ("the original source of every value you use (each dataset's `reference` / `bibtex`, "
                     "references.bib; index.csv lists the DOIs)"),
    "papers_bibliography": "the papers themselves when you use their data (doi / citation of each bibliography record)",
}

# name -> (title, BibTeX key, sources to cite as well); EXTRACTED: records holding values read from papers
EXTRACTED = {"orchid-data", "literature", "dft-benchmark", "phase-registry", "assessed-phase-assemblages", "espei",
             "bibliography"}
# per-record status text where the generic STATUS (literature readers) does not describe the record
STATUS_BY = {"assessed-phase-assemblages": (
    "Version 1, machine-extracted: phase fields of expert-drawn (assessed) diagrams named by deterministic rules, "
    "model readers, an adjudicator and partial expert review (basis and review state per record); one point deep "
    "inside each field; boundary and tie-line positions are drawing-based and of low credence."),
    "bibliography": (
    "Version 1, machine-extracted: roles and evidence counts derive from machine-extracted literature values and link "
    "tables; bibliographic fields are Crossref metadata as deposited (DOIs matched from printed references checked "
    "against Crossref: doi_match; mismatches dropped).")}
ABSTRACT = {
    "orchid-data": "Umbrella record of the ORCHID data release: literature values, calculated binary phase diagrams, "
                   "binary convex hulls, the DFT-vs-experiment formation-enthalpy benchmark, the phase registry, "
                   "assessed phase assemblages, ESPEI datasets and the bibliography.",
    "literature": "Thermodynamic and phase-equilibrium values read from published papers, each with its evidence "
                  "class, state, source paper, DOI, page and verbatim quote.",
    "phase-diagrams": "Binary phase diagrams calculated with pycalphad from the SGTE binary collection (Hallstedt "
                      "2025), each tied to its TDB file by SHA-256.",
    "hulls": "Binary formation-energy convex hulls from the Materials Project and OQMD (0 K, DFT) and isotherms of "
             "assessed CALPHAD descriptions (SGTE), as interval documents.",
    "dft-benchmark": "Literature formation enthalpies paired with the lowest-energy Materials Project and OQMD entries "
                     "at the same composition, with error statistics and a Kim et al. 2017 cross-check.",
    "phase-registry": "One entry per phase of each binary system: every name it goes by (TDB, corpus, paper labels, "
                      "Materials Project / OQMD ids), the CALPHAD models that describe it, where the calculated "
                      "diagrams make it stable, example literature references, and a structure type assigned only "
                      "when certain (standard name or unique DFT match), otherwise a review candidate list.",
    "assessed-phase-assemblages": "One (T, x, P) point deep inside each named phase field of the expert-drawn (assessed) "
                                  "binary phase diagrams compiled in Okamoto's Desk Handbook, with the phase assemblage "
                                  "present, phase-registry ids, coordinate uncertainty, naming basis, confidence, review "
                                  "state, a cross-check against calculated SGTE diagrams and ESPEI ZPF datasets; evidence "
                                  "class assessed_phase_diagram.",
    "espei": ("ESPEI-format datasets (thermochemical, phase-equilibrium ZPF and invariant-reaction ZPF) built from "
              "the ORCHID literature values, one folder per evidence class (experimental, model, DFT: literature, "
              "Materials Project, OQMD), every value traced to its published ORCHID observation id and source DOI, "
              "with its uncertainty in physical units (printed, ESPEI's experimental default, or calibrated for "
              "DFT against experiment) and the ESPEI weight derived from it; plus optional, tagged phase-constraint "
              "datasets (derived invariant-temperature brackets and liquid-stability points, assessed handbook "
              "assemblages) in the phase names of each assessment."),
    "bibliography": ("Every paper ORCHID knows about: Crossref bibliographic metadata (CC0), the role of each paper in "
                     "ORCHID, counts of the evidence released from it, and typed relations with provenance (citations of "
                     "values, TDB files, reassessments, handbook diagrams, Crossref reference lists)."),
}
# name -> (title, BibTeX key, sources to cite as well)
DATASETS = {
    "orchid-data": ("ORCHID data", "orchid_data_2026", ["papers", "sgte", "mp", "oqmd"]),
    "literature": ("ORCHID literature values", "orchid_literature_2026", ["papers"]),
    "phase-diagrams": ("ORCHID calculated binary phase diagrams", "orchid_phase_diagrams_2026", ["sgte"]),
    "hulls": ("ORCHID binary convex hulls (Materials Project, OQMD, assessed CALPHAD)", "orchid_hulls_2026",
              ["mp", "oqmd", "sgte"]),
    "dft-benchmark": ("ORCHID DFT-vs-experiment formation-enthalpy benchmark", "orchid_dft_benchmark_2026",
                      ["papers", "mp", "oqmd", "kim2017"]),
    "phase-registry": ("ORCHID phase registry (binary systems)", "orchid_phase_registry_2026",
                       ["sgte", "mp", "oqmd", "papers"]),
    "assessed-phase-assemblages": ("ORCHID assessed phase assemblages (binary systems)", "orchid_assessed_assemblages_2026",
                                   ["assessed", "sgte"]),
    "espei": ("ORCHID ESPEI datasets from the literature", "orchid_espei_2026", ["papers_espei", "mp", "oqmd"]),
    "bibliography": ("ORCHID bibliography", "orchid_bibliography_2026", ["papers_bibliography"]),
}
NAMES = tuple(DATASETS)


def placeholder(name: str) -> str:
    return PLACEHOLDER_PREFIX + name


def dois(data_root: str | None = None) -> dict:
    """name -> DOI from <orchid-data>/zenodo/dois.json (placeholders for anything missing)."""
    out = {n: placeholder(n) for n in NAMES}
    p = os.path.join(data_root or DEFAULT_DATA, "zenodo", "dois.json")
    if os.path.exists(p):
        with open(p, encoding="utf-8") as fh:
            d = json.load(fh)
        out.update({k: v for k, v in (d.get("dois") or {}).items() if k in out and v})
    return out


def is_placeholder(doi: str) -> bool:
    return doi.startswith(PLACEHOLDER_PREFIX)


def _text(name, doi):
    title = DATASETS[name][0]
    tag = "; Version 1, machine-extracted" if name in EXTRACTED else ""
    return "%s, %s (version %s%s), Zenodo (%s). https://doi.org/%s" % (AUTHORS_TEXT, title, VERSION, tag, YEAR, doi)


def _bibtex(name, doi):
    title, key, _ = DATASETS[name]
    return ("@misc{%s,\n"
            "  author    = {%s},\n"
            "  title     = {{%s}},\n"
            "  version   = {%s},\n"
            "  publisher = {Zenodo},\n"
            "  year      = {%s},\n"
            "  doi       = {%s},\n"
            "  url       = {https://doi.org/%s}%s\n"
            "}") % (key, AUTHORS_BIBTEX, title, VERSION, YEAR, doi, doi,
                    ",\n  note      = {Version 1, machine-extracted}" if name in EXTRACTED else "")


def entry(name: str, data_root: str | None = None) -> dict:
    """Citation of one record (text, BibTeX, DOI)."""
    doi = dois(data_root)[name]
    return {"title": DATASETS[name][0], "doi": doi, "url": "https://doi.org/" + doi, "text": _text(name, doi),
            "bibtex": _bibtex(name, doi)}


def cite_as(name: str, data_root: str | None = None) -> dict:
    """Top-level `cite_as` object embedded in JSON outputs: this dataset, the umbrella, and what else to cite."""
    d = {"dataset": name, "version": VERSION, "release_tag": RELEASE_TAG, "license": LICENSE,
         "status": STATUS_BY.get(name, STATUS) if name in EXTRACTED else "Version 1 (calculated / third-party derived data)",
         **entry(name, data_root),
         "also_cite": [THIRD_PARTY[k] for k in DATASETS[name][2]]}
    if name != UMBRELLA:
        d["part_of"] = entry(UMBRELLA, data_root)
    return d


def bibtex_block(names, data_root: str | None = None) -> str:
    """BibTeX of the named records (umbrella appended when missing), blank-line separated."""
    names = list(names) + ([UMBRELLA] if UMBRELLA not in names else [])
    return "\n\n".join(entry(n, data_root)["bibtex"] for n in names)


def cff(name: str, data_root: str | None = None, extra_message: str = "") -> str:
    """CITATION.cff (1.2.0) of one record."""
    e = entry(name, data_root)
    lines = ["cff-version: 1.2.0",
             'message: "If you use these data, please cite this Zenodo record and the original sources of the data%s."'
             % extra_message,
             'title: "%s"' % e["title"], "type: dataset", "authors:"]
    for fam, giv, aff, orcid in AUTHORS:
        lines += ['  - family-names: "%s"' % fam, '    given-names: "%s"' % giv,
                  '    orcid: "https://orcid.org/%s"' % orcid, '    affiliation: "%s"' % aff]
    lines += ['version: "%s"' % VERSION, 'doi: "%s"' % e["doi"], "license: %s" % LICENSE,
              'url: "%s"' % DATA_URL, 'repository-code: "https://github.com/orchid-calphad/orchid-data"',
              'abstract: "%s"' % (ABSTRACT[name] + (" " + STATUS_BY.get(name, STATUS) if name in EXTRACTED else "")).replace('"', "'")]
    if name != UMBRELLA:
        u = entry(UMBRELLA, data_root)
        lines += ["identifiers:", "  - type: doi", '    value: "%s"' % u["doi"],
                  '    description: "Umbrella record ORCHID data (all datasets of release %s)"' % RELEASE_TAG]
    lines.append("keywords: [CALPHAD, thermodynamics, phase diagrams, convex hull, DFT, formation enthalpy, AI-ready data]")
    return "\n".join(lines) + "\n"


def citation_txt(name: str, data_root: str | None = None, intro: str = "") -> str:
    """Plain-text CITATION.txt placed beside files that cannot carry metadata (CSV)."""
    c = cite_as(name, data_root)
    o = []
    if intro:
        o += [intro, ""]
    o += ["Please cite:", "", c["text"], "", "Part of: " + c["part_of"]["text"] if "part_of" in c else "", "",
          "Also cite " + "; ".join(c["also_cite"]) + ".", "", "Status: " + c["status"], "", "Licence: CC BY 4.0 "
          "(derived data; quoted sentences under fair use; source papers not redistributed).", "", "BibTeX:", "",
          bibtex_block([name], data_root), ""]
    return "\n".join(o)
