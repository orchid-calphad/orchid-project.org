# NOTICE — sources, attribution, and licenses

This repository catalogs thermodynamic database (`.tdb`) files from the public
sources below. Each file is governed by its own upstream terms. The MIT license
in `LICENSE` covers only the catalog, scripts, and documentation authored here,
except where a third-party-derived artifact is explicitly identified below.
`license_status` records whether redistribution evidence is `confirmed`,
`uncertain`, `restricted`, or `conflicting`. Public availability alone is not
permission to redistribute.

## Redistributed in this repository

### OpenCalphad diagnostic instrumentation patch — GPL-3.0

Source: https://github.com/sundmanbo/opencalphad at commit
`7233d49bbbab160711e7f0545d78855cac885522`.

`catalog/compatibility/opencalphad-6.121-phase-audit.patch` is a minimal patch
against OpenCalphad's GPL-3.0 source. The patch and included source context are
distributed under the same GPL-3.0 terms, not this repository's MIT license.
It adds two diagnostic writes for reproducible parser auditing and does not
alter parsing decisions. No OpenCalphad source tree or executable is vendored.
Upstream license: https://github.com/sundmanbo/opencalphad/blob/master/LICENSE

### SGTE binary collection — CC BY 4.0 (2,774 files)

Source and citation: Bengt Hallstedt, “The SGTE collection of binary datasets,”
*Calphad* 89 (2025) 102833.
https://doi.org/10.1016/j.calphad.2025.102833

The article and its supplementary release are available under the Creative
Commons Attribution 4.0 International license:
https://creativecommons.org/licenses/by/4.0/

The files are redistributed unchanged under `databases/sgte-binary/`. Retain
the headers and cite both Hallstedt and the relevant original assessment.
Archive hashes and additional provenance are in that directory's `SOURCE.md`.

### Two open-access assessment papers — CC BY (worked example)

`examples/thermo-extractor-worked-example/papers/` redistributes two articles
unchanged so the `thermo-extractor` walkthrough
(`notebooks/thermo_extractor_walkthrough.ipynb`) runs without downloading
anything. Full attribution and hashes: that directory's `ATTRIBUTION.md`.

- Dębski, Terlicka, Gąsior, Gierlotka, Pęska, Dworecka-Wójcik, Polański,
  “Calorimetric Studies of Magnesium-Rich Mg-Pd Alloys,” *Materials* 14 (2021)
  680, https://doi.org/10.3390/ma14030680 — Creative Commons Attribution 4.0,
  https://creativecommons.org/licenses/by/4.0/ (MDPI).
- Gao, Li, Guo, Du, “Thermodynamic re-assessment of the Hf-Si binary system,”
  *J. Phys.: Conf. Ser.* 1074 (2018) 012074,
  https://doi.org/10.1088/1742-6596/1074/1/012074 — Creative Commons
  Attribution 3.0, https://creativecommons.org/licenses/by/3.0/ (IOP Publishing).

The derived files beside them (paper manifests, hints, declared axes, a Camelot
report and a topology computed by this repository's tools from the papers'
printed parameters) are this repository's own work under its MIT licence; they
quote table cells and caption text only as locators.

### MatCalc databases via pycalphad/binder — ODbL/DbCL (3 files)

Mirror: https://github.com/pycalphad/binder
Current publisher information:
https://www.matcalc.at/index.php/databases/open-databases

The headers of `mc_al_v2.032.pycalphad.tdb`,
`mc_fe_v2.059.pycalphad.tdb`, and `mc_ni_v2.034.pycalphad.tdb` explicitly apply
the Open Database License 1.0 and Database Contents License 1.0. The files
retain their copyright, editor, attribution, reference, and license notices.
See `databases/binder/SOURCE.md`.

### pycalphad/scheil — MIT (1 file)

Source: https://github.com/pycalphad/scheil
Repository license: https://github.com/pycalphad/scheil/blob/main/LICENSE

`sl_chen.tdb` is redistributed under the source repository's MIT license. The
NIMS-copyrighted `alzn_mey.tdb` is excluded and manual-only.

### openiec/openiec — MIT repository (2 bundled, 2 source-fetched demos)

Source: https://github.com/openiec/openiec
Repository license: https://github.com/openiec/openiec/blob/master/LICENSE

The four demo TDB files are distributed with the MIT-licensed OpenIEC project.
`AlAgCuWitusiewicz2005.TDB` and `NiAlCrHuang1999.tdb` remain bundled.
`AlNiAnsara1997.TDB` and `NiAlHuang1999.tdb` remain cataloged and fetchable from
upstream but are not bundled following the SGTE binary consolidation decision.

### Mendeley Data — CC BY 4.0 / CC0 1.0 (6 bundled, 2 source-fetched files)

Eight versioned datasets cover Bi-Te, Bi-Te-Yb, Mg-Sc-Zn, Ag-Pd-Sn, Fe-Mg,
Ti-M, Ti-Nb-Ta, and Y2O3-Ta2O5. Six files remain bundled. `BITE.TDB` and
`TDB_Fe-Mg.TDB` remain cataloged and source-fetched under CC BY 4.0 following
the SGTE binary consolidation decision. Exact DOIs, upstream filenames,
historical hashes, authorship notes, and renaming details are in
`databases/mendeley/SOURCE.md`.

### Zenodo 20848011 — CC BY 4.0 (1 file)

Pranav P. S. Setty, “Thermodynamic database for the Mo–Si–Y–Hf quaternary
system,” https://doi.org/10.5281/zenodo.20848011. Zenodo records the dataset as
CC BY 4.0. The unchanged file retains detailed author and subsystem citations;
see `databases/zenodo/SOURCE.md`.

### pycalphad Fe-Co-Cr-Nb-Ti fixture — ODbL 1.0 (1 file)

Source: https://github.com/pycalphad/pycalphad

The header of `mc_fecocrnbti.tdb` states that the reduced test database is based
on MatCalc steel database 2.060 and used under ODbL. See
`databases/pycalphad/SOURCE.md`.

### Materials Project computed thermodynamics — CC BY 4.0 (derived data)

Source: https://next-gen.materialsproject.org/ (API: https://api.materialsproject.org)

The Materials Project terms of use place its data under the Creative Commons
Attribution 4.0 International license:
https://creativecommons.org/licenses/by/4.0/
License evidence: https://legacy.materialsproject.org/terms

Cite the Materials Project when using anything derived from it:

- A. Jain, S. P. Ong, G. Hautier, W. Chen, W. D. Richards, S. Dacek, S. Cholia,
  D. Gunter, D. Skinner, G. Ceder, K. A. Persson, "Commentary: The Materials
  Project: A materials genome approach to accelerating materials innovation,"
  *APL Materials* 1 (2013) 011002. https://doi.org/10.1063/1.4812323
- Plus the release identifier recorded in each artifact's `provenance.mp_release`
  (obtained from `MPRester.get_database_version()`), because Materials Project
  formation energies change between database releases.

**Changes made, as CC BY 4.0 requires stating.** No Materials Project file is
redistributed verbatim. `scripts/fetch_mp_hull_snapshot.py` retrieves selected
fields of the `thermo` and `summary` collections and normalizes them into a
snapshot; `scripts/build_binary_hulls.py` then derives per-system convex-hull
documents from that snapshot. The derived artifacts are transformations of CC BY
4.0 data and remain under CC BY 4.0 with attribution to the Materials Project,
not under this repository's MIT license. Energies, corrections, and stability
values originate upstream and are not recomputed; only the hull construction,
interval representation, and cross-check reporting are ours.

Materials Project data is also mirrored on AWS Open Data
(`s3://materialsproject-build`, `s3://materialsproject-parsed`). That registry
cites the "Materials Project Terms of Use" rather than naming a license, so the
terms page above is the license evidence of record.

### SGTE-derived CALPHAD hull corpus — CC BY 4.0 (derived data, 2,598 files)

`hulls/calphad-sgte/` holds hull documents computed from the SGTE binary
collection above, plus the pinned samples they were built from
(`samples/samples.tar.gz`). They are a transformation of CC BY 4.0 data and
remain CC BY 4.0 — **not** covered by this repository's MIT license.

Attribution runs to two places, and both are required:

- Bengt Hallstedt, "The SGTE collection of binary datasets," *Calphad* 89 (2025)
  102833. https://doi.org/10.1016/j.calphad.2025.102833
- The original assessment named in the header of the source file each document
  was built from, recorded in every document's `provenance.assessment` alongside
  `provenance.source_path` and the file's `provenance.source_sha256`.

Statement of changes, as CC BY 4.0 requires: equilibria were computed at
298.15 K and 101 325 Pa with pycalphad, converted to a formation basis against
the stable pure-element states of the same database, and fitted as piecewise
cubics. No source file was modified. The transformation, its tolerances, and its
known limits are documented in [`CALPHAD_HULLS.md`](CALPHAD_HULLS.md); the
per-document fit residual and the live-gate report travel with the data.

### Open Quantum Materials Database (OQMD) — CC BY 4.0 (derived data)

Source: https://oqmd.org/ (bulk download: https://oqmd.org/download/)

OQMD states that its data is available under the Creative Commons Attribution
4.0 International license: https://creativecommons.org/licenses/by/4.0/
License evidence: https://oqmd.org/download/

Cite OQMD when using anything derived from it:

- J. E. Saal, S. Kirklin, M. Aykol, B. Meredig, C. Wolverton, "Materials Design
  and Discovery with High-Throughput Density Functional Theory: The Open
  Quantum Materials Database (OQMD)," *JOM* 65 (2013) 1501-1509.
  https://doi.org/10.1007/s11837-013-0755-4
- S. Kirklin, J. E. Saal, B. Meredig, A. Thompson, J. W. Doak, M. Aykol,
  S. Rühl, C. Wolverton, "The Open Quantum Materials Database (OQMD): assessing
  the accuracy of DFT formation energies," *npj Computational Materials* 1
  (2015) 15010. https://doi.org/10.1038/npjcompumats.2015.10
- Plus the release identifier recorded in each artifact's
  `provenance.upstream_release`, because formation energies change between
  releases and between fits.

**Changes made, as CC BY 4.0 requires stating.** No OQMD file is redistributed
verbatim. A local copy of the `qmdb` database is queried and normalized into a
snapshot, from which convex-hull documents are derived. Formation energies are
OQMD's own stored `formation_energies.delta_e` under a single pinned `fit_id`
and are **not** recomputed; the hull construction, interval representation, and
reporting are ours. The elemental endpoints are emitted explicitly at zero
formation energy, because OQMD's standard fit adjusts elemental reference
energies to match experimental formation enthalpies and therefore no calculated
entry sits exactly at the reference.

Derived artifacts are a transformation of CC BY 4.0 data and remain CC BY 4.0
with attribution to OQMD, not under this repository's MIT license.

OQMD energies are **not comparable** with Materials Project energies: the two
apply different corrections to different quantities. Hull documents from the two
sources therefore carry different `reference.kind` values
(`dft_fitted_elemental_reference_0K` versus `dft_elemental_ground_state_0K`) and
live under separate layers, and the query module refuses to compare them without
an explicit declared conversion.

## Fetched from source, not redistributed here

`scripts/fetch_sources.py` downloads these files directly on the user's
machine. Fetching does not grant redistribution rights.

### pycalphad/binder — mixed/uncertain (9 files)

Source: https://github.com/pycalphad/binder

The current repository has no license; its former MIT license was removed in a
[2020 repository reset](https://github.com/pycalphad/binder/commit/0b59abbeb190223bd2f4fc7912acb3c0af440cde).
The `COST507.pycalphad.tdb` header permits free use but does not explicitly
grant redistribution. Eight other mirrored files have no current file-level
license. All nine are therefore fetch-only pending clearer permission.

### sundmanbo/opencalphad — mixed evidence (23 files)

Source: https://github.com/sundmanbo/opencalphad

The repository contains the GPL-3.0 text, but 22 example TDB files do not state
that the repository license applies to them. `iron4cd.TDB` separately declares
CC BY 4.0 in its header. All remain fetch-only to avoid implying that the
software repository's GPL necessarily controls third-party database content.

### richardotis/pycalphad-sandbox — mixed evidence (23 files)

Source: https://github.com/richardotis/pycalphad-sandbox

The repository has carried an
[MIT license](https://github.com/richardotis/pycalphad-sandbox/blob/master/LICENSE)
since its initial 2016 commit, so it
is not accurately described as an unlicensed repository. File-level provenance
is nevertheless unclear for 22 TDBs. `steel_database_fix.tdb` explicitly states
ODbL 1.0 / DbCL 1.0. All 23 remain fetch-only while the third-party provenance
is reviewed.

### pycalphad/pycalphad — uncertain test fixtures (24 files)

Source: https://github.com/pycalphad/pycalphad

These TDBs appear in pycalphad's MIT-licensed test suite, but their individual
headers and assessment provenance do not establish that the repository license
controls all database content. They remain fetch-only pending file-level review.
The hash-bound [source discovery review](catalog/discovery/source-review.md)
also evaluates newly observed pycalphad paths without treating repository
placement as permission or creating duplicate catalog records.

### Corrected COST507 gist — uncertain (1 file)

Source: https://gist.github.com/bocklund/c4714ddbc0500c78e6fe255a763e7550

The database header permits free use but does not explicitly address
redistribution. Brandon Bocklund's corrected variant is therefore fetch-only.

## Manual sources

These sources require a landing page, DOI, request, registration, or individual
rights review and are not automatically downloaded.

| Source | License evidence | Treatment |
|---|---|---|
| COST 531 lead-free solders | Project described as open; redistribution terms not explicit | Verify with provider |
| COST MP0602 / IPM Brno solders | Available on request; redistribution terms not explicit | Verify with provider |
| [NIST Materials Data Repository](https://materialsdata.nist.gov/page/tos) (12 catalog records) | Records are authored by third parties and display no provider-selected license | Manual; verify each record |
| [MatCalc native `mc_fe`, `mc_ni`, `mc_al`](https://www.matcalc.at/index.php/databases/open-databases) | ODbL 1.0 / DbCL 1.0 | Manual download; redistribution permitted subject to license terms |
| [Mendeley `mc_ni_L`](https://data.mendeley.com/datasets/mnwnpfdw4n/1) | Description says ODbL/DbCL; platform metadata says GPL-3.0 | Manual; resolve conflicting metadata |
| [NIMS CPDDB](https://cpddb.nims.go.jp/) | No reproduction, republication, or third-party distribution without written NIMS permission | Restricted/manual |
| NIMS-origin files mirrored by binder, scheil, pycalphad, and pycalphad-sandbox (14 records) | File headers state Copyright (C) NIMS | Restricted/manual; not fetched automatically |
| Mendeley CC BY-NC datasets (5 records) | CC BY-NC 3.0 clearly stated | Manual pending a repository policy decision on noncommercial licenses |
| TDBDB, Brown University | Index only | Check each linked database separately |
| [Yao Wang doctoral dissertation, University of Genoa](https://iris.unige.it/retrieve/e268c4c9-70ae-a6b7-e053-3a05fe0adea1/phdunige_4074852.pdf) | Author retains copyright; absent a separate license, reproduction is limited to study and research with citation | Manual local extraction; generated TDB is git-ignored and must not be redistributed |

Authoritative landing pages and per-record URLs are recorded in
`catalog/catalog.json` and shown by `scripts/fetch_sources.py --manual`.
