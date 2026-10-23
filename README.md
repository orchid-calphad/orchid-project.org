# orchid-project.org

Website and released data of **ORCHID** (Open Repository for CALPHAD Hypotheses, Intelligence, and Discovery),
served by GitHub Pages at https://orchid-project.org. Releases are archived on Zenodo (one DOI per version).

- `index.html`, `catalog/`, `best-practices.html` — static, self-contained pages (generated in the ORCHID tooling).
- `data/` — released data: catalog, redistributable TDB files, convex hulls. Licences: `data/LICENSES.md`, `data/NOTICE.md`.
- `llms.txt`, `data/croissant.json`, `CITATION.cff`, `.zenodo.json` — machine-readable metadata.

Only data with clear redistribution rights are published. Version 1, machine-extracted: values were read from papers by language-model readers constrained to quote their source and checked automatically (quote verification, field-placement rule, printed-page sign checks, DFT and Kim et al. 2017 cross-checks); a hand-checked precision sample is in progress.

## Zenodo release

Release v2026.10 is archived as five citable datasets plus an umbrella record (metadata and file manifests in the
`zenodo/` folder of `orchid-calphad/orchid-data`). `tools/zenodo_release.py plan | bundle | writeback` checks the
records, builds one deterministic zip per record into the ignored `dist/`, and writes the reserved DOIs back into
both repositories; every builder reads the DOIs from `zenodo/dois.json` through `tools/orchid_cite.py`.
