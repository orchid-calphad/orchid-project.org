# Contributing to orchid-project.org

This repository holds the ORCHID website (orchid-project.org), the database catalog and the release tools that
build the data pages and the data release (`tools/`).

- **A wrong value**: the **report** link next to it on the data pages, or the form
  [Report an error in a value](https://github.com/orchid-calphad/orchid-data/issues/new?template=report-value-error.yml).
  Papers, measurements, calculations and TDB files are contributed to
  [orchid-data](https://github.com/orchid-calphad/orchid-data/blob/main/CONTRIBUTING.md). The plain-English guide is
  [orchid-project.org/contribute/](https://orchid-project.org/contribute/).
- **A broken page, viewer or link**: open an issue with the bug form.
- **A TDB file for the catalog**: open an issue with its source, licence and the assessment DOI (or contribute it as
  type E in orchid-data).
- **Code** (`tools/`, pages): fork, branch, pull request. The tools are standard-library Python unless a file says
  otherwise, deterministic (same inputs, same bytes), and every builder has a `--check` mode; run the builders you
  touched with `--check` and the tests (`python3 -m pytest tools/tests`, which needs a local orchid-data checkout next
  to this repository). Generated pages are rebuilt by the maintainers; do not edit them by hand.

## Licence

Site code, tools and catalog: MIT (`LICENSE`). Data: CC BY 4.0 (see orchid-data). Third-party database files keep
their own licences (`data/NOTICE.md`, `data/LICENSES.md`). Sign off every commit (`git commit -s`): the `Signed-off-by` line certifies the
[Developer Certificate of Origin](https://developercertificate.org/), that you wrote the change or have the right to
submit it under the licence of this repository. CI checks the sign-off on every pull request.

## Maintainers

R. Arróyave (Texas A&M University, lead), J. Saal (Citrine Informatics) and D. Shin (Oak Ridge National Laboratory)
maintain all ORCHID repositories; how decisions are made is in
[GOVERNANCE.md](https://github.com/orchid-calphad/orchid-data/blob/main/GOVERNANCE.md). Everyone follows the
[Code of Conduct](CODE_OF_CONDUCT.md) (Contributor Covenant 2.1; reports to conduct@orchid-project.org).
