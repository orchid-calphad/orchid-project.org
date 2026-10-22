# SGTE binary dataset collection — source and license

This directory contains the 2,773 binary TDB datasets and the enriched SGTE
unary v5.0 file published as supplementary material to:

> Bengt Hallstedt, “The SGTE collection of binary datasets,” *Calphad* 89
> (2025) 102833. https://doi.org/10.1016/j.calphad.2025.102833

The publication and its supplementary release are open access under the
[Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).
The files are redistributed here unchanged. Each binary file's header cites
the original assessment from which its parameters were taken; users should
retain those headers and cite the relevant original assessment in scientific
work.

`MANIFEST.sha256` records the SHA-256 of every redistributed TDB file.
`scripts/import_sgte.py` verifies the two original archive hashes and all
members before replacing repository files, then regenerates that manifest.

Upstream supplementary files:

- `1-s2.0-S0364591625000367-mmc3.zip` — 2,773 binary datasets; SHA-256
  `6a0d53de8a1141fc77ae2a1636794bfe90d6c8bb74b144c70eae960ff976747d`
- `1-s2.0-S0364591625000367-mmc4.zip` — `unary50_bengt.tdb`; SHA-256
  `ebf30a5b8e4677a01a7e2780513e490d9166edbb1e2c72757b6690df88ea0adc`

The paper reports 2,773 datasets covering 1,238 distinct binary systems. The
files were standardized and tested with Thermo-Calc GES6; some may require
modification for other CALPHAD software.
