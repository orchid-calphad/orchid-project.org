# Fitting-input preflight

Run the preflight in an environment containing ESPEI and pycalphad, before starting a new fit:

```sh
python espei/check_fit.py --bundle espei --selection selection.json > preflight.json
```

Exit status 0 means that the selected inputs passed these checks and have no unresolved review flags. Exit status 1 blocks the fit. The checker reads the files without modifying values, weights, phase names or model parameters. A passing receipt does not establish model accuracy or guarantee that the optimizer will converge. Retain the receipt with the fitting run and use the exact files whose hashes it records. Changes made after this check, including `apply_tags` and residual-level multipliers, require another check of the actual fitting inputs or separate recorded run-level controls; this tool cannot intercept a running optimizer.

The selection is a JSON object with `target` and `datasets`. `target` requires `path` and `sha256` for the actual TDB. Each dataset requires `path`, relative to the bundle, and `sha256` for the released file. Paths for the target and optional `fit_file` are relative to the selection file, or absolute. Hashes are SHA-256 of the file bytes. Select files from `index.json`; the entire release is not a single fitting selection.

```json
{
  "target": {"path": "target.tdb", "sha256": "REPLACE_WITH_ACTUAL_SHA256"},
  "datasets": [
    {"path": "experimental/SYSTEM/FILE.json", "sha256": "REPLACE_WITH_ACTUAL_SHA256"}
  ]
}
```

For a weight override, provide a separate `fit_file` whose only change from the released dataset is `weight`, and add `weight_override` with `value` equal to that file's weight and a nonempty `reason`. The released weight is still checked against its effective uncertainty. The receipt records the old input hash, actual fitting-file hash and override. Do not replace the active weight with the superseded `suggested_weight` automatically.

The checker matches components and phases to the actual TDB, checks manual-solver sublattice ratios and species, and loads every selected file through ESPEI. Phase identities are compared using `phases/systems` beside the bundle. For a new model whose digest is not in that registry, `target.phase_identity` can explicitly map a phase name to an object containing `registry_id` and `reason`. Such assertions remain visible in the receipt; they are not inferred from the filename. Registered PhaseForge and PhaseForge+ provenance, including the potential, is attached only when its TDB digest matches.

Repeated observation identities, including known aliases and DFT mapping identities, block the selection. Shared derived observations, repeated numerical constraints and multiple encodings of one assessed diagram require review. A reviewed overlap may be recorded in `acknowledgements` as an object with the exact reported `code`, `paths` and a nonempty `reason`. This records a deliberate fitting choice; it does not establish statistical independence. Different cited sources with equal numbers are not automatically deleted or merged.

ZPF temperature and composition budgets remain distinct from the driving-force tolerance that ESPEI uses to weight these constraints. Shared-source covariance is not estimated. The historical Cu–Zr workflow remains pinned to its documented data and weight overrides; this preflight does not alter that workflow or rerun its fit.
