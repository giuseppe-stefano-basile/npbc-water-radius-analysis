# Reproducibility contract

## Exact production stages

The `npbc-water run` command executes four scripts in a fixed order.  The CLI
creates only symbolic links in a compatibility workspace; it does not copy or
rewrite the raw trajectories.

1. `analyze_radius_equivalence.py`
2. `compute_rdf_extensions.py`
3. `compute_radial_layer_analysis.py`
4. `make_stochkin_figure_suite.py`

The first stage must complete before the figure stage because the latter reads
its regional, RDF, geometry, and global tables.  The second and third stages
add the matched PBC and radial-layer tables used by Figures 17--28.

## What preflight guarantees

- all four system labels are present exactly once;
- radii, molecule counts, offsets, windows, and expected frames match the
  final validated protocol;
- every trajectory, frozen bias, layer-definition file, and experimental
  table exists and is nonempty;
- the selected Python executable is available.

The scientific kernels perform stricter frame-grid, topology, and finite-value
checks during execution.

## Expected outputs

The final result directory contains:

```text
data/                 numerical tables and metadata
figures/              28 PNG/PDF pairs
source_inputs/        small bias/experimental provenance files
figure_manifest.csv   caption and filename mapping
suite_metadata.json   conventions and source information
pipeline_run.json     exact commands and completion timestamp
```

## Verification levels

1. Unit tests: formulas and limiting cases.
2. Schema checks: required tables, rows, labels, and normalization.
3. Numerical checks: finite metrics and normalized histograms.
4. Figure checks: exactly 28 sequential vector figures and manifest rows.
5. Input checks: SHA-256 inventories written by the kernels.

Never validate a rerun solely by visual similarity of plots.
