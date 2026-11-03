# Learning path

The repository is designed for a master student who already knows basic
molecular dynamics and Python/NumPy but has not yet implemented finite-system
structural analysis.

## Phase 1: reproduce before modifying

1. Create the environment and run `pytest`.
2. Run `npbc-water check-reference reference_results`.
3. Rebuild all figures from the packaged CSV tables.
4. Compare regenerated file counts and `figure_manifest.csv` with the archive.
5. Read the SI methods while opening Figures 1, 14, 22, and 28.

The goal is to understand that a plot is the final layer of an estimator, not
the estimator itself.

## Phase 2: trace one observable end to end

Use normalized COM density as the first example:

```text
LAMMPS atoms -> water reconstruction -> molecular COM -> radial shell
-> per-frame count -> exact shell-volume normalization -> temporal mean/SD
-> CSV table -> Figure 1
```

Then repeat for tetrahedrality and for the edge-corrected O--O RDF.  At every
step identify the numerator, denominator, ensemble, time window, and unit.

## Phase 3: understand boundary corrections

Derive the sphere-intersection formula in `docs/02_estimators.md`.  Use the
tests in `tests/test_geometry.py` to check limiting cases.  Explain why an RDF
normalized by accessible volume and a direct coordination number answer
different physical questions near a wall.

## Phase 4: statistical interpretation

Distinguish:

- frame-to-frame temporal SD;
- SD across 100 ps blocks;
- SEM of block means;
- bootstrap confidence intervals;
- effect sizes and distribution distances.

Do not use these quantities interchangeably.  The plotted shaded bands are
explicitly identified in the figure captions and metadata.

## Phase 5: extend safely

Only after reproducing the reference outputs should the student change a
cutoff, shell width, RDF bin width, or time window.  Every extension should:

1. create a new configuration or branch;
2. add or update a unit test;
3. write the changed convention into metadata;
4. preserve the original reference outputs;
5. compare the new result against the baseline quantitatively.
