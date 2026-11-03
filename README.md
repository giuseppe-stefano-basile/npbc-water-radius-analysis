# NPBC water radius analysis

This repository reproduces the structural comparison of frozen-bias spherical
water cavities with radii 15, 20, and 30 Å against a matched 40 Å periodic
reference.  It is organized as a teaching repository: the exact validated
analysis kernels are preserved, while a typed configuration, preflight layer,
tests, reference outputs, and guided exercises make every estimator auditable.

## What the workflow computes

- molecular COM density in exact spherical shells;
- dipole, bond, and molecular-plane orientational order;
- O--H lengths and H--O--H angles;
- first-through-fifth O--O neighbors, direct coordination, tetrahedrality,
  and donated/accepted hydrogen bonds;
- per-shell temporal means and first temporal standard deviations;
- global COM displacement and angular-isotropy tensors;
- molecule-weighted and shell-volume-weighted distributions;
- the frozen Gaussian mean-field potential and analytic force;
- boundary-safe total and intermolecular O--O/O--H/H--H RDFs;
- exact spherical accessible-volume-corrected radial O--O RDFs;
- shell-resolved histogram comparisons against PBC;
- Jensen--Shannon, Wasserstein, standardized-effect, bootstrap, and practical-
  equivalence diagnostics;
- all 28 publication figures and their source tables.

The complete mathematical specification is in
[`paper/SI_R15_R20_R30_analysis_methods.pdf`](paper/SI_R15_R20_R30_analysis_methods.pdf).

## Repository layout

```text
analysis_scripts/     provenance-locked production kernels
configs/              explicit TOML configuration template
docs/                 guided theory, workflow, and student exercises
examples/             small parser/geometry demonstrations
paper/                SI LaTeX source and compiled document
reference_results/    final numerical tables and 28 PNG/PDF figure pairs
src/                   typed configuration, geometry, statistics, and CLI
tests/                 unit and reproducibility tests
```

Raw trajectories are deliberately not committed: they are large, cluster-
specific scientific records.  The final processed tables, figures, bias
profiles, experimental O--O table, metadata, and checks are included.  A local
configuration points the workflow to the raw dump files without copying them.

## Installation

Using `venv`:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Or with Conda/Mamba:

```bash
mamba env create -f environment.yml
mamba activate npbc-water-analysis
```

## First reproducibility checks

Run the unit tests and validate the packaged benchmark:

```bash
pytest
npbc-water check-reference reference_results
```

Regenerate all figures from the packaged tables, without raw trajectories:

```bash
npbc-water rebuild-figures reference_results
```

This is the recommended first exercise because it separates plotting from the
costly trajectory analysis.

## Running from raw trajectories

1. Copy the template and replace every `/path/to/...` entry:

   ```bash
   cp configs/r15_r20_r30.example.toml configs/r15_r20_r30.local.toml
   ```

2. Validate every file and every locked protocol value:

   ```bash
   npbc-water preflight configs/r15_r20_r30.local.toml
   ```

3. Inspect the exact commands without executing them:

   ```bash
   npbc-water run configs/r15_r20_r30.local.toml --dry-run
   ```

4. Run the full analysis:

   ```bash
   npbc-water run configs/r15_r20_r30.local.toml
   ```

The stage order is:

```text
whole-radius -> intermolecular RDF -> radial-layer analysis -> figures
```

Individual stages can be requested explicitly:

```bash
npbc-water run configs/r15_r20_r30.local.toml --stages whole rdf
npbc-water run configs/r15_r20_r30.local.toml --stages radial figures
```

## Why the protocol is locked

The final analysis uses two deliberate sampling conventions:

- Whole-radius profiles use 201 frames and include the endpoint.
- Radial-layer RDFs and distributions use 200 NPBC frames in half-open
  intervals, giving ten exactly equal 100 ps blocks.  PBC uses 160 frames and
  eight equal blocks.

Changing one convention silently would make the regenerated uncertainty bands
different from the archived figures.  The TOML file therefore records all
values, and the preflight refuses any mismatch.  A methodological extension
should be implemented on a branch with new tests and new output metadata.

## Scientific integrity rules

1. Never infer a missing frame count or molecule count.
2. Never replace a half-open interval with an inclusive interval implicitly.
3. Distinguish molecule weighting from shell-volume weighting.
4. Distinguish total site RDFs from intermolecular RDFs.
5. Near the wall, normalize by the exact accessible pair-shell volume.
6. Use `(N-1)/V` for distinct-particle O--O normalization.
7. Report temporal SD as temporal SD—not SEM or a confidence interval.
8. Treat the practical-equivalence screen as a multi-criterion diagnostic,
   not a mathematical proof of identical ensembles.

## Recommended learning path

Read the files in `docs/` in numerical order, reproduce Figures 1, 14, 22,
and 28, then attempt the exercises in `docs/05_student_project.md`.

## License and attribution

Repository code is MIT licensed.  The graphics convention is adapted from the
STOCHKIN project at the pinned commit recorded in `reference_results/`; its
license and attribution are preserved in `reference_results/style_reference/`.
Experimental data retain their original citation metadata.
