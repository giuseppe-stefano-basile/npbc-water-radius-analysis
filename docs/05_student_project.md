# Suggested master-student project

## Core reproduction

1. Reproduce the packaged figures from the archived tables.
2. Run the complete trajectory workflow on the original dumps.
3. Compare every regenerated input inventory and key output table with the
   reference result.
4. Explain Figures 1, 14, 22, and 28 from their estimators.

## Required written discussion

- Why does shell density require exact spherical volumes?
- Why is radial orientation undefined as a PBC comparison observable?
- Why do total O--H/H--H RDFs contain intramolecular peaks?
- Why is `(N-1)/V` the consistent density for distinct O--O pairs?
- Why may corrected RDFs look bulk-like while direct coordination decreases?
- Why is temporal SD not an uncertainty on the mean?

## Controlled extensions

Choose one extension and create it on a new branch:

1. Compare 0.5, 1.0, and 2.0 Å wall-layer widths while preserving exact
   accessible-volume correction.
2. Quantify sensitivity to RDF bin widths of 0.025, 0.05, and 0.10 Å.
3. Compute block-length sensitivity for 50, 100, and 200 ps blocks.
4. Add orientational distributions relative to the local wall normal for NPBC
   only, without inventing a PBC radial axis.
5. Add a finite-size ideal-gas Monte Carlo validation of the radial RDF
   denominator.

For every extension, report the changed estimator, computational cost,
expected bias/variance tradeoff, and whether the scientific conclusion changes.

## Assessment rubric

| Component | Weight |
|---|---:|
| Exact reproduction and provenance | 25% |
| Correct derivation and implementation | 30% |
| Tests and numerical validation | 20% |
| Physical interpretation | 15% |
| Code clarity and documentation | 10% |
