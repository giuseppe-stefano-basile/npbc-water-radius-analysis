# Provenance-locked analysis kernels

These four scripts are the validated kernels used to generate the archived
R15/R20/R30/PBC40 comparison:

1. `analyze_radius_equivalence.py` computes molecular observables, full-radius
   shell profiles, global diagnostics, regional histograms, practical-
   equivalence tests, and boundary-safe total RDFs.
2. `compute_rdf_extensions.py` computes intermolecular O--O, O--H, and H--H
   RDFs, matched PBC residuals, block fluctuations, and literature landmarks.
3. `compute_radial_layer_analysis.py` computes exact spherical edge-corrected
   O--O RDFs and shell-conditioned intrinsic-property distributions.
4. `make_stochkin_figure_suite.py` reconstructs frozen Gaussian fields and
   generates the 28 publication figures.

They are intentionally retained in the form used for the final analysis.  The
`npbc-water` command stages user-selected input files into their historical
layout using symbolic links and then executes these kernels without rewriting
scientific constants at runtime.  This makes the benchmark exactly
reproducible and makes protocol changes explicit code changes rather than
hidden configuration drift.
