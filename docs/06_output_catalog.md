# Output catalog

## Whole-radius and ensemble tables

- `temporal_shell_timeseries.csv.gz`: frame-by-shell observables.
- `temporal_shell_statistics.csv`: shell temporal mean, SD, SEM, and valid
  fraction.
- `shell_profiles.csv`: complete-radius profile summary.
- `regional_distributions.csv`: molecule- and volume-weighted histograms.
- `regional_summary.csv`: regional means and block uncertainty.
- `global_timeseries.csv`: COM shift and anisotropy by frame.
- `global_temporal_summary.csv`: global temporal means and SDs.
- `equivalence_tests.csv`: all practical-equivalence criteria.
- `shell_profile_comparisons.csv`: matched wall-distance profile metrics.
- `system_geometry.csv`: radii, water counts, densities, and outer fractions.

## Mean-field tables

- `potential_amplitudes_full_radius.csv`: physical/control Gaussian layers.
- `potential_fields_full_radius.csv`: exact `U(r)` and analytic `-dU/dr`.
- `potential_extrema.csv`: field ranges by system.
- `potential_input_inventory.csv`: frozen input checksums.

## RDF tables

- `boundary_safe_rdf.csv`: total O--O/O--H/H--H functions.
- `rdf_peak_metrics.csv`: total-RDF feature locations and heights.
- `intermolecular_rdf_with_block_std.csv`: intermolecular functions with block
  SD and SEM.
- `intermolecular_rdf_feature_metrics.csv`: intermolecular RDF landmarks.
- `npbc_pbc_rdf_comparisons.csv`: NPBC--PBC RMSE/MAE/max/tail comparisons.
- `pbc40_total_rdf.csv`: periodic total site RDFs.
- `rdf_extension_input_inventory.csv`: extension input checksums.
- `rdf_extension_metadata.json`: RDF conventions.

## Experimental comparison

- `experimental_gOO_skinner_2013_full.csv`: numerical experimental O--O curve.
- `experimental_gOO_skinner_2013_metadata.json`: provenance and units.
- `experimental_rdf_landmarks.csv`: published O--O/O--H/H--H landmarks.
- `gOO_experiment_overlay_aligned.csv`: experiment interpolated to simulation
  bins with residuals.
- `gOO_experiment_comparison_metrics.csv`: full/common-range error metrics.

## Radial-layer tables

- `radial_layer_gOO.csv`: coarse and fine edge-corrected RDFs.
- `radial_layer_gOO_metrics.csv`: errors, direct coordination, and features.
- `radial_layer_property_distributions_vs_pbc.csv`: 12 intrinsic histograms.
- `radial_layer_distribution_metrics_vs_pbc.csv`: JS/Wasserstein/mean changes.
- `radial_layer_input_inventory.csv`: trajectory checksums and windows.
- `radial_layer_analysis_metadata.json`: exact layer/normalization conventions.

## Provenance

- `dump_inventory.csv`: dump selection and checksums.
- `figure_manifest.csv`: maps all 28 figures to PNG, PDF, and caption.
- `suite_metadata.json`: global plotting and estimator conventions.
