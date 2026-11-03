# R15/R20/R30 STOCHKIN figure and data suite

This package redraws the frozen-validation comparison using the publication
graphics convention from
[`giuseppe-stefano-basile/stochkin`](https://github.com/giuseppe-stefano-basile/stochkin),
pinned to commit `3c99adb84f7fe5efd2a2b5e4b8ec6c094db8d067`.

The convention uses an Arial/Helvetica sans-serif fallback, inward major and
minor ticks, right-side y ticks, white axes, 0.75 pt spines, 1 pt curves,
compact publication dimensions, 300 dpi raster output, and vector PDF output.
The exact referenced `style.py`, commit metadata, and license are included in
`style_reference/`.

## Temporal-fluctuation convention

Every shaded radial band is

\[
\overline{A}(r)\pm\sigma_t[A(r)],
\]

where `A(r)` is first averaged within a shell for each saved trajectory frame
and `sigma_t` is the ordinary first temporal standard deviation across the 201
frames. It is **not** a standard error or confidence interval. The trajectories
cover matched 1 ns production windows at 5 ps resolution.

The central 0.25 Å shells have very small geometrical volumes and intermittent
occupancy. Their raw temporal fluctuations are preserved in the data and in
Figure 12. Figure 01 restricts the density y-axis to 0.5–1.5 so the physically
resolved full-radius profile remains visible; central uncertainty bands are
therefore clipped by the panel limits. Structural shells occupied in fewer
than half the frames are marked with open gray circles.

## Frozen mean field

The plotted conservative molecular mean field is reconstructed exactly as

\[
U(r)=0.0005\,{\rm eV}\sum_k
\left[g_k-g_0\right]
\exp\!\left[-\frac12\left(\frac{r-r_k}{\sigma_k}\right)^2\right],
\]

with `F_r=-dU/dr` in the LAMMPS sign convention used by the audited source.
Figures 06–08 cover the complete physical interval `0 <= r <= R`; Figure 09
aligns the same fields by distance from the wall. The R15, R20, and R30 frozen
profiles correspond to global updates 159, 212, and 188, respectively.

## Figures

Every numbered figure is supplied as a 300 dpi PNG and vector PDF.

1. Full-radius COM density, temporal mean ±1σ.
2. Full-radius dipole P1 and P2.
3. Full-radius tetrahedrality, O–O coordination, and hydrogen bonding.
4. Full-radius O–H length and H–O–H angle.
5. Full-radius first and fourth O–O neighbor distances.
6. Frozen Gaussian amplitudes across the complete absolute radius.
7. Gaussian amplitudes aligned by normalized radius.
8. Exact conservative potential and force across the complete radius.
9. Potential and force aligned by wall distance.
10. Boundary density mean ±1σ overlaid with the exact potential.
11. Global COM displacement and angular-anisotropy time traces, with mean ±1σ.
12. First temporal standard deviation of density and network observables.
13. Whole-cavity shell-volume-weighted distributions.
14. Boundary-safe OO, OH, and HH pair correlations.
15. Cross-radius equivalence effect-size matrix.
16. Surface-to-volume fractions and conditioned equivalence pass rates.
17. Boundary-safe and periodic `g_OO(r)` against published experimental
    first-maximum, first-minimum, and second-maximum landmarks.
18. Intermolecular-only `g_OH(r)` (same-molecule O--H pairs excluded), with
    published experimental landmarks.
19. Intermolecular-only `g_HH(r)` (same-molecule H--H pairs excluded), with
    published experimental landmarks.
20. Direct matched-state R15/R20/R30 NPBC versus 40 Å PBC MACE-OFF23 RDFs,
    including NPBC--PBC residuals.
21. Full common-range `g_OO(r)` overlay of R15/R20/R30/PBC40 with the
    2,000-point Skinner et al. experimental curve, its pointwise estimated
    error, and simulation-minus-experiment residuals.
22. O--O RDFs from the outer 0--2 Å layer inward, with matched PBC residuals.
23. Fine 1 Å radial map of `g_layer(r)-g_PBC(r)`.
24. Radial decay of RDF RMSE, direct coordination, and peak/minimum shifts.
25. Tetrahedrality, O--O coordination, and hydrogen-bond histograms by radial
    layer, overlaid with PBC.
26. O--H length and H--O--H angle histograms by radial layer and PBC.
27. First through fifth O--O-neighbor histograms by radial layer and PBC.
28. All-property Jensen--Shannon distance from PBC by radial layer.

## RDF extension and periodic reference

Figures 17--20 use the same 0.05 Å bins and 8 Å cutoff as Figure 14. NPBC
reference atoms are at least 8 Å from the spherical wall. Intermolecular O--H
and H--H functions exclude same-molecule pairs explicitly by molecule ID.

The periodic control is a 40.000 Å cubic MACE-OFF23 `small_2` simulation with
2383 waters at 300 K. Its number density is 0.037234375 molecule Å^-3, compared
with the NPBC target of 0.037242 molecule Å^-3. RDFs use minimum-image
distances over the stationary 200--1000 ps interval and are sampled every
25 ps, identically to the NPBC RDF sampling interval. Shading in Figures
17--20 is the first temporal standard deviation across 100 ps block RDF
estimates; it is not a confidence interval.

The O--O benchmark is the complete author-supplied numerical curve from
Skinner et al., *J. Chem. Phys.* **138**, 074506 (2013), DOI
`10.1063/1.4790861`: 2,000 points on a 0.01 Å grid from 0 to 19.99 Å with a
pointwise estimated error. Figure 21 uses the entire common simulation--
experiment interval, 1.1--7.975 Å. Figure 17 retains the independently
tabulated O--O landmarks as a compact feature comparison. O--H and H--H remain
landmark-only because no full numerical partial-RDF table was supplied; their
landmarks are from the joint X-ray/neutron EPSR benchmark of Soper and Benmore,
*Phys. Rev. Lett.* **101**, 065502 (2008), as tabulated by Zheng et al.,
*J. Chem. Phys.* **148**, 164505 (2018), DOI `10.1063/1.5023611`.

The direct matched-state comparison shows that the three boundary-safe NPBC
RDFs reproduce the PBC MACE-OFF23 structure. Across the selected structural
ranges, NPBC--PBC RDF RMSE values are 0.0179--0.0458 (OO), 0.0120--0.0302
(OH), and 0.0106--0.0265 (HH), decreasing systematically from R15 to R30.
Consequently, differences from ambient experimental landmarks should be
interpreted primarily as state/model effects rather than spherical-boundary
effects.

The continuous experimental comparison reaches the same conclusion. Over
2.2--7.975 Å, the `g_OO` RMSE values against experiment are 0.1312 (R15),
0.1194 (R20), 0.1180 (R30), and 0.1088 (PBC40). The PBC reference is therefore
the closest of the four curves, while the NPBC residual decreases with cavity
radius. The remaining common deviation is much larger than the pointwise
experimental error and is shared by PBC40, so it cannot be attributed to the
spherical boundary. This is a descriptive cross-state comparison: the
experiment is at 22 °C and the simulation and experiment are not matched in
number density.

## Radially resolved interface analysis

Figures 22--28 condition each NPBC reference water on its molecular-COM
distance from the spherical wall. The O--O RDFs use 0.05 Å bins, the
distinct-particle finite-N density `(N_O-1)/V`, and the analytic intersection
volume between every O-centered pair-distance shell and the spherical cavity.
This removes the geometrical loss of available solid angle near the wall.
Fine RDF maps use 1 Å wall layers; distribution panels use 2 Å wall layers.

The analysis uses half-open intervals and equal 100 ps blocks: 200 frames for
each NPBC system and 160 frames for PBC40. The plotted fluctuation bands are
the first temporal SD across these blocks. The strongest restructuring is
confined to 0--2 Å, 2--4 Å is transitional, and beyond about 4 Å the
network, neighbor, and RDF distributions are close to PBC. The complete
method, checks, numerical effect sizes, and interpretation are in
`RADIAL_LAYER_RDF_AND_DISTRIBUTION_ASSESSMENT.md`.

`figure_manifest.csv` maps every image to its caption.

## Relevant data

- `data/temporal_shell_timeseries.csv.gz`: every frame × shell observation.
- `data/temporal_shell_statistics.csv`: temporal mean, SD, SEM, and valid-frame
  fraction for every shell and observable.
- `data/potential_amplitudes_full_radius.csv`: physical and control Gaussian
  coefficients.
- `data/potential_fields_full_radius.csv`: exact U(r) and F(r) on a 0.01 Å grid.
- `data/potential_extrema.csv`: potential/force ranges by radius.
- `data/global_timeseries.csv` and `data/global_temporal_summary.csv`: global
  COM and anisotropy fluctuations.
- `data/regional_distributions.csv`: molecule- and shell-volume-weighted
  probability distributions.
- `data/equivalence_tests.csv`: all pairwise practical-equivalence tests.
- `data/boundary_safe_rdf.csv`: OO/OH/HH RDF curves.
- `data/intermolecular_rdf_with_block_std.csv`: intermolecular OO/OH/HH RDFs
  for R15/R20/R30 and PBC40, with 100 ps block temporal SD and SEM.
- `data/intermolecular_rdf_feature_metrics.csv`: peak/minimum locations and
  heights for every system and pair.
- `data/npbc_pbc_rdf_comparisons.csv`: RMSE, MAE, maximum absolute deviation,
  and 7--8 Å tail mean against PBC40.
- `data/experimental_rdf_landmarks.csv`: traceable literature benchmark
  positions, heights, available uncertainties, and DOIs.
- `data/experimental_gOO_skinner_2013_full.csv`: complete author-supplied
  experimental O--O RDF and pointwise estimated error.
- `data/gOO_experiment_overlay_aligned.csv`: experimental curve interpolated
  to every simulation RDF bin center, together with simulation temporal SD,
  SEM, and residuals.
- `data/gOO_experiment_comparison_metrics.csv`: full-common-range and
  structural-range RMSE, MAE, bias, maximum deviation, and error-normalized
  descriptive diagnostics.
- `data/experimental_gOO_skinner_2013_metadata.json`: state point, DOI,
  source hash, grids, and comparison caveat.
- `data/radial_layer_gOO.csv`: coarse- and fine-layer edge-corrected O--O
  RDFs, counts, denominators, and block SD/SEM, plus the PBC curve.
- `data/radial_layer_gOO_metrics.csv`: layer-level RDF residual, coordination,
  peak, minimum, and temporal-fluctuation metrics against PBC.
- `data/radial_layer_property_distributions_vs_pbc.csv`: normalized
  shell-conditioned histograms for 12 intrinsic properties with block SD/SEM
  and pointwise PBC differences.
- `data/radial_layer_distribution_metrics_vs_pbc.csv`: Jensen--Shannon,
  Wasserstein, and mean-shift metrics, including block temporal SD.
- `data/radial_layer_input_inventory.csv` and
  `data/radial_layer_analysis_metadata.json`: exact trajectory hashes,
  selected intervals, definitions, and sampling conventions.
- `data/pbc40_total_rdf.csv`: total periodic site-site RDFs, retaining
  intramolecular O--H/H--H peaks for comparison with Figure 14 conventions.
- `data/rdf_extension_metadata.json` and
  `data/rdf_extension_input_inventory.csv`: RDF definitions, state point,
  model identity, sampling, source hashes, and uncertainty convention.
- `data/shell_profiles.csv`: original 0.25 Å profile statistics.
- `data/system_geometry.csv`: molecule counts, number densities, and boundary
  volume fractions.
- `source_inputs/`: the exact frozen bias and layer-definition files used for
  potential reconstruction, plus the PBC production log, final configuration,
  and stationary-window extraction script. The 247 MB raw PBC dump is omitted
  because all RDF bins, counts, block fluctuations, hashes, and state metadata
  needed for the figure suite are already packaged.

## Reproduction

Regenerate all 28 PNG/PDF figure pairs directly from the data packaged here:

```bash
python3 make_stochkin_figure_suite.py --figures-only
```

For a complete upstream rebuild from the raw dumps, run from the workspace
containing the three-radius dump analysis and bias inputs:

```bash
python3 make_stochkin_figure_suite.py \
  --workspace /path/to/workspace \
  --output /path/to/stochkin_radius_figures_20260921
```

The script rebuilds the temporal shell table directly from the dumps and uses
atomic, validated PNG/PDF export to prevent incomplete figure files.

To recompute the extended RDF tables from the source trajectories before
redrawing the figures:

```bash
python3 compute_rdf_extensions.py \
  --workspace /path/to/workspace \
  --pbc-dump /path/to/traj_pbc40_prod_stationary_200_1000ps.dump \
  --output /path/to/stochkin_radius_figures_20260921

python3 make_stochkin_figure_suite.py --figures-only
```

To recompute the radial-layer RDFs and histograms from the source trajectories:

```bash
python3 compute_radial_layer_analysis.py \
  --workspace /path/to/workspace \
  --pbc-dump /path/to/traj_pbc40_prod_stationary_200_1000ps.dump \
  --output /path/to/stochkin_radius_figures_20260921

python3 make_stochkin_figure_suite.py --figures-only
```
