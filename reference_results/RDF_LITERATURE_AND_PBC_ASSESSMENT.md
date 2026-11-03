# RDF literature and matched-PBC assessment

## Conclusion

The boundary-safe R15, R20, and R30 radial distribution functions reproduce
the matched 40 Å periodic MACE-OFF23 reference. The residual differences are
small, decrease with cavity radius, and are compatible with the larger block
fluctuations expected from the smaller number of boundary-safe reference atoms
in R15. This supports the conclusion that the observed liquid structure is not
created by the spherical NPBC boundary.

Agreement with ambient-water experiments is feature-dependent. The first O--O
peak and intermolecular H--H maximum are close to the published benchmarks.
The O--O first minimum is shifted outward and the broad second O--O maximum is
flattened/shifted inward. Intermolecular O--H peaks are also higher than the
tabulated ambient benchmark. Because the matched PBC simulation shows the same
features, these deviations are state/model effects (including the simulated
number density) rather than radius-dependent NPBC artifacts.

The complete Skinner experimental `g_OO(r)` curve confirms this interpretation,
rather than changing it. The first peak is reproduced closely, but every
simulation retains too much probability between roughly 3.1 and 3.8 Å,
underestimates the experimental second-shell maximum near 4.5 Å, and shows
smaller phase/amplitude differences through the subsequent oscillations. The
same pattern is present in PBC40.

## Matched NPBC--PBC differences

RDF RMSE over the structural range:

| Pair | R15 | R20 | R30 |
|---|---:|---:|---:|
| OO, 2.3--8 Å | 0.0458 | 0.0248 | 0.0179 |
| OH intermolecular, 1.3--8 Å | 0.0302 | 0.0160 | 0.0120 |
| HH intermolecular, 1.8--8 Å | 0.0265 | 0.0138 | 0.0106 |

The 7--8 Å tails are close to unity for every pair and radius. The direct
curves and residuals are plotted in Figure 20; complete RMSE, MAE, maximum
deviation, and tail statistics are in
`data/npbc_pbc_rdf_comparisons.csv`.

## Full experimental O--O curve

Figure 21 overlays the complete numerical Skinner et al. curve over the full
common interval, 1.1--7.975 Å. The source has 2,000 points at 0.01 Å spacing
from 0 to 19.99 Å and includes a pointwise estimated error. No curve was
digitized from an image.

Descriptive errors over the structurally informative 2.2--7.975 Å interval:

| System | RMSE | MAE | Mean signed difference | Maximum absolute difference |
|---|---:|---:|---:|---:|
| R15 | 0.1312 | 0.0978 | +0.0268 | 0.3873 |
| R20 | 0.1194 | 0.0886 | +0.0161 | 0.3581 |
| R30 | 0.1180 | 0.0897 | +0.0157 | 0.3480 |
| PBC40 | 0.1088 | 0.0840 | +0.0056 | 0.3302 |

PBC40 is closest to experiment, and the NPBC error decreases from R15 to R30.
More importantly, all four curves share the same structured residual. The
remaining disagreement is therefore not an NPBC-radius effect. The pointwise
experimental error is substantially narrower than the model residual over
most of the interval, but the corresponding error-normalized values are not a
formal chi-squared test: the 22 °C experiment and the 300 K simulations are not
matched in number density, and the simulation uncertainty shown is a temporal
fluctuation statistic rather than experimental measurement uncertainty.

The raw curve, aligned per-bin values, residuals, and metrics are in
`data/experimental_gOO_skinner_2013_full.csv`,
`data/gOO_experiment_overlay_aligned.csv`, and
`data/gOO_experiment_comparison_metrics.csv`.

## Experimental landmarks

| Feature | Experiment | R15 | R20 | R30 | PBC40 |
|---|---:|---:|---:|---:|---:|
| OO first max: r (Å) | 2.80 | 2.775 | 2.775 | 2.775 | 2.775 |
| OO first max: g | 2.57 | 2.672 | 2.696 | 2.623 | 2.613 |
| OO first min: r (Å) | 3.45 | 3.925 | 3.875 | 3.875 | 3.775 |
| OO first min: g | 0.84 | 0.884 | 0.911 | 0.946 | 0.953 |
| OO second max: r (Å) | 4.50 | 4.425 | 4.375 | 4.175 | 4.325 |
| OO second max: g | 1.12 | 1.047 | 1.032 | 1.001 | 0.993 |
| OH H-bond max: r (Å) | 1.73 | 1.825 | 1.825 | 1.825 | 1.825 |
| OH H-bond max: g | 1.01 | 1.155 | 1.187 | 1.188 | 1.182 |
| HH intermolecular max: r (Å) | 2.38 | 2.375 | 2.375 | 2.375 | 2.325 |
| HH intermolecular max: g | 1.29 | 1.371 | 1.345 | 1.325 | 1.320 |

Peak/minimum positions are resolved on a 0.05 Å RDF grid, so comparisons finer
than one bin are not meaningful. Experimental symbols are published tabulated
landmarks, not values digitized from plotted curves.

## Definitions and uncertainty

- RDF bin width: 0.05 Å; maximum distance: 8 Å.
- NPBC: reference atoms satisfy `r <= R - 8 Å`.
- PBC: minimum-image distances in a 40.000 Å cubic cell.
- O--H and H--H: atom pairs sharing a molecule ID are excluded.
- Sampling: every 25 ps.
- Shading: first temporal standard deviation across 100 ps RDF blocks, not SEM
  or a confidence interval.
- Experimental O--O band: author-supplied pointwise estimated error, not the
  simulation temporal SD.
- PBC state: 2383 waters, 64,000 Å³, 0.037234375 molecule Å^-3, 300 K NVT,
  stationary 200--1000 ps interval.
- NPBC target density: 0.037242 molecule Å^-3.

## Literature sources

1. L. B. Skinner et al., *J. Chem. Phys.* **138**, 074506 (2013), DOI
   `10.1063/1.4790861` — complete numerical benchmark ambient-water O--O pair
   distribution and pointwise estimated error used in Figure 21.
2. A. K. Soper and C. J. Benmore, *Phys. Rev. Lett.* **101**, 065502 (2008),
   DOI `10.1103/PhysRevLett.101.065502` — joint X-ray/neutron EPSR water
   structure used for the O--H and H--H landmarks.
3. L. Zheng et al., *J. Chem. Phys.* **148**, 164505 (2018), DOI
   `10.1063/1.5023611` — Tables I--III compiling the RDF landmark values used
   here.
