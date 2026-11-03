# Radially resolved RDF and distribution assessment

## Scope

This extension resolves the validated R15, R20, and R30 NPBC trajectories by
the reference-water molecular-COM distance from the spherical wall and
compares every layer directly with the matched 40 Å periodic MACE-OFF23
reference. It is designed to distinguish a localized interfacial response
from a cavity-wide change of the liquid.

The analysis uses half-open production intervals and equal 100 ps temporal
blocks:

| System | Interval (ps) | Frames | Reference-water samples |
|---|---:|---:|---:|
| R15 | [100, 1100) | 200 | 105,200 |
| R20 | [0, 1000) | 200 | 249,600 |
| R30 | [100, 1100) | 200 | 842,200 |
| PBC40 | [200, 1000) | 160 | 381,280 |

All shaded bands are the first temporal standard deviation across the equal
100 ps blocks. They are not SEMs or confidence intervals.

## RDF definition and validation

The shell-conditioned RDF is an O--O RDF. The radial layer is assigned from
the reference water COM, while the pair coordinate is the O--O separation.
NPBC denominators use the analytic intersection of each O-centered
pair-distance shell with the spherical cavity. This removes the geometrical
loss of solid angle at the wall. Distinct-particle normalization uses
`(N_O-1)/V`, including the finite-N correction. The bin width is 0.05 Å and
the maximum pair distance is 8 Å.

The following internal checks passed:

- every expected frame and reference-water sample was recovered;
- coarse- and fine-layer RDFs contain identical total pair counts;
- all 192 property histograms normalize to unity;
- all 54 layer-level RDF metric rows are finite;
- equal, half-open 100 ps blocks are used throughout.

The small number of configurations in which an oxygen lies just outside the
nominal COM cavity is retained and reported in
`data/radial_layer_analysis_metadata.json`. The fractions are 0.277% (R15),
0.185% (R20), and 0.123% (R30) of reference samples and do not control the
conclusions.

## Main result: the boundary response is localized

The same spatial crossover is observed for every cavity. The 0--2 Å layer is
strongly restructured, the 2--4 Å layer is a transition region, and layers
beyond about 4 Å are close to the PBC distribution. The effect is local and
nearly radius-independent at fixed wall distance.

| System | Layer | RDF RMSE vs PBC, 2.2--6 Å | Direct coordination change, <=3.5 Å |
|---|---|---:|---:|
| R15 | 0--2 Å | 0.1571 | -2.174 |
| R15 | 2--4 Å | 0.0281 | -0.385 |
| R15 | >=8 Å | 0.0406 | +0.202 |
| R20 | 0--2 Å | 0.1537 | -2.080 |
| R20 | 2--4 Å | 0.0265 | -0.381 |
| R20 | >=8 Å | 0.0206 | +0.102 |
| R30 | 0--2 Å | 0.1469 | -2.018 |
| R30 | 2--4 Å | 0.0248 | -0.384 |
| R30 | >=8 Å | 0.0150 | +0.090 |

The fine 1 Å analysis resolves the crossover more sharply. RDF RMSE drops
from 0.184--0.192 in the outermost 0--1 Å to 0.120--0.131 in 1--2 Å,
0.046--0.050 in 2--3 Å, and 0.017--0.026 in 3--4 Å. At 3--4 Å the direct
coordination difference is already +0.034, +0.000, and -0.030 for R15, R20,
and R30, respectively.

R15 has visibly larger interior scatter because the same 1 Å layer contains
far fewer molecules near the center: its `>=12 Å` category has only 870
reference samples, compared with 16,113 for R20 and 183,639 for R30. The
larger R15 interior RMSE should therefore not be interpreted as a renewed
long-range wall effect.

## What changes in the outer layer

The direct molecular histograms independently confirm the RDF result. Relative
to PBC, the 0--2 Å layer has lower tetrahedrality and coordination, fewer
hydrogen bonds, and progressively displaced third-to-fifth neighbors. The
intramolecular O--H and H--O--H distributions change only slightly.

| Outer-layer mean minus PBC | R15 | R20 | R30 |
|---|---:|---:|---:|
| Tetrahedrality q | -0.256 | -0.232 | -0.214 |
| O--O coordination <=3.5 Å | -2.174 | -2.080 | -2.018 |
| Total H bonds | -0.679 | -0.641 | -0.617 |
| First-neighbor distance (Å) | +0.047 | +0.045 | +0.044 |
| Fourth-neighbor distance (Å) | +0.332 | +0.303 | +0.284 |
| Fifth-neighbor distance (Å) | +0.469 | +0.434 | +0.411 |
| Mean O--H length (Å) | -0.00361 | -0.00360 | -0.00343 |
| H--O--H angle (deg) | -0.271 | -0.231 | -0.254 |

Across all 12 intrinsic distributions, the median Jensen--Shannon distance
from PBC is 0.277, 0.264, and 0.256 in the R15, R20, and R30 outer 0--2 Å
layers. It falls to 0.037, 0.034, and 0.033 in 2--4 Å, then to 0.029, 0.021,
and 0.015 in 4--6 Å. Thus the network observables, nearest-neighbor
distributions, and edge-corrected RDF all locate the perturbation in the same
interfacial region.

## Interpretation

The three cavities do not contain three different bulk liquids. They share a
common local interfacial response and increasingly large PBC-like interiors.
The whole-cavity average can nevertheless vary with radius because the
fraction of molecules in the 0--4 Å surface region scales with
surface-to-volume ratio. R15 therefore mixes a much larger fraction of
interfacial water into a whole-cavity observable than R30.

Radial orientation observables are intentionally absent from the PBC
histogram comparison: a translationally invariant periodic liquid has no
physical cavity center or radial axis. Their NPBC radial profiles remain in
Figures 02 and 10 of the complete suite.

## New figures and data

- Figure 22: coarse radial-layer O--O RDFs and residuals against PBC.
- Figure 23: fine 1 Å map of `g_layer - g_PBC`.
- Figure 24: radial decay of RDF RMSE, coordination, and peak/minimum shifts.
- Figure 25: tetrahedrality, coordination, and H-bond histograms by layer.
- Figure 26: O--H length and H--O--H angle histograms by layer.
- Figure 27: first through fifth O--O neighbor histograms by layer.
- Figure 28: all-property Jensen--Shannon distance heatmap.

The complete numerical results are in the six `data/radial_layer_*` files,
and the full calculation is reproducible with
`compute_radial_layer_analysis.py`.
