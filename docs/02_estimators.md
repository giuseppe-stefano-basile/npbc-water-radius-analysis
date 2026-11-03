# Estimators

## Radial density

For shell `s = [r_-, r_+)`,

\[
V_s=\frac{4\pi}{3}(r_+^3-r_-^3),\qquad
\widetilde\rho_s(t)=\frac{n_s(t)}{V_s\rho_0}.
\]

The mean is taken over per-frame normalized densities.  An empty shell is a
valid zero for density but is missing for a molecular property whose shell
mean would otherwise divide by zero.

## Molecule and volume weighting

Molecule weighting pools all molecules before normalization.  Volume
weighting first normalizes the distribution in each shell and then averages
shell distributions with equal shell weight.  They answer different questions:

- molecule weighting: what is observed for a randomly selected molecule?
- volume weighting: what is observed for a randomly selected equal-volume
  shell element?

## Boundary-safe RDF

For maximum pair distance `r_max`, NPBC reference atoms satisfy

\[
|\mathbf r_{\mathrm{ref}}|\le R-r_{\max}.
\]

Their complete pair-distance spheres lie within the cavity; ordinary shell
normalization is then geometrically valid.  Neighbors are not restricted.

## Exact radial-layer RDF

For a reference oxygen at distance `d` from the cavity center, the available
volume of pair bin `[r_-,r_+)` is

\[
A=\mathcal V(R,r_+,d)-\mathcal V(R,r_-,d),
\]

where `V` is the analytic two-sphere intersection volume.  The layer RDF is

\[
g_L(r)=\frac{C_L(r)}{[(N_O-1)/V_{\mathrm{cav}}]\sum_{i\in L}A_i(r)}.
\]

Counts are directed because every molecule in the selected layer acts as a
reference.  The `(N_O-1)/V` density is essential: an oxygen cannot pair with
itself.

## Direct coordination versus corrected RDF

The direct coordination count retains the real absence of molecules outside
the wall.  The corrected RDF asks whether correlations in the geometrically
available volume differ from bulk.  A coordination deficit and a nearly
bulk-like corrected RDF can therefore occur simultaneously.

## Frozen Gaussian field

The molecular mean field is

\[
U(r)=s\sum_k(a_k-a_0)
\exp[-(r-r_k)^2/(2\sigma_k^2)],\qquad s=0.0005\ \mathrm{eV},
\]

and the plotted outward force is evaluated analytically as `-dU/dr`.  It is
not obtained by numerically differentiating a plotted curve.
