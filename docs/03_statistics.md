# Statistical conventions

## First temporal standard deviation

For a shell observable evaluated at each frame, the displayed band is

\[
\overline A_s\pm\sqrt{\frac{1}{n-1}\sum_t[A_s(t)-\overline A_s]^2}.
\]

It is a temporal standard deviation, not an SEM and not a confidence interval.

## Block uncertainty

RDF and radial-distribution uncertainties use equal 100 ps blocks.  The RDF is
recomputed from each block's numerator and denominator before calculating the
SD across block estimates.  Averaging pre-normalized frame RDFs would generally
be a different estimator.

## Distribution distances

The workflow reports:

- base-2 Jensen--Shannon **distance** (the square root of the divergence);
- one-dimensional Wasserstein distance on histogram bin midpoints;
- standardized mean difference using the pooled distribution SD.

Fixed bin edges are shared by every radius and PBC.  Comparing histograms with
different bins is not allowed.

## Practical-equivalence screen

An observable passes only if all conditions hold:

```text
|standardized mean difference| <= 0.2
Jensen--Shannon distance       <= 0.1
normalized Wasserstein         <= 0.2
95% bootstrap CI lies inside +/- 0.2 pooled SD
```

This is a conservative descriptive screen.  It does not prove equality of
probability measures and should not be reported as a null-hypothesis test.

## Experimental comparison

The Skinner O--O curve is interpolated to simulation bin centers.  Continuous
experimental O--H and H--H curves are not fabricated; those comparisons use
published landmarks only.  The experimental and simulated state points are
not exactly density matched, so residuals are descriptive.
