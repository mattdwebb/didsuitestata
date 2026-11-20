# Choosing inference

## Native and jackknife uncertainty

Native inference preserves the estimator's own variance calculation and reference
distribution. It is unavailable for the suite UN-DID and DID-INT aggregates.
The wrapper does not label their jackknife results as native analytic results.

With G inference clusters, every cluster is deleted once. If the deleted estimates
are theta[-1], ..., theta[-G], their mean is m and the mean-centered jackknife
variance is `(G-1)/G * sum((theta[-g]-m)^2)`. The test uses Student t(G-1).

Each deletion refits the supported estimator. When a cohort disappears, the
aggregate is recomputed over surviving estimable effects, including weights and
denominators. Every required scalar deletion must exist. A failed deletion or
zero/invalid variance makes that inference unavailable; a valid full-sample
point estimate is retained. Failed deletions are never dropped to reduce G.

## A runnable RI example

The wrapper accepts supplied assignments rather than generating them itself.
The native TWFE dependency can generate a count-preserving reference set on the
same data. After the installation steps in the main README, run:

```stata
clear
set type double
insheet using "examples/school_policy.csv", comma names clear
fasttwfe test_score baseline_score, time(year) cluster(school) ///
    gvar(first_treat) inference(ri_beta) rireps(99) seed(123)
matrix A=e(assignments)
didsuite test_score baseline_score, time(year) cluster(school) ///
    gvar(first_treat) estimators(twfe did2s bjs) datamode(rcs) ///
    inference(full) assignments(A)
```

Here 99 alternatives keep the teaching example small; the p-value resolution
is 0.01. Choose an adequate reference count for the actual research design.
Use the saved matrix to compare estimators on identical assignments.

Rows of A are unique, nonobserved alternative assignments. Columns follow sorted
cluster identifiers, matching `e(observed_assignment)` and `e(cluster_order)`.
Each row must preserve the observed cohort counts. Stata matrix column names
are not currently used to validate alignment, so preserve the numeric ordering.

Only a zero sharp null is supported by these RI routes. The finite-sample
randomization justification requires the observed and reference assignments to
follow the appropriate same assignment mechanism under that null. Studentizing
does not make an arbitrary observational assignment mechanism exact.

See [the runnable RI example](examples/randomization.do) and [the main user guide](README.md).
