> Historical development notes retained for provenance. For current behavior, installation and examples, use the [main user guide](../README.md). Old audit statements below may describe superseded behavior.

# didsuitestata

Native Stata 16+ harness, MIT copyright 2026 Matt Webb. Development milestone,
not full nine-estimator certification. Install with:

```stata
net install didsuite, from("absolute/path/to/didsuitestata")
didsuite_discover
didsuite y, time(year) cluster(state) gvar(g) estimators(all) ///
  inference(full) assignments(A) stackedwindow(-1 1) ///
  stackedcontrol(both) lpdidwindow(0 1) csv(results.csv)
```

Independently install optional native commands: fasttwfe, fastdid2s, fastbjs,
fastflex, fastcsdid, suite_undid, suite_didint, faststacked, fastlpdid. Every route
uses its native package, including TWFE. The current fasttwfe milestone explicitly
rejects individual-panel FE; the harness does not silently substitute areg for
that unsupported case. No production cross-language calls or Java-dependent CSV routines.
Original code is MIT; estimator/oracle licenses remain separate and unchanged.

## Contract

Common gvar uses 0 never-treated. Alternatively supply absorbing treatment();
inconsistent or reversing treatment fails. ivar() and datamode(panel) validate
unique panel units and nesting. Assignment and inference clusters must coincide.
Required model variables cannot be missing; analytic weights must be positive.
Weights are unsupported for UN-DID/DID-INT/LP except the default unit weights.

The supplied matrix has columns in Stata sorted cluster-label order, as in
e(observed_assignment). Column names are not validated. Generate and align the
common matrix outside Stata before passing it in. Automatic assignment generation
and canonical SHA256 are not implemented; hash_status is explicit, never invented.

Native controls are metadata: TWFE all observations; DID2S/BJS/FLEX untreated
never and not-yet; CSDID/UN-DID/DID-INT never; Stacked both; LP clean through
horizon. stackedwindow() routes only to Stacked and lpdidwindow() only to LP.
Stacked translates cohort0 to missing and uses true point calls for refits.
csdidmethod(none|reg|ipw|dripw) and csdidcontrol(never|not_yet) route through
all inference refits. CSDID defaults dripw, RCS only. LP defaults reweighted with
lag-one covariates, panel only. Other LP option selection is not yet exposed.

## Results And Validation

e(results) numeric matrix, e(estimators), e(observed_assignment); optional CSV/DTA
with status, failure reasons, native-control labels and target_compatible=0.
Runtime fields and assignment space are currently uninstrumented missing values;
estimator versions are adapter milestone labels. nobs currently means input N,
not estimator-specific trimmed/stacked N. No common-target equivalence claimed.

tests/smoke.do and tests/oracle.do compare native TWFE, physical deletions and
all RI statistics with direct Stata regress. validation/master.do uses the shared
sibling spec fixture/matrix and eight native RCS families. Tests use explicit
results logs; the /e auto log must have a different filename. Successful exit
alone is not evidence: inspect DIDSUITE_*_PASS markers and comparison CSV.
Stata tests were actually run under StataNow SE19; hosted CI needs a licensed
self-hosted runner. No broad performance claim or large Monte Carlo run.
