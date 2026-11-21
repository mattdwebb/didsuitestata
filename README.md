# DiD Suite for Stata

**Compare difference-in-differences estimates and their uncertainty in one table.**

Suppose some schools adopt a policy earlier than others. You may want to see how
the estimated effect changes across TWFE, imputation, group-time comparisons,
and other DiD approaches. `didsuite` runs the supported native Stata commands,
organizes their results, and exports a comparison file. It also applies a common
delete-cluster jackknife without requiring you to write the deletion loop.

The suite does not make the estimators target the same weighted treatment effect.
Their treatment windows, comparison groups, and aggregation rules still matter.

## Install

Requires Stata 16 or later. No R or Python runtime is used for estimation.
The estimator packages are separate: installing this wrapper alone does not
install them. Install the wrapper and the three backends used below from GitHub:

```stata
net install didsuite, from("https://raw.githubusercontent.com/mattdwebb/didsuitestata/main/") replace
net install fasttwfe, from("https://raw.githubusercontent.com/mattdwebb/fasttwfestata/main/") replace
net install fastdid2s, from("https://raw.githubusercontent.com/mattdwebb/fastdid2sstata/main/") replace
net install fastbjs, from("https://raw.githubusercontent.com/mattdwebb/fastbjsstata/main/") replace
discard
didsuite_discover
```

For CSDID, use `net install fastcsdid, from("https://raw.githubusercontent.com/mattdwebb/fastcsdid/v1.0.1/stata/") replace`.
The public Fast CSDID distribution is the canonical CSDID dependency. For
UN-DID, install the separately named `suite_undid` adapter; it does not replace
the upstream SSC `undid` command. Source checkouts can also be added with
`adopath ++` for a session-local trial.

The commands found by `didsuite_discover` determine what can run. Do not substitute
unrelated installed commands for missing suite dependencies.

## A first comparison

The bundled fictional example has 960 pupil observations in twelve schools over
eight years. Four schools adopt in 2013, four in 2015, and four remain untreated.
The outcome is a test score; the adjustment variable is a predetermined baseline
score. Pupils are repeated cross sections, and uncertainty is clustered by school.
See the [data guide](examples/DATA_GUIDE.md) for definitions and the generator.

```stata
clear
set type double
insheet using "examples/school_policy.csv", comma names clear
didsuite test_score baseline_score, time(year) cluster(school) ///
    gvar(first_treat) estimators(twfe did2s bjs) datamode(rcs) ///
    inference(standard) csv(comparison.csv)
matrix list e(results)
```

The executed comparison includes:

| Estimator | Estimate (score points) | Native SE | Jackknife SE |
|---|---:|---:|---:|
| twfe | 2.012599 | 0.367339 | 0.401294 |
| did2s | 2.316869 | 0.285375 | 0.336814 |
| bjs | 2.316869 | 0.290170 | 0.336814 |

TWFE estimates about 2.01 score points; did2s and BJS estimate about 2.32 points.
The latter use the same untreated-outcome model and aggregate here, so their
point estimates and deletion estimates agree. Their analytic variance formulas
differ, which explains their different native standard errors. The jackknife
uses all twelve school deletions and Student t(11) inference.

The artificial data generator adds a two-point effect, but these numbers are a
worked software example, not evidence about an actual school policy or a test
of the estimators' general performance. `comparison.csv` contains estimates,
standard errors, p-values, targets, and status information at full precision.
The command restores the input dataset after building its results table.

`do examples/school_policy.do` reproduces the comparison and verifies that the
point estimates and jackknife standard errors are available.

## Add estimators and choose a specification

| Estimator | Stata dependency | Principal scope |
|---|---|---|
| TWFE | fasttwfe | Pooled coefficient, cluster/time fixed effects |
| Gardner two-stage DiD | fastdid2s | Static aggregate from untreated first stage |
| BJS imputation | fastbjs | Static imputation aggregate |
| FLEX | fastflex | Historical state-group flexible specification |
| CSDID | fastcsdid | Stationary repeated cross sections; none/reg/ipw/dripw |
| UN-DID | suite_undid | Never-treated comparisons; no native analytic SE |
| DID-INT | suite_didint | Fully interacted fitted-value specification; no native analytic SE |
| Corrective-weighted stacked DiD | faststacked | Explicit event window containing -1 |
| LP-DiD | fastlpdid | Genuine panel; explicit horizon window |

For the eight supported RCS estimators, after making their repositories available:

```stata
didsuite test_score baseline_score, time(year) cluster(school) ///
    gvar(first_treat) datamode(rcs) ///
    estimators(twfe did2s bjs flex csdid undid didint stacked) ///
    csdidmethod(reg) stackedwindow(-1 1) inference(standard)
```

This explicitly selects regression-adjusted CSDID. `csdidmethod(dripw)` selects
the nonlinear doubly robust method; its exact-refit jackknife may take longer.
The default CSDID method is dripw. Missing native inference for UN-DID and
DID-INT is expected and is not replaced by a jackknife under a native label.

For LP-DiD, use the separate genuine school panel:

```stata
insheet using "examples/school_panel.csv", comma names clear
didsuite test_score baseline_score, time(year) cluster(school) ///
    ivar(school) gvar(first_treat) datamode(panel) estimators(lpdid) ///
    lpdidwindow(0 1) inference(standard)
```

Do not request all nine on RCS data and expect nine supported rows: LP-DiD uses
panel equations. Likewise, the supported CSDID implementation uses RCS equations.
Omit `baseline_score` for an unadjusted model. Never-treated cohorts are zero in
the common interface; the wrapper translates that convention for stacked DiD.

## Inference and unavailable results

`inference(point)` requests estimates only; `inference(standard)` adds native
clustered and jackknife inference. The jackknife refits each deleted sample,
reaggregating surviving estimable effects and weights. It requires one estimate
for every input inference cluster. A failed deletion does not erase a valid
full-sample estimate, and is not dropped to manufacture a standard error.

Read each row's status and failure reason. `PARTIAL` can mean that a point
estimate and jackknife are valid but native inference is unavailable.
`UNSUPPORTED` identifies an unavailable package or specification. Missing
values are not zero effects or zero p-values. The exported `nobs` is input
sample size, not the size of every estimator's internally trimmed/stacked sample.

## Help and reproducibility

- [Command help](didsuite.sthlp)
- [Method and targets](METHOD_SPECIFICATION.md)
- [Inference definitions](INFERENCE.md)
- [Technical development reference](docs/DEVELOPMENT_REFERENCE.md)
- [License](LICENSE); estimator packages retain their own licenses and citations.

Record the release commit manifest along with your analysis. The package uses
native targets; harmonized common-target/common-sample modes are not advertised
as implemented. Developers can run `tests/run_release_tests.ps1` with a licensed
Stata installation. Report reproducible problems through this repository's Issues
page using artificial or otherwise shareable data.


## Inference scope and deletion support

The suite provides coefficient RI (`ri_beta`) and native-studentized RI
(`ri_t_native`) in addition to native and standalone jackknife inference.
These tests are implemented by the suite. Standalone Fast CSDID intentionally
has a narrower public inference API, so its API checker rejects those option
names; that restriction does not remove the suite's separately tested methods.
See [jackknife support and identification](JACKKNIFE_SUPPORT.md) for surviving
ATT aggregation, full-sample failures, overlap limits, and interface differences.
