# Estimators and targets

`didsuite` routes each requested estimator to its native Stata package. It does
not change the estimator to force a numerical comparison. Inputs name the outcome,
time, inference cluster and first treatment date; zero denotes never treated.
Treatment is absorbing and constant within cluster and period. Required model
variables must be observed. Panel mode additionally requires a genuine unit ID.

| Family | Calculation | Default comparison/aggregation |
|---|---|---|
| TWFE | Policy coefficient with cluster and time fixed effects | All included observations; a coefficient, not automatically ATT |
| did2s | Untreated first stage followed by adjusted-outcome estimation | Never and not-yet-treated first-stage observations; treated-observation aggregate |
| BJS | Untreated-outcome prediction and imputation | Treated-observation aggregate |
| FLEX | Flexible state-group regression and fitted effects | Supported treated state-time effects |
| CSDID | Stationary RCS group-time comparisons | Never-treated controls; full-window cohort masses |
| UN-DID | Treated-state/time contrasts with controls | Never-treated controls; treated post-period counts |
| DID-INT | State-year fitted models standardized at treated covariates | Fully interacted slopes; eligible controls and supported effects |
| Stacked | Cohort stacks and corrective-weighted event regression | Clean future-treated and never-treated controls; explicit window |
| LP-DiD | Local-projection contrasts around adoption | Reweighted pooled horizons; lag-one covariates by default |

The wrapper translates zero never-treated codes to the standalone stacked
convention. Its default CSDID backend is the suite implementation; selecting
the separately supported original backend changes the documented comparison-mass
convention. It never silently substitutes a different command when a dependency
is absent. `didsuite_discover` shows the available commands.

The native definitions need not target the same weighted ATT under heterogeneous
effects. Harmonized/common-sample modes and distinct assignment/inference units
are not available. An explicit event window applies only where the estimator
supports it; `stackedwindow()` and `lpdidwindow()` select separate windows.

See each dependency's method specification and original-paper references:
[TWFE](https://github.com/mattdwebb/fasttwfestata),
[did2s](https://github.com/mattdwebb/fastdid2sstata),
[BJS](https://github.com/mattdwebb/fastbjsstata),
[FLEX](https://github.com/mattdwebb/fastflexstata),
[CSDID](https://github.com/mattdwebb/fastcsdid/tree/main/stata),
[UN-DID](https://github.com/mattdwebb/suiteundidstata),
[DID-INT](https://github.com/mattdwebb/suitedidintstata),
[stacked DiD](https://github.com/mattdwebb/faststackedstata), and
[LP-DiD](https://github.com/mattdwebb/fastlpdidstata).
Development repositories remain access-controlled until their visibility is
changed separately.
