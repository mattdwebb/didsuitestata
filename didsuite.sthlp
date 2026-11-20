{smcl}
{title:Native DiD Suite}
{p 4 4 2}
Development harness, version 0.1.0; not universal estimator certification.

{title:Syntax}
{p 4 4 2}
{cmd:didsuite} outcome [covariates] [if] [in] [aw=weight],
{cmd:time(}time{cmd:)} {cmd:cluster(}cluster{cmd:)}
[{cmd:gvar(}cohort{cmd:)} {cmd:treatment(}binary{cmd:)}
{cmd:estimators(}names{cmd:)} {cmd:inference(}methods{cmd:)}
{cmd:assignments(}matrix{cmd:)} {cmd:ivar(}unit{cmd:)}
{cmd:datamode(}auto|rcs|panel{cmd:)}
{cmd:stackedwindow(}lo hi{cmd:)} {cmd:stackedcontrol(}both|never_treated|not_yet_treated{cmd:)}
{cmd:lpdidwindow(}lo hi{cmd:)} {cmd:csdidmethod(}method{cmd:)}
{cmd:csdidcontrol(}never|not_yet{cmd:)} {cmd:csv(}path{cmd:)} {cmd:saving(}path{cmd:)}]

{p 4 4 2}
Supply cohort or absorbing treatment. Cohort zero means never. Default native
controls are preserved. Stacked uses both and translates zero to missing.
Native windows are estimator-specific; global window() remains available for
single-estimator calls. Harmonized/common sample modes are unsupported.

{title:Inference And Results}
{p 4 4 2}
point, standard, full, or native jackknife ri_beta ri_t_native .
 aliases . All deletions and assigned denominators are refitted.
Undefined methods produce missing values and PARTIAL, retaining point estimates.
RI requires a unique supplied matrix excluding observed assignment, preserving
cohort counts, columns in sorted cluster-label order, and null zero.

{p 4 4 2}
Returns e(results), e(estimators), e(observed_assignment); optional CSV/DTA.
See README for exact adapter availability, unsupported methods, schema limitations,
validation provenance and independent installation of native estimators.

{title:Estimable support and weights}

{p 4 4 2}The jackknife refits after every cluster deletion, allowing adoption
cohorts to disappear and reweighting the surviving estimable effects. A valid
variance requires all deletion estimates. {cmd:e(deletions)} stores them in
sorted cluster order, with columns in estimator order, for standalone jackknife requests.

{p 4 4 2}UN-DID tests whether the post-period coefficient is estimable in each
local regression. DID-INT compares predicted changes only where both predictions
at the target covariate mean are estimable. Within each treated state-year the
supported controls are reweighted, followed by reweighting across supported
treated state-years. The exported cells_used, min_controls and max_controls
describe this realized support. LP-DiD requires a longitudinal individual panel.

{p 4 4 2}{cmd:undidweights(post)} is the native default. For the merit paper's
simple UN-DID, {cmd:undidweights(comparison)} weights treated-state changes inside
each cohort-year by both periods and cohort-years by their post counts.
{cmd:csdidweights(fullwindow)} is the native CSDID default;
{cmd:csdidweights(comparison)} weights each supported cohort-year by treated
counts or sampling-weight mass in its two comparison years. The same choices
apply to all deletions and reassignments. Unidentified CSDID cells are excluded,
with their counts reported as cells_used and cells_eligible.

{p 4 4 2}{cmd:csdidbackend(original)} explicitly calls the separately installed
Stata {cmd:csdid} command with {cmd:agg(simple)} and its comparison-period
weights. Use {cmd:csdidmethod(dripw)} for the original adjusted merit example.
The default backend is the independent suite implementation. The CSV identifies
the selected backend; reference and suite support conventions may differ.
