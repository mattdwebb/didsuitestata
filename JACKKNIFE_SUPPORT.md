# Jackknife support and identification

The default delete-cluster rule is the same in the suite and its estimator
packages: retain a deletion when at least one sub-aggregate ATT supported by
the estimator remains estimable. Recompute the surviving aggregate, including
its weights and denominator. A missing cohort or changed weights alone is not
a failure. The point estimator, control eligibility, and event window remain
the estimator's documented definition.

Every input inference cluster contributes one deletion estimate. A deletion
with no supported ATT makes that jackknife unavailable; it is never omitted
from the variance calculation. The variance is centered on the mean of all
deletion estimates, and the standalone test uses Student t(G-1). At least two
clusters are required by the formula; each estimator must also identify each
deleted-sample ATT. Two clusters do not by themselves guarantee feasibility.

CSDID full-sample fits require every requested ATT(g,t) cell. An unidentified
cell returns or propagates IDENTIFICATION_FAILURE, rather than a reweighted
full-sample estimate labeled PASS. During a changing-support deletion, only
the surviving identified cells enter the aggregation. Python and R record per-deletion support changes. The Stata suite
retains all deletion estimates; it does not yet return a separate per-deletion
cell-support table. No successful-deletion-only denominator is used.

Python and R expose an optional fixed-support CSDID policy for users who
explicitly require the original cells. The default remains changing support.
The Stata suite currently implements changing support only; the standalone
Fast CSDID Stata command also offers jksupport(fixed). Thus the default policy
agrees across languages, while optional interface coverage differs.

The R reference-compatible propensity prefit rejects a cell when its maximum
fitted propensity reaches 0.999. This overlap diagnostic is distinct from the
0.995 control-observation trimming rule in the RCS IPW/DR cell calculations.
An overlap or rank failure is an identification failure, not an instruction
to drop a full-sample ATT cell. Do not interpret numerical agreement on tested
designs as certification of every overlap configuration.
