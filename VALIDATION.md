# Validation Milestone

2026-09-29, actual StataNow SE19, sequential hidden batch processes, native
thread caps one. Explicit markers verified, not exit-code assumptions:

- DIDSUITE_STATA_SMOKE_PASS
- DIDSUITE_STATA_ORACLE_PASS
- DIDSUITE_STATA_ROUTES_PASS
- DIDSUITE_STATA_MASTER_PASS
- DIDSUITE_STATA_EXPORT_COMPARISON_PASS (export comparison)

Smoke/oracle scripts test direct regress versus areg point/CR1, physical
mean-centered deletion SE, native and JK reference probabilities, all supplied RI
statistics, weights, string IDs, panel/treatment modes, aliases, invalid matrices
and method failure isolation. Routes test weighted CSDID none/reg/ipw/dripw and
complete JK against direct native commands, not-yet controls, and LP panel full
inference against its native command.

The shared 549-row fixture, 14 clusters and nine times, uses exactly three
supplied alternatives from sibling didsuite_spec/fixtures/shared_assignments.csv.
Expected externally supplied hash:
bea08f08f138e2343a8e1c618a7542169a9c56f966dd174dbac3d240abc3b4ec
Stata does not compute this hash itself; table honestly reports UNSUPPORTED_SHA256.

Eight RCS point estimates and all available native/JK/RI outputs match the
supplied Python CSV at 1e-8; all valid RI p-values .25. UN-DID/DID-INT native
VCE/RI-t remain missing and rows PARTIAL. LP correctly rejects RCS. Stacked
uses full window(-1,1), both controls; all other native targets unrestricted.
See validation/master_python_comparison.csv for exact final numeric gaps.

Precision correction: result postfile columns are explicitly double to avoid
float rounding; exports use 17-digit formatting. No estimator/oracle source
or license was copied, modified or relicensed.

Remaining: native SHA256/generator, timings, exact version metadata and trimmed
sample sizes, cell-support-signature comparison, harmonized/common sample
contracts, more LP options, broader external-oracle evidence. Stata install
package includes all native helpers/help/license; remote CI requires a separately
licensed self-hosted runner and was not remotely executed at this checkpoint.
