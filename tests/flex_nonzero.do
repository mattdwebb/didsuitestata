version 16
clear all
set more off
adopath ++ ".."
adopath ++ "../../fastflexstata"
insheet using "../examples/school_policy.csv", comma names clear double
fastflex test_score baseline_score, time(year) cluster(school) gvar(first_treat) inference(native)
scalar null=e(estimate)-.5*e(native_se)
scalar expecteddf=e(native_df)
didsuite test_score baseline_score, time(year) cluster(school) gvar(first_treat) estimators(flex) inference(native) null(`=null')
matrix r=e(results)
assert abs(r[1,3]-.5)<1e-10
assert r[1,4]==expecteddf
assert abs(r[1,5]-2*ttail(expecteddf,.5))<1e-12
assert abs(r[1,5]-2*normal(-.5))>1e-3
display "FLEX_NONZERO_STUDENT_WRAPPER_PASS"
