version 16
clear all
set more off
set type double
adopath ++ "."
foreach package in fasttwfestata fastdid2sstata fastbjsstata {
    adopath ++ "../`package'"
}
insheet using "examples/school_policy.csv", comma names clear
didsuite test_score baseline_score, time(year) cluster(school) gvar(first_treat) estimators(twfe did2s bjs) datamode(rcs) inference(standard) csv(examples/school_results.csv)
matrix school_results=e(results)
forvalues j=1/3 {
    assert !missing(school_results[`j',1],school_results[`j',6])
}
display "SCHOOL_EXAMPLE_PASS"
