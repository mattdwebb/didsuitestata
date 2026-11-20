version 16
clear all
set more off
set type double
adopath ++ "."
foreach package in fasttwfestata fastdid2sstata fastbjsstata {
    adopath ++ "../`package'"
}
insheet using "examples/school_policy.csv", comma names clear
fasttwfe test_score baseline_score, time(year) cluster(school) ///
    gvar(first_treat) inference(ri_beta) rireps(99) seed(123)
matrix A=e(assignments)
assert rowsof(A)==99 & colsof(A)==12
didsuite test_score baseline_score, time(year) cluster(school) ///
    gvar(first_treat) estimators(twfe did2s bjs) datamode(rcs) ///
    inference(full) assignments(A)
matrix results=e(results)
forvalues j=1/3 {
    forvalues k=10/12 {
        assert !missing(results[`j',`k'])
        assert results[`j',`k']>=.01 & results[`j',`k']<=1
    }
}
display "SCHOOL_RANDOMIZATION_PASS"
