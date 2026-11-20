clear all
set more off
set type double
capture log close
log using master_results.log, text replace
adopath ++ ".."
foreach id in did2s bjs flex csdid stacked {
    adopath ++ "../../fast`id'stata"
}
adopath ++ "../../../fastlpdidstata"
adopath ++ "../../suiteundidstata"
adopath ++ "../../suitedidintstata"
insheet using "../../didsuite_spec/fixtures/shared_assignments.csv", comma names clear
mkmat s00-s13, matrix(A)
insheet using "../../didsuite_spec/fixtures/stochastic_rcs.csv", comma names clear
didsuite y, time(year) cluster(state) gvar(g) datamode(rcs) inference(full) ///
    assignments(A) stackedwindow(-1 1) stackedcontrol(both) lpdidwindow(0 1) ///
    csv(master_stata.csv)
matrix R=e(results)
assert abs(R[1,1]-2.042422508994194)<1e-8
assert abs(R[2,1]-2.0550223209969785)<1e-8
assert abs(R[3,1]-2.055022320996979)<1e-8
assert abs(R[4,1]-2.0550223209969767)<1e-8
assert abs(R[5,1]-2.014521456826598)<1e-8
assert abs(R[6,1]-2.014521456826598)<1e-8
assert abs(R[7,1]-2.014521456826598)<1e-8
assert abs(R[8,1]-2.0202851572330545)<1e-8
assert abs(R[5,2]-.04707464786604753)<1e-8
assert abs(R[5,6]-.052713775270101536)<1e-8
foreach i in 1 2 3 4 5 8 {
    forvalues j=10/12 {
        assert R[`i',`j']==.25
    }
}
display "DIDSUITE_STATA_MASTER_PASS"
log close
