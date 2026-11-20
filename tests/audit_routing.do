clear all
set more off
set type double
adopath ++ ".."
adopath ++ "../../fastlpdidstata"
adopath ++ "../../fasttwfestata"
log using audit_routing_results.log, text replace
set obs 225
gen int state=ceil(_n/25)
bysort state: gen int time=ceil(_n/5)
bysort state time: gen int person=_n
gen double g=cond(state<=3,3,cond(state<=6,4,0))
gen byte D=g>0 & time>=g
gen double x=sin(person*2+state*.3+time*.7)
gen double y=state*.3+time*.2+.4*x+2*D+sin(state*time)*.3+cos(person+state)*.2
gen long unit=state*10+person
foreach cov in none x {
    local xopt
    local nativeopt
    if "`cov'"=="x" {
        local xopt x
        local nativeopt covariates(x)
    }
    fastlpdid y, unit(unit) time(time) cluster(state) gvar(g) window(0 1) `nativeopt'
    scalar refb=e(estimate)
    scalar refs=e(native_se)
    scalar refdf=e(native_df)
    scalar refp=e(native_p)
    if "`cov'"=="none" assert refdf==8
    else assert missing(refdf)
    didsuite y `xopt', time(time) cluster(state) gvar(g) ivar(unit) datamode(panel) ///
        estimators(lpdid) lpdidwindow(0 1) inference(native)
    matrix R=e(results)
    assert abs(R[1,1]-refb)<1e-10
    assert abs(R[1,2]-refs)<1e-10
    assert R[1,4]==refdf
    assert abs(R[1,5]-refp)<1e-12
}
fasttwfe y x, time(time) cluster(state) gvar(g) inference(native)
scalar refb=e(estimate)
scalar refs=e(native_se)
scalar refdf=e(native_df)
scalar refp=e(native_p)
didsuite y x, time(time) cluster(state) gvar(g) estimators(twfe) inference(native)
matrix R=e(results)
assert abs(R[1,1]-refb)<1e-10
assert abs(R[1,2]-refs)<1e-10
assert R[1,4]==refdf
assert abs(R[1,5]-refp)<1e-12
preserve
drop if inlist(state,5,6)
fastlpdid y, unit(unit) time(time) cluster(state) gvar(g) window(0 1) inference(point)
assert e(estimate)<.
* The late cohort has only state 4. Its deletion leaves the early-cohort ATT
* estimable, so retain that deletion and recompute the supported pooled effect.
fastlpdid y, unit(unit) time(time) cluster(state) gvar(g) window(0 1) inference(jackknife)
matrix support_deletions=e(jk_deletions)
scalar reported_jk=e(jk_se)
assert rowsof(support_deletions)==7
assert e(jk_df)==6
levelsof state, local(surviving_states)
tempfile support_sample
quietly save "`support_sample'"
local row=0
foreach state_id of local surviving_states {
    local ++row
    quietly use "`support_sample'", clear
    drop if state==`state_id'
    quietly fastlpdid y, unit(unit) time(time) cluster(state) gvar(g) window(0 1) inference(point)
    assert e(estimate)<.
    assert abs(e(estimate)-support_deletions[`row',1])<1e-10
}
mata: st_numscalar("reference_jk",sqrt(6/7*sum((st_matrix("support_deletions"):-mean(st_matrix("support_deletions"))):^2)))
assert abs(reported_jk-reference_jk)<1e-10
restore
display "DIDSUITE_AUDIT_ROUTING_PASS"
log close
