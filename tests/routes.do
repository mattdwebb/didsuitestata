clear all
set more off
set type double
adopath ++ ".."
adopath ++ "../../fastcsdid/stata"
adopath ++ "../../fastlpdidstata"
adopath ++ "../../fasttwfestata"
log using routes_results.log, text replace
set obs 225
generate int state=ceil(_n/25)
bysort state: generate int time=ceil(_n/5)
bysort state time: generate int person=_n
generate double g=cond(state<=3,3,cond(state<=6,4,0))
generate byte D=g>0 & time>=g
generate double x=sin(person*2+state*.3+time*.7)
generate double y=state*.3+time*.2+.4*x+2*D+sin(state*time)*.3+cos(person+state)*.2
generate long unit=state*10+person
generate double w=1+mod(person,3)/4
foreach m in none reg ipw dripw {
    local xopt
    local covopt
    if "`m'"!="none" {
        local xopt x
        local covopt covariates(x)
    }
    fastcsdid y, time(time) cluster(state) gvar(g) method(`m') weights(w) `covopt'
    scalar refb=_b[ATT]
    scalar refs=e(native_se)
    didsuite y `xopt' [aw=w], time(time) cluster(state) gvar(g) ///
        estimators(csdid) csdidmethod(`m') inference(standard)
    matrix R=e(results)
    assert abs(R[1,1]-refb)<1e-10
    assert abs(R[1,2]-refs)<1e-10
    assert R[1,6]<.
}
display "PASS all four weighted CSDID methods and complete JK routing"
fastcsdid y, time(time) cluster(state) gvar(g) control(not_yet)
scalar refb=_b[ATT]
didsuite y, time(time) cluster(state) gvar(g) estimators(csdid) ///
    csdidcontrol(not_yet) inference(point)
matrix R=e(results)
assert abs(R[1,1]-refb)<1e-10
display "PASS CSDID not-yet control routing"
matrix A=(4,4,4,0,0,0,3,3,3 \ 0,0,0,3,3,3,4,4,4 \ 3,3,4,4,4,0,0,0,3)
fastlpdid y, unit(unit) time(time) cluster(state) gvar(g) window(0 1)
scalar refb=e(estimate)
scalar refs=e(native_se)
scalar refdf=e(native_df)
scalar refp=e(native_p)
didsuite y, time(time) cluster(state) gvar(g) ivar(unit) datamode(panel) ///
    estimators(lpdid) lpdidwindow(0 1) inference(full) assignments(A)
matrix R=e(results)
assert abs(R[1,1]-refb)<1e-10
assert abs(R[1,2]-refs)<1e-10
assert R[1,4]==refdf
assert abs(R[1,5]-refp)<1e-12
forvalues j=6/12 {
    assert R[1,`j']<.
}
display "PASS LP panel full inference routing"
fasttwfe y x [aw=w], time(time) cluster(state) gvar(g) inference(native)
scalar refb=e(estimate)
scalar refs=e(native_se)
scalar refdf=e(native_df)
scalar refp=e(native_p)
didsuite y x [aw=w], time(time) cluster(state) gvar(g) estimators(twfe) inference(native)
matrix R=e(results)
assert abs(R[1,1]-refb)<1e-10
assert abs(R[1,2]-refs)<1e-10
assert R[1,4]==refdf
assert abs(R[1,5]-refp)<1e-12
display "PASS TWFE native package routing"
display "DIDSUITE_STATA_ROUTES_PASS"
log close
