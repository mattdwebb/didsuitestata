version 16
clear all
set more off
adopath ++ ".."
adopath ++ "../../fasttwfestata"
log using "smoke_results.log", text replace
set obs 225
generate int state=ceil(_n/25)
bysort state: generate int time=ceil(_n/5)
bysort state time: generate int person=_n
generate double g=cond(state<=3,3,cond(state<=6,4,0))
generate byte D=g>0 & time>=g
generate double x=sin(person*2+state*.3+time*.7)
generate double y=state*.3+time*.2+.4*x+2*D+sin(state*time)*.3+cos(person+state)*.2
generate long unit=state*10+person
quietly regress y D x i.state i.time, vce(cluster state)
scalar refb=_b[D]
scalar refse=_se[D]
didsuite y x, time(time) cluster(state) gvar(g) estimators(twfe) inference(standard)
matrix result=e(results)
assert abs(result[1,1]-refb)<1e-9
assert abs(result[1,2]-refse)<1e-9
assert result[1,8]==8
assert result[1,9]<.
display "PASS native reference and complete JK"
matrix A=(4,4,4,0,0,0,3,3,3 \ 0,0,0,3,3,3,4,4,4 \ 3,3,4,4,4,0,0,0,3)
didsuite y x, time(time) cluster(state) gvar(g) estimators(twfe lpdid) inference(full) assignments(A) csv(results.csv) saving(results.dta)
matrix result=e(results)
assert result[1,10]<.
assert result[1,11]<.
assert result[1,12]<.
assert missing(result[2,1])
assert _N==225
display "PASS full RI and failure isolation with caller data restored"
didsuite y x, time(time) cluster(state) treatment(D) estimators(twfe) inference(point)
matrix alt=e(results)
assert abs(alt[1,1]-refb)<1e-9
display "PASS treatment-only alternative"
preserve
replace D=0 if state==1 & time==5
capture didsuite y x, time(time) cluster(state) treatment(D) estimators(twfe) inference(point)
assert _rc!=0
restore
display "PASS nonabsorbing rejection"
matrix bad=A \ A[1,1...]
capture didsuite y x, time(time) cluster(state) gvar(g) estimators(twfe) inference(ri_beta) assignments(bad)
assert _rc!=0
display "PASS duplicate assignment rejection"
didsuite y x, time(time) cluster(state) gvar(g) estimators(twfe) inference(point) targetmode(harmonized) window(0 2)
matrix alt=e(results)
assert missing(alt[1,1])
display "PASS explicit unsupported harmonization"
didsuite y x, time(time) cluster(state) gvar(g) ivar(unit) datamode(panel) estimators(twfe) inference(standard)
matrix alt=e(results)
assert missing(alt[1,1])
display "PASS native TWFE panel limitation is explicit"
display "DIDSUITE_STATA_SMOKE_PASS"
log close

