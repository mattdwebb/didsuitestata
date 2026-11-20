version 16
clear all
set more off
adopath ++ ".."
foreach p in fasttwfestata fastflexstata suiteundidstata suitedidintstata fastcsdid/stata {
    adopath ++ "../../`p'"
}
set obs 192
gen int state=ceil(_n/32)
bys state: gen int year=ceil(_n/8)
bys state year: gen byte person=_n
gen int g=cond(state==1,2,cond(state==2,3,0))
gen byte D=g>0 & year>=g
gen double x=cond(state==3,0,mod(person,2))
gen double y=cos(person+year*state)+state*.1+year*.15+x*.2+D*.3
* One treated state survives deletion of either singleton. All states must be deleted.
didsuite y x, time(year) cluster(state) gvar(g) estimators(undid didint flex csdid) ///
    inference(standard) csdidmethod(reg)
matrix d=e(deletions)
forvalues s=1/6 {
    quietly didsuite y x if state!=`s', time(year) cluster(state) gvar(g) ///
        estimators(undid didint flex csdid) inference(point) csdidmethod(reg)
    matrix b=e(results)
    forvalues j=1/4 {
        assert !missing(d[`s',`j'], b[`j',1])
        assert abs(d[`s',`j']-b[`j',1])<1e-10
    }
}
* With the late adoption in the final year, deleting the early adopter leaves exactly one ATT.
preserve
replace g=4 if state==2
quietly didsuite y x, time(year) cluster(state) gvar(g) estimators(undid didint csdid) inference(standard) csdidmethod(reg)
matrix single=e(deletions)
quietly didsuite y x if state!=1, time(year) cluster(state) gvar(g) estimators(undid didint csdid) inference(point) csdidmethod(reg)
matrix remaining=e(results)
forvalues j=1/3 {
    assert !missing(single[1,`j'],remaining[`j',1])
    assert abs(single[1,`j']-remaining[`j',1])<1e-10
}
restore
* State 3's constant x cannot identify the target x=.5 for DID-INT.
quietly suite_didint y x, time(year) cluster(state) gvar(g)
matrix c=e(att_st)
assert c[1,8]==3 & c[1,9]==4
quietly suite_undid y x, time(year) cluster(state) gvar(g)
matrix c=e(att_st)
assert c[1,8]==4
* A dependence among nuisance covariates leaves an identified prediction unchanged.
gen double duplicate=2*x
scalar theta=e(estimate)
quietly suite_undid y duplicate x, time(year) cluster(state) gvar(g)
assert abs(e(estimate)-theta)<1e-10
* If the last treated cluster is deleted, preserve the other valid deletion estimates.
replace g=0 if state==2
quietly didsuite y, time(year) cluster(state) gvar(g) estimators(twfe) inference(standard)
matrix b=e(results)
matrix d=e(deletions)
assert b[1,1]<. & missing(b[1,6])
assert missing(d[1,1]) & d[2,1]<. & d[6,1]<.
display "SURVIVING_ATT_TEST_PASS"
