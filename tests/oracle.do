version 16

clear all

set more off

adopath ++ ".."

adopath ++ "../../fasttwfestata"

log using "oracle_results.log", text replace

set obs 225

generate int state=ceil(_n/25)

bysort state: generate int time=ceil(_n/5)

bysort state time: generate int person=_n

generate double g=cond(state<=3,3,cond(state<=6,4,0))

generate byte D=g>0 & time>=g

generate double x=sin(person*2+state*.3+time*.7)

generate double y=state*.3+time*.2+.4*x+2*D+sin(state*time)*.3+cos(person+state)*.2

matrix A=(4,4,4,0,0,0,3,3,3 \ 0,0,0,3,3,3,4,4,4 \ 3,3,4,4,4,0,0,0,3)

quietly didsuite y x, time(time) cluster(state) gvar(g) estimators(twfe) inference(full) assignments(A)

matrix got=e(results)

matrix stats=J(4,2,.)

forvalues draw=0/3 {

    generate double gd=g

    if `draw'>0 {

        forvalues c=1/9 {

            quietly replace gd=A[`draw',`c'] if state==`c'

        }

    }

    generate byte dt=gd>0 & time>=gd

    quietly regress y dt x i.state i.time, vce(cluster state)

    local row=`draw'+1

    matrix stats[`row',1]=_b[dt]

    matrix stats[`row',2]=_b[dt]/_se[dt]

    if `draw'==0 {
        matrix deletions=J(9,1,.)
        forvalues c=1/9 {
            quietly regress y dt x i.state i.time if state!=`c'
            matrix deletions[`c',1]=_b[dt]
        }
        mata: st_numscalar("jkse",sqrt(8/9*sum((st_matrix("deletions"):-mean(st_matrix("deletions"))):^2)))
        assert abs(got[1,6]-jkse)<1e-9
        assert abs(got[1,7]-stats[1,1]/jkse)<1e-9
        assert abs(got[1,9]-2*ttail(8,abs(stats[1,1]/jkse)))<1e-12
    }
    drop gd dt

}

forvalues method=1/2 {

    scalar hits=0

    forvalues i=2/4 {

        scalar hits=hits+(abs(stats[`i',`method'])>=abs(stats[1,`method'])-1e-12)

    }

    assert abs(got[1,`=`method'+9']-(1+hits)/4)<1e-12

}

display "PASS independent dummy-regression native and all three RI statistics"

generate double w=1+mod(person,3)/4

quietly regress y D x i.state i.time [aw=w], vce(cluster state)

scalar wb=_b[D]

scalar ws=_se[D]

quietly didsuite y x [aw=w], time(time) cluster(state) gvar(g) estimators(twfe) inference(standard)

matrix got=e(results)

assert abs(got[1,1]-wb)<1e-9

assert abs(got[1,2]-ws)<1e-9

display "PASS weighted native oracle"

generate str10 sid="state"+string(state)

quietly didsuite y x, time(time) cluster(sid) gvar(g) estimators(twfe) inference(point)

matrix got=e(results)

assert got[1,1]<.

display "PASS string cluster labels"

preserve

drop if inlist(state,2,3)

quietly didsuite y x, time(time) cluster(state) gvar(g) estimators(twfe) inference(standard)

matrix got=e(results)

assert got[1,1]<. & got[1,6]<.

matrix dels=e(deletions)

quietly regress y D x i.state i.time if state!=1

assert abs(dels[1,1]-_b[D])<1e-9

replace g=0 if state!=1

quietly didsuite y x, time(time) cluster(state) gvar(g) estimators(twfe) inference(standard)

matrix got=e(results)

assert got[1,1]<. & missing(got[1,6])

restore

display "PASS singleton deletion does not erase point"

quietly didsuite y x, time(time) cluster(state) gvar(g) estimators(all) inference(point)

matrix got=e(results)

assert rowsof(got)==9

display "PASS all estimator rows retained"

matrix Bad=(3,3,3,4,4,4,0,0,0)

capture didsuite y x, time(time) cluster(state) gvar(g) estimators(twfe) inference(ri_beta) assignments(Bad)

assert _rc!=0

matrix Bad=A

matrix Bad[1,1]=5

capture didsuite y x, time(time) cluster(state) gvar(g) estimators(twfe) inference(ri_beta) assignments(Bad)

assert _rc!=0

display "PASS observed and cohort-count rejection"



display "DIDSUITE_STATA_ORACLE_PASS"

log close

exit, clear



