*! didsuitestata 0.1.0, MIT 2026 Matt Webb
program define didsuite, eclass
    version 16
    syntax varlist(min=1 numeric) [if] [in] [aw], TIME(varname numeric) CLuster(varname) ///
        [GVAR(varname numeric) TREATment(varname numeric) ESTimators(string) ///
        INference(string) ASSignments(name) IVAR(varname) ASSIGNUnit(varname) ///
        DATAMode(string) TARGETMode(string) SAMPLEMode(string) ///
        WINDOW(numlist integer min=2 max=2) STACKEDWindow(numlist integer min=2 max=2) ///
        STACKEDControl(string) CSDIDMethod(string) CSDIDControl(string) CSDIDWeights(string) CSDIDBackend(string) UNDIDWeights(string) ///
        LPDIDWindow(numlist integer min=2 max=2) ///
        NULL(real 0) SAVing(string) CSV(string)]
    if "`gvar'`treatment'"=="" exit 198
    if missing(`null') exit 198
    if "`assignunit'"!="" & "`assignunit'"!="`cluster'" exit 198
    if "`estimators'"=="" | "`estimators'"=="all" local estimators "twfe did2s bjs flex csdid undid didint stacked lpdid"
    if !inlist("`undidweights'","","post","comparison") exit 198
    if !inlist("`csdidweights'","","fullwindow","comparison") exit 198
    if !inlist("`csdidbackend'","","suite","original") exit 198
    if "`csdidbackend'"=="original" & "`csdidweights'"=="fullwindow" exit 198
    if "`csdidbackend'"=="original" local csdidweights comparison
    local allowed twfe did2s bjs flex csdid undid didint stacked lpdid
    local seen
    foreach id of local estimators {
        if !`: list id in allowed' | `: list id in seen' exit 198
        local seen `seen' `id'
    }
    if "`inference'"=="" local inference "standard"
    if "`inference'"=="point" local inference ""
    if "`inference'"=="standard" local inference "native jackknife"
    if "`inference'"=="full" local inference "native jackknife ri_beta ri_t_native"
    local inference : subinstr local inference "ribeta" "ri_beta", all
    local inference : subinstr local inference "ritnative" "ri_t_native", all
    local inference : list uniq inference
    foreach method of local inference {
        if !inlist("`method'","native","jackknife","ri_beta","ri_t_native") exit 198
    }
    local hasri= strpos("`inference'","ri_")>0
    if `hasri' & ("`assignments'"=="" | `null'!=0) {
        di as error "UNSUPPORTED: zero-null RI requires supplied assignments(matrix)"
        exit 198
    }
    if !`hasri' & "`assignments'"!="" exit 198
    if "`datamode'"=="" | "`datamode'"=="auto" {
        local datamode "rcs"
        if "`ivar'"!="" {
            capture isid `ivar'
            if _rc local datamode "panel"
        }
    }
    if !inlist("`datamode'","panel","rcs") exit 198
    if "`datamode'"=="panel" & "`ivar'"=="" exit 198
    if "`targetmode'"=="" local targetmode "native"
    if "`samplemode'"=="" local samplemode "native"
    if !inlist("`targetmode'","native","harmonized") | !inlist("`samplemode'","native","common") exit 198
    marksample touse, novarlist
    preserve
    quietly keep if `touse'
    foreach v of varlist `varlist' `time' `cluster' `gvar' `treatment' `ivar' {
        quietly count if missing(`v')
        if r(N)>0 exit 459
    }
    if _N==0 exit 2000
    tempvar cc gg dd w
    quietly egen long `cc'=group(`cluster')
    quietly gen double `w'=1
    if "`weight'"!="" quietly replace `w' `exp'
    quietly count if missing(`w') | `w'<=0
    if r(N)>0 exit 459
    local wopt "[aw=`w']"
    if "`datamode'"=="panel" {
        quietly isid `ivar' `time'
        bysort `ivar': assert `cc'==`cc'[1]
    }
    if "`gvar'"!="" {
        quietly gen double `gg'=`gvar'
        bysort `cc': assert `gg'==`gg'[1] & `gg'>=0
    }
    else {
        quietly assert inlist(`treatment',0,1)
        quietly bysort `cc': egen double `gg'=min(cond(`treatment'==1,`time',.))
        quietly replace `gg'=0 if missing(`gg')
    }
    quietly gen byte `dd'=(`gg'>0 & `time'>=`gg')
    if "`treatment'"!="" quietly assert `dd'==`treatment'
    quietly levelsof `cluster', local(cluster_order)
    quietly levelsof `cc', local(ids)
    local G : word count `ids'
    quietly levelsof `gg', local(cohorts)
    tempname observed
    matrix `observed'=J(1,`G',.)
    foreach c of local ids {
        quietly summarize `gg' if `cc'==`c', meanonly
        local g=r(min)
        matrix `observed'[1,`c']=`g'
        if `g'>0 {
            quietly count if `time'==`g'
            if r(N)==0 exit 459
            quietly count if `cc'==`c' & `time'<`g'
            if r(N)==0 exit 459
            quietly count if `cc'==`c' & `time'>=`g'
            if r(N)==0 exit 459
        }
    }
    quietly count if `gg'>0
    if r(N)==0 exit 459
    local B=0
    if `hasri' {
        confirm matrix `assignments'
        local B=rowsof(`assignments')
        if `B'==0 | colsof(`assignments')!=`G' exit 503
        forvalues i=1/`B' {
            mata: st_local("valid",strofreal(all(sort(st_matrix("`assignments'")[`i',]',1):==sort(st_matrix("`observed'")',1)) & !all(st_matrix("`assignments'")[`i',]:==st_matrix("`observed'"))))
            if "`valid'"!="1" exit 459
            if `i'>1 {
                forvalues j=1/`=`i'-1' {
                    mata: st_local("same",strofreal(all(st_matrix("`assignments'")[`i',]:==st_matrix("`assignments'")[`j',])))
                    if "`same'"=="1" exit 459
                }
            }
        }
    }
    local n : word count `estimators'
    tempname results posts deletions delcol
    matrix `deletions'=J(`G',`n',.)
    matrix colnames `deletions'=`estimators'
    matrix `results'=J(`n',13,.)
    matrix colnames `results'=estimate native_se native_stat native_df native_p jackknife_se jackknife_t jackknife_df jackknife_p ri_beta_p ri_t_native_p N N_clust
    matrix rownames `results'=`estimators'
    tempfile table
    postfile `posts' str12 estimator str8 language str24 status str120 failure_reason ///
        str100 target_label double estimate double native_se double native_p ///
        double jackknife_se double jackknife_p double ri_beta_p ///
        double ri_t_native_p double nobs double nclusters ///
        double ri_reps cells_used cells_eligible min_controls max_controls using "`table'", replace
    local index=0
    foreach id of local estimators {
        local ++index
        local reason ""
        local used .
        local eligible .
        local mincontrols .
        local maxcontrols .
        local status "PASS"
        local target ""
        local opts time(`time') cluster(`cc') gvar(`gg') estimator(`id') datamode(`datamode')
        if "`ivar'"!="" local opts `opts' ivar(`ivar')
        if "`window'"!="" local opts `opts' window(`window')
        if "`id'"=="stacked" & "`stackedwindow'"!="" {
            if "`window'"!="" exit 198
            local opts `opts' window(`stackedwindow')
        }
        if "`id'"=="stacked" & "`stackedcontrol'"!="" local opts `opts' controlgroup(`stackedcontrol')
        if "`id'"=="lpdid" & "`lpdidwindow'"!="" {
            if "`window'"!="" exit 198
            local opts `opts' window(`lpdidwindow')
        }
        if "`id'"=="csdid" {
            if "`csdidmethod'"!="" local opts `opts' method(`csdidmethod')
            if "`csdidcontrol'"!="" local opts `opts' controlgroup(`csdidcontrol')
        }
        if "`id'"=="undid" & "`undidweights'"=="comparison" local opts `opts' cohortfirst
        if "`id'"=="csdid" & "`csdidbackend'"=="original" local opts `opts' originalcs
        if "`id'"=="csdid" & "`csdidweights'"!="" local opts `opts' aggweights(`csdidweights')
        capture quietly _didsuite_eval `varlist' `wopt', `opts'
        local rc=_rc
        if !`rc' & missing(r(estimate)) local rc=498
        if "`targetmode'"!="native" | "`samplemode'"!="native" local rc=198
        if "`window'"!="" & !inlist("`id'","stacked","lpdid") local rc=198
        if `rc' {
            local status "UNSUPPORTED"
            if !inlist(`rc',198,199) local status "FAILURE"
            if `rc'==498 local status "IDENTIFICATION_FAILURE"
            local reason "point: rc=`rc'; no estimate substituted"
            if "`id'"=="lpdid" & "`datamode'"!="panel" local reason "LP-DiD requires an individual panel with repeated outcomes" 
        }
        else {
            local used=r(cells_used)
            local eligible=r(cells_eligible)
            local mincontrols=r(min_controls)
            local maxcontrols=r(max_controls)
            matrix `results'[`index',1]=r(estimate)
            local target "`r(target_label)'"
            foreach method of local inference {
                local rc=0
                if "`method'"=="native" {
                    capture quietly _didsuite_eval `varlist' `wopt', `opts' native
                    local rc=_rc
                    if !`rc' {
                        tempname se df stat
                        scalar `se'=r(native_se)
                        scalar `df'=r(native_df)
                        if missing(scalar(`se')) | scalar(`se')<=0 local rc=498
                        else {
                            scalar `stat'=(`results'[`index',1]-`null')/scalar(`se')
                            matrix `results'[`index',2]=scalar(`se')
                            matrix `results'[`index',3]=scalar(`stat')
                            matrix `results'[`index',4]=scalar(`df')
                            if missing(scalar(`df')) matrix `results'[`index',5]=2*normal(-abs(scalar(`stat')))
                            else matrix `results'[`index',5]=2*ttail(scalar(`df'),abs(scalar(`stat')))
                        }
                    }
                }
                else if "`method'"=="jackknife" {
                    capture quietly _didsuite_jk `varlist' `wopt', `opts'
                    local rc=_rc
                    if !`rc' {
                        matrix `delcol'=r(deletions)
                        matrix `deletions'[1,`index']=`delcol'
                        matrix `results'[`index',6]=r(se)
                        matrix `results'[`index',7]=(`results'[`index',1]-`null')/r(se)
                        matrix `results'[`index',8]=r(df)
                        matrix `results'[`index',9]=2*ttail(r(df),abs(`results'[`index',7]))
                        if missing(`results'[`index',6]) local rc=498
                    }
                }
                else {
                    local column=10
                    if "`method'"=="ri_t_native" local column=11
                    tempname obs den vals
                    scalar `obs'=`results'[`index',1]
                    if "`method'"=="ri_t_native" {
                        capture quietly _didsuite_eval `varlist' `wopt', `opts' native
                        local rc=_rc
                        if !`rc' {
                            scalar `den'=r(native_se)
                            if missing(scalar(`den')) | scalar(`den')<=0 local rc=498
                            else scalar `obs'=scalar(`obs')/scalar(`den')
                        }
                    }
                    matrix `vals'=J(`B',1,.)
                    if !`rc' {
                        forvalues i=1/`B' {
                            foreach c of local ids {
                                quietly replace `gg'=`assignments'[`i',`c'] if `cc'==`c'
                            }
                            local nat
                            if "`method'"=="ri_t_native" local nat "native"
                            capture quietly _didsuite_eval `varlist' `wopt', `opts' `nat'
                            local drawrc=_rc
                            if !`drawrc' {
                                matrix `vals'[`i',1]=r(estimate)
                                if "`method'"=="ri_t_native" {
                                    scalar `den'=r(native_se)
                                    if missing(scalar(`den')) | scalar(`den')<=0 local drawrc=498
                                    else matrix `vals'[`i',1]=`vals'[`i',1]/scalar(`den')
                                }
                            }
                            if `drawrc' local rc=`drawrc'
                        }
                    }
                    foreach c of local ids {
                        quietly replace `gg'=`observed'[1,`c'] if `cc'==`c'
                    }
                    if !`rc' {
                        tempname hits
                        mata: st_numscalar("`hits'",sum(abs(st_matrix("`vals'")):>=abs(st_numscalar("`obs'")):-64*epsilon(1)*rowmax((J(`B',1,1),abs(st_matrix("`vals'")),J(`B',1,abs(st_numscalar("`obs'")))))))
                        matrix `results'[`index',`column']=(1+scalar(`hits'))/(`B'+1)
                    }
                }
                if `rc' {
                    local status "PARTIAL"
                    local reason "`reason' `method':rc=`rc';"
                }
            }
        }
        matrix `results'[`index',12]=_N
        matrix `results'[`index',13]=`G'
        post `posts' ("`id'") ("Stata") ("`status'") ("`reason'") ("`target'") ///
            (`results'[`index',1]) (`results'[`index',2]) (`results'[`index',5]) ///
            (`results'[`index',6]) (`results'[`index',9]) (`results'[`index',10]) ///
            (`results'[`index',11]) (_N) (`G') (`B') (`used') (`eligible') (`mincontrols') (`maxcontrols')
    }
    postclose `posts'
    quietly use "`table'", clear
    generate str12 target_mode="`targetmode'"
    generate byte target_compatible=0
    generate str12 sample_mode="`samplemode'"
    generate str8 data_mode="`datamode'"
    generate str12 csdid_backend=cond(estimator=="csdid","suite","")
    if "`csdidbackend'"=="original" quietly replace csdid_backend="original" if estimator=="csdid"
    generate str64 assignment_hash=""
    generate str40 hash_status="UNSUPPORTED_SHA256"
    generate str80 event_window="`window'"
    if "`stackedwindow'"!="" quietly replace event_window="`stackedwindow'" if estimator=="stacked"
    if "`lpdidwindow'"!="" quietly replace event_window="`lpdidwindow'" if estimator=="lpdid"
    generate str32 inference_cluster="`cluster'"
    generate str32 assignment_unit="`cluster'"
    generate str24 estimator_version="adapter-development-0.1"
    generate double native_stat=el(`results',_n,3)
    generate double native_df=el(`results',_n,4)
    generate double jackknife_t=el(`results',_n,7)
    generate double jackknife_df=el(`results',_n,8)
    generate str80 cohorts_used="`cohorts'"
    generate str64 control_group="native"
    quietly replace control_group="all observations; TWFE coefficient" if estimator=="twfe"
    quietly replace control_group="untreated observations: never and not-yet" if inlist(estimator,"did2s","bjs","flex")
    quietly replace control_group="never" if estimator=="csdid"
    quietly replace control_group="never" if inlist(estimator,"undid","didint")
    if "`csdidcontrol'"!="" quietly replace control_group="`csdidcontrol'" if estimator=="csdid"
    quietly replace control_group="both" if estimator=="stacked"
    if "`stackedcontrol'"!="" quietly replace control_group="`stackedcontrol'" if estimator=="stacked"
    quietly replace control_group="native clean through pooled horizon" if estimator=="lpdid"
    generate str80 covariate_mode="native supplied covariates; no freezing"
    generate str80 weighting_scheme="native; positive analytic weights"
    quietly replace weighting_scheme="post-year treated state counts" if estimator=="undid"
    if "`undidweights'"=="comparison" quietly replace weighting_scheme="two-period weights within cohort, post counts across cells" if estimator=="undid"
    quietly replace weighting_scheme="full-window cohort mass" if estimator=="csdid"
    if "`csdidweights'"=="comparison" quietly replace weighting_scheme="comparison-period cohort mass" if estimator=="csdid"
    generate str64 native_method="estimator-specific; see adapter"
    foreach part in setup point native jk ri total {
        generate double runtime_`part'=.
    }
    generate str40 runtime_status="NOT_INSTRUMENTED"
    generate double assignment_space=.
    generate double n_treated_clusters=.
    if "`saving'"!="" quietly save "`saving'", replace
    if "`csv'"!="" {
        format estimate native_se native_p jackknife_se jackknife_p ri_beta_p ri_t_native_p native_stat native_df jackknife_t jackknife_df %24.17g
        quietly outsheet using "`csv'", comma replace
    }
    di as text _n "DiD Suite: native estimates and clustered inference"
    di as text "{hline 78}"
    di as text %12s "Estimator" " " %9s "Estimate" " " %9s "Native SE" " " %9s "JK SE" " " %9s "Native p" " " %9s "JK p" " " %9s "Status"
    di as text "{hline 78}"
    forvalues row=1/`n' {
        di as text %12s estimator[`row'] " " as result %9.4f estimate[`row'] " " %9.4f native_se[`row'] " " %9.4f jackknife_se[`row'] " " %9.3g native_p[`row'] " " %9.3g jackknife_p[`row'] " " as text %9s status[`row']
    }
    di as text "{hline 78}"
    if `hasri' {
        di as text _n "Randomization p-values (`B' alternative assignments)"
        di as text %12s "Estimator" " " %12s "Coefficient" " " %12s "Native t" " "
        forvalues row=1/`n' {
            di as text %12s estimator[`row'] " " as result %12.4f ri_beta_p[`row'] " " %12.4f ri_t_native_p[`row']
        }
    }
    di as text "JK: delete-cluster jackknife. A dot means unavailable or not requested."
    di as text "Targets can differ across estimators; inspect the exported target labels."
    quietly count if status!="PASS"
    if r(N)>0 {
        di as text _n "Rows needing attention (valid point estimates are retained):"
        forvalues row=1/`n' {
            if status[`row']!="PASS" di as text "  " estimator[`row'] ": " failure_reason[`row']
        }
    }
    restore
    ereturn clear
    if strpos(" `inference' "," jackknife ") ereturn matrix deletions=`deletions'
    ereturn matrix results=`results'
    ereturn matrix observed_assignment=`observed'
    ereturn local cluster_order `"`cluster_order'"'
    ereturn local estimators "`estimators'"
    ereturn local cmd "didsuite"
    ereturn local certification "DEVELOPMENT_LOCAL_TESTS_ONLY"
    ereturn local hash_status "UNSUPPORTED_SHA256"
end

