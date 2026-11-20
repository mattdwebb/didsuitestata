program define _didsuite_eval, rclass
    version 16
    syntax varlist(min=1 numeric) [if] [in] [aw], TIME(varname numeric) ///
        CLuster(varname numeric) GVAR(varname numeric) ESTimator(string) ///
        [IVAR(varname) DATAMode(string) WINDOW(numlist integer min=2 max=2) ///
        METHOD(string) CONTROLGroup(string) COHORTFirst AGGWeights(string) ORIGINALCS NATIVE DELETION]
    marksample use, novarlist
    gettoken y x : varlist
    local wopt
    if "`weight'"!="" local wopt "[aw`exp']"
    if inlist("`estimator'","twfe","did2s","bjs","flex") {
        capture which fast`estimator'
        if _rc exit 199
        local mode "point"
        if "`native'"!="" local mode "native"
        local unitopt
        if "`ivar'"!="" local unitopt "ivar(`ivar')"
        quietly fast`estimator' `varlist' if `use' `wopt', time(`time') ///
            cluster(`cluster') gvar(`gvar') inference(`mode') datamode(`datamode') `unitopt'
        return scalar estimate=e(estimate)
        return scalar native_se=e(native_se)
        return scalar native_df=e(native_df)
        return local target_label "`e(target_label)'"
        return local native_method "fast`estimator' native port"
        exit
    }
    if inlist("`estimator'","undid","didint") {
        if "`native'"!="" exit 198
        capture which suite_`estimator'
        if _rc exit 199
        preserve
        quietly keep if `use'
        tempvar ww
        quietly gen double `ww'=1
        if "`weight'"!="" quietly replace `ww' `exp'
        quietly count if `ww'!=1
        if r(N)>0 exit 198
        local unitopt
        if "`ivar'"!="" local unitopt "ivar(`ivar')"
        quietly suite_`estimator' `varlist', time(`time') cluster(`cluster') ///
            gvar(`gvar') datamode(`datamode') inference(point) `unitopt' `cohortfirst'
        tempname b cells minc maxc
        matrix `cells'=e(att_st)
        local used=rowsof(`cells')
        local eligible=e(n_cells_eligible)
        mata: st_numscalar("`minc'",min(st_matrix("`cells'")[,8]))
        mata: st_numscalar("`maxc'",max(st_matrix("`cells'")[,8]))
        scalar `b'=e(estimate)
        local target "`e(target_label)'"
        restore
        return scalar estimate=scalar(`b')
        return scalar native_se=.
        return scalar native_df=.
        return local target_label "`target'"
        return scalar cells_used=`used'
        return scalar cells_eligible=`eligible'
        return scalar min_controls=scalar(`minc')
        return scalar max_controls=scalar(`maxc')
        return local native_method "UNSUPPORTED native analytic VCE"
        exit
    }
    if "`estimator'"=="csdid" {
        if "`datamode'"!="rcs" exit 198
        if "`originalcs'"!="" {
            if !inlist("`controlgroup'","","never","not_yet") exit 198
            if "`method'"=="" local method dripw
            if "`method'"=="none" & "`x'"!="" exit 198
            if inlist("`method'","none","ra") local method reg
            if "`method'"=="dr" local method dripw
            if "`method'"=="ipw" local method stdipw
            local controls nevertreated
            if "`controlgroup'"=="not_yet" local controls notyet
            local iw
            if "`weight'"!="" local iw "[iw`exp']"
            quietly csdid `varlist' if `use' `iw', time(`time') gvar(`gvar') ///
                cluster(`cluster') method(`method') `controls' rcs base_period(universal) analytical pointwise
            quietly estat simple, post
            if missing(_b[ATT]) exit 498
            return scalar estimate=_b[ATT]
            return scalar native_se=_se[ATT]
            return scalar native_df=.
            return local target_label "original csdid simple ATT (comparison-period weights)"
            exit
        }
        capture which fastcsdid
        if _rc exit 199
        if "`method'"=="" local method "dripw"
        if "`controlgroup'"=="" local controlgroup "never"
        tempvar w
        quietly gen double `w'=1 if `use'
        if "`weight'"!="" quietly replace `w' `exp' if `use'
        local covopt
        if "`x'"!="" local covopt "covariates(`x')"
        if "`deletion'"!="" {
            local delweight fullwindow
            if "`aggweights'"=="comparison" local delweight comparison
            quietly _didsuite_csdid_comparison `y' `x' if `use', time(`time') gvar(`gvar') cluster(`cluster') weights(`w') method(`method') control(`controlgroup') aggweights(`delweight') deletion
            return add
            exit
        }
        quietly fastcsdid `y' if `use', time(`time') cluster(`cluster') ///
            gvar(`gvar') weights(`w') method(`method') control(`controlgroup') `covopt'
        tempname cells
        matrix `cells'=e(attgt)
        if "`aggweights'"=="comparison" {
            quietly _didsuite_csdid_comparison `y' `x' if `use', time(`time') gvar(`gvar') cluster(`cluster') weights(`w') method(`method') control(`controlgroup') cells(`cells')
            return add
            exit
        }
        return scalar cells_used=rowsof(`cells')
        return scalar cells_eligible=rowsof(`cells')
        return scalar estimate=_b[ATT]
        return scalar native_se=e(native_se)
        return scalar native_df=.
        return local target_label "CSDID RCS `method' cohort-mass ATT"
        return local native_method "CSDID cluster IF normal"
        exit
    }
    if "`estimator'"=="lpdid" {
        if "`datamode'"!="panel" | "`ivar'"=="" | "`window'"=="" exit 198
        capture which fastlpdid
        if _rc exit 199
        preserve
        quietly keep if `use'
        tempvar uu ww
        quietly egen long `uu'=group(`ivar')
        quietly gen double `ww'=1
        if "`weight'"!="" quietly replace `ww' `exp'
        quietly count if `ww'!=1
        if r(N)>0 exit 198
        local covopt
        if "`x'"!="" local covopt "covariates(`x')"
        local mode "point"
        if "`native'"!="" local mode "native"
        quietly fastlpdid `y', unit(`uu') time(`time') cluster(`cluster') ///
            gvar(`gvar') window(`window') datamode(panel) inference(`mode') `covopt'
        tempname b se df
        scalar `b'=e(estimate)
        scalar `se'=e(native_se)
        scalar `df'=e(native_df)
        local target "`e(target_label)'"
        restore
        return scalar estimate=scalar(`b')
        return scalar native_se=scalar(`se')
        return scalar native_df=scalar(`df')
        return local target_label "`target'"
        return local native_method "LP-DiD native CR1 or RA CR0"
        exit
    }
    if "`estimator'"=="stacked" {
        capture which faststacked
        if _rc exit 199
        if "`window'"=="" exit 198
        preserve
        quietly keep if `use'
        tempvar cohort w
        quietly gen double `cohort'=cond(`gvar'==0,.,`gvar')
        local unitopt
        if "`ivar'"!="" local unitopt "unit(`ivar')"
        local covopt
        if "`x'"!="" local covopt "covariates(`x')"
        local weightopt
        if "`weight'"!="" {
            quietly gen double `w' `exp'
            local weightopt "weightvar(`w')"
        }
        local mode "point"
        if "`native'"!="" local mode "native"
        if "`controlgroup'"=="" local controlgroup "both"
        quietly faststacked `y', time(`time') cluster(`cluster') gvar(`cohort') ///
            window(`window') datamode(`datamode') `unitopt' `covopt' `weightopt' ///
            controlgroup(`controlgroup') inference(`mode')
        tempname b se df
        scalar `b'=r(estimate)
        scalar `se'=r(native_se)
        scalar `df'=r(df)
        local target "`r(target_label)'"
        restore
        return scalar estimate=scalar(`b')
        return scalar native_se=scalar(`se')
        return scalar native_df=scalar(`df')
        return local target_label "`target'"
        return local native_method "Stacked original-cluster CR1"
        exit
    }
    exit 199
end

