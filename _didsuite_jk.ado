program define _didsuite_jk, rclass
    version 16
    syntax varlist(min=1 numeric) [aw], TIME(varname numeric) CLuster(varname numeric) ///
        GVAR(varname numeric) ESTimator(string) [IVAR(varname) DATAMode(string) ///
        WINDOW(numlist integer min=2 max=2) METHOD(string) CONTROLGroup(string) AGGWeights(string) COHORTFirst ORIGINALCS]
    local wopt
    if "`weight'"!="" local wopt "[aw`exp']"
    local opts time(`time') cluster(`cluster') gvar(`gvar') estimator(`estimator') datamode(`datamode')
    if "`ivar'"!="" local opts `opts' ivar(`ivar')
    if "`window'"!="" local opts `opts' window(`window')
    if "`method'"!="" local opts `opts' method(`method')
    if "`controlgroup'"!="" local opts `opts' controlgroup(`controlgroup')
    if "`aggweights'"!="" local opts `opts' aggweights(`aggweights')
    local opts `opts' `cohortfirst' `originalcs'
    quietly levelsof `cluster', local(ids)
    local G : word count `ids'
    if `G'<2 exit 498
    tempname dels full scale
    quietly _didsuite_eval `varlist' `wopt', `opts'
    scalar `full'=r(estimate)
    matrix `dels'=J(`G',1,.)
    local i=0
    foreach c of local ids {
        local ++i
        capture quietly _didsuite_eval `varlist' if `cluster'!=`c' `wopt', `opts' deletion
        if !_rc matrix `dels'[`i',1]=r(estimate)
    }
    tempname se
    mata: st_numscalar("`se'",sqrt((rows(st_matrix("`dels'"))-1)/rows(st_matrix("`dels'"))*sum((st_matrix("`dels'"):-mean(st_matrix("`dels'"))):^2)))
    mata: st_numscalar("`scale'",max((1,abs(st_numscalar("`full'")),max(abs(st_matrix("`dels'"))))))
    mata: st_local("incomplete",strofreal(any(st_matrix("`dels'"):>=.)))
    if `incomplete' | missing(scalar(`se')) | scalar(`se')<=64*2.220446049250313e-16*scalar(`scale') scalar `se'=.
    return scalar se=scalar(`se')
    return scalar df=`G'-1
    return matrix deletions=`dels'
end

