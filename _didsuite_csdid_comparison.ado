program define _didsuite_csdid_comparison, rclass
    version 16
    syntax varlist(numeric min=1) [if] [in], TIME(varname numeric) GVAR(varname numeric) CLuster(varname numeric) WEIGHTS(varname numeric) METHOD(string) CONTROL(string) [CELLS(name) AGGWeights(string) DELETION]
    marksample use
    gettoken y x : varlist
    if "`method'"=="dripw" local method dr
    if "`method'"=="ra" local method reg
    if "`method'"=="stdipw" local method ipw
    tempname result
    mata: _didsuite_comparison_if("`y' `time' `gvar' `cluster' `weights' `x'","`use'","`method'","`control'","`cells'","`result'","`aggweights'","`deletion'")
    return scalar estimate=`result'[1,1]
    return scalar native_se=`result'[1,2]
    return scalar native_df=.
    return scalar cells_used=`result'[1,3]
    return scalar cells_eligible=`result'[1,4]
    return local target_label "CSDID RCS comparison-period treated-mass ATT"
    return local native_method "CSDID cluster IF normal; estimated aggregation weights"
end
mata:
// Exact Fast CSDID 1.0.1 cell kernel, MIT, Matthew D. Webb.
// The standalone ado scopes its Mata helpers internally; this adapter needs
// cell influence vectors to reaggregate with comparison-period mass weights.
// Upstream source SHA256: fe2060d9ba998336a40a1c9e089cea2d0b27b79abb130906ff7e0057fc8334f7
// Only the function prefix changes. Matched R/Python/Stata tests cover native
// variances and every deletion, including deletion of singleton cohorts.
real scalar _ds_cs_ratio(real colvector w, real colvector y)
{
    if (sum(w)<=0) return(.)
    return(sum(w:*y)/sum(w))
}
real colvector _ds_cs_rif(real colvector w, real colvector y)
{
    return(w:*(y:-_ds_cs_ratio(w,y))/mean(w))
}
real scalar _ds_cs_wls(real matrix X, real colvector y, real colvector w,
    real colvector m, real matrix IF)
{
    real matrix H
    real colvector b
    H = quadcross(X,w,X)
    if (rank(H)<cols(X)) return(0)
    b = qrsolve(H,quadcross(X,w,y))
    m = X*b
    IF = qrsolve(H/rows(X),(X:*(w:*(y-m)))')'
    return(1)
}
real scalar _ds_cs_logit(real matrix X, real colvector D, real colvector w,
    real colvector p, real matrix IF)
{
    real colvector b, step, eta, bn
    real matrix H
    real scalar iter, scale, loss, newer
    b = J(cols(X),1,0)
    for(iter=1;iter<=100;iter++) {
        p = invlogit(X*b)
        H = quadcross(X,w:*p:*(1:-p),X)
        if(rank(H)<cols(X)) return(0)
        step = qrsolve(H,quadcross(X,w,D-p))
        eta = X*b
        loss = sum(w:*(ln(1:+exp(eta)):-D:*eta))
        scale = 1
        while(scale>2^-25) {
            bn = b+scale*step
            eta = X*bn
            newer = sum(w:*(ln(1:+exp(eta)):-D:*eta))
            // The summed objective grows with N. Allow only floating-point
            // roundoff in its comparison; keep the Newton-step convergence
            // threshold unchanged.
            if(newer<=loss+64*2^-52*max((1,abs(loss),abs(newer)))) break
            scale = scale/2
        }
        b = bn
        if(max(abs(step))<1e-10) break
        if(max(abs(b))>40) return(0)
    }
    if(iter>100) return(0)
    p = rowmin((invlogit(X*b),J(rows(X),1,1-1e-6)))
    H = quadcross(X,w:*p:*(1:-p),X)
    if(rank(H)<cols(X)) return(0)
    IF = qrsolve(H/rows(X),(X:*(w:*(D-p)))')'
    return(1)
}
real colvector _ds_cs_cell(real colvector y, real colvector P, real colvector D,
    real matrix X, real colvector w, string scalar method)
{
    real colvector w10,w11,w00,w01,wd,p,trim,z,m00,m01,m10,m11,inf,mp,mb
    real matrix o00,o01,o10,o11,pi
    real scalar att,r10,r11,r00,r01
    if(rows(y)==0 | sum(w)<=0) return(J(rows(y)+1,1,.))
    w=w/mean(w)
    w10=w:*D:*(1:-P); w11=w:*D:*P
    w00=w:*(1:-D):*(1:-P); w01=w:*(1:-D):*P
    if(min((sum(w10),sum(w11),sum(w00),sum(w01)))<=0) return(J(rows(y)+1,1,.))
    if(method=="none") {
        att=_ds_cs_ratio(w11,y)-_ds_cs_ratio(w10,y)-_ds_cs_ratio(w01,y)+_ds_cs_ratio(w00,y)
        inf=_ds_cs_rif(w11,y)-_ds_cs_rif(w10,y)-_ds_cs_rif(w01,y)+_ds_cs_rif(w00,y)
        return(att\inf)
    }
    m00=m01=m10=m11=J(rows(y),1,.)
    o00=o01=o10=o11=pi=J(rows(y),cols(X),.)
    if(method=="reg" | method=="dr") {
        if(!_ds_cs_wls(X,y,w00,m00,o00) | !_ds_cs_wls(X,y,w01,m01,o01)) return(J(rows(y)+1,1,.))
    }
    if(method=="reg") {
        wd=w:*D
        att=_ds_cs_ratio(w11,y)-_ds_cs_ratio(w10,y)-_ds_cs_ratio(wd,m01-m00)
        inf=_ds_cs_rif(w11,y)-_ds_cs_rif(w10,y)-_ds_cs_rif(wd,m01-m00)-(o01-o00)*(X'*wd/sum(wd))
        return(att\inf)
    }
    p=J(rows(y),1,.)
    if(!_ds_cs_logit(X,D,w,p,pi)) return(J(rows(y)+1,1,.))
    trim=(D:==1):|(p:<.995)
    w10=trim:*w10; w11=trim:*w11
    w00=trim:*p:/(1:-p):*w00; w01=trim:*p:/(1:-p):*w01
    if(min((sum(w00),sum(w01)))<=0) return(J(rows(y)+1,1,.))
    z=y
    if(method=="dr") z=y-P:*m01-(1:-P):*m00
    r10=_ds_cs_ratio(w10,z); r11=_ds_cs_ratio(w11,z)
    r00=_ds_cs_ratio(w00,z); r01=_ds_cs_ratio(w01,z)
    att=r11-r10-r01+r00
    inf=_ds_cs_rif(w11,z)-_ds_cs_rif(w10,z)-_ds_cs_rif(w01,z)+_ds_cs_rif(w00,z) -
        pi*(X'*(w01:*(z:-r01))/sum(w01)-X'*(w00:*(z:-r00))/sum(w00))
    if(method=="dr") {
        if(!_ds_cs_wls(X,y,w:*D:*(1:-P),m10,o10) | !_ds_cs_wls(X,y,w:*D:*P,m11,o11)) return(J(rows(y)+1,1,.))
        wd=trim:*w:*D
        att=att+_ds_cs_ratio(wd,m11-m01)-_ds_cs_ratio(w11,m11-m01)-_ds_cs_ratio(wd,m10-m00)+_ds_cs_ratio(w10,m10-m00)
        mp=X'*(wd/sum(wd)-w11/sum(w11))
        mb=X'*(wd/sum(wd)-w10/sum(w10))
        // beta_00 loading is treated-pre minus weighted-control-pre.
        inf=inf+o01*(X'*(w01/sum(w01)-w11/sum(w11)))+o00*(X'*(w10/sum(w10)-w00/sum(w00))) +
            _ds_cs_rif(wd,m11-m01)-_ds_cs_rif(w11,m11-m01)-_ds_cs_rif(wd,m10-m00)+_ds_cs_rif(w10,m10-m00) +
            (o11-o01)*mp-(o10-o00)*mb
    }
    return(att\inf)
}

void _didsuite_comparison_if(string scalar vars,string scalar use,string scalar method,string scalar control,string scalar cellname,string scalar resultname,string scalar weighting,string scalar deletion)
{
    real matrix Z,M,X,I,mass,CI
    real colvector y,t,g,c,w,idx,v,raw,active,alpha,score,ids,times,groups,pre,post
    real scalar N,K,k,theta,se,h,j,q
    Z=st_data(.,tokens(vars),use)
    y=Z[,1];t=Z[,2];g=Z[,3];c=Z[,4];w=Z[,5]
    if(cellname!="") M=st_matrix(cellname)
    else {
        M=J(0,5,.)
        times=uniqrows(sort(t,1));groups=uniqrows(sort(select(g,g:>0),1))
        for(j=1;j<=rows(groups);j++) {
            pre=select(times,times:<groups[j]);post=select(times,times:>=groups[j])
            if(!rows(pre)) continue
            for(q=1;q<=rows(post);q++) M=M\(groups[j],max(pre),post[q],post[q]-groups[j],.)
        }
    }
    N=rows(Z);K=rows(M);X=J(N,1,1)
    if(!K) _error(498)
    if(cols(Z)>5) X=(X,Z[,6..cols(Z)])
    I=J(N,K,.);mass=J(N,K,0)
    for(k=1;k<=K;k++) {
        idx=selectindex(((g:==M[k,1]):|(g:==0):|((control=="not_yet"):&(g:>M[k,3]))):&((t:==M[k,2]):|(t:==M[k,3])))
        v=_ds_cs_cell(y[idx],t[idx]:==M[k,3],g[idx]:==M[k,1],X[idx,],w[idx],method)
        if(cellname=="") M[k,5]=v[1]
        if(v[1]<.) {
            I[,k]=J(N,1,0);I[idx,k]=N/rows(idx)*v[2..rows(v)]
        }
        mass[,k]=w:*(g:==M[k,1])
        if(weighting!="fullwindow" & !(weighting=="" & deletion!="")) mass[,k]=mass[,k]:*((t:==M[k,2]):|(t:==M[k,3]))
    }
    raw=colsum(mass)';active=selectindex((raw:>0):&(M[,5]:<.))
    if(!rows(active) | (deletion=="" & rows(active)!=K)) _error(498)
    alpha=raw[active]/sum(raw[active]);theta=alpha'*M[active,5]
    score=I[,active]*alpha+N*mass[,active]*(M[active,5]:-theta)/sum(raw[active])
    ids=uniqrows(sort(c,1));CI=J(rows(ids),1,.)
    for(h=1;h<=rows(ids);h++) CI[h]=sum(select(score,c:==ids[h]))
    se=sqrt(sum(CI:^2))/N
    st_matrix(resultname,(theta,se,rows(active),K))
}
end
