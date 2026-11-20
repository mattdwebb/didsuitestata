version 16
clear all
set more off
adopath ++ ".."
adopath ++ "../../fastcsdid/stata"
foreach method in reg dr {
 foreach control in never not_yet {
  foreach weighting in cohort comparison {
   import delimited using csdid_weighting_reference.csv, clear varnames(1) asdouble
   keep if method=="`method'" & control=="`control'" & weighting=="`weighting'"
   sort cluster
   mkmat estimate, matrix(ref)
   scalar refse=native_se[1]
   scalar refjk=jk_se[1]
   import delimited using csdid_weighting_data.csv, clear varnames(1) asdouble
   local option fullwindow
   if "`weighting'"=="comparison" local option comparison
   didsuite y x z [aw=w], time(time) cluster(cluster) gvar(g) datamode(rcs) estimators(csdid) csdidmethod(`method') csdidcontrol(`control') csdidweights(`option')
   matrix res=e(results)
   matrix del=e(deletions)
   assert abs(res[1,1]-ref[1,1])<1e-8
   assert abs(res[1,2]-refse)<1e-8
   assert abs(res[1,6]-refjk)<1e-8
   assert res[1,8]==11
   forvalues i=1/12 {
    assert abs(del[`i',1]-ref[`i'+1,1])<1e-8
   }
  }
 }
}
display "CSDID_ALL_WEIGHTS_CONTROLS_SURVIVING_COHORTS_PARITY_PASS"
