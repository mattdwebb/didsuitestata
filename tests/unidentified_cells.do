version 16
clear all
set more off
adopath ++ ".."
adopath ++ "../../fastcsdid/stata"
tempfile summary
tempname records
postfile `records' str3 method str10 weighting double cluster double deletion double att double se using "`summary'", replace
foreach method in reg dr {
 foreach weighting in fullwindow comparison {
  import delimited using unidentified_cell.csv, clear varnames(1) asdouble
  didsuite y x z [aw=w], time(time) cluster(cluster) gvar(g) datamode(rcs) estimators(csdid) csdidmethod(`method') csdidweights(`weighting') inference(point)
  matrix bad=e(results)
  assert missing(bad[1,1])
  tempfile table
  didsuite y x z [aw=w], time(time) cluster(cluster) gvar(g) datamode(rcs) estimators(csdid) csdidmethod(`method') csdidweights(`weighting') inference(point) saving("`table'")
  preserve
  use "`table'", clear
  assert status=="IDENTIFICATION_FAILURE"
  restore
  import delimited using partial_deletion.csv, clear varnames(1) asdouble
  didsuite y x z [aw=w], time(time) cluster(cluster) gvar(g) datamode(rcs) estimators(csdid) csdidmethod(`method') csdidweights(`weighting')
  matrix good=e(results)
  matrix del=e(deletions)
  assert !missing(good[1,1],good[1,6])
  assert good[1,8]==12
  assert rowsof(del)==13
  forvalues j=1/13 {
   assert !missing(del[`j',1])
   post `records' ("`method'") ("`weighting'") (`j') (del[`j',1]) (good[1,1]) (good[1,6])
  }
  if "`weighting'"=="fullwindow" {
   fastcsdid y, time(time) cluster(cluster) gvar(g) weights(w) covariates(x z) method(`method') jackknife
   matrix standalone=e(deleted)
   forvalues j=1/13 {
    assert abs(standalone[`j',1]-del[`j',1])<1e-8
   }
  }
 }
}
postclose `records'
use "`summary'", clear
format cluster deletion att se %24.17g
export delimited using unidentified_results.csv, replace datafmt
display "FULL_SAMPLE_FAILURE_AND_SURVIVING_DELETIONS_PASS"
