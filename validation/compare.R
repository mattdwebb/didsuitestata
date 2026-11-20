# Development-only cross-language export comparison, never called by ado code.
spec <- if(length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "../../didsuite_spec"
s <- read.csv("master_stata.csv",na.strings=c("","."))
p <- read.csv(file.path(spec,"validation","results","master_python.csv"))
fields <- c("estimate","native_se","jackknife_se","ri_beta_p","ri_t_native_p")
out <- do.call(rbind,lapply(fields,function(v) {
  rhs <- p[[v]][match(s$estimator,p$estimator)]
  lhs <- s[[v]]
  data.frame(estimator=s$estimator,field=v,Stata=lhs,Python=rhs,
    absolute_gap=abs(lhs-rhs),status=ifelse(is.na(lhs)&is.na(rhs),"BOTH_NA",
      ifelse(is.na(lhs)|is.na(rhs),"AVAILABILITY_DIFF",ifelse(abs(lhs-rhs)<1e-8,"PASS","DIFFERENCE"))))
}))
write.csv(out,"master_python_comparison.csv",row.names=FALSE,na="")
print(out[out$status %in% c("DIFFERENCE","AVAILABILITY_DIFF"),])
cat("MAX ABSOLUTE GAP",max(out$absolute_gap,na.rm=TRUE),"\n")
stopifnot(!any(out$status=="DIFFERENCE"))
cat("DIDSUITE_STATA_EXPORT_COMPARISON_PASS\n")
