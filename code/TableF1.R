library(fBasics)
################################################################################
#
#   Filename: TableF1.R
#   Purpose: Generate and export a table of parameter estimates with their  
#            standard errors (in parentheses) for GMM and FM-TMVN models  
#            fitted to the HSCT dataset.
#   Input data files: function/SE.FMTMVN.na.R 
#                     data/hsct_result.RData
#   Output data files: results/TableF1.csv 
#   R Version: R-4.4.1
#   Required R packages: fBasics
#
################################################################################

source(paste0(SPATH,"/function/SE.FMTMVN.na.R"))
load(paste0(SPATH,"/data/hsct_result.RData"))

sd.fit1=SE.MTMVN.na(Y.na, a=c(-Inf,-Inf,-Inf,-Inf), b=c(Inf,Inf,Inf,Inf), para.est=fit_GMM)
sd.fit2=SE.MTMVN.na(Y.na, a=a, b=b, para.est=fit_FMTMVN)
row_idx = c(1:7, 18:21, 32:35, 46:49, 8:17, 22:31, 36:45, 50:59)
table_result = cbind.data.frame(sd.fit1[row_idx,],sd.fit2[row_idx,])
colnames(table_result) = c("est_GMM", "se_GMM", "est_FMTMVN", "se_FMTMVN")
write.csv(table_result, file = paste0(SPATH, "/results/TableS1.csv"), row.names = TRUE)
