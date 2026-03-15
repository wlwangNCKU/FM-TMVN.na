################################################################################
#
#   Filename: Table3.R
#   Purpose: Generate and export a cross-tabulation of true versus estimated 
#            cluster labels for GMM and FM-TMVN models on the HSCT dataset, 
#            and compute ARI and CCR for each method.
#   Input data files: data/hsct_result.RData
#   Output data files: results/Table3.csv 
#   R Version: R-4.4.1
#
################################################################################

load(paste0(SPATH,"/data/hsct_result.RData"))

table1 = table(true.clus, fit_GMM$post.clus)
table2 = table(true.clus, fit_FMTMVN$post.clus)
table_result = cbind(table1, table2)
ac = matrix(c(round(fit_GMM$ARI,3), round(fit_GMM$CCR,3), rep("",6), round(fit_FMTMVN$ARI,3), round(fit_FMTMVN$CCR,3), rep("",6)), nrow = 2)
table_result = rbind(table_result,ac)
rownames(table_result) = c(1,2,3,4,"ARI","CCR")
write.csv(table_result, file = paste0(SPATH, "/results/Table3.csv"), row.names = TRUE)
