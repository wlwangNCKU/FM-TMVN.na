################################################################################
#
#   Filename: Table1.R
#   Purpose:    Compute and tabulate sample standard deviations (STD) and 
#               integrated mean squared errors (IMSE) of MLE parameter estimates 
#               from simulation under varying missing rates (0%, 10%, 20%, 30%) 
#               and sample sizes (n = 300, 600, 1200, 1800).
#   Input data files: results/simulation/para/*.csv (estimated parameters)
#                     results/simulation/sd/*.csv   (estimated standard errors)
#   Output data files: results/Table1.csv 
#   R Version: R-4.4.1
#
################################################################################

rates = c("0%", "10%", "20%", "30%")
params = c(
  "pi_1",
  "mu_11","mu_12","mu_21","mu_22",
  "Sigma_111","Sigma_121","Sigma_122",
  "Sigma_211","Sigma_221","Sigma_222"
)
table_result = data.frame(
  missing_rate = rep(rates, each = length(params)),
  parameter    = rep(params, times = length(rates)),
  stringsAsFactors = FALSE
)
std_re = imse_re = NA
for (i in 1:16){
  position = switch(i,"(300_0)","(300_0.1)","(300_0.2)","(300_0.3)","(600_0)","(600_0.1)","(600_0.2)","(600_0.3)","(1200_0)","(1200_0.1)","(1200_0.2)","(1200_0.3)","(1800_0)","(1800_0.1)","(1800_0.2)","(1800_0.3)")
  col_idx = c(1:3, 7:8, 4:6, 9:11)
  sd = read.csv(paste0(SPATH,"/results/simulation/sd/",position,".csv"),header = F)[ , col_idx]
  paraest = read.csv(paste0(SPATH,"/results/simulation/para/",position,".csv"),header = F)[ , col_idx]
  std_re = cbind(std_re, t(round(apply(paraest, 2, sd),3)))
  imse_re = cbind(imse_re, t(round(colMeans(sd),3)))
}
table_result$STD_300 = std_re[2:45]; table_result$IMSE_300 = imse_re[2:45]
table_result$STD_600 = std_re[46:89]; table_result$IMSE_600 = imse_re[46:89]
table_result$STD_1200 = std_re[90:133]; table_result$IMSE_1200 = imse_re[90:133]
table_result$STD_1800 = std_re[134:177]; table_result$IMSE_1800 = imse_re[134:177]
write.csv(table_result, file = paste0(SPATH, "/results/Table1.csv"), row.names = TRUE)
