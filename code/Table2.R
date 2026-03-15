################################################################################
#
#   Filename: Table2.R
#   Purpose:  Generate summary statistics (sample size, mean, standard deviation, 
#             and number of missing values) for each of the four defined cell 
#             clusters in the HSCT dataset.
#   Input data files: data/source/hsct.csv 
#   Output data files: results/Table2.csv 
#   R Version: R-4.4.1
#
################################################################################

hsct = read.csv(paste0(SPATH,"/data/source/hsct.csv"), header = T, stringsAsFactors = FALSE)

# Replace 0 with NA
hsct[hsct==0] = NA
vars    = c("FL1.H","FL2.H","FL3.H","FL4.H")
clusts  = sort(unique(hsct$cls))
# Number of samples per group
ns = setNames(as.integer(table(hsct$cls)), clusts)

# For each cluster and each variable, compute the mean, standard deviation, and count of missing values
res_list = lapply(clusts, function(cl) {
  sub = hsct[hsct$cls == cl, vars]
  t( sapply(sub, function(x) c(
    `sample mean`            = round(mean(x, na.rm=TRUE), 3),
    `sample deviation`       = round(sd(x,   na.rm=TRUE), 3),
    `number of missing values` = sum(is.na(x))
  )) )
})
names(res_list) = paste0("Class ", clusts, " (n=", ns, ")")
res_mat = do.call(cbind, res_list)
write.csv(res_mat, file = paste(SPATH, "/results/Table2.csv",sep=""), row.names = TRUE)
