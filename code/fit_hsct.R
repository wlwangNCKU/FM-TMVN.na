library(mvtnorm)
library(MomTrunc)
library(mclust)
library(mnormt)
library(VIM)
################################################################################
#
#   Filename: fit_hsct.R
#   Purpose: Fit and compare GMM and FM-TMVN mixture models on the HSCT dataset, 
#            using KNN imputation for missing values and ECM algorithms with g = 4.
#   Input data files: function/GMIX.na.EM.R; 
#                     function/FMTMVN.na.EM.R;
#                     data/source/hsct.csv
#   Output data files: data/hsct.result.Rdata 
#   R Version: R-4.4.1
#   Required R packages: mvtnorm, MomTrunc, mclust, mnormt, VIM 
#
################################################################################

source(paste0(SPATH,"/function/GMIX.na.EM.R"))
source(paste0(SPATH,"/function/FMTMVN.na.EM.R"))

data = read.csv(paste0(SPATH,"/data/source/hsct.csv",T))[,-1]

true.clus = data$cls
Y.na = data[,1:4]
Y.na[Y.na==0]=NA
a = rep(0,4)
b = rep(Inf,4)
p=ncol(Y.na)
g = 4
set.seed(20250319)
Y.knn = kNN(Y.na, k = 5)[,1:p]
fit_GMM=GMIXM.VVV.EM(Y.na, Y.knn, g, init.clus=true.clus, true.clus=true.clus, tol=1e-6, max.iter.EM=3000, per=100)
fit_FMTMVN=MTMVN.na.EM(Y.na, Y.knn, a=a, b=b, g, init.clus=true.clus, true.clus=true.clus, tol=1e-6, max.iter=3000, per=100)
