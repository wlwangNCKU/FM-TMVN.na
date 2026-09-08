library(mvtnorm)
library(MomTrunc)
library(mclust)
library(mnormt)
library(VIM)
################################################################################
#
#   Filename: fit_stone.R
#   Purpose: Fit and compare GMM and FM-TMVN mixture models on the Stone Flakes 
#            dataset for cluster counts g = 1~6, using KNN imputation for missing  
#            values and ECM algorithms.
#   Input data files: function/gener.na.R
#                     function/GMIX.na.EM.R; 
#                     function/FMTMVN.na.EM.R;
#                     data/source/StoneFlakes.csv
#   Output data files: data/stone.result.Rdata 
#   R Version: R-4.4.1
#   Required R packages: mvtnorm, MomTrunc, mclust, mnormt, VIM
#
################################################################################
SPATH = "C:/研究/FM-TMVN/Data_and_code_FMTMVN.na"
source(paste0(SPATH,"/function/gener.na.R"))
source(paste0(SPATH,"/function/GMIX.na.EM.R"))
source(paste0(SPATH,"/function/FM-TMVN.na.EM.R"))

Y.na = as.matrix(read.csv(paste0(SPATH,"/data/source/StoneFlakes.csv"))[,c(1,4,5,6,7)])
sc.y = scale(Y.na)
stonefit = function(Y, g, seed){
  p = ncol(Y)
  a=apply(Y,2,min,na.rm=T)
  b=apply(Y,2,max,na.rm=T)
  set.seed(seed)
  Y.knn = kNN(Y, k = 5)[,1:p]
  set.seed(seed)
  km.clus = kmeans(Y.knn, g)$cluster
  fit_GMM = GMIXM.VVV.EM(Y.na = Y, Y.knn=Y.knn, g=g, init.clus=km.clus , true.clus=F , tol=1e-6 , max.iter.EM=3000,per=100)
  fit_FMTMVN = MTMVN.na.EM(Y.na = Y, Y.knn=Y.knn, a=a, b=b, g=g, init.clus=km.clus, true.clus=F, tol=1e-6, max.iter=3000, per=100)
  return(list(fit_GMM = fit_GMM, fit_FMTMVN = fit_FMTMVN))
}
s_1 = stonefit(Y = sc.y, g=1, seed=20250526)
s_2 = stonefit(Y = sc.y, g=2, seed=20250526)
s_3 = stonefit(Y = sc.y, g=3, seed=20250526)
s_4 = stonefit(Y = sc.y, g=4, seed=20250526)
s_5 = stonefit(Y = sc.y, g=5, seed=20250526)
s_6 = stonefit(Y = sc.y, g=6, seed=20250526)
total = list()
for (g in 1:6) {
  total[[paste0(g,"_1")]] = get(paste0("s_",g))$fit_GMM
  total[[paste0(g,"_2")]] = get(paste0("s_",g))$fit_FMTMVN
}