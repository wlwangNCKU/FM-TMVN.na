library(mvtnorm)
library(MomTrunc)
library(mclust)
library(relliptical)
library(mnormt)
library(VIM)
library(fBasics)
################################################################################
#
#   Filename: simulation.R
#   Purpose: Conduct simulation study for the FM-TMVN model with missing data,
#            estimate parameters using ECM algorithm, and store estimates and
#            standard errors across various sample sizes and missing rates
#   Input data files: function/gener.na.R;
#                     function/Sim_GMIX.na.EM.r.R;
#                     function/Sim_FMTMVN.na.EM.r.R;
#                     function/SE.MTMVN.na.R
#   Output data files: results/simulation/para/(n_dropout).csv (parameter estimates)
#                      results/simulation/sd/(n_dropout).csv   (standard errors)
#                      results/simulation/rdata/n_dropout.RData (rdatas)
#   R Version: R-4.4.1
#   Required R packages: mvtnorm, MomTrunc, mclust, relliptical, mnormt, VIM, 
#                        fBasics
#
################################################################################

source(paste(SPATH, "/function/gener.na.R",sep=""))
source(paste(SPATH, "/function/Sim_GMIX.na.EM.r",sep=""))
source(paste(SPATH, "/function/Sim_FMTMVN.na.EM.r",sep=""))
source(paste(SPATH, "/function/SE.FMTMVN.na.R",sep=""))

simu = function(M, n, na.rate)
{
  ##----true parameters----##
  p=2
  g=2
  a=c(0,0)
  b=c(5,5)
  w = c(0.4, 0.6)
  ni = as.vector(rmultinom(1, size = n, prob = w))
  true.clus = rep(1:g, times = ni)
  mu = matrix(c(1.2, 1.8, 3.8, 3.2), p, g)
  Sigma = array(NA, dim = c(p, p, g))
  Sigma[,,1] = matrix(c(1.4, 0.9, 0.9, 0.8), 2, 2)
  Sigma[,,2] = matrix(c(0.7, -0.6, -0.6, 1.3), 2, 2)
  for(d in 1:M){
    repeat{
      ##----generate FM-TMVN data----##
      Y=NULL
      #seed.num = 500+d
      #set.seed(seed.num)
      for(i in 1:g) Y = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a, upper=b, dist="Normal")))
      #set.seed(seed.num)
      Y.na = gener.na(Y, na.rate)
      na.posi = is.na(Y.na)
      #----------------------------------------SIMULATION------------------------------------------------------------
      fit1=try(GMIXM.VVV.EM(Y.na, Y.knn=Y, g=2, init.clus=true.clus, true.clus=true.clus, tol=1e-6, max.iter.EM=1000, per=100), silent=T)
      fit2=try(MTMVN.na.EM(Y.na, Y.knn=Y, a=a, b=b, g=2, init.clus=true.clus, true.clus=true.clus, tol=1e-6, max.iter=1000, per=100), silent=T)
      sd=try(SE.MTMVN.na(Y.na, a, b, para.est=fit2), silent=T)
      if (!inherits(fit1, "try-error") && !inherits(fit2, "try-error") && !inherits(sd, "try-error")) break
    }
    save(Y, Y.na, fit1, fit2, sd, file = paste0(SPATH,"/results/simulation/rdata/",n,"_",na.rate,"/",d,".RData"))
    paraest=t(sd[,1]);sd1=t(sd[,2])
    write.table(paraest,paste0(SPATH,"/results/simulation/para/(",n,"_",na.rate,").csv"),sep=",",append=T,row.names=F,col.names=F, na = "NA")
    write.table(sd1,paste0(SPATH,"/results/simulation/sd/(",n,"_",na.rate,").csv"),sep=",",append=T,row.names=F,col.names=F, na = "NA")
  }
}
#-------------------------------------------------------------------------
simu(M=500, n=300, na.rate=0)
simu(M=500, n=600, na.rate=0)
simu(M=500, n=1200, na.rate=0)
simu(M=500, n=1800, na.rate=0)


simu(M=500, n=300, na.rate=0.1)
simu(M=500, n=600, na.rate=0.1)
simu(M=500, n=1200, na.rate=0.1)
simu(M=500, n=1800, na.rate=0.1)


simu(M=500, n=300, na.rate=0.2)
simu(M=500, n=600, na.rate=0.2)
simu(M=500, n=1200, na.rate=0.2)
simu(M=500, n=1800, na.rate=0.2)

simu(M=500, n=300, na.rate=0.3)
simu(M=500, n=600, na.rate=0.3)
simu(M=500, n=1200, na.rate=0.3)
simu(M=500, n=1800, na.rate=0.3)
