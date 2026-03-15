library(mvtnorm)
library(MomTrunc)
library(mclust)
library(relliptical)
library(mnormt)
library(VIM)
library(fBasics)
SPATH = "C:/研究/FM-TMVN/Data_and_code_FMTMVN.na"
source(paste(SPATH, "/function/gener.na.R",sep=""))
source(paste(SPATH, "/function/GMIX.na.EM.R",sep=""))
source(paste(SPATH, "/function/Sim2_FMTMVN.na.EM.R",sep=""))
source(paste(SPATH, "/function/Sim3_FMTMVN.na.EM.R",sep=""))

simu_g=function(n)
{
  cc=sqrt(2/pi)
  p = 5; q = 2; g = 3
  mu=matrix(c(1.47,-0.04,0.96,1.01,-0.14,0.34,0.43,-0.70,-0.36,0.11,-1.65, 1.65, 0.59, 0.02, -0.57),p)
  B =array(NA,dim=c(p,q,g))
  B[,,1] = matrix(c(0.05,0.87,0.12,0.69,0.01,0.46,0.07,0.02,0.59,0.66),p,q)
  B[,,2] = matrix(c(0.47,0.35,0.03,0.67,0.50,0.04,0.95,0.48,0.84,0.47),p,q)
  B[,,3] = matrix(c(0.84,0.89,0.35,0.05,0.99,0.05,0.94,0.68,0.39,0.52),p,q)
  D = cbind(rep(0.25,p),rep(0.5,p),rep(0.75,p))
  Sigma =array(NA,dim=c(p,p,g))
  ni = c(rmultinom(n=1, size=n, prob=rep(1/g,g)))
  a=c(0,0,0,0,0)
  b=c(3,3,3,3,3)
  Y=NULL
  for(i in 1:g)
  {
    Sigma[,,i]= B[,,1] %*% t(B[,,1])+ diag(D[,i])
    Y = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a, upper=b, dist="Normal")))
    true.clus=rep(1:g, times=ni)
  }
  list(Y=Y, true.clus=true.clus, mu=mu, B=B, D=D)
}



simu = function(M, n, na.rate)
{
  for(d in 1:M){
    repeat{
      ##----generate FM-TMVN data----##
      Data=simu_g(n)
      Y=Data$Y
      true.clus=Data$true.clus
      a=c(0,0,0,0,0)
      b=c(3,3,3,3,3)
      Y.na = gener.na(Y, na.rate)
      na.posi = is.na(Y.na)
      #----------------------------------------SIMULATION------------------------------------------------------------
      fit1=try(GMIXM.VVV.EM(Y.na, Y.knn=Y, g=3, init.clus=true.clus, true.clus=true.clus, tol=1e-6, max.iter.EM=3000, per=100), silent=T)
      fit2=try(MTMVN.na.EMup(Y.na, Y.knn=Y, a=a, b=b, g=3, init.clus=true.clus, true.clus=true.clus, tol=1e-6, max.iter=3000, per=100), silent=T)
      fit3=try(MTMVN.na.EM(Y.na, Y.knn=Y, a=a, b=b, g=3, init.clus=true.clus, true.clus=true.clus, tol=1e-6, max.iter=3000, per=100), silent=T)
      if (!inherits(fit1, "try-error") && !inherits(fit2, "try-error") && !inherits(fit3, "try-error")) break
    }
    save(Y, Y.na, na.posi, fit1, fit2, file = paste0(SPATH,"/test_input4/up/rdata/",n,"_",na.rate,"/",d,".RData"))
    save(Y, Y.na, na.posi, fit1, fit3, file = paste0(SPATH,"/test_input4/rdata/",n,"_",na.rate,"/",d,".RData"))
    y_true = Y[na.posi]
    GMM_input = mean((y_true - fit1$post.pred[na.posi])^2)
    FMTMVN_input = mean((y_true - fit2$post.pred[na.posi])^2)
    result = c(fit1$BIC, fit2$BIC, fit1$ARI, fit2$ARI, fit1$CCR, fit2$CCR, GMM_input, FMTMVN_input, fit1$sec.time, fit2$sec.time)
    write.table(t(result),paste0(SPATH,"/test_input4/up/vs/(",n,"_",na.rate,").csv"),sep=",",append=T,row.names=F,col.names=F, na = "NA")
    FMTMVN1_input = mean((y_true - fit3$post.pred[na.posi])^2)
    result1 = c(fit1$BIC, fit3$BIC, fit1$ARI, fit3$ARI, fit1$CCR, fit3$CCR, GMM_input, FMTMVN1_input, fit1$sec.time, fit3$sec.time)
    write.table(t(result1),paste0(SPATH,"/test_input4/vs/(",n,"_",na.rate,").csv"),sep=",",append=T,row.names=F,col.names=F, na = "NA")
  }
}
simu(M=100, n=300, na.rate=0.1)
#simu(M=10, n=600, na.rate=0.1)
simu(M=100, n=900, na.rate=0.1)
#simu(M=10, n=1200, na.rate=0.1)
#simu(M=10, n=1500, na.rate=0.1)
simu(M=100, n=1800, na.rate=0.1)


simu(M=100, n=300, na.rate=0.2)
#simu(M=10, n=600, na.rate=0.2)
simu(M=100, n=900, na.rate=0.2)
#simu(M=10, n=1200, na.rate=0.2)
#simu(M=10, n=1500, na.rate=0.2)
simu(M=100, n=1800, na.rate=0.2)

simu(M=100, n=300, na.rate=0.3)
#simu(M=10, n=600, na.rate=0.3)
simu(M=100, n=900, na.rate=0.3)
#simu(M=10, n=1200, na.rate=0.3)
#simu(M=10, n=1500, na.rate=0.3)
simu(M=100, n=1800, na.rate=0.3)

for (i in 1:3){
  position = switch(i,"300_0.1","600_0.1","900_0.1")
  for (j in 1:10) {
    load(paste0(SPATH,"/test_input3/rdata/",position,"/",j,".RData"))
    result = c(fit1$BIC, fit2$BIC, fit1$ARI, fit2$ARI, fit1$ICL, fit2$ICL, fit1$CCR, fit2$CCR)
    write.table(t(result),paste0(SPATH,"/test_input3/(",position,")1.csv"),sep=",",append=T,row.names=F,col.names=F, na = "NA")
  }
}
