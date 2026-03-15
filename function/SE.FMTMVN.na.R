# Standard deviation calculation for the FM-TMVN model
SE.MTMVN.na=function(Y.na, a, b, para.est)
{
  n = nrow(Y.na)
  p = ncol(Y.na)
  w = para.est$para$w
  mu = para.est$para$mu
  S = para.est$para$Sigma
  g = length(w)
  na.posi=is.na(Y.na) # Positions of missing values
  po=p-rowSums(na.posi) # Number of observed variables in each row
  ind.na=colSums(t(na.posi)*2^(1:p-1)) # Pattern of missing data
  num.class.na=length(unique(ind.na)) # Number of missingness patterns
  row.posi=O.list=M.list=as.list(numeric(num.class.na))
  uni.ind=unique(ind.na)
  for(i in 1:num.class.na)
  {
    row.posi[[i]]=which(ind.na==uni.ind[i]) # Get the row indices for each missing data pattern
    O.list[[i]]=matrix(diag(p)[!na.posi[row.posi[[i]][1],],],ncol=p) # O matrix for each missing data pattern
    M.list[[i]]=matrix(diag(p)[na.posi[row.posi[[i]][1],],],ncol=p) # M matrix for each missing data pattern
  }
  Y=Y.na
  Y[is.na(Y)]=99999
  Ip=diag(p)
  Y.hat=array(NA,dim=c(n,p,g))
  log.w.den=matrix(NA,n,g)
  Sig.mo=Om=list()
  cov.o=array(NA,dim=c(p,p,n))
  den.mo=numeric(n)
  for(i in 1:g)
  {
    #Y.hat[,,i]=as.matrix(Y)
    cent=t(Y)-mu[,i]
    for(j in 1:num.class.na)
    {
      ind=row.posi[[j]] # Row indices of observations with the j-th missing data pattern
      O=O.list[[j]] # O matrix for the j-th missing data pattern
      M=M.list[[j]] # M matrix for the j-th missing data pattern
      a.m=a*M #A.m={Y.m|a.m<y.m<b.m}
      b.m=b*M #A.m={Y.m|a.m<y.m<b.m}
      Y.o=as.matrix(Y[ind,])%*%t(O) # Observed data for the j-th missing data pattern
      mu.o= O%*%mu[,i] 
      OSO=O%*%S[,,i]%*%t(O) 
      Soo=t(O)%*%solve(OSO)%*%O
      SiSoo=S[,,i]%*%Soo
      cent.ind=cent[,ind]
      if(!all(M==0)) # Case with missing data
      {
        mu.m=M%*%mu[,i]
        MSM=M%*%S[,,i]%*%t(M)
        mu.mo=M%*%(mu[,i]+SiSoo%*%cent.ind)
        Sig.mo[[j]]=M%*%(Ip-SiSoo)%*%S[,,i]%*%t(M)
        eta=matrix(NA,nrow(mu.mo),length(ind))
        Psi=array(NA,dim = c(nrow(M),nrow(M),length(ind)))
        for(s in 1:length(ind))
        {
          tmp=meanvarTMD(lower = a.m[which(M!=0)],upper = b.m[which(M!=0)], mu =mu.mo[,s], Sigma = Sig.mo[[j]],  dist = 'normal')
          eta[,s]=tmp$mean # Y.m|Y.o~TN_pjm(mu.mo,Sig.mo;A.m) compute the conditional expectation of Y.m|Y.o denote as eta.mo
          Psi[,,s] = tmp$varcov # Y.m|Y.o~TN_pjm(mu.mo,Sig.mo;A.m) compute the conditional covariance matrix of Y.m|Y.o denote as Psi.mmo
          den.mo[ind]=abs(sadmvn(lower=a.m[which(M!=0)], upper=b.m[which(M!=0)], mean=mu.mo[,s], varcov=Sig.mo[[j]])[1]) # Compute the integral of N_pjm(Y.m;mu.mo,Sig.mo) for the j-th missing data pattern
          cov.o[,,ind[s]]=t(M)%*%Psi[,,s]%*%M # Cov(Y|Y.o)
        }
        Y.hat[ind,,i]=t(t(O)%*%t(Y.o))+t(t(M)%*%eta) #E-step: compute Y.hat
        log.Den=log(sadmvn(lower=a, upper=b, mean=mu[,i], varcov=S[,,i])[1])
        log.den = log(dmvnorm(Y.o,mu.o,OSO))+log(den.mo[ind]) # Compute the numerator of the likelihood f(Yj.o), then take the log
        log.w.den[ind,i]=log(w[i])+log.den-log.Den #obsrved log likelihood
      }
      else{ # all(M == 0) indicates no missing data; the observed data Y.o equals the full data Y
        den.mo[ind]=0
        cov.o[,,ind]=matrix(0,p,p)
        Y.hat[ind,,i]=Y.o # If there is no missing data, then Y.hat equals Y
        log.den = log(dmvnorm(Y,mu[,i],S[,,i]))-log(abs(pmvnorm(lower=a, upper=b, mean=mu[,i], sigma=S[,,i],seed = 1)[1]))
        log.w.den[ind,i]=log(w[i])+log.den[ind] #f(Yj.o)=f(Yj)
      }
    }
    Om[[i]]=cov.o
  }
  log.max.wden = apply(log.w.den, 1, max)
  w.den = exp(log.w.den - log.max.wden)
  indv.den = rowSums(w.den)
  Z=w.den/indv.den
  EX=matrix(NA,p,g)
  EXX=array(NA,dim = c(p,p,g))
  for(i in 1:g)
  {
    mom = meanvarTMD(lower = a-mu[,i],upper = b-mu[,i], mu = rep(0,p), Sigma = S[,,i],  dist = 'normal')
    EX[,i] = mom$mean
    EXX[,,i] = mom$EYY
  }
  # Calculate the score vector
  S.w = rep(0,g-1)
  if(g > 1){
    S.w = matrix(0,n,g-1)
    for(i in 1:(g-1)){ 
      S.w[,i] = Z[,i]/w[i]-Z[,g]/w[g]
    }
  }else{S.w = 0}
  S.mu=0
  for(i in 1:g){ 
    S.mui = matrix(0,n,p)
    for(j in 1:n){
      S.mui[j,] = Z[j,i]*solve(S[,,i])%*%t((t(Y.hat[j,,i])-mu[,i]-EX[,i]))
    } 
    S.mu=cbind(S.mu,S.mui)
  }
  S.mu = S.mu[,-1]
  S.omega = 0
  for(i in 1:g){
    y.cent = t(Y.hat[,,i]) - mu[,i]
    S.oi = matrix(0,n,p*(p+1)/2)
    for(j in 1:n){
      A = (solve(S[,,i])%*%(y.cent[,j] %*% t(y.cent[,j]) + cov.o[,,j])%*%solve(S[,,i])-solve(S[,,i])%*%EXX[,,i]%*%solve(S[,,i]))
      S.sigma = Z[j,i]*(A-1/2*diag(diag(A)))
      S.oi[j,]=vech(S.sigma)
    }
    S.omega=cbind(S.omega,S.oi)
  }
  S.omega=S.omega[,-1]
  # Calculate the Fisher information matrix
  Q=cbind(S.w,S.mu,S.omega)
  sd=sqrt(diag(solve(t(Q)%*%Q)))
  para = w[1:(g-1)]
  para.name = paste("w",1:(g-1),sep="")
  r=p*(p+1)/2
  se1=0
  for(i in 1:g)
  {
    para = round(c(para,mu[,i],vech(S[,,i])),3)
    para.name = c(para.name,paste("mu_",i,",",1:p,sep=""))
    para.name = c(para.name,paste("S_",i,",",vech(outer(1:p,1:p,paste,sep=",")),sep=""))
    se1=round(c(se1,sd[((i-1)*p+g):((i-1)*p+g+(p-1))],sd[((i-1)*r+p*g+g):(i*r+p*g+g-1)]),3)
  }
  se1=se1[-1];se1=round(c(sd[1:(g-1)],se1),3)
  SE1 = cbind(est=para,sd=se1)
  rownames(SE1)=para.name
  return(SE1=SE1)
}
