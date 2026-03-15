# EM algorithm for the FM-TMVN model
MTMVN.na.EM=function(Y.na, Y.knn, a, b, g, init.clus, true.clus, tol, max.iter, per)
{
  begin=proc.time()[1]
  n=nrow(Y.na)
  p=ncol(Y.na)
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
  w=rep(NA,g)
  mu=matrix(NA,p,g)
  S=array(NA,dim=c(p,p,g))
  #initial value
  for(i in 1:g){ 
    ind.Y = as.matrix(Y.knn[init.clus==i,])#Find the positions of observations in group i
    w[i]=nrow(ind.Y)/n #Calculate the weight for component i
  mu[,i]=colMeans(ind.Y,na.rm = T) #Calculate group i's mean vector (NA removed)
  S[,,i]=var(ind.Y, na.rm = T) #Calculate group i's covariance matrix (NA removed)
  }
  Y.hat=array(NA,dim=c(n,p,g))
  log.w.den=matrix(NA,n,g)
  Sig.mo=Om=list()
  cov.o=array(NA,dim=c(p,p,n))
  den.mo=numeric(n)
  for(i in 1:g)
  {
    cent=t(Y)-mu[,i]
    for(j in 1:num.class.na)
    {
      ind=row.posi[[j]] # Row indices of observations with the j-th missing data pattern
      O=O.list[[j]] # O matrix for the j-th missing data pattern
      M=M.list[[j]] # M matrix for the j-th missing data pattern
      a.m=a*M #A.m={Y.m|a.m<y.m<b.m}
      b.m=b*M #A.m={Y.m|a.m<y.m<b.m}
      Y.o=matrix(as.numeric(unlist(Y[ind, ])), ncol = p, byrow = FALSE)%*%t(O) # Observed data for the j-th missing data pattern
      mu.o= O%*%mu[,i] 
      OSO=O%*%S[,,i]%*%t(O) 
      Soo=t(O)%*%solve(OSO)%*%O
      SiSoo=S[,,i]%*%Soo
      cent.ind=cent[,ind]
      if(dim(M)[1]!=0) # Case with missing data
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
          eta[,s] = tmp$mean #Y.m|Y.o~TN_pjm(mu.mo,Sig.mo;A.m) compute the conditional expectation of Y.m|Y.o denote as eta.mo
          Psi[,,s] = tmp$varcov #Y.m|Y.o~TN_pjm(mu.mo,Sig.mo;A.m) compute the conditional covariance matrix of Y.m|Y.o denote as Psi.mmo
          den.mo[ind] = abs(sadmvn(lower=a.m[which(M!=0)], upper=b.m[which(M!=0)], mean=mu.mo[,s], varcov=Sig.mo[[j]])[1]) #Compute the integral of N_pjm(Y.m;mu.mo,Sig.mo) for the j-th missing data pattern
          cov.o[,,ind[s]]=t(M)%*%Psi[,,s]%*%M #Cov(Y|Y.o)
        }
        Y.hat[ind,,i]=t(t(O)%*%t(Y.o))+t(t(M)%*%eta) #E-step: compute Y.hat
        log.Den=log(abs(sadmvn(lower=a, upper=b, mean=mu[,i], varcov=S[,,i])[1]))
        log.den = log(dmvnorm(Y.o,mu.o,OSO))+log(den.mo[ind]) # Compute the numerator of the likelihood f(Yj.o), then take the log
        log.w.den[ind,i]=log(w[i])+log.den-log.Den #obsrved log likelihood
      }
      else{ # all(M == 0) indicates no missing data; the observed data Y.o equals the full data Y
        den.mo[ind]=0
        cov.o[,,ind]=matrix(0,p,p)
        Y.hat[ind,,i]=Y.o # If there is no missing data, then Y.hat equals Y
        log.den = log(dmvnorm(Y,mu[,i],S[,,i]))-log(abs(sadmvn(lower=a, upper=b, mean=mu[,i], varcov=S[,,i])[1]))
        log.w.den[ind,i]=log(w[i])+log.den[ind] #f(Yj.o)=f(Yj)
      }
    }
    Om[[i]]=cov.o
  }
  log.max.wden = apply(log.w.den, 1, max)
  w.den = exp(log.w.den - log.max.wden)
  indv.den = rowSums(w.den)
  log.indv.den = log(indv.den) + log.max.wden
  loglike.old = sum(log.indv.den) 
  iter = 0
  cat(paste(rep('-', 20), sep = '', collapse = ''), 'Truncated Normal EM: ',  paste(rep('-', 20), sep = '', collapse = ''), '\n')
  cat("iter =", iter, "\t loglik =", loglike.old, "\n")
  #iter.loglik=loglik.old
  #E-step
  repeat{
    iter=iter+1
    Z=w.den/indv.den
    which(is.nan(Z), arr.ind = TRUE)
    ni=colSums(Z)
    w=ni/n
    for(i in 1:g)
    {
      mom = mom.mtruncnorm(2, mean=rep(0,p), varcov=S[,,i], lower=a-mu[,i], upper=b-mu[,i])
      EX = mom$cum1
      EXX = mom$order2$m2
      #M-step
      mu[,i] = colSums(Z[,i]*Y.hat[,,i])/ni[i]-EX
      y.cent = t(Y.hat[,,i]) - mu[,i]
      yy.hat_2 = matrix(0,p,p)
      for(j in 1:n){
        yy.hat_2 = yy.hat_2 +  Z[j,i]*Om[[i]][,,j]
      }
      yy.hat_1 = t(sqrt(Z[,i])*t(y.cent))%*%(sqrt(Z[,i])*t(y.cent))
      yy.hat = yy.hat_1+yy.hat_2
      S[,,i] = yy.hat/ni[i] + S[,,i] - EXX
      cent = t(Y) - mu[,i]
      for(j in 1:num.class.na)
      {
        ind=row.posi[[j]] # Row indices of observations corresponding to the i-th missing data pattern
        O=O.list[[j]] # O matrix for the i-th missing data pattern
        M=M.list[[j]] # M matrix for the i-th missing data pattern
        a.m=a*M #A.m={Y.m|a.m<y.m<b.m}
        b.m=b*M
        Y.o=matrix(as.numeric(unlist(Y[ind, ])), ncol = p, byrow = FALSE)%*%t(O)
        mu.o=O%*%mu[,i]
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
          Psi=array(NA,dim = c(nrow(mu.mo),nrow(mu.mo),length(ind)))
          for(s in 1:length(ind))
          {
            tmp=meanvarTMD(lower = a.m[which(M!=0)],upper = b.m[which(M!=0)], mu =mu.mo[,s], Sigma = Sig.mo[[j]],  dist = 'normal')
            eta[,s] = tmp$mean #Y.m|Y.o~TN_pjm(mu.mo,Sig.mo;A.m) compute the conditional expectation of Y.m|Y.o denote as eta.mo
            Psi[,,s] = tmp$varcov #Y.m|Y.o~TN_pjm(mu.mo,Sig.mo;A.m) compute the conditional covariance matrix of Y.m|Y.o denote as Psi.mmo
            den.mo[ind]=abs(sadmvn(lower=a.m[which(M!=0)], upper=b.m[which(M!=0)], mean=mu.mo[,s], varcov=Sig.mo[[j]])[1]) #計算第j個遺失pattern中N_pjm(Y.m;mu.mo,Sig.mo)的積分
            cov.o[,,ind[s]]=t(M)%*%Psi[,,s]%*%M #Cov(Y|Y.o)
          } 
          Y.hat[ind,,i]=t(t(O)%*%t(Y.o))+t(t(M)%*%eta) #E-step: compute Y.hat
          log.Den=log(sadmvn(lower=a, upper=b, mean=mu[,i], varcov=S[,,i])[1])
          log.den = log(dmvnorm(Y.o,mu.o,OSO))+log(den.mo[ind]) # Compute the numerator of the likelihood f(Yj.o), then take the log
          log.w.den[ind,i]=log(w[i])+log.den-log.Den
        }
        else{ # all(M == 0) indicates no missing data; the observed data Y.o equals the full data Y
          den.mo[ind]=0
          cov.o[,,ind]=matrix(0,p,p)
          Y.hat[ind,,i]=Y.o # If there is no missing data, then Y.hat equals Y
          log.den = log(dmvnorm(Y,mu[,i],S[,,i]))-log(abs(sadmvn(lower=a, upper=b, mean=mu[,i], varcov=S[,,i])[1]))
          log.w.den[ind,i]=log(w[i])+log.den[ind] #f(Yj.o)=f(Yj)
        }
      }
      Om[[i]]=cov.o
    }
    log.max.wden = apply(log.w.den, 1, max)
    w.den = exp(log.w.den - log.max.wden)
    indv.den = rowSums(w.den)
    log.indv.den = log(indv.den) + log.max.wden
    loglike.new = sum(log.indv.den)
    #iter.loglik=c(iter.loglik,loglik.new)
    diff=loglike.new-loglike.old
    if(iter%%per==0) cat("iter=",iter,",\t obs.logli=",loglike.new,"\t diff=",diff,"\n")
    if(diff<tol|iter==max.iter) break
    loglike.old=loglike.new
  }
  Z=w.den/indv.den
  
  post.clus=matrix(apply(Z,1,order),nrow=g)[g,]
  if(is.logical(true.clus)==F){
    if(length(unique(post.clus))==length(unique(true.clus)))
    {
      CCR=1-classError(true.clus,post.clus)$errorRate
    }
    else {
      CCR=sum(apply(table(post.clus,true.clus),1,max))/n
    }
    ARI=adjustedRandIndex(true.clus,post.clus)
  }
  else{
    CCR = ARI = NULL
  }
  m = (g-1)+g*p+g*(p*(p+1)/2)
  AIC=2*m-2*loglike.new
  BIC=m*log(n)-2*loglike.new
  ICL=BIC-2*sum(Z*log(Z+(1e-300)))
  # predicting missing value 
  post.pred = apply(rep(Z, p) * aperm(Y.hat, perm = c(1,3,2)), 3, rowSums)
  cat(paste(rep("-",60),sep="",collapse=""),"\n")
  cat("iter =", iter, "\t logli =", loglike.new, "\t diff =", diff, "\n")
  para=list(w=w, mu=mu, Sigma=S)
  end=proc.time()[1]
  sec.time=end-begin
  cat("It took",sec.time,"seconds. \n")
  return(list(iter=iter, sec.time=sec.time, m=m, loglik=loglike.new, para=para, post.clus=post.clus, AIC=AIC, BIC=BIC, ICL=ICL, ARI=ARI, CCR=CCR, post.pred=post.pred, diff=diff))
}
