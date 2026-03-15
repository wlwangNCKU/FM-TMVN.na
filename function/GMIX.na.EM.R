# EM algorithm for the GMM model
GMIXM.VVV.EM = function(Y.na, Y.knn, g, init.clus , true.clus , tol , max.iter.EM,per)
{
  begin = proc.time()[1]
  n = nrow(Y.na) 
  p = ncol(Y.na) 
  na.posi = is.na(Y.na) 
  po = p - rowSums(na.posi) 
  ind.na = colSums(t(na.posi) * 2 ^ (1:p - 1))
  num.class.na = length(unique(ind.na))
  row.posi = O.list = as.list(numeric(num.class.na))
  uni.ind = unique(ind.na)
  for(i in 1:num.class.na)
  {
    row.posi[[i]] = which(ind.na == uni.ind[i])
    O.list[[i]] = matrix(diag(p)[!na.posi[row.posi[[i]][1],],], ncol = p)
  }
  Ip = diag(p)
  #observed pdf
  dmvnorm.o = function (po, det.o, delta.o)
  {
    (2 * pi) ^ (- po / 2) * det.o ^ (- 1 / 2) * exp(- delta.o / 2)
  }
  Y.na[na.posi] = 99999
  #initial values
  w = rep(NA,g);mu = matrix(NA, p, g)
  Sigma = array(NA, dim = c(p,p,g))
  for(i in 1:g)
  {
  ind.Y = as.matrix(Y.knn[init.clus==i,]) #Find the positions of observations in group i
  w[i]=nrow(ind.Y)/n #Calculate the weight for component i
  mu[,i]=colMeans(ind.Y,na.rm = T) #Calculate group i's mean vector (NA removed)
  Sigma[,,i]=var(ind.Y, na.rm = T) #Calculate group i's covariance matrix (NA removed)
  }
  #calculate log-likelihood
  det.o = rep(NA, n)
  w.den = delta.o = matrix(NA, n, g)
  Y.hat = array(NA, dim = c(n, p, g))
  na.cov = array(NA, dim = c(p, p, num.class.na, g))
  for(i in 1:g)
  {
    cent = t(Y.na) - mu[,i]
    for(j in 1:num.class.na)
    {
      O = O.list[[j]]
      OSO = O %*% Sigma[,,i] %*% t(O)
      ind = row.posi[[j]]
      cent.ind = cent[,ind]
      det.o[ind] = det(OSO)
      Soo = t(O) %*% solve(OSO) %*% O
      delta.o[ind, i] = colSums(cent.ind * (Soo %*% cent.ind))
      SiSoo = Sigma[,,i] %*% Soo
      Y.hat[ind,,i] = t(mu[,i] + SiSoo %*% cent.ind)
      na.cov[,,j,i] = (Ip - SiSoo) %*% Sigma[,,i]
    }
    w.den[,i] = w[i] * dmvnorm.o(po = po, det.o = det.o, delta.o = delta.o[,i])
  }
  indv.den = rowSums(w.den)
  loglike.old = sum(log(indv.den))
  iter = 0
  cat(paste(rep('-', 20), sep = '', collapse = ''), ' Normal EM: ',  paste(rep('-', 20), sep = '', collapse = ''), '\n')
  cat("iter =", iter, "\t loglik =", loglike.old, "\n")
  repeat
  {
    iter = iter + 1
    #E-step
    Z.mt = w.den / indv.den
    #M-step
    n.z = colSums(Z.mt)
    w = n.z / n
    for(i in 1:g)
    {
      mu[,i] = colSums(Z.mt[,i] * Y.hat[,,i]) / n.z[i]
      cent = t(Y.hat[,,i]) - mu[,i]
      sum.na.cov = matrix(0, p, p)
      for(j in 1:num.class.na)
      {
        sum.na.cov = sum.na.cov + na.cov[,,j,i] * sum(Z.mt[row.posi[[j]],i])
      }
      Sigma[,,i] = (cent %*% (Z.mt[,i] * t(cent)) + sum.na.cov) / n.z[i]
    }
    #calculate log-likelihood
    for(i in 1:g)
    {
      cent = t(Y.na) - mu[,i]
      for(j in 1:num.class.na)
      {
        O = O.list[[j]]
        OSO = O %*% Sigma[,,i] %*% t(O)
        ind = row.posi[[j]]
        cent.ind = cent[,ind]
        det.o[ind] = det(OSO)
        Soo = t(O) %*% solve(OSO) %*% O
        delta.o[ind, i] = colSums(cent.ind * (Soo %*% cent.ind))
        SiSoo = Sigma[,,i] %*% Soo
        Y.hat[ind,,i] = t(mu[,i] + SiSoo %*% cent.ind)
        na.cov[,,j,i] = (Ip - SiSoo) %*% Sigma[,,i]
      }
      w.den[,i] = w[i] * dmvnorm.o(po = po, det.o = det.o, delta.o = delta.o[,i])
    }
    indv.den = rowSums(w.den)
    loglike.new = sum(log(indv.den))
    diff=loglike.new-loglike.old 
    if(iter%%per==0) cat("iter=",iter,",\t obs.logli=",loglike.new,"\t diff=",diff,"\n")
    if(diff < tol | iter == max.iter.EM) break
    loglike.old = loglike.new
  }
  #posterior classification
  Z.mt = w.den / indv.den 
  post.clus = matrix(apply(Z.mt, 1, order),nrow=g)[g,] 
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
  ICL=BIC-2*sum(Z.mt*log(Z.mt+(1e-300)))
  #predicting missing value
  post.pred = apply(rep(Z.mt, p) * aperm(Y.hat, perm = c(1,3,2)), 3, rowSums)
  cat(paste(rep("-",60),sep="",collapse=""),"\n")
  cat("iter =", iter, "\t logli =", loglike.new, "\t diff =", diff, "\n")
  para=list(w=w,mu=mu,Sigma=Sigma)
  end = proc.time()[1]
  sec.time=end-begin
  cat("It took",sec.time,"seconds. \n")
  return(list(iter=iter, sec.time=sec.time, m=m, loglik=loglike.new, para=para, post.clus=post.clus, AIC=AIC, BIC=BIC, ICL=ICL, ARI=ARI, CCR=CCR, post.pred=post.pred, diff=diff))
}
