library(mvtnorm)
library(relliptical)
################################################################################
#
#   Filename: figF1.R
#   Purpose: Generate scatter–histogram plots with marginal histograms and contour 
#            overlays for one simulation replicate (n = 600).
#   Input data files: function/gener.na.R
#   Output data files: results/figF1.eps 
#   R Version: R-4.4.1
#   Required R packages: mvtnorm, relliptical
#
################################################################################

source(paste0(SPATH, "/function/gener.na.R"))
n=600
p=2
g=2
a = c(0, 0)
b = c(5, 5)
w = c(0.4, 0.6)
ni = as.vector(rmultinom(1, size = n, prob = w))
true.clus = rep(1:g, times = ni)
mu = matrix(c(1.2, 1.8, 3.8, 3.2), p, g)
Sigma = array(NA, dim = c(p, p, g))
Sigma[,,1] = matrix(c(1.4, 0.9, 0.9, 0.8), 2, 2)
Sigma[,,2] = matrix(c(0.7, -0.6, -0.6, 1.3), 2, 2)
set.seed(325)
Y=NULL
for(i in 1:g){
  Y = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a, upper=b, dist="Normal")))
}
Y.na = gener.na(Y, 0.1)
na.posi = is.na(Y.na)
na.row = which(rowSums(is.na(Y.na)) != 0)

postscript(paste0(SPATH, "results/figS1.eps", width=8, height=10))
layout (matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE),widths = c(3, 1.5), heights = c(1.5, 3))
par(mar=c(4, 5, 0, 0), cex.lab=1.5, bty='o')
col_0=c("pink", "lightgreen")
col_1=c( "red", "darkgreen")
m=50
tmp1=range(Y[,1])
tmp2=range(Y[,2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+w[i]*(dmvnorm(xx, mean = mu[,i], sigma = Sigma[,,i])/pmvnorm(lower=a, upper=b, mean=mu[,i], sigma=Sigma[,,i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col=gray(.6), lwd=0.005, drawlabels = F, nlevels=10,xlim=c(min(Y[,1]),b[1]), ylim=c(min(Y[,2]),b[2]), xlab='Variable 1', ylab='Variable 2')
points(Y[-na.row,], col = col_0[true.clus[-na.row]], pch=true.clus[-na.row]+14, cex=0.7)
text(Y[na.row,], labels = "x", col = col_1[true.clus[na.row]], cex = 1.2, font = 0.5)
y1 = list(y11 = Y[true.clus==1, 1], y12 = Y[true.clus==2, 1])
h1 = lapply(y1, hist, breaks = seq(min(Y[,1]), max(Y[,1]), length=40), plot=F)
t1 = rbind(h1[[1]]$counts, h1[[2]]$counts)
rownames(t1)=names(h1)
colnames(t1)=h1[[1]]$mids
par(mar=c(0, 4, 0.5, 0), cex.lab=1.5, bty='n')
a1 = which(h1$y11$breaks >= 1 & h1$y11$breaks<=4)
barplot(t1[1,a1], ylim=c(0,max(t1)),col=0, border = 2, xaxt='n', las=1, space=0, xaxt='n', yaxt='n')
barplot(t1[2,a1], ylim=c(0,max(t1)),col=0, border = 3, xaxt='n', las=1, space=0, add=T, xaxt='n', yaxt='n')
y2 = list(y21 = Y[true.clus==1, 2], y22 = Y[true.clus==2, 2])
h2 = lapply(y2, hist, breaks = seq(min(Y[,2]), max(Y[,2]), length=80), plot=F)
t2 = rbind(h2[[1]]$counts, h2[[2]]$counts)
rownames(t2)=names(h2)
colnames(t2)=h2[[1]]$mids
par(mar=c(4, 0, 0, 0.5), cex.lab=1.5, bty='n')
a2 = which(h2$y21$breaks >= 2 & h2$y21$breaks<=3.5)
barplot(t2[2,a2], xlim=c(0,max(t2)), horiz=T, col=0, border = 3, xaxt='n', las=1, space=0, xaxt='n', yaxt='n')
barplot(t2[1,a2], xlim=c(0,max(t2)), horiz=T, col=0, border = 2, xaxt='n', las=1, space=0, add=T, xaxt='n', yaxt='n')
dev.off()
