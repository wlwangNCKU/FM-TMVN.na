################################################################################
#
#   Filename: fig5.R
#   Purpose: Generate pairwise scatter plots of the Stone flakes data with observed  
#            points, GMM and FM-TMVN imputed values, and overlay component  
#            contours for each mixture model.
#   Input data files: data/stone_result.RData
#   Output data files: results/fig5.eps 
#   R Version: R-4.4.1
#
################################################################################

load(paste(SPATH, "/data/stone_result.RData", sep=''))

postscript(paste0(SPATH, "results/fig5.eps", width = 8, height = 14))
par(mfrow = c(4, 3), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))
g = 2
v = c("LBI", "FLA", "PSF", "FSF", "ZDF1")
m = 500
#--1--
v1 = 1; v2 = 2
tmp1=range(total$`2_1`$post.pred[,v1])
tmp2=range(total$`2_1`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+total$`2_1`$para$w[i]*(dmvnorm(xx, mean = total$`2_1`$para$mu[c(v1,v2),i], sigma = total$`2_1`$para$Sigma[c(v1,v2),c(v1,v2),i]))
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col="#FFD700",lty = 5, lwd=2, drawlabels = F, nlevels=8,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2])
tmp1=range(total$`2_2`$post.pred[,v1])
tmp2=range(total$`2_2`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
g = 2
den=0
for(i in 1:g) {
  den=den+total$`2_2`$para$w[i]*(dmvnorm(xx, mean = total$`2_2`$para$mu[c(v1,v2),i], sigma = total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])/pmvnorm(lower=a[c(v1,v2)], upper=b[c(v1,v2)], mean=total$`2_2`$para$mu[c(v1,v2),i], sigma=total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col= "black",lty = 1, lwd=0.9, drawlabels = F, nlevels=8,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2], add = TRUE)
na.posi = is.na(Y.na[,c(v1,v2)])
points(sc.y[,c(v1,v2)], col = "#0000FF", pch = 1, cex = 1)
points(total$`2_2`$post.pred[,c(v1,v2)][na.posi], col = "green",pch = 17, cex = 1)
points(total$`2_1`$post.pred[,c(v1,v2)][na.posi], col = "red", pch = 8 , cex = 1)
axis(1, at = seq(floor(a[v1] - 1), ceiling(b[v1] + 1), by = 1))
axis(2, at = seq(floor(a[v2] - 0.5), ceiling(b[v2] + 0.5), by = 1))
#--2--
v1 = 1; v2 = 3
tmp1=range(total$`2_1`$post.pred[,v1])
tmp2=range(total$`2_1`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+total$`2_1`$para$w[i]*(dmvnorm(xx, mean = total$`2_1`$para$mu[c(v1,v2),i], sigma = total$`2_1`$para$Sigma[c(v1,v2),c(v1,v2),i]))
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col="#FFD700",lty = 5, lwd=2, drawlabels = F, nlevels=10,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2])
tmp1=range(total$`2_2`$post.pred[,v1])
tmp2=range(total$`2_2`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
g = 2
den=0
for(i in 1:g) {
  den=den+total$`2_2`$para$w[i]*(dmvnorm(xx, mean = total$`2_2`$para$mu[c(v1,v2),i], sigma = total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])/pmvnorm(lower=a[c(v1,v2)], upper=b[c(v1,v2)], mean=total$`2_2`$para$mu[c(v1,v2),i], sigma=total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col= "black",lty = 1, lwd=0.9, drawlabels = F, nlevels=6,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2], add = TRUE)
na.posi = is.na(Y.na[,c(v1,v2)])
points(sc.y[,c(v1,v2)], col = "#0000FF", pch = 1, cex = 1)
points(total$`2_2`$post.pred[,c(v1,v2)][na.posi], col = "green",pch = 17, cex = 1)
points(total$`2_1`$post.pred[,c(v1,v2)][na.posi], col = "red", pch = 8 , cex = 1)
axis(1, at = seq(floor(a[v1] - 1), ceiling(b[v1] + 1), by = 1))
axis(2, at = seq(floor(a[v2] - 0.5), ceiling(b[v2] + 0.5), by = 1))
#--3--
v1 = 1; v2 = 4
tmp1=range(total$`2_1`$post.pred[,v1])
tmp2=range(total$`2_1`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+total$`2_1`$para$w[i]*(dmvnorm(xx, mean = total$`2_1`$para$mu[c(v1,v2),i], sigma = total$`2_1`$para$Sigma[c(v1,v2),c(v1,v2),i]))
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col="#FFD700",lty = 5, lwd=2, drawlabels = F, nlevels=7,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2])
tmp1=range(total$`2_2`$post.pred[,v1])
tmp2=range(total$`2_2`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
g = 2
den=0
for(i in 1:g) {
  den=den+total$`2_2`$para$w[i]*(dmvnorm(xx, mean = total$`2_2`$para$mu[c(v1,v2),i], sigma = total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])/pmvnorm(lower=a[c(v1,v2)], upper=b[c(v1,v2)], mean=total$`2_2`$para$mu[c(v1,v2),i], sigma=total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col= "black",lty = 1, lwd=0.9, drawlabels = F, nlevels=7,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2], add = TRUE)
na.posi = is.na(Y.na[,c(v1,v2)])
points(sc.y[,c(v1,v2)], col = "#0000FF", pch = 1, cex = 1)
points(total$`2_2`$post.pred[,c(v1,v2)][na.posi], col = "green",pch = 17, cex = 1)
points(total$`2_1`$post.pred[,c(v1,v2)][na.posi], col = "red", pch = 8 , cex = 1)
axis(1, at = seq(floor(a[v1] - 1), ceiling(b[v1] + 1), by = 1))
axis(2, at = seq(floor(a[v2] - 0.5), ceiling(b[v2] + 0.5), by = 1))
#--4--
v1 = 1; v2 = 5
tmp1=range(total$`2_1`$post.pred[,v1])
tmp2=range(total$`2_1`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+total$`2_1`$para$w[i]*(dmvnorm(xx, mean = total$`2_1`$para$mu[c(v1,v2),i], sigma = total$`2_1`$para$Sigma[c(v1,v2),c(v1,v2),i]))
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col="#FFD700",lty = 5, lwd=2, drawlabels = F, nlevels=7,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2])
tmp1=range(total$`2_2`$post.pred[,v1])
tmp2=range(total$`2_2`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
g = 2
den=0
for(i in 1:g) {
  den=den+total$`2_2`$para$w[i]*(dmvnorm(xx, mean = total$`2_2`$para$mu[c(v1,v2),i], sigma = total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])/pmvnorm(lower=a[c(v1,v2)], upper=b[c(v1,v2)], mean=total$`2_2`$para$mu[c(v1,v2),i], sigma=total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col= "black",lty = 1, lwd=0.9, drawlabels = F, nlevels=7,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2], add = TRUE)
na.posi = is.na(Y.na[,c(v1,v2)])
points(sc.y[,c(v1,v2)], col = "#0000FF", pch = 1, cex = 1)
points(total$`2_2`$post.pred[,c(v1,v2)][na.posi], col = "green",pch = 17, cex = 1)
points(total$`2_1`$post.pred[,c(v1,v2)][na.posi], col = "red", pch = 8 , cex = 1)
axis(1, at = seq(floor(a[v1] - 1), ceiling(b[v1] + 1), by = 1))
axis(2, at = seq(floor(a[v2] - 0.5), ceiling(b[v2] + 0.5), by = 1))
#--5--
v1 = 2; v2 = 3
tmp1=range(total$`2_1`$post.pred[,v1])
tmp2=range(total$`2_1`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+total$`2_1`$para$w[i]*(dmvnorm(xx, mean = total$`2_1`$para$mu[c(v1,v2),i], sigma = total$`2_1`$para$Sigma[c(v1,v2),c(v1,v2),i]))
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col="#FFD700",lty = 5, lwd=2, drawlabels = F, nlevels=10,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2])
tmp1=range(total$`2_2`$post.pred[,v1])
tmp2=range(total$`2_2`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
g = 2
den=0
for(i in 1:g) {
  den=den+total$`2_2`$para$w[i]*(dmvnorm(xx, mean = total$`2_2`$para$mu[c(v1,v2),i], sigma = total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])/pmvnorm(lower=a[c(v1,v2)], upper=b[c(v1,v2)], mean=total$`2_2`$para$mu[c(v1,v2),i], sigma=total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col= "black",lty = 1, lwd=0.9, drawlabels = F, nlevels=12,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2], add = TRUE)
na.posi = is.na(Y.na[,c(v1,v2)])
points(sc.y[,c(v1,v2)], col = "#0000FF", pch = 1, cex = 1)
points(total$`2_2`$post.pred[,c(v1,v2)][na.posi], col = "green",pch = 17, cex = 1)
points(total$`2_1`$post.pred[,c(v1,v2)][na.posi], col = "red", pch = 8 , cex = 1)
axis(1, at = seq(floor(a[v1] - 1), ceiling(b[v1] + 1), by = 1))
axis(2, at = seq(floor(a[v2] - 0.5), ceiling(b[v2] + 0.5), by = 1))
#--6--
v1 = 2; v2 = 4
tmp1=range(total$`2_1`$post.pred[,v1])
tmp2=range(total$`2_1`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+total$`2_1`$para$w[i]*(dmvnorm(xx, mean = total$`2_1`$para$mu[c(v1,v2),i], sigma = total$`2_1`$para$Sigma[c(v1,v2),c(v1,v2),i]))
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col="#FFD700",lty = 5, lwd=2, drawlabels = F, nlevels=6,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2])
tmp1=range(total$`2_2`$post.pred[,v1])
tmp2=range(total$`2_2`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
g = 2
den=0
for(i in 1:g) {
  den=den+total$`2_2`$para$w[i]*(dmvnorm(xx, mean = total$`2_2`$para$mu[c(v1,v2),i], sigma = total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])/pmvnorm(lower=a[c(v1,v2)], upper=b[c(v1,v2)], mean=total$`2_2`$para$mu[c(v1,v2),i], sigma=total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col= "black",lty = 1, lwd=0.9, drawlabels = F, nlevels=6,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2], add = TRUE)
na.posi = is.na(Y.na[,c(v1,v2)])
points(sc.y[,c(v1,v2)], col = "#0000FF", pch = 1, cex = 1)
points(total$`2_2`$post.pred[,c(v1,v2)][na.posi], col = "green",pch = 17, cex = 1)
points(total$`2_1`$post.pred[,c(v1,v2)][na.posi], col = "red", pch = 8 , cex = 1)
axis(1, at = seq(floor(a[v1] - 1), ceiling(b[v1] + 1), by = 1))
axis(2, at = seq(floor(a[v2] - 0.5), ceiling(b[v2] + 0.5), by = 1))
#--7--
v1 = 2; v2 = 5
tmp1=range(total$`2_1`$post.pred[,v1])
tmp2=range(total$`2_1`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+total$`2_1`$para$w[i]*(dmvnorm(xx, mean = total$`2_1`$para$mu[c(v1,v2),i], sigma = total$`2_1`$para$Sigma[c(v1,v2),c(v1,v2),i]))
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col="#FFD700",lty = 5, lwd=2, drawlabels = F, nlevels=8,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2])
tmp1=range(total$`2_2`$post.pred[,v1])
tmp2=range(total$`2_2`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
g = 2
den=0
for(i in 1:g) {
  den=den+total$`2_2`$para$w[i]*(dmvnorm(xx, mean = total$`2_2`$para$mu[c(v1,v2),i], sigma = total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])/pmvnorm(lower=a[c(v1,v2)], upper=b[c(v1,v2)], mean=total$`2_2`$para$mu[c(v1,v2),i], sigma=total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col= "black",lty = 1, lwd=0.9, drawlabels = F, nlevels=8,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2], add = TRUE)
na.posi = is.na(Y.na[,c(v1,v2)])
points(sc.y[,c(v1,v2)], col = "#0000FF", pch = 1, cex = 1)
points(total$`2_2`$post.pred[,c(v1,v2)][na.posi], col = "green",pch = 17, cex = 1)
points(total$`2_1`$post.pred[,c(v1,v2)][na.posi], col = "red", pch = 8 , cex = 1)
axis(1, at = seq(floor(a[v1] - 1), ceiling(b[v1] + 1), by = 1))
axis(2, at = seq(floor(a[v2] - 0.5), ceiling(b[v2] + 0.5), by = 1))
#--8--
v1 = 3; v2 = 4
tmp1=range(total$`2_1`$post.pred[,v1])
tmp2=range(total$`2_1`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+total$`2_1`$para$w[i]*(dmvnorm(xx, mean = total$`2_1`$para$mu[c(v1,v2),i], sigma = total$`2_1`$para$Sigma[c(v1,v2),c(v1,v2),i]))
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col="#FFD700",lty = 5, lwd=2, drawlabels = F, nlevels=6,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2])
tmp1=range(total$`2_2`$post.pred[,v1])
tmp2=range(total$`2_2`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
g = 2
den=0
for(i in 1:g) {
  den=den+total$`2_2`$para$w[i]*(dmvnorm(xx, mean = total$`2_2`$para$mu[c(v1,v2),i], sigma = total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])/pmvnorm(lower=a[c(v1,v2)], upper=b[c(v1,v2)], mean=total$`2_2`$para$mu[c(v1,v2),i], sigma=total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col= "black",lty = 1, lwd=0.9, drawlabels = F, nlevels=6,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2], add = TRUE)
na.posi = is.na(Y.na[,c(v1,v2)])
points(sc.y[,c(v1,v2)], col = "#0000FF", pch = 1, cex = 1)
points(total$`2_2`$post.pred[,c(v1,v2)][na.posi], col = "green",pch = 17, cex = 1)
points(total$`2_1`$post.pred[,c(v1,v2)][na.posi], col = "red", pch = 8 , cex = 1)
axis(1, at = seq(floor(a[v1] - 1), ceiling(b[v1] + 1), by = 1))
axis(2, at = seq(floor(a[v2] - 0.5), ceiling(b[v2] + 0.5), by = 1))
#--9--
v1 = 3; v2 = 5
tmp1=range(total$`2_1`$post.pred[,v1])
tmp2=range(total$`2_1`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+total$`2_1`$para$w[i]*(dmvnorm(xx, mean = total$`2_1`$para$mu[c(v1,v2),i], sigma = total$`2_1`$para$Sigma[c(v1,v2),c(v1,v2),i]))
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col="#FFD700",lty = 5, lwd=2, drawlabels = F, nlevels=7,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2])
tmp1=range(total$`2_2`$post.pred[,v1])
tmp2=range(total$`2_2`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
g = 2
den=0
for(i in 1:g) {
  den=den+total$`2_2`$para$w[i]*(dmvnorm(xx, mean = total$`2_2`$para$mu[c(v1,v2),i], sigma = total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])/pmvnorm(lower=a[c(v1,v2)], upper=b[c(v1,v2)], mean=total$`2_2`$para$mu[c(v1,v2),i], sigma=total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col= "black",lty = 1, lwd=0.9, drawlabels = F, nlevels=7,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2], add = TRUE)
na.posi = is.na(Y.na[,c(v1,v2)])
points(sc.y[,c(v1,v2)], col = "#0000FF", pch = 1, cex = 1)
points(total$`2_2`$post.pred[,c(v1,v2)][na.posi], col = "green",pch = 17, cex = 1)
points(total$`2_1`$post.pred[,c(v1,v2)][na.posi], col = "red", pch = 8 , cex = 1)
axis(1, at = seq(floor(a[v1] - 1), ceiling(b[v1] + 1), by = 1))
axis(2, at = seq(floor(a[v2] - 0.5), ceiling(b[v2] + 0.5), by = 1))
#--10--
v1 = 4; v2 = 5
tmp1=range(total$`2_1`$post.pred[,v1])
tmp2=range(total$`2_1`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
den=0
for(i in 1:g) {
  den=den+total$`2_1`$para$w[i]*(dmvnorm(xx, mean = total$`2_1`$para$mu[c(v1,v2),i], sigma = total$`2_1`$para$Sigma[c(v1,v2),c(v1,v2),i]))
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col="#FFD700",lty = 5, lwd=2, drawlabels = F, nlevels=5,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2])
tmp1=range(total$`2_2`$post.pred[,v1])
tmp2=range(total$`2_2`$post.pred[,v2])
x1 = seq(tmp1[1],tmp1[2],length=m)
x2 = seq(tmp2[1],tmp2[2],length=m)
xx=expand.grid(x1,x2)
g = 2
den=0
for(i in 1:g) {
  den=den+total$`2_2`$para$w[i]*(dmvnorm(xx, mean = total$`2_2`$para$mu[c(v1,v2),i], sigma = total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])/pmvnorm(lower=a[c(v1,v2)], upper=b[c(v1,v2)], mean=total$`2_2`$para$mu[c(v1,v2),i], sigma=total$`2_2`$para$Sigma[c(v1,v2),c(v1,v2),i])[1])
}
den=matrix(den,m,m)
contour(x1, x2, den, labcex = 0.3, col= "black",lty = 1, lwd=0.9, drawlabels = F, nlevels=10,xlim=c(a[v1]-1,b[v1]+1), ylim=c(a[v2]-0.5,b[v2]+0.5), xlab=v[v1], ylab=v[v2], add = TRUE)
na.posi = is.na(Y.na[,c(v1,v2)])
points(sc.y[,c(v1,v2)], col = "#0000FF", pch = 1, cex = 1)
points(total$`2_2`$post.pred[,c(v1,v2)][na.posi], col = "green",pch = 17, cex = 1)
points(total$`2_1`$post.pred[,c(v1,v2)][na.posi], col = "red", pch = 8 , cex = 1)
axis(1, at = seq(floor(a[v1] - 1), ceiling(b[v1] + 1), by = 1))
axis(2, at = seq(floor(a[v2] - 0.5), ceiling(b[v2] + 0.5), by = 1))

plot.new()
legend("center",
       legend = c("Observed data", "GMM Imputed value", "FM-TMVN Imputed value", "GMM contours", "FM-TMVN contours"),
       col = c("blue", "red", "green", "#FFD700", "black"),
       lty = c(NA, NA, NA, 1, 1),
       lwd = c(NA, NA, NA, 2, 2),
       pch = c(1, 8, 17, NA, NA),
       text.col = "black",
       cex = 1.3,
       pt.cex = 1.4,
       bty = "n")
dev.off()
