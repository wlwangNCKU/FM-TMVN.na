library(mvtnorm)
library(MASS)
library(relliptical)
################################################################################
#
#   Filename: fig1.R
#   Purpose: Generate 3D perspective and contour plots of finite mixtures of 
#            truncated bivariate normal distributions (FM-TMVN) for varying 
#            numbers of components (g = 1, 2, 3) and truncation schemes (no, 
#            right, left, double).
#   Output data files: results/fig1.eps 
#   R Version: R-4.4.1
#   Required R packages: mvtnorm, MASS, relliptical
#
################################################################################
# Bivariate Normal PDF
bivN = function(x1, x2, mu, Sigma)
{
  p = 2
  sig21 = Sigma[2,1]; sig11 = Sigma[1,1]; sig22 = Sigma[2,2]
  rho = sig21 / sqrt(sig11*sig22)
  del = ((x1-mu[1])^2/sig11 + (x2-mu[2])^2/sig22 - 2*rho*(x1-mu[1])*(x2-mu[2])/sqrt(sig11*sig22)) / (1-rho^2)
  den = 1/((2*pi)^(p/2)*sqrt(det(Sigma))) *  exp(-del/2)
  return(den)
}
# Bivariate Mixture of Truncated Normal PDF
bivMTN = function(x1, x2, w, mu, Sigma, lower, upper)
{
  g = length(w)
  f_x = 0
  for (i in 1:g) {
    f_i = bivN(x1, x2, mu=mu[,i], Sigma=Sigma[,,i])
    cdf = pmvnorm(lower=lower, upper=upper, mean=mu[,i], sigma=Sigma[,,i])[1]
    f_x = f_x + w[i] * f_i/cdf
  }
  Tden = f_x
  return(Tden)
}

# Truncate region
ubd = 0.15; ubd2 = 1.5; ubd3 = 3
a.low = c(-Inf,-Inf); a.upp = rep(Inf, 2)
a.low1 = c(-Inf,-Inf); a.upp1 = c(ubd,ubd)
a.low2 = c(-ubd,-ubd); a.upp2 = c(Inf,Inf)
a.low3 = c(-ubd,-ubd); a.upp3 = c(ubd,ubd)
a.low21 = c(-Inf,-Inf); a.upp21 = c(ubd2,ubd2)
a.low22 = c(-ubd2,-ubd2); a.upp22 = c(Inf,Inf)
a.low23 = c(-ubd2,-ubd2); a.upp23 = c(ubd2,ubd2)
a.low31 = c(-Inf,-Inf); a.upp31 = c(ubd3,ubd3)
a.low32 = c(-ubd3,-ubd3); a.upp32 = c(Inf,Inf)
a.low33 = c(-ubd3,-ubd3); a.upp33 = c(ubd3,ubd3)

# Plot 
postscript(paste0(SPATH, "results/fig1.eps", width = 7, height = 10)) 
layout(matrix(c(0,1:19), 5, 4), c(1,rep(5,3)), c(1,rep(5,4)))
par(mar=c(0,0,0,0))
plot(0:1,0:1, type='n',xlab='',ylab='',xaxt='n',yaxt='n', bty='n')
mtext('No Truncation', 2, line=-2, cex=1, font=2)
plot(0:1,0:1, type='n',xlab='',ylab='',xaxt='n',yaxt='n', bty='n')
mtext('Right Truncation', 2, line=-2, cex=1, font=2)
plot(0:1,0:1, type='n',xlab='',ylab='',xaxt='n',yaxt='n', bty='n')
mtext('Left Truncation', 2, line=-2, cex=1, font=2)
plot(0:1,0:1, type='n',xlab='',ylab='',xaxt='n',yaxt='n', bty='n')
mtext('Double Truncation', 2, line=-2, cex=1, font=2)
set.seed(20250311)
# Case 1
n = 900
p = 2
g = 1
ni=c(n/g)
w=c(1)
mu = matrix(c(0,0), p, g)
Sigma = array(NA, dim=c(p,p,g))
Sigma[,,1] = matrix(c(3, 0, 0, 3), 2, 2)
# No Truncation
Y = NULL
for(i in 1:g) Y = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low, upper=a.upp, dist="Normal")))
x = seq(min(Y[,1]), max(Y[,1]), length=50)
y = seq(min(Y[,2]), max(Y[,2]), length=50)
den = outer(x, y, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low, upper=a.upp)

# Right Truncation
Yt1 = NULL
for(i in 1:g) Yt1 = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low1, upper=a.upp1, dist="Normal")))
x1 = seq(min(Yt1[,1]), max(Yt1[,1]), length=50)
y1 = seq(min(Yt1[,2]), max(Yt1[,2]), length=50)
Tden1 = outer(x1, y1, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low1, upper=a.upp1)

# Left Truncation
Yt2 = NULL
for(i in 1:g) Yt2 = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low2, upper=a.upp2, dist="Normal")))
x2 = seq(min(Yt2[,1]), max(Yt2[,1]), length=50)
y2 = seq(min(Yt2[,2]), max(Yt2[,2]), length=50)
Tden2 = outer(x2, y2, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low2, upper=a.upp2)

# Double Truncation
Yt3 = NULL
for(i in 1:g) Yt3 = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low3, upper=a.upp3, dist="Normal")))
x3 = seq(min(Yt3[,1]), max(Yt3[,1]), length=50)
y3 = seq(min(Yt3[,2]), max(Yt3[,2]), length=50)
Tden3 = outer(x3, y3, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low3, upper=a.upp3)

plot(0:1,0:1, type='n',xlab='',ylab='',xaxt='n',yaxt='n', bty='n')
mtext("g=1", 1, line=-2, cex=1, font=1)

# No Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran = persp(x,y,den, xlim = range(x), ylim=range(y), zlim=c(-0.035, max(den, na.rm=T)), box=T, theta=40, lwd=0.01,
             axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clines1 = contourLines(x, y, den, levels=seq(0, 0.01, length=3))
contourline = lapply(clines1, function(contour){ lines(trans3d(contour$x,contour$y, -0.035, tran), col='coral2') })
clines2 = contourLines(x, y, den, levels=seq(0.02, 0.05, length=4))
contourline = lapply(clines2, function(contour){ lines(trans3d(contour$x,contour$y, -0.035, tran), col='pink1') })
clines3 = contourLines(x, y, den, levels=seq(0.06, 0.1, length=5))
contourline = lapply(clines3, function(contour){ lines(trans3d(contour$x,contour$y, -0.035, tran), col='tomato') })
clines3 = contourLines(x, y, den, levels=seq(0.1, max(den), length=6))
contourline = lapply(clines3, function(contour){ lines(trans3d(contour$x,contour$y, -0.035, tran), col='red3') })
# Right Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran1 = persp(x1,y1,Tden1, xlim = range(x1)+c(-1,1), ylim=range(y1)+c(-1,1), zlim=c(-0.075, max(Tden1, na.rm=T)), box=T, theta=40, lwd=.01,
              axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clinesT11 = contourLines(x1, y1, Tden1, levels=seq(0.001, 0.02, length=3))
contourlineT1 = lapply(clinesT11, function(contour){ lines(trans3d(contour$x,contour$y, -0.075, tran1), col='coral2') })
clinesT12 = contourLines(x1, y1, Tden1, levels=seq(0.02, 0.05, length=5))
contourlineT1 = lapply(clinesT12, function(contour){ lines(trans3d(contour$x,contour$y, -0.075, tran1), col='pink1') })
clinesT13 = contourLines(x1, y1, Tden1, levels=seq(0.06, 0.1, length=4))
contourlineT1 = lapply(clinesT13, function(contour){ lines(trans3d(contour$x,contour$y, -0.075, tran1), col='tomato') })
clinesT13 = contourLines(x1, y1, Tden1, levels=seq(0.1, max(Tden1), length=5))
contourlineT1 = lapply(clinesT13, function(contour){ lines(trans3d(contour$x,contour$y, -0.075, tran1), col='red3') })

# Left Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran2 = persp(x2,y2,Tden2, xlim = range(x2)+c(-1,1), ylim=range(y2)+c(-1,1), zlim=c(-0.15, max(Tden2, na.rm=T)), box=T, theta=40, lwd=.01,
              axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clinesT21 = contourLines(x2, y2, Tden2, levels=seq(0.001, 0.02, length=3))
contourlineT2 = lapply(clinesT21, function(contour){ lines(trans3d(contour$x,contour$y, -0.15, tran2), col='coral2') })
clinesT22 = contourLines(x2, y2, Tden2, levels=seq(0.02, 0.05, length=5))
contourlineT2 = lapply(clinesT22, function(contour){ lines(trans3d(contour$x,contour$y, -0.15, tran2), col='pink1') })
clinesT23 = contourLines(x2, y2, Tden2, levels=seq(0.06, 0.1, length=4))
contourlineT2 = lapply(clinesT23, function(contour){ lines(trans3d(contour$x,contour$y, -0.15, tran2), col='tomato') })
clinesT23 = contourLines(x2, y2, Tden2, levels=seq(0.1, max(Tden2), length=5))
contourlineT2 = lapply(clinesT23, function(contour){ lines(trans3d(contour$x,contour$y, -0.15, tran2), col='red3') })

# Double Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran3 = persp(x3,y3,Tden3, xlim = range(x3)+c(-1,1), ylim=range(y3)+c(-1,1), zlim=c(-4.5, max(Tden3, na.rm=T)), box=T, theta=40, lwd=.01,
              axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clinesT31 = contourLines(x3, y3, Tden3, levels=seq(0, 0.01, length=2))
contourlineT3 = lapply(clinesT31, function(contour){ lines(trans3d(contour$x,contour$y, -4.5, tran3), col='coral2') })
clinesT32 = contourLines(x3, y3, Tden3, levels=seq(0.02, 0.05, length=3))
contourlineT3 = lapply(clinesT32, function(contour){ lines(trans3d(contour$x,contour$y, -4.5, tran3), col='pink1') })
clinesT33 = contourLines(x3, y3, Tden3, levels=seq(0.06, 0.1, length=4))
contourlineT3 = lapply(clinesT33, function(contour){ lines(trans3d(contour$x,contour$y, -4.5, tran3), col='tomato') })
clinesT33 = contourLines(x3, y3, Tden3, levels=seq(0.1, max(Tden3), length=6))
contourlineT3 = lapply(clinesT33, function(contour){ lines(trans3d(contour$x,contour$y, -4.5, tran3), col='red3') })

# Case 2
n = 900
p = 2
g = 2
ni=c(n/g, n/g)
w=c(0.5,0.5)
mu = matrix(c(-2,-2,2,2), p, g)
Sigma = array(NA, dim=c(p,p,g))
Sigma[,,1] = matrix(c(0.8, 0.7, 0.7, 1), 2, 2)
Sigma[,,2] = matrix(c(1.5, -0.3, -0.3, 0.8), 2, 2)


# No Truncation
Y = NULL
for(i in 1:g) Y = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low, upper=a.upp, dist="Normal")))
x = seq(min(Y[,1]), max(Y[,1]), length=50)
y = seq(min(Y[,2]), max(Y[,2]), length=50)
den2 = outer(x, y, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low, upper=a.upp)

# Right Truncation
Yt1 = NULL
for(i in 1:g) Yt1 = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low21, upper=a.upp21, dist="Normal")))
x1 = seq(min(Yt1[,1]), max(Yt1[,1]), length=50)
y1 = seq(min(Yt1[,2]), max(Yt1[,2]), length=50)
Tden21 = outer(x1, y1, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low21, upper=a.upp21)

# Left Truncation
Yt2 = NULL
for(i in 1:g) Yt2 = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low22, upper=a.upp22, dist="Normal")))
x2 = seq(min(Yt2[,1]), max(Yt2[,1]), length=50)
y2 = seq(min(Yt2[,2]), max(Yt2[,2]), length=50)
Tden22 = outer(x2, y2, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low22, upper=a.upp22)

# Double Truncation
Yt3 = NULL
for(i in 1:g) Yt3 = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low23, upper=a.upp23, dist="Normal")))
x3 = seq(min(Yt3[,1]), max(Yt3[,1]), length=50)
y3 = seq(min(Yt3[,2]), max(Yt3[,2]), length=50)
Tden23 = outer(x3, y3, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low23, upper=a.upp23)

plot(0:1,0:1, type='n',xlab='',ylab='',xaxt='n',yaxt='n', bty='n')
mtext("g=2", 1, line=-2, cex=1, font=1)
# No Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran20 = persp(x,y,den2, xlim = range(x), ylim=range(y), zlim=c(-0.105, max(den2, na.rm=T)), box=T, theta=40, lwd=.01,
               axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clines21 = contourLines(x, y, den2, levels=5e-7)
contourline = lapply(clines21, function(contour){ lines(trans3d(contour$x,contour$y, -0.105, tran20), col='coral2') })
clines21 = contourLines(x, y, den2, levels=seq(0.0001, 0.01, length=2))
contourline = lapply(clines21, function(contour){ lines(trans3d(contour$x,contour$y, -0.105, tran20), col='pink1') })
clines22 = contourLines(x, y, den2, levels=seq(0.02, 0.05, length=3))
contourline = lapply(clines22, function(contour){ lines(trans3d(contour$x,contour$y, -0.105, tran20), col='pink') })
clines23 = contourLines(x, y, den2, levels=seq(0.06, 0.1, length=4))
contourline = lapply(clines23, function(contour){ lines(trans3d(contour$x,contour$y, -0.105, tran20), col='tomato') })
clines23 = contourLines(x, y, den2, levels=seq(0.1, max(den2), length=6))
contourline = lapply(clines23, function(contour){ lines(trans3d(contour$x,contour$y, -0.105, tran20), col='red3') })

# Right Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran21 = persp(x1,y1,Tden21, xlim = range(x1)+c(-1,1), ylim=range(y1)+c(-1,1), zlim=c(-.5, max(Tden21, na.rm=T)), box=T, theta=40, lwd=.01,
               axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clines2T11 = contourLines(x1, y1, Tden21, levels=5e-7)
contourlineT1 = lapply(clines2T11, function(contour){ lines(trans3d(contour$x,contour$y, -.5, tran21), col='coral2') })
clines2T11 = contourLines(x1, y1, Tden21, levels=seq(0.0001, 0.01, length=2))
contourlineT1 = lapply(clines2T11, function(contour){ lines(trans3d(contour$x,contour$y, -.5, tran21), col='pink1') })
clines2T12 = contourLines(x1, y1, Tden21, levels=seq(0.02, 0.05, length=3))
contourlineT1 = lapply(clines2T12, function(contour){ lines(trans3d(contour$x,contour$y, -.5, tran21), col='pink') })
clines2T13 = contourLines(x1, y1, Tden21, levels=seq(0.06, 0.1, length=4))
contourlineT1 = lapply(clines2T13, function(contour){ lines(trans3d(contour$x,contour$y, -.5, tran21), col='tomato') })
clines2T13 = contourLines(x1, y1, Tden21, levels=seq(0.1, max(Tden21), length=6))
contourlineT1 = lapply(clines2T13, function(contour){ lines(trans3d(contour$x,contour$y, -.5, tran21), col='red3') })

# Left Truncation
par(mar=c(0,0,0,0))
tran22 = persp(x2,y2,Tden22, xlim = range(x2)+c(-1,1), ylim=range(y2)+c(-1,1), zlim=c(-0.5, max(Tden22, na.rm=T)), box=T, theta=40, lwd=.01,
               axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clines2T21 = contourLines(x2, y2, Tden22, levels=5e-7)
contourlineT2 = lapply(clines2T21, function(contour){ lines(trans3d(contour$x,contour$y, -0.5, tran22), col='coral2') })
clines2T21 = contourLines(x2, y2, Tden22, levels=seq(0.0001, 0.01, length=2))
contourlineT2 = lapply(clines2T21, function(contour){ lines(trans3d(contour$x,contour$y, -0.5, tran22), col='pink1') })
clines2T22 = contourLines(x2, y2, Tden22, levels=seq(0.02, 0.05, length=3))
contourlineT2 = lapply(clines2T22, function(contour){ lines(trans3d(contour$x,contour$y, -0.5, tran22), col='pink') })
clines2T23 = contourLines(x2, y2, Tden22, levels=seq(0.06, 0.1, length=4))
contourlineT2 = lapply(clines2T23, function(contour){ lines(trans3d(contour$x,contour$y, -0.5, tran22), col='tomato') })
clines2T23 = contourLines(x2, y2, Tden22, levels=seq(0.1, max(Tden22), length=6))
contourlineT2 = lapply(clines2T23, function(contour){ lines(trans3d(contour$x,contour$y, -0.5, tran22), col='red3') })

# Double Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran23 = persp(x3,y3,Tden23, xlim = range(x3)+c(-1,1), ylim=range(y3)+c(-1,1), zlim=c(-.45, max(Tden23, na.rm=T)), box=T, theta=40, lwd=.01,
               axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clines2T31 = contourLines(x3, y3, Tden23, levels=seq(0, 0.01, length=2))
contourlineT3 = lapply(clines2T31, function(contour){ lines(trans3d(contour$x,contour$y, -.45, tran23), col='coral2') })
clines2T32 = contourLines(x3, y3, Tden23, levels=seq(0.02, 0.05, length=3))
contourlineT3 = lapply(clines2T32, function(contour){ lines(trans3d(contour$x,contour$y, -.45, tran23), col='pink1') })
clines2T33 = contourLines(x3, y3, Tden23, levels=seq(0.06, 0.1, length=4))
contourlineT3 = lapply(clines2T33, function(contour){ lines(trans3d(contour$x,contour$y, -.45, tran23), col='tomato') })
clines2T33 = contourLines(x3, y3, Tden23, levels=seq(0.1, max(Tden23), length=6))
contourlineT3 = lapply(clines2T33, function(contour){ lines(trans3d(contour$x,contour$y, -.45, tran23), col='red3') })

# Case 3
n = 900
p = 2
g = 3
ni=c(n/g, n/g, n/g)
w=c(1/3,1/3,1/3)
mu = matrix(c(-3,-3,3,3,0,0), p, g)
Sigma = array(NA, dim=c(p,p,g))
Sigma[,,1] = matrix(c(2, 0.8, 0.8, 1.5), 2, 2)
Sigma[,,2] = matrix(c(1.2, -0.7, -0.7, 2), 2, 2)
Sigma[,,3] = matrix(c(1.5, 0.3, 0.3, 1), 2, 2)
Y = NULL
for(i in 1:g) Y = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low, upper=a.upp, dist="Normal")))
x = seq(min(Y[,1]), max(Y[,1]), length=50)
y = seq(min(Y[,2]), max(Y[,2]), length=50)
den3 = outer(x, y, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low, upper=a.upp)

# Right Truncation
Yt1 = NULL
for(i in 1:g) Yt1 = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low31, upper=a.upp31, dist="Normal")))
x1 = seq(min(Yt1[,1]), max(Yt1[,1]), length=50)
y1 = seq(min(Yt1[,2]), max(Yt1[,2]), length=50)
Tden31 = outer(x1, y1, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low31, upper=a.upp31)

# Left Truncation
Yt2 = NULL
for(i in 1:g) Yt2 = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low32, upper=a.upp32, dist="Normal")))
x2 = seq(min(Yt2[,1]), max(Yt2[,1]), length=50)
y2 = seq(min(Yt2[,2]), max(Yt2[,2]), length=50)
Tden32 = outer(x2, y2, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low32, upper=a.upp32)

# Double Truncation
Yt3 = NULL
for(i in 1:g) Yt3 = as.data.frame(rbind(Y, rtelliptical(ni[i] ,mu=mu[,i], Sigma=Sigma[,,i], lower=a.low33, upper=a.upp33, dist="Normal")))
x3 = seq(min(Yt3[,1]), max(Yt3[,1]), length=50)
y3 = seq(min(Yt3[,2]), max(Yt3[,2]), length=50)
Tden33 = outer(x3, y3, FUN = bivMTN, w, mu=mu, Sigma=Sigma, lower=a.low33, upper=a.upp33)

plot(0:1,0:1, type='n',xlab='',ylab='',xaxt='n',yaxt='n', bty='n')
mtext("g=3", 1, line=-2, cex=1, font=1)
# No Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran30 = persp(x,y,den3, xlim = range(x), ylim=range(y), zlim=c(-0.025, max(den3, na.rm=T)), box=T, theta=40, lwd=.01,
               axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clines31 = contourLines(x, y, den3, levels=1e-18)
contourline = lapply(clines31, function(contour){ lines(trans3d(contour$x,contour$y, -0.025, tran30), col='coral2') })
clines31 = contourLines(x, y, den3, levels=seq(0.0001, 0.01, length=2))
contourline = lapply(clines31, function(contour){ lines(trans3d(contour$x,contour$y, -0.025, tran30), col='pink1') })
clines32 = contourLines(x, y, den3, levels=seq(0.02, 0.005, length=3))
contourline = lapply(clines32, function(contour){ lines(trans3d(contour$x,contour$y, -0.025, tran30), col='pink') })
clines33 = contourLines(x, y, den3, levels=seq(0.06, 0.07, length=6))
contourline = lapply(clines33, function(contour){ lines(trans3d(contour$x,contour$y, -0.025, tran30), col='tomato') })
clines33 = contourLines(x, y, den3, levels=seq(0.07, max(den3), length=6))
contourline = lapply(clines33, function(contour){ lines(trans3d(contour$x,contour$y, -0.025, tran30), col='red3') })

# Right Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran31 = persp(x1,y1,Tden31, xlim = range(x1)+c(-1,1), ylim=range(y1)+c(-1,1), zlim=c(-0.1, max(Tden31, na.rm=T)), box=T, theta=40, lwd=.01,
               axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clines3T11 = contourLines(x1, y1, Tden31, levels=1e-18)
contourlineT1 = lapply(clines3T11, function(contour){ lines(trans3d(contour$x,contour$y, -0.1, tran31), col='coral2') })
clines3T11 = contourLines(x1, y1, Tden31, levels=seq(0.0001, 0.01, length=2))
contourlineT1 = lapply(clines3T11, function(contour){ lines(trans3d(contour$x,contour$y, -0.1, tran31), col='pink1') })
clines3T12 = contourLines(x1, y1, Tden31, levels=seq(0.02, 0.05, length=3))
contourlineT1 = lapply(clines3T12, function(contour){ lines(trans3d(contour$x,contour$y, -0.1, tran31), col='pink') })
clines3T13 = contourLines(x1, y1, Tden31, levels=seq(0.06, 0.1, length=4))
contourlineT1 = lapply(clines3T13, function(contour){ lines(trans3d(contour$x,contour$y, -0.1, tran31), col='tomato') })
clines3T13 = contourLines(x1, y1, Tden31, levels=seq(0.1, max(Tden31), length=6))
contourlineT1 = lapply(clines3T13, function(contour){ lines(trans3d(contour$x,contour$y, -0.1, tran31), col='red3') })

# Left Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran32 = persp(x2,y2,Tden32, xlim = range(x2)+c(-1,1), ylim=range(y2)+c(-1,1), zlim=c(-0.08, max(Tden32, na.rm=T)), box=T, theta=40, lwd=.01,
               axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clines3T21 = contourLines(x2, y2, Tden32, levels=1e-18)
contourlineT2 = lapply(clines3T21, function(contour){ lines(trans3d(contour$x,contour$y, -0.08, tran32), col='coral2') })
clines3T21 = contourLines(x2, y2, Tden32, levels=seq(0.0001, 0.01, length=2))
contourlineT2 = lapply(clines3T21, function(contour){ lines(trans3d(contour$x,contour$y, -0.08, tran32), col='pink1') })
clines3T22 = contourLines(x2, y2, Tden32, levels=seq(0.02, 0.05, length=3))
contourlineT2 = lapply(clines3T22, function(contour){ lines(trans3d(contour$x,contour$y, -0.08, tran32), col='pink') })
clines3T23 = contourLines(x2, y2, Tden32, levels=seq(0.06, 0.1, length=4))
contourlineT2 = lapply(clines3T23, function(contour){ lines(trans3d(contour$x,contour$y, -0.08, tran32), col='tomato') })
clines3T23 = contourLines(x2, y2, Tden32, levels=seq(0.1, max(Tden32), length=6))
contourlineT2 = lapply(clines3T23, function(contour){ lines(trans3d(contour$x,contour$y, -0.08, tran32), col='red3') })

# Double Truncation
par(mar=c(0,0,0,0), lwd=0.2)
tran33 = persp(x3,y3,Tden33, xlim = range(x3)+c(-1,1), ylim=range(y3)+c(-1,1), zlim=c(-0.1, max(Tden33, na.rm=T)), box=T, theta=40, lwd=.01,
               axes=T, expand=1.2, col='salmon', phi=7, border="pink1", xlab='y1', ylab='y2', zlab='Density', main = '')
clines3T31 = contourLines(x3, y3, Tden33, levels=seq(0.001, 0.02, length=2))
contourlineT3 = lapply(clines3T31, function(contour){ lines(trans3d(contour$x,contour$y, -0.1, tran33), col='coral2') })
clines3T32 = contourLines(x3, y3, Tden33, levels=seq(0.02, 0.05, length=3))
contourlineT3 = lapply(clines3T32, function(contour){ lines(trans3d(contour$x,contour$y, -0.1, tran33), col='pink1') })
clines3T33 = contourLines(x3, y3, Tden33, levels=seq(0.06, 0.1, length=4))
contourlineT3 = lapply(clines3T33, function(contour){ lines(trans3d(contour$x,contour$y, -0.1, tran33), col='tomato') })
clines3T33 = contourLines(x3, y3, Tden33, levels=seq(0.1, max(Tden23), length=6))
contourlineT3 = lapply(clines3T33, function(contour){ lines(trans3d(contour$x,contour$y, -0.1, tran33), col='red3') })
dev.off()
