#-----------------------------------------------------------------------------------------------------------

haspackage = require("EPSY905R")
if (haspackage==FALSE){
  if (!require("devtools")) install.packages("devtools")
  devtools::install_github("jonathantemplin/EPSY905R")
}
library("EPSY905R")

haspackage = require("fields")
if (haspackage==FALSE){
  install.packages("fields")
}
library("fields")

haspackage = require("mnormt")
if (haspackage==FALSE){
  install.packages("mnormt")
}
library("mnormt")


haspackage = require("nlme")
if (haspackage==FALSE){
  install.packages("nlme")
}
library("nlme")

#-----------------------------------------------------------------------------------------------------------

#building normal distribution for single observation
mu_x = 114.4
sigma2_x = 5.29

x=seq(100,130,.1)
y=dnorm(x,mu_x,sqrt(sigma2_x))

plot(x,y,type="l",xlab="IQ",ylab="f(x)",cex=5,lwd=5)

xloc = 114.4
yloc = dnorm(xloc,mu_x,sqrt(sigma2_x))
lines(c(xloc,xloc),c(-1,yloc),cex=5,lwd=5,lty=2)
lines(c(00,xloc),c(yloc,yloc),cex=5,lwd=5,lty=2)

xloc=110
yloc = dnorm(xloc,mu_x,sqrt(sigma2_x))
lines(c(xloc,xloc),c(-1,yloc),cex=5,lwd=5,lty=3)
lines(c(00,xloc),c(yloc,yloc),cex=5,lwd=5,lty=3)

######################################################################
#plot for likelihood for mean of one observation
mu = seq(100,124,.01)
x = 112
sigma2 = 5.29
y = dnorm(x,mu,sqrt(sigma2))

plot(mu,y,type="l",xlab=expression(mu[x]),ylab=expression(paste("L(",mu[x],"|x,",sigma[x]^2,")")),cex=3,lwd=5,main=expression(paste("Likelihood of ",mu[x])))
muloc=112
yloc=dnorm(x,muloc,sqrt(sigma2))
lines(c(muloc,muloc),c(-1,yloc),cex=5,lwd=5,lty=2)
lines(c(00,muloc),c(yloc,yloc),cex=5,lwd=5,lty=2)


######################################################################
#plot for sample likelihood for mean
perf = matrix(c(10,12,14,16,12))
iq = c(112,113,115,118,114)

#for the likelihood (without a log)
mu = seq(100,124,.1)
mudim = dim(matrix(mu))[1]
y = mu
sigma2_x = 5.29

for (i in 1:mudim) {
  likelihood = 1
  for (j in 1:5) {
    likelihood=likelihood*dnorm(iq[j],mu[i],sqrt(sigma2_x))
  }
  y[i]=likelihood
}

#maximum likelihood:
maxL = max(y)

#at location:
MLEmu = mu[which.max(y)]

plot(mu,y,type="l",xlab=expression(mu[x]),ylab=expression(paste("L(",mu[x],"|x,",sigma[x]^2,")")),lwd=5)
lines(c(MLEmu,MLEmu),c(-1,maxL),lwd=5,lty=2)
lines(c(90,MLEmu),c(maxL,maxL),lwd=5,lty=2)

#for the log-likelihood (with a log)
mu = seq(100,124,.1)
mudim = dim(matrix(mu))[1]
y = mu

for (i in 1:mudim) {
  likelihood = 0
  for (j in 1:5) {
    likelihood=likelihood+dnorm(iq[j],mu[i],sqrt(sigma2_x),log=TRUE)
  }
  y[i]=likelihood
}
#maximum likelihood:
maxLogL = max(y)

#at location:
MLEmu = mu[which.max(y)]

plot(mu,y,type="l",xlab=expression(mu),ylab="logL",lwd=5)
lines(c(MLEmu,MLEmu),c(-125,maxLogL),lwd=5,lty=2)
lines(c(90,MLEmu),c(maxLogL,maxLogL),lwd=5,lty=2)

#plotting tangent line at maximum
plot(mu,y,type="l",xlab=expression(mu),ylab="logL",lwd=5)
lines(c(MLEmu-5,MLEmu+5),c(maxLogL,maxLogL),lwd=5,lty=3)

##############################################################################################################
#section for plotting the two dimensional likelihood Surface

iq = matrix(c(112,113,115,118,114))
n = dim(matrix(iq))[1]
mu = matrix(seq(114.3,114.5,.01))
nstep = dim(matrix(mu))[1]

sigma = matrix(seq(2.0581,2.0601,.0001))
dim(sigma)
dim(mu)

x=matrix(0,nstep,1)
y=matrix(0,nstep,1)
z=matrix(0,nstep,nstep)

maxval = -999999
maxmu = 0
maxsigma = 0
for (i in 1:nstep){
  for (j in 1:nstep){
    x[i]=mu[i]
    y[j]=sigma[j]^2
    z[i,j]=0
    
    for (k in 1:n) {
      z[i,j] = z[i,j]+dnorm(iq[k],mu[i],sigma[j],log=TRUE)
    }
    
    if (z[i,j]>maxval) {
      maxval = z[i,j]
      maxmu=mu[i]
      maxsigma=sigma[j]
    }
    
  }
  
}

maxmu
maxsigma

grid.list=list(x=x,y=y)
mygrid=make.surface.grid(grid.list)
out=list(x=grid.list$x,y=grid.list$y,z=z)
par(mfrow=c(1,2))
plot.surface(out,type="p",xlab="Mean",ylab="Variance",zlab="Log L",main="Log Likelihood Function")
plot.surface(out,type="c",xlab="test",ylab="Slope",zlab="SSE",main="Log Likelihood Function")
max(z)


#ML Analyses =================================


data01 = data.frame(iq = iq, perf = perf)

#empty model predicting IQ
model01 = gls(iq~1,data=data01,method="ML")
summary(model01)
names(model01)
anova(model01)

#emtpy model predicting performance
model02 = gls(perf~1,data=data01,method="ML")
summary(model02)

#centering IQ at mean of 114.4
data01$iq114 = data01$iq-114.4

#Regression with ML:
model03a = gls(perf~iq114,data=data01,method="ML")
summary(model03a)

lrt02v03 = -2*(model02$logLik-model03a$logLik)

lrt02v03
pchisq(lrt02v03,df=1,lower.tail=FALSE)

anova(model02,model03a)

(model02$sigma^2-model03a$sigma^2)/model02$sigma^2

#note the residual standard error is off
model03b = lm(perf~iq114,data=data01) 
summary(model03b)

#Regression with ML:
model03c = gls(perf~iq114,data=data01,method="REML")
summary(model03c)

