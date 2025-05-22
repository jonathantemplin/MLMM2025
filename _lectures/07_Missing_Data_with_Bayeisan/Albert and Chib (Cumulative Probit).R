# File: Albert and Chib (Ordered Probit).r
# Generate Ordered Categorical Data and run Gibbs Sampler Using A&C Section 4.1
# CUMULATIVE PROBIT MODEL: pr(class j)=phi(int_j - x*beta) 
# 		Note: higher beta --> more likely to be in higher class
# NOTE: May require to converge depending on choice of intercept values (more extreme take longer to converge)
# August 28, 2015
##################################
library(mvtnorm)
library(MASS)  				# For polr function

# Simulate data
 set.seed(250)
 n<-1000				      # Sample Size
 int1<--1				      # Intercept 1 has to be less than int2
 int2<-.5             # Intecept 2 has to be greater than int1
 beta<--.75		        # Regression Coefs 
 x<-sample(0:4,n,replace=T)	# Covariate
 eta<-x*beta				 
 
 z<-rnorm(n,eta,1)
 y<-rep(0,n)
 y[z<=int1]<-1
 y[int1<z& z<=int2]<-2
 y[int2<z]<-3
 ns<-table(y)	 	        # Sample size in each category

# Priors on Beta=c(alpha, beta)
  beta0<-0			    # Prior Mean for Beta
  vbeta0<-100		    # Prior Cov of Beta (vague)

#Initial Values
  int1<-int2<-0
  beta<-0

# Create vectors to store results
  nsim<-10000
  z<-rep(0,n)				# Latent Normal Variable
  Beta<-rep(0,nsim)
  Int1<-Int2<-rep(0,nsim)

# Posterior Variance of Beta (Can update outside of sampler)
  prec0<-1/vbeta0
  vbeta<-1/(prec0+sum(x^2))

###################
# GIBBS SAMPLER	  #
###################
tmp<-proc.time()
for (i in 1:nsim) {

 # Draw Latent Variable, z, from its full conditional, given y
   muz<-x*beta							# Update Mean of Z
   
 # Update z using inverse CDF method (could also use tnorm function in msm package)
   z[y==1]<-qnorm(runif(n,0,pnorm(int1,muz,1)),muz,1)[y==1]
   z[y==2]<-qnorm(runif(n,pnorm(int1,muz,1),pnorm(int2,muz,1)),muz,1)[y==2]
   z[y==3]<-qnorm(runif(n,pnorm(int2,muz,1),1),muz,1)[y==3]
   
 # Update intercepts
   Int1[i]<-int1<-runif(1,max(z[y==1]),min(min(z[y==2]),int2))  
   Int2[i]<-int2<-runif(1,max(max(z[y==2]),int1),min(z[y==3]))  

 # Posterior mean of Beta|Z
   mbeta <- vbeta%*%(prec0%*%beta0+crossprod(x,z))

 # Draw Beta=c(alpha,beta) from its full conditional
   Beta[i]<-beta<-rnorm(1,mbeta,sqrt(vbeta))
  
   if(i%%100==0){
     print(i)
    # print(int1)
    # print(int2)
    # print(beta)
   }
 }

proc.time()-tmp

# Get Summaries
mean.beta<-mean(Beta[5001:nsim])
mean.int1<-mean(Int1[5001:nsim])
mean.int2<-mean(Int2[5001:nsim])

# Compare Estimates
polr(as.factor(y)~x,method="probit")
c(mean.int1,mean.int2,mean.beta)

# Trace Plots
par(mfrow=c(3,1))
plot(1:nsim,Int1[1:nsim],type="l",col="lightgreen")
abline(h=mean.int1,col="blue4")
plot(1:nsim,Int2[1:nsim],type="l",col="lightgreen")
abline(h=mean.int2,col="blue4")
plot(1:nsim,Beta[1:nsim],type="l",col="lightgreen")
abline(h=mean.beta,col="blue4")
