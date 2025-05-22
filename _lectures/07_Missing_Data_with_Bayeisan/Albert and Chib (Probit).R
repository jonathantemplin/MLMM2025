# File: Albert and Chib (Probit).r
# Generate Binary Data and run Gibbs Sampler
# Assumes conjugate normal prior on Beta
# Albert and Chib (1993) section 3
# Nov 15, 2007
##################################
library(mvtnorm)
library(msm)          			# For rtnorm function

# Simulate data
set.seed(250)
n<-1000					# Sample Size
x<-sample(0:4,n,replace=T)		# Covariates
X<-matrix(c(rep(1,n),x), ncol=2)   	# Design matrix
tbeta<-c(-.75,.75)          		# True beta
p<-pnorm(X%*%tbeta) 	        	# Prob. Vector
y<-rbinom(n,1,p)				# Observations
k<-2 						# Number of parms
n1<-sum(y)                  		# Number of successes
n0<-n-n1

# Priors on Beta=c(alpha, beta)
  beta0<-c(0,0)				# Prior Mean for Beta
  vbeta0<-diag(10,2)			# Prior Cov of Beta (vague)

# Initial Values
  fit<-glm(y~x,family=binomial(link=probit))
  Beta<-rep(0,k)

# Create vectors to store results
  nsim<-1000  				# Number of Iterations of Gibbs Sampler
  z<-rep(0,n)				# Latent Normal Variable
  Betamat<-matrix(0,nrow=nsim,ncol=k) 	# Store Results

###################
# GIBBS SAMPLER	#
###################
 # Posterior Variance of Beta
	prec0<-solve(vbeta0)
	vbeta<-solve(prec0+crossprod(X,X))

for (i in 2:nsim) {

 # Draw Latent Variable, z, from its full conditional, given y
  muz<-X%*%Beta			# Update Mean of Z

  # z[y==0]<-qnorm(runif(n,0,pnorm(0,muz,1)),muz,1)[y==0]  # using inverse-CDF method
  # z[y==1]<-qnorm(runif(n,pnorm(0,muz,1),1),muz,1)[y==1]

  z[y==0]<-rtnorm(n0,muz[y==0],1,-Inf,0)  # Using truncated normal function -- seems to give similar results
  z[y==1]<-rtnorm(n1,muz[y==1],1,0,Inf)

 # Posterior mean of Beta|Z
   mbeta <- vbeta%*%(prec0%*%beta0+crossprod(X,z))
   Betamat[i,]<-Beta<-c(rmvnorm(1,mbeta,vbeta))

  if (i%%100==0) print(i)
 }

# Get Summaries
apply(Betamat[501:nsim,],2,mean)

plot(1:nsim,Betamat[,2], type="l", col="lightgreen")
abline(h=mean(Betamat[501:nsim,2]),col="blue4")

### Modifying to accommodate missing data


# create a set of missing cases
y[1:10] = NA

# Initial Values
fit<-glm(y~x,family=binomial(link=probit))
Beta<-rep(0,k)

# Create vectors to store results
nsim<-1000  				# Number of Iterations of Gibbs Sampler
z<-rep(0,n)				# Latent Normal Variable
Betamat<-matrix(0,nrow=nsim,ncol=k) 	# Store Results
zmat = matrix(0, nrow=nsim, ncol=10)  # Store missing imputations

###################
# GIBBS SAMPLER	#
###################
# Posterior Variance of Beta
prec0<-solve(vbeta0)
vbeta<-solve(prec0+crossprod(X,X))


n1<-length(which(y==1))                 		# Number of successes
n0<-length(which(y==0))
nNA = length(which(is.na(y)))


for (i in 2:nsim) {
  
  # Draw Latent Variable, z, from its full conditional, given y
  muz<-X%*%Beta			# Update Mean of Z
  
  z[which(y==0)]<-rtnorm(n0,muz[which(y==0)],1,-Inf,0)  # Using truncated normal function -- seems to give similar results
  z[which(y==1)]<-rtnorm(n1,muz[which(y==1)],1,0,Inf)
  z[which(is.na(y))]<-rnorm(nNA,muz[which(is.na(y))],1)
  
  # Posterior mean of Beta|Z
  mbeta <- vbeta%*%(prec0%*%beta0+crossprod(X,z))
  Betamat[i,]<-Beta<-c(rmvnorm(1,mbeta,vbeta))
  zmat[i,] = z[which(is.na(y))]
  
  if (i%%100==0) print(i)
}

# Get Summaries
apply(Betamat[501:nsim,],2,mean)

plot(1:nsim,Betamat[,2], type="l", col="lightgreen")
abline(h=mean(Betamat[501:nsim,2]),col="blue4")

matplot(1:nsim,zmat, type="l", col=1:10)
