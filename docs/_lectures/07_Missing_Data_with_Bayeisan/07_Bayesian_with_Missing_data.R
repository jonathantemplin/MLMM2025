rm(list = ls())
# =====================================================================================================================
if (!require(mvtnorm)){
  install.packages("mvtnorm")
}
library(mvtnorm)

haspackage = require("lavaan")
if (haspackage==FALSE){
  install.packages("lavaan")
}
library(lavaan)

haspackage = require("R2jags")
if (haspackage==FALSE){
  install.packages("R2jags")
}
library(R2jags)

# simulated data for today's example =================================================================================================
set.seed(20250223)

Nmen = 121
Nwomen = 229

MeanMEN = c(5.081333347, 9.723030258, 51.94404048, 52.66452548, 34.04606078, 77.06181845, 14.75389904)
MeanWOMEN = c(4.80418631, 10.60486174, 50.34834542, 48.13359134, 30.61321679, 71.77082955, 13.75449003)

Corr = matrix(c(
  1.0,  .15, .12, .48, .44, .47, .44,
  .15, 1.0, .06, .25, .20, .23, .23,
  .12, .06, 1.0, .40, .32, .19, .14,
  .48, .25, .40, 1.0, .87, .61, .54,
  .44, .20, .32, .87, 1.0, .56, .51,
  .47, .23, .19, .61, .56, 1.0, .70,
  .44, .23, .14, .54, .51, .70, 1.0
), nrow = 7, byrow = TRUE)

SD = c(1.2, 6.0, 15.2, 16.6, 10.9, 10.5, 2.8)
SDdiag = diag(SD)

# Create the covariance matrix
Cov = SDdiag %*% Corr %*% SDdiag

# Generate multivariate normal data
xmen = rmvnorm(Nmen, MeanMEN, Cov)
xwomen = rmvnorm(Nwomen, MeanWOMEN, Cov)

# Add group identifiers
xmen = cbind(0, xmen)  # 1 for men
xwomen = cbind(1, xwomen)  # 0 for women

# Combine datasets
allx = as.data.frame(rbind(xmen, xwomen))

# Convert to dataframe and rename columns
names(allx) = c("female", "hsl", "cc", "use", "msc", "mas", "mse", "perf")


# Introduce missingness in CC, PERF, and USE based on remaining variables
M_CC_intercept = -3
M_CC_female = 1
M_CC_msc = .01
M_CC_mas = .01
M_CC_mse = .01

missingCCmodelLogit = 
  M_CC_intercept +
  M_CC_female*allx[,"female"] +
  M_CC_msc*allx[,"msc"] +
  M_CC_mas*allx[,"mas"] +
  M_CC_mse*allx[,"mse"]

missingCCmodelProb = exp(missingCCmodelLogit)/(1+exp(missingCCmodelLogit))

# create missing data
makeMissingCC = which(runif(Nmen+Nwomen) < missingCCmodelProb)
allx$cc[makeMissingCC] = NA
allx$cc

# next for perf
M_perf_intercept = -3
M_perf_female = .5
M_perf_msc = .001
M_perf_mas = .02
M_perf_mse = .01

missingPERFmodelLogit = 
  M_perf_intercept +
  M_perf_female*allx[,"female"] +
  M_perf_msc*allx[,"msc"] +
  M_perf_mas*allx[,"mas"] +
  M_perf_mse*allx[,"mse"]

missingPERFmodelProb = exp(missingPERFmodelLogit)/(1+exp(missingPERFmodelLogit))

# create missing data
makeMissingPerf = which(runif(Nmen+Nwomen) < missingPERFmodelProb)
allx$perf[makeMissingPerf] = NA
allx$perf

# next for use
M_use_intercept = -3
M_use_female = .5
M_use_msc = .001
M_use_mas = .02
M_use_mse = .01

missingUSEmodelLogit = 
  M_use_intercept +
  M_use_female*allx[,"female"] +
  M_use_msc*allx[,"msc"] +
  M_use_mas*allx[,"mas"] +
  M_use_mse*allx[,"mse"]

missingUSEmodelProb = exp(missingUSEmodelLogit)/(1+exp(missingUSEmodelLogit))

# create missing data
makeMissingUse = which(runif(Nmen+Nwomen) < missingUSEmodelProb)
allx$use[makeMissingUse] = NA
allx$use



data02 = as.data.frame(allx)

#centering CC at 10: ===================================================================================================
par(mfrow = c(1,2))
boxplot(data02$cc, main="Boxplot of College Experience (CC)")
hist(data02$cc, main = "Histogram of College Experience (CC)", xlab = "CC")
par(mfrow = c(1,1))

data02$cc10 = data02$cc - 10

#creating interation between female and cc: ============================================================================
data02$femXcc10 = data02$female*data02$cc10


newData = data02[c("perf", "use", "female", "cc10")]

# determine if any observations are missing on all variables: 
allNA = apply(newData, 1, function(x) all(is.na(x)))
which(allNA)

# none are missing, we proceed

#Building Analysis Model #1: an empty model ============================================================================
# here, female and cc10 both predict perf and use in separate equations, estimated simultaneously

# let's first get a count of how many observations are missing both perf and use
sum(is.na(newData$perf) & is.na(newData$use))

# that leaves this many observations that should be in the model:
sum(!is.na(newData$perf) | !is.na(newData$use))

#analysis syntax
model01.lavaan.syntax = "

#Means:
perf ~ 1
use  ~ 1

#Variances:
perf ~~ perf
use  ~~ use

#Covariance:
perf ~~ use

"
  

#analysis estimation
model01.lavaan.fit = lavaan(model01.lavaan.syntax, data=data02, estimator = "MLR", mimic = "MPLUS")
summary(model01.lavaan.fit, standardized=TRUE, fit.measures=TRUE)

# R2jags function with Wishart prior for the inverse covariance matrix using the R2jags library

model01a.jags.syntax <- function(){

    # Priors for means
    mu.perf ~ dnorm(0, 0.001)
    mu.use ~ dnorm(0, 0.001)

    # Prior for inverse covariance matrix (precision matrix)
    Omega[1:2, 1:2] ~ dwish(R, df)
    Sigma <- inverse(Omega) # Covariance matrix

    # Likelihood
    for (i in 1:n) {
      Y[i, 1:2] ~ dmnorm(mu[1:2], Omega[1:2,1:2])
    }

    # Means vector
    mu[1] <- mu.perf
    mu[2] <- mu.use

    # Derived parameters (variances and covariance)
    variance.perf <- Sigma[1, 1]
    variance.use <- Sigma[2, 2]
    covariance <- Sigma[1, 2]
    correlation <- covariance / (sqrt(variance.perf) * sqrt(variance.use))
}


# Prepare data for JAGS
model01a.jags.data = list(
  Y = data02[c("perf","use")],
  n = nrow(data02),
  df = 2, # Degrees of freedom for Wishart (should be >= dimension of matrix)
  R = diag(2) # Scale matrix for Wishart (identity matrix)
) 
  
# Parameters to monitor
model01a.jags.params = c("mu.perf", "mu.use", "variance.perf", "variance.use", "covariance", "correlation")

# Run the model using the R2jags package
model01.jags.samples = jags(
  data = model01a.jags.data,
  parameters.to.save = model01a.jags.params,
  model.file = model01a.jags.syntax,
  n.chains = 5,
  n.iter = 6000, # R2jags includes burn-in in n.iter
  n.burnin = 1000,
  n.thin = 1
)

# JAGS problem: Multivariate distributions cannot contain missing data 

## What things would look like if there were complete cases:
data03 = data02[c("perf", "use")]
data03 = data03[complete.cases(data03),]

# Prepare data for JAGS
model01b.jags.data = list(
  Y = data03,
  n = nrow(data03),
  df = 2, # Degrees of freedom for Wishart (should be >= dimension of matrix)
  R = diag(2) # Scale matrix for Wishart (identity matrix)
) 

# Run the model using the R2jags package
model01b.jags.samples = jags(
  data = model01b.jags.data,
  parameters.to.save = model01a.jags.params,
  model.file = model01a.jags.syntax,
  n.chains = 4,
  n.iter = 6000, # R2jags includes burn-in in n.iter
  n.burnin = 1000,
  n.thin = 1,
)

# examine convergence, first with Rhat:
model01b.jags.samples

# next with a plot of the chain:
plot(as.mcmc(model01b.jags.samples))

# extracting MCMC chain from R2jags object
model01b.chainList = mcmc.list(
  lapply(X = 1:dim(model01b.jags.samples$BUGSoutput$sims.array)[2], FUN = function(x) {
    return(as.mcmc(model01b.jags.samples$BUGSoutput$sims.array[,x,]))
  })
)

# gathering autocorrelation information (from coda package, which loads with R2jags)
coda::autocorr.diag(model01b.chainList)

# plotting autocorrelation information
coda::autocorr.plot(model01b.chainList)

#analysis estimation
model01.lavaan.fit = lavaan(model01.lavaan.syntax, data=data03, estimator = "MLR", mimic = "MPLUS")
summary(model01.lavaan.fit, standardized=TRUE, fit.measures=TRUE)


# now back to missing data...there are two ways to handle missing data in JAGS:
# (1) omit missing cases from the likelihood function (same likelihood as maximum likelihood)
# (2) have JAGS impute data using Bayes' Theorem
# both have to use a factored regression approach as JAGS won't be able to run multivariate distributions with partially missing data



# f(perf,use) = f(perf|use)f(use)
# compare with factored regression specification in lavaan: perf ~ use and use ~ 1 --> same log likelihood
model01c.lavaan.syntax = "

# exogenous mean:
use ~ 1

# exogenous variance:
use ~~ use

# regression:
perf ~ 1 + use

# residual variance:
perf ~~ perf
"

model01c.lavaan.fit = lavaan(model01c.lavaan.syntax, data=data02, estimator = "MLR", mimic = "MPLUS")
summary(model01c.lavaan.fit, standardized=TRUE, fit.measures=TRUE)

# examining gamma priors on precision parameters
plot(
  x = seq(.001,20,.001),
  y = dgamma(x = seq(.001,20,.001), shape = .0001, rate = .0001),
  type = "l"
)



# Factored regression; Here we can model perf|use and use:
model01c.jags.syntax <- function(){
  
  # model for use: -----------------
  
  for (obs in 1:N){
    # use mean model: 
    muUse[obs] <- use.intercept
    
    # use model data distribution
    use[obs] ~ dnorm(muUse[obs], tauUse)  
  }
  
  # use prior distributions
  use.intercept ~ dnorm(use.intercept.mean0, use.intercept.precision0)
  tauUse ~ dgamma(use.tau.shape0, use.tau.rate0)
  
  # use variance
  sigma2Use <- 1/tauUse
  
  # model for perf: -----------------
  
  # perf model data distribution
  for (obs in 1:N){
    # perf mean model
    
    muPerf[obs] <- perf.intercept + perf.use*(use[obs])
    perf[obs] ~ dnorm(muPerf[obs], tauPerf)    
  }
  

  # perf prior distributions
  perf.intercept ~ dnorm(perf.intercept.mean0, perf.intercept.precision0)
  perf.use ~ dnorm(perf.use.mean0, perf.use.precision0)
  tauPerf ~ dgamma(perf.tau.shape0, perf.tau.rate0)
  
  # perf residual variance
  sigma2Perf <- 1/tauPerf
  
  
}

# Prepare data for JAGS -- here perf and use are separate variables
model01c.jags.data = list(
  use = data02$use,
  perf = data02$perf,
  N = nrow(data02),
  use.intercept.mean0 = 0,
  use.intercept.precision0 = 1/1000,
  use.tau.shape0 = 1,
  use.tau.rate0 = 1,
  perf.intercept.mean0 = 0,
  perf.intercept.precision0 = 1/1000,
  perf.use.mean0 = 0,
  perf.use.precision0 = 1/1000,
  perf.tau.shape0 = 1,
  perf.tau.rate0 = 1
) 

model01c.jags.params = c(
  "use.intercept", 
  "sigma2Use", 
  "perf.intercept", 
  "perf.use", 
  "sigma2Perf",
  "perf[3]",
  "use[8]"
)

# Run the model using the R2jags package
model01c.jags.samples = jags(
  data = model01c.jags.data,
  parameters.to.save = model01c.jags.params,
  model.file = model01c.jags.syntax,
  n.chains = 4,
  n.iter = 6000, # R2jags includes burn-in in n.iter
  n.burnin = 1000,
  n.thin = 1
)


model01c.jags.samples

plot(as.mcmc(model01c.jags.samples))

# extracting MCMC chain from R2jags object
model01c.chainList = mcmc.list(
  lapply(X = 1:dim(model01c.jags.samples$BUGSoutput$sims.array)[2], FUN = function(x) {
    return(as.mcmc(model01c.jags.samples$BUGSoutput$sims.array[,x,]))
  })
)

# gathering autocorrelation information (from coda package, which loads with R2jags)
coda::autocorr.diag(model01c.chainList)

# plotting autocorrelation information
coda::autocorr.plot(model01c.chainList)


# second model--add predictors as main effects only: ===================================================================

# we will create missing data on female (missing completely at random for now)
data02 = data02
data02$female[sample(x = 1:nrow(data02), replace = TRUE, size = .1*nrow(data02) )] = NA
data02$female

#analysis syntax
model02a.lavaan.syntax = "

#Means:
perf ~ 1 + female + cc10
use  ~ 1 + female + cc10

#Variances:
perf ~~ perf
use  ~~ use

#Covariance:
perf ~~ use

"

#analysis estimation
model02.fit = lavaan(model02a.lavaan.syntax, data=data02, estimator = "MLR", mimic = "MPLUS")

#analysis summary (note the additional terms: standardized = TRUE for standardized estimates and fit.measures=TRUE for model fit indices)
summary(model02.fit, standardized=TRUE, fit.measures=TRUE)

# here, in JAGS, we will have the same issue--that missing data will cause a problem

model02a.jags.syntax = function(){
  
  # model for use: -----------------
  
  for (obs in 1:N){
    # use mean model: 
    muUse[obs] <- use.intercept + use.female*female[obs] + use.cc10*cc10[obs]
    
    # use model data distribution
    use[obs] ~ dnorm(muUse[obs], tauUse)  
  }
  
  # use prior distributions -- empirical prior for regression parameters
  use.intercept ~ dnorm(useMean0, usePrecision0)
  use.female ~ dnorm(useMean0, usePrecision0)
  use.cc10 ~ dnorm(useMean0, usePrecision0)
  
  useMean0 ~ dnorm(0, 1/1000)
  usePrecision0 ~ dgamma(1, 1)
  
  tauUse ~ dgamma(use.tau.shape0, use.tau.rate0)
  
  # use variance
  sigma2Use <- 1/tauUse
  
  # model for perf: -----------------
  
  # perf model data distribution
  for (obs in 1:N){
    # perf mean model
    
    muPerf[obs] <- perf.intercept + perf.female*female[obs] + perf.cc10*cc10[obs]
    perf[obs] ~ dnorm(muPerf[obs], tauPerf)    
  }
  

  # perf prior distributions
  perf.intercept ~ dnorm(perfMean0, perfPrecision0)
  perf.female ~ dnorm(perfMean0, perfPrecision0)
  perf.cc10 ~ dnorm(perfMean0, perfPrecision0)
  
  perfMean0 ~ dnorm(0, 1/1000)
  perfPrecision0 ~ dgamma(1, 1)
  
  tauPerf ~ dgamma(perf.tau.shape0, perf.tau.rate0)
  
  # perf residual variance
  sigma2Perf <- 1/tauPerf
  

}

# Prepare data for JAGS -- here perf and use are separate variables
model02a.jags.data = list(
  use = data02$use,
  perf = data02$perf,
  female = data02$female,
  cc10 = data02$cc10,
  N = nrow(data02),
  use.tau.shape0 = 1,
  use.tau.rate0 = 1,
  perf.tau.shape0 = 1,
  perf.tau.rate0 = 1
) 

model02a.jags.params = c(
  "use.intercept", 
  "use.female",
  "use.cc10",
  "sigma2Use", 
  "perf.intercept", 
  "perf.female",
  "perf.cc10",
  "sigma2Perf",
  "perf[3]",
  "use[8]"
)

# Run the model using the R2jags package
model02a.jags.samples = jags(
  data = model02a.jags.data,
  parameters.to.save = model02a.jags.params,
  model.file = model02a.jags.syntax,
  n.chains = 4,
  n.iter = 6000, # R2jags includes burn-in in n.iter
  n.burnin = 1000,
  n.thin = 1
)

# problem...predictors have missing data...solution: enter into the likelihood
# f(perf, use, female, cc10) = f(perf|female, cc10)f(use|female, cc10)f(female)f(cc10)\
# f(perf, use, female, cc10) = f(perf, use | female, cc10)f(female)f(cc10)
# f(perf, use, female, cc10) = f(perf|female, cc10, use)f(use|female, cc10)f(female|cc10)f(cc10)
# f(female, cc10) = 
model02b.jags.syntax = function(){
  
  # model for female -------------
  for (obs in 1:N){
    female[obs] ~ dbern(pFemale)  
  }
  
  # priors for female
  pFemale ~ dbeta(1, 1)
  varFemale = pFemale*(1-pFemale) # for standardized coefficents
  
  # model for cc10
  for (obs in 1:N){
    cc10[obs] ~ dnorm(muCC10, tauCC10)  
  }
  
  # priors for cc10
  muCC10 ~ dnorm(0, 1/1000)
  tauCC10 ~ dgamma(1, 1)
  
  varCC10 = 1/tauCC10 # for standardized coefficients
  
  # model for use: -----------------
  
  for (obs in 1:N){
    # use mean model: 
    muUse[obs] <- use.intercept + use.female*female[obs] + use.cc10*cc10[obs]
    
    # use model data distribution
    use[obs] ~ dnorm(muUse[obs], tauUse)  
  }
  
  # use prior distributions -- empirical prior for regression parameters
  use.intercept ~ dnorm(useMean0, usePrecision0)
  use.female ~ dnorm(useMean0, usePrecision0)
  use.cc10 ~ dnorm(useMean0, usePrecision0)
  
  useMean0 ~ dnorm(0, 1/1000)
  usePrecision0 ~ dgamma(1, 1)
  useSigma2_0 = 1/usePrecision0
  
  tauUse ~ dgamma(use.tau.shape0, use.tau.rate0)
  
  # use variance
  sigma2Use <- 1/tauUse
  
  # create standardized coefficients
  varUse = use.female^2*varFemale + use.cc10^2*varCC10 + sigma2Use
  
  use.femaleS = use.female*sqrt(varFemale)/sqrt(varUse)
  use.cc10S = use.cc10*sqrt(varCC10)/sqrt(varUse)
  
  # model for perf: -----------------
  
  # perf model data distribution
  for (obs in 1:N){
    # perf mean model
    muPerf[obs] <- perf.intercept + perf.female*female[obs] + perf.cc10*cc10[obs]
    perf[obs] ~ dnorm(muPerf[obs], tauPerf)    
  }
  
  
  # perf prior distributions
  perf.intercept ~ dnorm(perfMean0, perfPrecision0)
  perf.female ~ dnorm(perfMean0, perfPrecision0)
  perf.cc10 ~ dnorm(perfMean0, perfPrecision0)
  
  perfMean0 ~ dnorm(0, 1/1000)
  perfPrecision0 ~ dgamma(1, 1)
  perfSigma2_0 = 1/perfPrecision0
  
  tauPerf ~ dgamma(perf.tau.shape0, perf.tau.rate0)
  
  # perf residual variance
  sigma2Perf <- 1/tauPerf
  
  # create standardized coefficients
  varPerf = use.female^2*varFemale + use.cc10^2*varCC10 + sigma2Perf
  
  perf.femaleS = perf.female*sqrt(varFemale)/sqrt(varUse)
  perf.cc10S = perf.cc10*sqrt(varCC10)/sqrt(varUse)
  
  
}

# Prepare data for JAGS -- here perf and use are separate variables
model02b.jags.data = list(
  use = data02$use,
  perf = data02$perf,
  female = data02$female,
  cc10 = data02$cc10,
  N = nrow(data02),
  use.tau.shape0 = 1,
  use.tau.rate0 = 1,
  perf.tau.shape0 = 1,
  perf.tau.rate0 = 1
) 

model02b.jags.params = c(
  "pFemale",
  "muCC10",
  "tauCC10",
  "use.intercept", 
  "use.female",
  "use.cc10",
  "useMean0",
  "useSigma2_0",
  "sigma2Use", 
  "perf.intercept", 
  "perf.female",
  "perf.cc10",
  "sigma2Perf",
  "perfMean0",
  "perfSigma2_0",
  "perf[3]",
  "use[8]",
  "female[5]",
  "cc10[6]",
  "varFemale",
  "varCC10",
  "varUse",
  "use.femaleS",
  "use.cc10S",
  "varPerf",
  "perf.femaleS",
  "perf.cc10S"
)

# Run the model using the R2jags package
model02b.jags.samples = jags(
  data = model02b.jags.data,
  parameters.to.save = model02b.jags.params,
  model.file = model02b.jags.syntax,
  n.chains = 4,
  n.iter = 6000, # R2jags includes burn-in in n.iter
  n.burnin = 1000,
  n.thin = 1
)

model02b.jags.samples 

# extracting MCMC chain from R2jags object
model02b.chainList = mcmc.list(
  lapply(X = 1:dim(model02b.jags.samples$BUGSoutput$sims.array)[2], FUN = function(x) {
    return(as.mcmc(model02b.jags.samples$BUGSoutput$sims.array[,x,]))
  })
)

# gathering autocorrelation information (from coda package, which loads with R2jags)
coda::autocorr.diag(model02b.chainList)

# plotting autocorrelation information
coda::autocorr.plot(model02b.chainList)


# problem--regression isn't saturated...need residual covariance...no easy solution to this with missing data

# adding interactions ------------------------------------

model03a.jags.syntax = function(){
  
  # model for female -------------
  for (obs in 1:N){
    female[obs] ~ dbern(pFemale)  
  }
  
  # priors for female
  pFemale ~ dbeta(1, 1)
  
  # model for cc10
  for (obs in 1:N){
    cc10[obs] ~ dnorm(muCC10, tauCC10)  
  }
  
  
  # priors for cc10
  muCC10 ~ dnorm(0, 1/1000)
  tauCC10 ~ dgamma(1, 1)
  
  # model for use: -----------------
  
  for (obs in 1:N){
    # use mean model: 
    muUse[obs] <- use.intercept + use.female*female[obs] + use.cc10*cc10[obs] + 
      use.cc10xfemale*cc10[obs]*female[obs]
    
    # use model data distribution
    use[obs] ~ dnorm(muUse[obs], tauUse)  
  }
  
  # use prior distributions -- empirical prior for regression parameters
  use.intercept ~ dnorm(useMean0, usePrecision0)
  use.female ~ dnorm(useMean0, usePrecision0)
  use.cc10 ~ dnorm(useMean0, usePrecision0)
  use.cc10xfemale ~ dnorm(useMean0, usePrecision0)
  
  useMean0 ~ dnorm(0, 1/1000)
  usePrecision0 ~ dgamma(1, 1)
  useSigma2_0 = 1/usePrecision0
  
  tauUse ~ dgamma(use.tau.shape0, use.tau.rate0)
  
  # use variance
  sigma2Use <- 1/tauUse
  
  # model for perf: -----------------
  
  # perf model data distribution
  for (obs in 1:N){
    # perf mean model
    muPerf[obs] <- perf.intercept + perf.female*female[obs] + perf.cc10*cc10[obs] + 
      perf.cc10xfemale*cc10[obs]*female[obs]
    perf[obs] ~ dnorm(muPerf[obs], tauPerf)    
  }
  
  
  # perf prior distributions
  perf.intercept ~ dnorm(perfMean0, perfPrecision0)
  perf.female ~ dnorm(perfMean0, perfPrecision0)
  perf.cc10 ~ dnorm(perfMean0, perfPrecision0)
  perf.cc10xfemale ~ dnorm(perfMean0, perfPrecision0)
  
  perfMean0 ~ dnorm(0, 1/1000)
  perfPrecision0 ~ dgamma(1, 1)
  perfSigma2_0 = 1/perfPrecision0
  
  tauPerf ~ dgamma(perf.tau.shape0, perf.tau.rate0)
  
  # perf residual variance
  sigma2Perf <- 1/tauPerf
  
  
}

# Prepare data for JAGS -- here perf and use are separate variables
model03a.jags.data = list(
  use = data02$use,
  perf = data02$perf,
  female = data02$female,
  cc10 = data02$cc10,
  N = nrow(data02),
  use.tau.shape0 = 1,
  use.tau.rate0 = 1,
  perf.tau.shape0 = 1,
  perf.tau.rate0 = 1
) 

model03a.jags.params = c(
  "pFemale",
  "muCC10",
  "tauCC10",
  "use.intercept", 
  "use.female",
  "use.cc10",
  "use.cc10xfemale",
  "useMean0",
  "useSigma2_0",
  "sigma2Use", 
  "perf.intercept", 
  "perf.female",
  "perf.cc10",
  "perf.cc10xfemale",
  "sigma2Perf",
  "perfMean0",
  "perfSigma2_0",
  "perf[3]",
  "use[8]"
)

# Run the model using the R2jags package
model03a.jags.samples = jags(
  data = model03a.jags.data,
  parameters.to.save = model03a.jags.params,
  model.file = model03a.jags.syntax,
  n.chains = 4,
  n.iter = 6000, # R2jags includes burn-in in n.iter
  n.burnin = 1000,
  n.thin = 1
)

model03a.jags.samples 

# extracting MCMC chain from R2jags object
model03a.chainList = mcmc.list(
  lapply(X = 1:dim(model03a.jags.samples$BUGSoutput$sims.array)[2], FUN = function(x) {
    return(as.mcmc(model03a.jags.samples$BUGSoutput$sims.array[,x,]))
  })
)

# gathering autocorrelation information (from coda package, which loads with R2jags)
coda::autocorr.diag(model03a.chainList)

# plotting autocorrelation information
coda::autocorr.plot(model03a.chainList)

# auxiliary variables ----------------------------------------------------------------------------------

# we will add the remaining variables to the model as auxiliary variables
# but, as these also have missing data, we need to form their joint model as a series of factored regressions

model06.syntax = "

# Exogenous Variables (now with auxiliary variables)
#Means:
cc10 ~ 1
female ~ 1

# Auxiliary variables
hsl ~ 1
msc ~ 1
mas ~ 1
mse ~ 1

#Variances:
cc10 ~~ cc10
female ~~ female

# Auxiliary variables:
hsl ~~ hsl
msc ~~ msc
mas ~~ mas
mse ~~ mse

#Covariances:
cc10 ~~ female + hsl + msc + mas + mse
female ~~ hsl + msc + mas + mse
hsl ~~ msc + mas + mse
msc ~~ mas + mse
mas ~~ mse


# Endogenous Variables 
#Regressions:
perf ~ 1 + female + cc10
use  ~ 1 + female + cc10

#Residual Variances:
perf ~~ perf
use  ~~ use

#Residual Covariance:
perf ~~ use

#Endogenous/exogenous covariances
perf ~~ hsl + msc + mas + mse
use  ~~ hsl + msc + mas + mse

"


model04a.jags.syntax = function(){
  
  # model for hsl ---------------

  for (obs in 1:N){
    hsl[obs] ~ dnorm(muHSL, tauHSL)  
  }
  
  # priors for hsl
  muHSL ~ dnorm(0, 1/1000)
  tauHSL ~ dgamma(1, 1)
  
  # model for msc ---------------
  
  for (obs in 1:N){
    mscMean[obs] <- msc.intercept + msc.hsl*hsl[obs]
    msc[obs] ~ dnorm(mscMean[obs], tauMSC)
  }
  
  # priors for msc
  msc.intercept ~ dnorm(muMSC0, tauMSC0)
  msc.hsl ~ dnorm(muMSC0, tauMSC0)
  muMSC0 ~ dnorm(0, 1/1000)
  tauMSC0 ~ dgamma(1, 1)
  tauMSC ~ dgamma(1, 1)
  
  # model for mas ---------------
  
  for (obs in 1:N){
    masMean[obs] <- mas.intercept + mas.hsl*hsl[obs] + mas.msc*msc[obs]
    mas[obs] ~ dnorm(masMean[obs], tauMAS)
  }
  
  # priors for mas
  mas.intercept ~ dnorm(muMAS0, taoMAS0)
  mas.hsl ~ dnorm(muMAS0, taoMAS0)
  mas.msc ~ dnorm(muMAS0, taoMAS0)
  muMAS0 ~ dnorm(0, 1/1000)
  taoMAS0 ~ dgamma(1, 1)
  tauMAS ~ dgamma(1, 1)
  
  # model for mse ---------------
  
  for (obs in 1:N){
    mseMean[obs] <- mse.intercept + mse.hsl*hsl[obs] + mse.msc*msc[obs] + mse.mas*mas[obs]
    mse[obs] ~ dnorm(mseMean[obs], tauMSE)
  }
  
  # priors for mse
  mse.intercept ~ dnorm(muMSE0, tauMSE0)
  mse.hsl ~ dnorm(muMSE0, tauMSE0)
  mse.msc ~ dnorm(muMSE0, tauMSE0)
  mse.mas ~ dnorm(muMSE0, tauMSE0)
  muMSE0 ~ dnorm(0, 1/1000)
  tauMSE0 ~ dgamma(1, 1)
  tauMSE ~ dgamma(1, 1)
  
    
  # model for female -------------
  for (obs in 1:N){
    pFemale[obs] = ilogit(female.intercept + female.hsl*hsl[obs] + 
                            female.msc*msc[obs] + female.mas*mas[obs] + female.mse*mse[obs])
    female[obs] ~ dbern(pFemale[obs])
  }
  
  # priors for female
  female.intercept ~ dnorm(muFemale0, tauFemale0)
  female.hsl ~ dnorm(muFemale0, tauFemale0)
  female.msc ~ dnorm(muFemale0, tauFemale0)
  female.mas ~ dnorm(muFemale0, tauFemale0)
  female.mse ~ dnorm(muFemale0, tauFemale0)
  muFemale0 ~ dnorm(0, 1/1000)
  tauFemale0 ~ dgamma(1, 1)
  
  
  # model for cc10
  for (obs in 1:N){
    muCC10[obs] = cc10.intercept + cc10.hsl*hsl[obs] + cc10.msc*msc[obs] + 
      cc10.mas*mas[obs] + cc10.mse*mse[obs]
    cc10[obs] ~ dnorm(muCC10[obs], tauCC10)  
  }
  
  
  # priors for cc10
  cc10.intercept ~ dnorm(muCC10_0, tauCC10_0)
  cc10.hsl ~ dnorm(muCC10_0, tauCC10_0)
  cc10.msc ~ dnorm(muCC10_0, tauCC10_0)
  cc10.mas ~ dnorm(muCC10_0, tauCC10_0)
  cc10.mse ~ dnorm(muCC10_0, tauCC10_0)
  muCC10_0 ~ dnorm(0, 1/1000)
  tauCC10_0 ~ dgamma(1, 1)
  
  tauCC10 ~ dgamma(1, 1)
  
  # model for use: -----------------
  
  for (obs in 1:N){
    # use mean model: 
    muUse[obs] <- use.intercept + use.female*female[obs] + use.cc10*cc10[obs] + 
      use.cc10xfemale*cc10[obs]*female[obs] + use.hsl*hsl[obs] + use.msc*msc[obs] + 
      use.mas*mas[obs] + use.mse*mse[obs]
    
    # use model data distribution
    use[obs] ~ dnorm(muUse[obs], tauUse)  
  }
  
  # use prior distributions -- empirical prior for regression parameters
  use.intercept ~ dnorm(useMean0, usePrecision0)
  use.female ~ dnorm(useMean0, usePrecision0)
  use.cc10 ~ dnorm(useMean0, usePrecision0)
  use.cc10xfemale ~ dnorm(useMean0, usePrecision0)
  use.hsl ~ dnorm(useMean0, usePrecision0)
  use.msc ~ dnorm(useMean0, usePrecision0)
  use.mas ~ dnorm(useMean0, usePrecision0)
  use.mse ~ dnorm(useMean0, usePrecision0)
  
  useMean0 ~ dnorm(0, 1/1000)
  usePrecision0 ~ dgamma(1, 1)
  useSigma2_0 = 1/usePrecision0
  
  tauUse ~ dgamma(use.tau.shape0, use.tau.rate0)
  
  # use variance
  sigma2Use <- 1/tauUse
  
  # model for perf: -----------------
  
  # perf model data distribution
  for (obs in 1:N){
    # perf mean model
    muPerf[obs] <- perf.intercept + perf.female*female[obs] + perf.cc10*cc10[obs] + 
      perf.cc10xfemale*cc10[obs]*female[obs] + perf.hsl*hsl[obs] + perf.msc*msc[obs] +
      perf.mas*mas[obs] + perf.mse*mse[obs]
    perf[obs] ~ dnorm(muPerf[obs], tauPerf)    
  }
  
  
  # perf prior distributions
  perf.intercept ~ dnorm(perfMean0, perfPrecision0)
  perf.female ~ dnorm(perfMean0, perfPrecision0)
  perf.cc10 ~ dnorm(perfMean0, perfPrecision0)
  perf.cc10xfemale ~ dnorm(perfMean0, perfPrecision0)
  perf.hsl ~ dnorm(perfMean0, perfPrecision0)
  perf.msc ~ dnorm(perfMean0, perfPrecision0)
  perf.mas ~ dnorm(perfMean0, perfPrecision0)
  perf.mse ~ dnorm(perfMean0, perfPrecision0)
  
  perfMean0 ~ dnorm(0, 1/1000)
  perfPrecision0 ~ dgamma(1, 1)
  perfSigma2_0 = 1/perfPrecision0
  
  tauPerf ~ dgamma(perf.tau.shape0, perf.tau.rate0)
  
  # perf residual variance
  sigma2Perf <- 1/tauPerf
  
  
}

# Prepare data for JAGS -- here perf and use are separate variables
model04a.jags.data = list(
  use = data02$use,
  perf = data02$perf,
  female = data02$female,
  cc10 = data02$cc10,
  hsl = data02$hsl,
  msc = data02$msc,
  mas = data02$mas,
  mse = data02$mse,
  N = nrow(data02),
  use.tau.shape0 = 1,
  use.tau.rate0 = 1,
  perf.tau.shape0 = 1,
  perf.tau.rate0 = 1
) 

model04a.jags.params = c(
  "muHSL",
  "tauHSL",
  "msc.intercept",
  "msc.hsl",
  "tauMSC",
  "mas.intercept",
  "mas.hsl",
  "mas.msc",
  "tauMAS",
  "mse.intercept",
  "mse.hsl",
  "mse.msc",
  "mse.mas",
  "tauMSE",
  "female.intercept",
  "female.hsl",
  "female.msc",
  "female.mas",
  "female.mse",
  "cc10.intercept",
  "cc10.hsl",
  "cc10.msc",
  "cc10.mas",
  "cc10.mse",
  "muCC10_0",
  "tauCC10_0",
  "use.intercept", 
  "use.female",
  "use.cc10",
  "use.cc10xfemale",
  "use.hsl",
  "use.msc",
  "use.mas",
  "use.mse",
  "useMean0",
  "useSigma2_0",
  "sigma2Use", 
  "perf.intercept", 
  "perf.female",
  "perf.cc10",
  "perf.cc10xfemale",
  "perf.hsl",
  "perf.msc",
  "perf.mas",
  "perf.mse",
  "sigma2Perf",
  "perfMean0",
  "perfSigma2_0",
  "perf[3]",
  "use[8]"
)

# Run the model using the R2jags package
model04a.jags.samples = jags(
  data = model04a.jags.data,
  parameters.to.save = model04a.jags.params,
  model.file = model04a.jags.syntax,
  n.chains = 4,
  n.iter = 6000, # R2jags includes burn-in in n.iter
  n.burnin = 1000,
  n.thin = 1
)

model04a.jags.samples 

# extracting MCMC chain from R2jags object
model04a.chainList = mcmc.list(
  lapply(X = 1:dim(model04a.jags.samples$BUGSoutput$sims.array)[2], FUN = function(x) {
    return(as.mcmc(model04a.jags.samples$BUGSoutput$sims.array[,x,]))
  })
)

# gathering autocorrelation information (from coda package, which loads with R2jags)
coda::autocorr.diag(model04a.chainList)

# plotting autocorrelation information
coda::autocorr.plot(model04a.chainList)

# from Chapter 6: using data augmentation for Gibbs sampling for models with outcomes that are binary variables
# we can use the same approach for missing data in JAGS
