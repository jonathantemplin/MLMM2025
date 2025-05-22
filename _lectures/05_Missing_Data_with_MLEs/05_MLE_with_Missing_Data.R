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

haspackage = require("semPlot")
if (haspackage==FALSE){
  install.packages("semPlot")
}
library(semPlot)

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
model01.syntax = "

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
model01.fit = lavaan(model01.syntax, data=data02, estimator = "MLR", mimic = "MPLUS")
summary(model01.fit, standardized=TRUE, fit.measures=TRUE)

# mapping analysis onto multivariate normal distribution

#to see the model-implied mean vector and covariance matrix, use the fitted() function
fitted(model01.fit)

# mean vector in data:
meanVec = colMeans(data02[c("perf", "use")], na.rm = TRUE)
meanVec

# compare to covariance matrix in data:
cov(data02[c("perf", "use")], use = "complete.obs")

# difference is because ML divides by N and not N-1
N = sum(!is.na(newData$perf) | !is.na(newData$use))
covMat = cov(data02[c("perf", "use")], use = "complete.obs")*(N-1)/N
covMat

#to see the saturated model mean vector and covariance matrix
inspect(model01.fit, what="sampstat.h1")

#to see the discrepancy between the model-implied and saturated mean vector and covariance matrix, use the residuals() function
residuals(model01.fit, type = "raw")

#to see the normalized residuals:
residuals(model01.fit, type = "normalized")

#to see R-squared values for DVs == nothing is a DV in this model
inspect(model01.fit, what="r2") #r-squared values for DVs

#plot path diagram: initial version without numbers
semPaths(model01.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, what ="path")

semPaths(model01.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, whatLabels ="est")

# next, let's show that the loglikelihood in the model is the same as what we would obtain in R
completeData = data02[, c("perf", "use")]

# remove cases that are missing on all variables
completeData = completeData[!is.na(completeData$perf) | !is.na(completeData$use),]
logLike = rep(NA, nrow(completeData))
for (obs in 1:nrow(completeData)){
  obsData = completeData[obs,]
  obsVars = names(completeData)[which(!is.na(obsData))]
  logLike[obs] = dmvnorm(
    x = obsData[obsVars], 
    mean = meanVec[obsVars], 
    sigma = covMat[obsVars, obsVars, drop=FALSE], 
    log = TRUE
  )
}
# our calculation
sum(logLike)

# lavaan's calculation (from each case)
sum(inspect(model01.fit, what="loglik.casewise"), na.rm = TRUE)

# compare with factored regression specification: perf ~ use and use ~ 1 --> same log likelihood
model01b.syntax = "

# exogenous mean:
use ~ 1

# exogenous variance:
use ~~ use

# regression:
perf ~ 1 + use

# residual variance:
perf ~~ perf
"

model01b.fit = lavaan(model01b.syntax, data=data02, estimator = "MLR", mimic = "MPLUS")
summary(model01b.fit, standardized=TRUE, fit.measures=TRUE)

# compare other factored regression specification: use ~ perf and perf ~ 1 --> same log likelihood
model01c.syntax = "

# exogenous mean:
perf ~ 1

# exogenous variance:
perf ~~ perf

# regression:
use ~ 1 + perf

# residual variance:
use ~~ use

"

model01c.fit = lavaan(model01c.syntax, data=data02, estimator = "MLR", mimic = "MPLUS")
summary(model01c.fit, standardized=TRUE, fit.measures=TRUE)


# second model--add predictors as main effects only: ===================================================================

#analysis syntax
model02.syntax = "

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
model02.fit = lavaan(model02.syntax, data=data02, estimator = "MLR", mimic = "MPLUS")

#analysis summary (note the additional terms: standardized = TRUE for standardized estimates and fit.measures=TRUE for model fit indices)
summary(model02.fit, standardized=TRUE, fit.measures=TRUE)

#to see the model-implied mean vector and covariance matrix, use the fitted() function
fitted(model02.fit)

#to see the saturated model mean vector and covariance matrix
inspect(model02.fit, what="sampstat.h1")

#to see the discrepancy between the model-implied and saturated mean vector and covariance matrix, use the residuals() function
residuals(model02.fit, type = "raw")

#to see the normalized residuals:
residuals(model02.fit, type = "normalized")

#to see modification indices:
modindices(model02.fit)

#to see R-squared values for DVs
inspect(model02.fit, what="r2") #r-squared values for DVs

#plot path diagram: initial version without numbers
semPaths(model02.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, what ="path")

semPaths(model02.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, whatLabels ="est")

# next, bringing more variables into the likelihood function 
model03.syntax = "
# Exogenous Variables 
#Mean:
cc10 ~ 1

#Variance:
cc10 ~~ cc10

# Endogenous Variables 
#Regressions:
perf ~ 1 + female + cc10
use  ~ 1 + female + cc10

#Residual Variances:
perf ~~ perf
use  ~~ use

#Residual Covariance:
perf ~~ use

"

#analysis estimation -- note: no warnings about missing data
model03.fit = lavaan(model03.syntax, data=data02, estimator = "MLR", mimic = "MPLUS")

#analysis summary (note the additional terms: standardized = TRUE for standardized estimates and fit.measures=TRUE for model fit indices)
summary(model03.fit, standardized=TRUE, fit.measures=TRUE)

#to see the model-implied mean vector and covariance matrix, use the fitted() function
fitted(model03.fit)

#to see the saturated model mean vector and covariance matrix
inspect(model03.fit, what="sampstat.h1")

#to see the discrepancy between the model-implied and saturated mean vector and covariance matrix, use the residuals() function
residuals(model03.fit, type = "raw")

#to see the normalized residuals:
residuals(model03.fit, type = "normalized")

#to see modification indices:
modindices(model03.fit)

#to see R-squared values for DVs
inspect(model03.fit, what="r2") #r-squared values for DVs

#plot path diagram: initial version without numbers
semPaths(model03.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, what ="path")

semPaths(model03.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, whatLabels ="est")


# problem: what if we have missing on the categorical variable female?

# example 2: factored regression in lavaan but also EM algorithm

# we will create missing data on female (missing completely at random for now)
data03 = data02
data03$female[sample(x = 1:nrow(data03), replace = TRUE, size = .1*nrow(data03) )] = NA
data03$female

# if we run the previous model (not having female specified in the likelihood): missing data warning
model03b.fit = lavaan(model03.syntax, data=data03, estimator = "MLR", mimic = "MPLUS")


# add female to the likelihood function:
model04.syntax = "

# Exogenous Variables 
#Mean:
cc10 ~ 1
female ~ 1

#Variance:
cc10 ~~ cc10
female ~~ female

#Covariance:
cc10 ~~ female

# Endogenous Variables 
#Regressions:
perf ~ 1 + female + cc10
use  ~ 1 + female + cc10

#Residual Variances:
perf ~~ perf
use  ~~ use

#Residual Covariance:
perf ~~ use


"

model04.fit = lavaan(model04.syntax, data=data03, estimator = "MLR", mimic = "MPLUS")

#analysis summary (note the additional terms: standardized = TRUE for standardized estimates and fit.measures=TRUE for model fit indices)
summary(model04.fit, standardized=TRUE, fit.measures=TRUE)

#to see the model-implied mean vector and covariance matrix, use the fitted() function
fitted(model04.fit)

#to see the saturated model mean vector and covariance matrix
inspect(model04.fit, what="sampstat.h1")

#to see the discrepancy between the model-implied and saturated mean vector and covariance matrix, use the residuals() function
residuals(model04.fit, type = "raw")

#to see the normalized residuals:
residuals(model04.fit, type = "normalized")

#to see modification indices:
modindices(model04.fit)

#to see R-squared values for DVs
inspect(model04.fit, what="r2") #r-squared values for DVs

#plot path diagram: initial version without numbers
semPaths(model04.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, what ="path")

semPaths(model04.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, whatLabels ="est")

# problem: female is not normally distributed....must fix using factored regression and EM algorithm
# regression will be this: f(use, perf, cc10, female) = 
# f(use|perf, cc10, female) * f(perf|cc10, female) * f(cc10|female) * f(female)

# first create random parameter values
# f(female):
female_intercept = 0

# f(cc10|female) 
cc10_intercept = 0
cc10_female = 1

# f(perf|cc10, female):
perf_intercept = 0
perf_cc10 = 1
perf_female = 1

# f(use|perf, cc10, female):
use_intercept = 0
use_perf = 1
use_cc10 = 1
use_female = 1

# next we impute data in the E-step (we need to note which values are missing):
femaleMissing = which(is.na(data03$female))
cc10Missing = which(is.na(data03$cc10))
perfMissing = which(is.na(data03$perf))
useMissing = which(is.na(data03$use))

# here is our working copy of the data:
algorithmData = data03

maxIterations = 10

iteration = 1
while(iteration < maxIterations){
  
  # we will use the EM algorithm to impute the missing data
  algorithmData$female[femaleMissing] = 
    round(exp(female_intercept)/(1+exp(female_intercept)), 0)
  
  algorithmData$cc10[cc10Missing] = 
    cc10_intercept + 
    cc10_female*algorithmData$female[cc10Missing]
  
  algorithmData$perf[perfMissing] =
    perf_intercept + 
    perf_cc10*algorithmData$cc10[perfMissing] + 
    perf_female*algorithmData$female[perfMissing]
  
  algorithmData$use[useMissing] =
    use_intercept + 
    use_perf*algorithmData$perf[useMissing] + 
    use_cc10*algorithmData$cc10[useMissing] +
    use_female*algorithmData$female[useMissing]
  
  # now the M-step, we (1) calculate new regression coefficients for each factored model
  # then we (2) calculate the loglikelihood for each observation
  femaleModel = glm(female ~ 1, data = algorithmData, family = binomial(link="logit"))
  female_intercept = femaleModel$coefficients[1]
  logLikeFemale = dbinom(
    x = algorithmData$female, 
    size = 1, 
    prob = exp(female_intercept)/(1+exp(female_intercept)), 
    log = TRUE
  )
  
  cc10Model = lm(cc10 ~ female, data = algorithmData)
  
  cc10_intercept = cc10Model$coefficients["(Intercept)"]
  cc10_female = cc10Model$coefficients["female"]
  cc10_residVar = anova(cc10Model)$`Mean Sq`[2]
  
  logLikeCC10 = dnorm(
    x = algorithmData$cc10, 
    mean = cc10_intercept + cc10_female*algorithmData$female,
    sd = sqrt(cc10_residVar),
    log = TRUE
  )
  
  perfModel = lm(perf ~ cc10 + female, data = algorithmData)
  
  perf_intercept = perfModel$coefficients["(Intercept)"]
  perf_cc10 = perfModel$coefficients["cc10"]
  perf_female = perfModel$coefficients["female"]
  perf_residVar = anova(perfModel)$`Mean Sq`[2]
  
  logLikePerf = dnorm(
    x = algorithmData$perf, 
    mean = perf_intercept + perf_cc10*algorithmData$cc10 + perf_female*algorithmData$female,
    sd = sqrt(perf_residVar),
    log = TRUE
  )
  
  useModel = lm(use ~ perf + cc10 + female, data = algorithmData)
  use_intercept = useModel$coefficients["(Intercept)"]
  use_perf = useModel$coefficients["perf"]
  use_cc10 = useModel$coefficients["cc10"]
  use_female = useModel$coefficients["female"]
  
  use_residVar = anova(useModel)$`Mean Sq`[2]
  
  logLikeUse = dnorm(
    x = algorithmData$use, 
    mean = use_intercept + use_perf*algorithmData$perf + use_cc10*algorithmData$cc10 + use_female*algorithmData$female,
    sd = sqrt(use_residVar),
    log = TRUE
  )
  
  # next compile all loglikelihoods
  logLikeMatrix = cbind(
    logLikeFemale,
    logLikeCC10,
    logLikePerf,
    logLikeUse
  )
  
  loglik = sum(rowSums(logLikeMatrix))
  
  # write iteration number and loglikelihood
  cat("Iteration: ", iteration, "Loglikelihood: ", loglik, "\n")
  iteration = iteration + 1
}
summary(femaleModel)
summary(cc10Model)
summary(perfModel)
summary(useModel)

# adding interactions using the just-another-variable approach (SEM framework only) ------------------------------------

# create interaction variable between female and cc10
data03$femXcc10 = data03$female*data03$cc10

# add variable to model (and explicitly have lavaan estimate its mean/variance)

model05.syntax = "

# Exogenous Variables 
#Mean:
cc10 ~ 1
female ~ 1
femXcc10 ~ 1

#Variances:
cc10 ~~ cc10
female ~~ female
femXcc10 ~~ femXcc10

#Covariances:
cc10 ~~ female
cc10 ~~ femXcc10
female ~~ femXcc10

# Endogenous Variables 
#Regressions:
perf ~ 1 + female + cc10 + femXcc10
use  ~ 1 + female + cc10 + femXcc10

#Residual Variances:
perf ~~ perf
use  ~~ use

#Residual Covariance:
perf ~~ use


"

model05.fit = lavaan(model05.syntax, data=data03, estimator = "MLR", mimic = "MPLUS")

#analysis summary (note the additional terms: standardized = TRUE for standardized estimates and fit.measures=TRUE for model fit indices)
summary(model05.fit, standardized=TRUE, fit.measures=TRUE)

#to see the model-implied mean vector and covariance matrix, use the fitted() function
fitted(model05.fit)

#to see the saturated model mean vector and covariance matrix
inspect(model05.fit, what="sampstat.h1")

#to see the discrepancy between the model-implied and saturated mean vector and covariance matrix, use the residuals() function
residuals(model05.fit, type = "raw")

#to see the normalized residuals:
residuals(model05.fit, type = "normalized")

#to see modification indices:
modindices(model05.fit)

#to see R-squared values for DVs
inspect(model05.fit, what="r2") #r-squared values for DVs

#plot path diagram: initial version without numbers
semPaths(model05.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, what ="path")

semPaths(model05.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, whatLabels ="est")

# now with factored regression:
# regression will be this: f(use, perf, cc10, female, femXcc10) = 
# f(use|perf, cc10, female, femXcc10) * f(perf|cc10, female, femXcc10) * f(cc10|female) * f(female)

# first create random parameter values
# f(female):
female_intercept = 0

# f(cc10|female) 
cc10_intercept = 0
cc10_female = 1

# f(perf|cc10, female):
perf_intercept = 0
perf_cc10 = 1
perf_female = 1
perf_femXcc10 = 1

# f(use|perf, cc10, female):
use_intercept = 0
use_perf = 1
use_cc10 = 1
use_female = 1
use_femXcc10 = 1

# next we impute data in the E-step (we need to note which values are missing):
femaleMissing = which(is.na(data03$female))
cc10Missing = which(is.na(data03$cc10))
perfMissing = which(is.na(data03$perf))
useMissing = which(is.na(data03$use))

# here is our working copy of the data:
algorithmData = data03

maxIterations = 10

iteration = 1
while(iteration < maxIterations){
  
  # we will use the EM algorithm to impute the missing data
  algorithmData$female[femaleMissing] = 
    round(exp(female_intercept)/(1+exp(female_intercept)), 0)
  
  algorithmData$cc10[cc10Missing] = 
    cc10_intercept + 
    cc10_female*algorithmData$female[cc10Missing]
  
  algorithmData$perf[perfMissing] =
    perf_intercept + 
    perf_cc10*algorithmData$cc10[perfMissing] + 
    perf_female*algorithmData$female[perfMissing] +
    perf_femXcc10*algorithmData$cc10[perfMissing]*algorithmData$female[perfMissing]
  
  algorithmData$use[useMissing] =
    use_intercept + 
    use_perf*algorithmData$perf[useMissing] + 
    use_cc10*algorithmData$cc10[useMissing] +
    use_female*algorithmData$female[useMissing] +
    use_femXcc10*algorithmData$cc10[useMissing]*algorithmData$female[useMissing]
  
  # now the M-step, we (1) calculate new regression coefficients for each factored model
  # then we (2) calculate the loglikelihood for each observation
  femaleModel = glm(female ~ 1, data = algorithmData, family = binomial(link="logit"))
  female_intercept = femaleModel$coefficients[1]
  logLikeFemale = dbinom(
    x = algorithmData$female, 
    size = 1, 
    prob = exp(female_intercept)/(1+exp(female_intercept)), 
    log = TRUE
  )
  
  cc10Model = lm(cc10 ~ female, data = algorithmData)
  
  cc10_intercept = cc10Model$coefficients["(Intercept)"]
  cc10_female = cc10Model$coefficients["female"]
  cc10_residVar = anova(cc10Model)$`Mean Sq`[2]
  
  logLikeCC10 = dnorm(
    x = algorithmData$cc10, 
    mean = cc10_intercept + cc10_female*algorithmData$female,
    sd = sqrt(cc10_residVar),
    log = TRUE
  )
  
  perfModel = lm(perf ~ cc10 + female + cc10:female, data = algorithmData)
  
  perf_intercept = perfModel$coefficients["(Intercept)"]
  perf_cc10 = perfModel$coefficients["cc10"]
  perf_female = perfModel$coefficients["female"]
  perf_femXcc10 = perfModel$coefficients["cc10:female"]
  perf_residVar = anova(perfModel)$`Mean Sq`[2]
  
  logLikePerf = dnorm(
    x = algorithmData$perf, 
    mean = perf_intercept + 
      perf_cc10*algorithmData$cc10 + 
      perf_female*algorithmData$female +
      perf_femXcc10*algorithmData$cc10*algorithmData$female,
    sd = sqrt(perf_residVar),
    log = TRUE
  )
  
  useModel = lm(use ~ perf + cc10 + female + cc10:female, data = algorithmData)
  use_intercept = useModel$coefficients["(Intercept)"]
  use_perf = useModel$coefficients["perf"]
  use_cc10 = useModel$coefficients["cc10"]
  use_female = useModel$coefficients["female"]
  use_femXcc10 = useModel$coefficients["cc10:female"]
  
  use_residVar = anova(useModel)$`Mean Sq`[2]
  
  logLikeUse = dnorm(
    x = algorithmData$use, 
    mean = use_intercept + 
      use_perf*algorithmData$perf + 
      use_cc10*algorithmData$cc10 + 
      use_female*algorithmData$female +
      use_femXcc10*algorithmData$cc10*algorithmData$female,
    sd = sqrt(use_residVar),
    log = TRUE
  )
  
  # next compile all loglikelihoods
  logLikeMatrix = cbind(
    logLikeFemale,
    logLikeCC10,
    logLikePerf,
    logLikeUse
  )
  
  loglik = sum(rowSums(logLikeMatrix))
  
  # write iteration number and loglikelihood
  cat("Iteration: ", iteration, "Loglikelihood: ", loglik, "\n")
  iteration = iteration + 1
}
summary(femaleModel)
summary(cc10Model)
summary(perfModel)
summary(useModel)


# auxiliary variables (SEM framework) ----------------------------------------------------------------------------------

# we will add the remaining variables to the model as auxiliary variables

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

model06.fit = lavaan(model06.syntax, data=data03, estimator = "MLR", mimic = "MPLUS")

#analysis summary (note the additional terms: standardized = TRUE for standardized estimates and fit.measures=TRUE for model fit indices)
summary(model06.fit, standardized=TRUE, fit.measures=TRUE)

#to see the model-implied mean vector and covariance matrix, use the fitted() function
fitted(model06.fit)

#to see the saturated model mean vector and covariance matrix
inspect(model06.fit, what="sampstat.h1")

#to see the discrepancy between the model-implied and saturated mean vector and covariance matrix, use the residuals() function
residuals(model06.fit, type = "raw")

#to see the normalized residuals:
residuals(model06.fit, type = "normalized")

#to see modification indices:
modindices(model06.fit)

#to see R-squared values for DVs
inspect(model06.fit, what="r2") #r-squared values for DVs

#plot path diagram: initial version without numbers
semPaths(model06.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, what ="path")

semPaths(model06.fit, intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=2,
         optimizeLatRes = TRUE, sizeMan = 8, whatLabels ="est")
