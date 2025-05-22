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

#Playing with Matrices =================================== ==============================================================
#IMPORT DATA AND PUT INTO DATASET
data01 = read.csv("jobperf.csv")

# descriptive statistics

# means:
apply(data01[c("IQ", "perfC")], MARGIN=2, FUN=mean)

# mean vector:
t(t(apply(data01[c("IQ", "perfC")], MARGIN=2, FUN=mean)))

# covariance matrix:
cov(data01[c("IQ", "perfC")])

# correlation matrix:
cor(data01[c("IQ", "perfC")])

# creating a scatterplot of our data but have the points be filled in with black
plot(data01$IQ, data01$perfC, xlab = "IQ", ylab = "Performance (complete)", cex = 2, pch = 19)

# add means to plot 
points(x = 100, y = 10.35, col = "blue", pch = 19, cex = 2, type = "p")

# add regression line to plot
abline(lm(data01$perfC ~ data01$IQ), col = "red")

# creating correlation matrix from covariance
S = cov(data01[c("IQ", "perfC")])
S

# get diagonal matrix of standard deviations
D = diag(sqrt(diag(S)))
D

# create correlation matrix
R = solve(D) %*% S %*% solve(D)
R

# similarly, convert correlation matrix to convariance matrix
S_new = D %*% R %*% D
S_new

# generalized variance (determinant of S)
det(S)

# generalized variance if no covariance
S2 = S
S2[1,2] = S2[2,1] = 0
det(S2)

# total sample variance (trace of S)
sum(diag(S))


#use dmvnorm() function to return likelihood and log-likelihood from MVN:
XBAR = t(t(apply(data01[c("IQ", "perfC")], MARGIN=2, FUN=mean)))

data01[2,c("IQ", "perfC")]
likelihood_case2 = dmvnorm(x=data01[2,c("IQ", "perfC")],mean = XBAR,sigma = S, log=FALSE)
likelihood_case2

data01[2,c("IQ", "perfC")]
loglikelihood_case2 = dmvnorm(x=data01[2,c("IQ", "perfC")],mean = XBAR,sigma = S, log=TRUE)
loglikelihood_case2

data01[5, c("IQ", "perfC")]
loglikelihood_case5 = dmvnorm(x=data01[5, c("IQ", "perfC")],mean = XBAR,sigma = S, log=TRUE)
loglikelihood_case5

XBAR
loglikelihood_caseXBAR = dmvnorm(x=t(XBAR), mean = XBAR, sigma = S, log=TRUE)
loglikelihood_caseXBAR

# using lavaan for ML estimation of a model where IQ predicts perf =====================================================

model01.syntax = "

  # regressions are indicated by a ~
  # here, perfMAR is predicted by an intercept (the 1; not usually needed in syntax) and IQ
  
  perfMAR ~ 1 + IQ  
  
  # variances are indicated by ~~
  # here, we are estimating the residual variance of perfMar (also usually not needed but included for demonstration)
  
  perfMAR ~~ perfMAR 
  
"


#analysis estimation
model01.fit = sem(model01.syntax, data=data01, mimic = "MPLUS", estimator = "MLR")

#analysis summary (note the additional terms: standardized = TRUE for standardized estimates and fit.measures=TRUE for model fit indices)
summary(model01.fit, standardized=TRUE, fit.measures=TRUE, rsquare = TRUE)

# compare results to OLS regression:
summary(lm(data01$perfMAR ~ data01$IQ))

# development of EM algorithm ==========================================================================================

# Load necessary library
library(MASS)

# Set random seed for reproducibility
set.seed(123)

# Generate synthetic bivariate normal data
n = 100
mu_X = 5
mu_Y = 10
sigma_X = 2
sigma_Y = 3
rho_XY = 0.7
Sigma = matrix(c(sigma_X^2, rho_XY * sigma_X * sigma_Y, 
                 rho_XY * sigma_X * sigma_Y, sigma_Y^2), 2, 2)

data = MASS::mvrnorm(n, mu = c(mu_X, mu_Y), Sigma = Sigma)
colnames(data) = c("X", "Y")

# Introduce missing values in Y
missing_idx = sample(1:n, size = 30, replace = FALSE)  # 30% missing in Y
data[missing_idx, "Y"] = NA

data = data01[c("IQ", "perfMAR")]
missing_idx = which(is.na(data$perfMAR))
colnames(data) = c("X", "Y")


# Function to compute log-likelihood
log_likelihood = function(data, mu_X, mu_Y, sigma_X, sigma_Y, rho) {
  observed_idx = !is.na(data[, "Y"])
  missing_idx = is.na(data[, "Y"])
  
  X_obs = data[observed_idx, "X"]
  Y_obs = data[observed_idx, "Y"]
  X_mis = data[missing_idx, "X"]  # X values corresponding to missing Y
  
  # Compute joint covariance matrix
  Sigma = matrix(c(sigma_X^2, rho * sigma_X * sigma_Y, 
                   rho * sigma_X * sigma_Y, sigma_Y^2), 2, 2)
  mu_vec = c(mu_X, mu_Y)
  
  # Log-likelihood for observed cases
  ll_obs = sum(dmvnorm(cbind(X_obs, Y_obs), mean = mu_vec, sigma = Sigma, log = TRUE))
  
  # Log-likelihood for missing Y cases: Integrate over P(Y | X)
  beta = Sigma[1,2] / Sigma[1,1]  # Regression coefficient of Y on X
  sigma_Y_given_X = sqrt(Sigma[2,2] - beta^2 * Sigma[1,1])
  
  ll_missing = sum(dnorm(X_mis, mean = mu_X, sd = sigma_X, log = TRUE)) + 
    sum(dnorm(mu_Y + beta * (X_mis - mu_X), mean = mu_Y, sd = sigma_Y_given_X, log = TRUE))
  
  return(ll_obs + ll_missing)
}

# EM Algorithm for parameter estimation
em_algorithm = function(data, tol = 1e-6, max_iter = 1000) {
  observed_idx = !is.na(data[, "Y"])
  missing_idx = is.na(data[, "Y"])
  
  X_all = data[, "X"]
  Y_obs = data[observed_idx, "Y"]
  X_obs = data[observed_idx, "X"]
  
  # Initialize parameters
  mu_X_hat = mean(X_all, na.rm = TRUE)
  mu_Y_hat = mean(Y_obs, na.rm = TRUE)
  sigma_X_hat = sd(X_all, na.rm = TRUE)
  sigma_Y_hat = sd(Y_obs, na.rm = TRUE)
  rho_hat = cor(X_obs, Y_obs, use = "complete.obs")
  
  Sigma_hat = matrix(c(sigma_X_hat^2, rho_hat * sigma_X_hat * sigma_Y_hat,
                       rho_hat * sigma_X_hat * sigma_Y_hat, sigma_Y_hat^2), 2, 2)
  
  iter = 0
  diff = Inf
  
  # Initialize history storage
  iter_history = matrix(ncol = 5, nrow = 0)
  colnames(iter_history) = c("Iteration", "mu_Y", "sigma_Y", "rho", "LogLikelihood")
  
  while (diff > tol && iter < max_iter) {
    iter = iter + 1
    mu_Y_old = mu_Y_hat
    Sigma_old = Sigma_hat
    
    # E-step: Compute expected values for missing Y
    beta = Sigma_hat[1,2] / Sigma_hat[1,1]  # Regression coefficient of Y on X
    sigma_Y_given_X = sqrt(Sigma_hat[2,2] - beta^2 * Sigma_hat[1,1])
    
    Y_missing_est = mu_Y_hat + beta * (data[missing_idx, "X"] - mu_X_hat)
    
    # M-step: Update estimates
    Y_filled = data[, "Y"]
    Y_filled[missing_idx] = Y_missing_est
    
    mu_Y_hat = mean(Y_filled, na.rm = TRUE)
    sigma_Y_hat = sqrt(var(Y_filled, na.rm = TRUE))
    rho_hat = cor(X_all, Y_filled, use = "complete.obs")
    
    Sigma_hat = matrix(c(sigma_X_hat^2, rho_hat * sigma_X_hat * sigma_Y_hat,
                         rho_hat * sigma_X_hat * sigma_Y_hat, sigma_Y_hat^2), 2, 2)
    
    # Compute log-likelihood
    log_lik = log_likelihood(data, mu_X_hat, mu_Y_hat, sigma_X_hat, sigma_Y_hat, rho_hat)
    
    # Store iteration results
    iter_history = rbind(iter_history, c(iter, mu_Y_hat, sigma_Y_hat, rho_hat, log_lik))
    
    diff = sum(abs(mu_Y_hat - mu_Y_old)) + sum(abs(Sigma_hat - Sigma_old))
  }
  
  # Convert iteration history to a data frame and print
  iter_history = as.data.frame(iter_history)
  print(iter_history)
  
  list(mu_Y = mu_Y_hat, sigma_Y = sigma_Y_hat, rho = rho_hat, Sigma = Sigma_hat, iterations = iter)
}

# Run the EM algorithm
em_results = em_algorithm(data)


# simulated data for today's second example =================================================================================================
set.seed(0)

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
allx = rbind(xmen, xwomen)

# Convert to dataframe and rename columns
colnames(allx) = c("female", "hsl", "cc", "use", "msc", "mas", "mse", "perf")

# Introduce missingness in CC, PERF, and USE based on remaining variables
missing_prob = function(hsl, msc, mas, mse) {
  return(0.3 * (hsl > median(allx[, "hsl"])) +
           0.3 * (msc < median(allx[, "msc"])) +
           0.2 * (mas > median(allx[, "mas"])) +
           0.2 * (mse < median(allx[, "mse"])))
}

set.seed(123)
missing_indices = runif(nrow(allx)) < missing_prob(allx[, "hsl"], allx[, "msc"], allx[, "mas"], allx[, "mse"])
allx[missing_indices, c("cc", "perf", "use")] = NA
data02 = as.data.frame(allx)

#centering CC at 10: ===================================================================================================
par(mfrow = c(1,2))
boxplot(data02$cc, main="Boxplot of College Experience (CC)")
hist(data02$cc, main = "Histogram of College Experience (CC)", xlab = "CC")
par(mfrow = c(1,1))

data02$cc10 = data02$cc - 10

#creating interation between female and cc: ============================================================================
data02$femXcc10 = data02$female*data02$cc10


#Building Analysis Model #1: an empty model ============================================================================


#analysis syntax
model01.syntax = "

#Means:
perf ~ 1 + 0*female + 0*cc10 + 0*femXcc10
use  ~ 1 + 0*female + 0*cc10 + 0*femXcc10

#Variances:
perf ~~ perf
use  ~~ use

#Covariance:
perf ~~ use

"

#analysis estimation
model01.fit = sem(model01.syntax, data=data02, mimic = "MPLUS", estimator = "MLR")

#analysis syntax
model02.syntax = "

#Exogenous variables into likelihood function--estimate parameters about them

# means
cc10 ~ 1
femXcc10 ~ 1
female ~ 1

# covariances
cc10 ~~ femXcc10 + female
femXcc10 ~~ female

#Means:
perf ~ 1 + 0*female + 0*cc10 + 0*femXcc10
use  ~ 1 + 0*female + 0*cc10 + 0*femXcc10

#Variances:
perf ~~ perf
use  ~~ use

#Covariance:
perf ~~ use

"

#analysis estimation
model02.fit = sem(model02.syntax, data=data02, mimic = "MPLUS", estimator = "MLR")


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


#model 03: all parameters included -----------------------------------------------------------------------
model03.syntax = "

# means
cc10 ~ 1
femXcc10 ~ 1
female ~ 1

# covariances
cc10 ~~ femXcc10 + female
femXcc10 ~~ female


#Means:
perf ~ 1 + (p_f)*female + (p_cc)*cc10 + (intPerf)*femXcc10
use  ~ 1 + (u_f)*female + (u_cc)*cc10 + (intUse)*femXcc10

#Variances:
perf ~~ perf
use  ~~ use

#Covariance:
perf ~~ use


"
#analysis estimation
model03.fit = sem(model03.syntax, data=data02, mimic = "MPLUS", estimator = "MLR")
#analysis summary (note the additional terms: standardized = TRUE for standardized estimates and fit.measures=TRUE for model fit indices)
summary(model03.fit, standardized=TRUE, fit.measures=TRUE, rsquare=TRUE)


#to compare model fit vs. model01 with scaled likelihood ratio test
anova(model02.fit, model03.fit)

#to see the model-implied mean vector and covariance matrix, use the fitted() function
fitted(model03.fit)
#to see the saturated model mean vector and covariance matrix
inspect(model03.fit, what="sampstat.h1")

#to see the discrepancy between the model-implied and saturated mean vector and covariance matrix, use the residuals() function
residuals(model03.fit)



