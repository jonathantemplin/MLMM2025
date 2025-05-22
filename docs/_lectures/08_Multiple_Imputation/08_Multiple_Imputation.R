rm(list = ls()) # Clear workspace

# --- 1. Load necessary libraries ---
# install.packages("mvtnorm") # Uncomment if not installed
# install.packages("mice")    # Uncomment if not installed
# install.packages("MASS")    # Needed for mvrnorm (for Bayesian sim)
# install.packages("lavaan")  # Uncomment if not installed
# install.packages("semTools") # Uncomment if not installed (provides lavaan.mi)

library(mvtnorm)
library(mice)
library(lavaan)
library(semTools) # For lavaan.mi
library(MASS) # For simulating Bayesian draws (alternative approach shown below)


# --- 2. Replicate Data Simulation (from your provided script) ---
set.seed(20250223) # Use the original seed for data generation

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

Cov = SDdiag %*% Corr %*% SDdiag

xmen = rmvnorm(Nmen, MeanMEN, Cov)
xwomen = rmvnorm(Nwomen, MeanWOMEN, Cov)

xmen = cbind(0, xmen)
xwomen = cbind(1, xwomen)

allx = as.data.frame(rbind(xmen, xwomen))
names(allx) = c("female", "hsl", "cc", "use", "msc", "mas", "mse", "perf")

# Introduce Missingness in CC, PERF, USE (as per your script)
M_CC_intercept = -3; M_CC_female = 1; M_CC_msc = .01; M_CC_mas = .01; M_CC_mse = .01
missingCCmodelLogit = M_CC_intercept + M_CC_female*allx[,"female"] + M_CC_msc*allx[,"msc"] + M_CC_mas*allx[,"mas"] + M_CC_mse*allx[,"mse"]
missingCCmodelProb = exp(missingCCmodelLogit)/(1+exp(missingCCmodelLogit))
makeMissingCC = which(runif(Nmen+Nwomen) < missingCCmodelProb)
allx$cc[makeMissingCC] = NA

M_perf_intercept = -3; M_perf_female = .5; M_perf_msc = .001; M_perf_mas = .02; M_perf_mse = .01
missingPERFmodelLogit = M_perf_intercept + M_perf_female*allx[,"female"] + M_perf_msc*allx[,"msc"] + M_perf_mas*allx[,"mas"] + M_perf_mse*allx[,"mse"]
missingPERFmodelProb = exp(missingPERFmodelLogit)/(1+exp(missingPERFmodelLogit))
makeMissingPerf = which(runif(Nmen+Nwomen) < missingPERFmodelProb)
allx$perf[makeMissingPerf] = NA

M_use_intercept = -3; M_use_female = .5; M_use_msc = .001; M_use_mas = .02; M_use_mse = .01
missingUSEmodelLogit = M_use_intercept + M_use_female*allx[,"female"] + M_use_msc*allx[,"msc"] + M_use_mas*allx[,"mas"] + M_use_mse*allx[,"mse"]
missingUSEmodelProb = exp(missingUSEmodelLogit)/(1+exp(missingUSEmodelLogit))
makeMissingUse = which(runif(Nmen+Nwomen) < missingUSEmodelProb)
allx$use[makeMissingUse] = NA

data02 = as.data.frame(allx)
data02$cc10 = data02$cc - 10

# Select the relevant variables for imputation
data_to_impute_orig <- data02[, c("perf", "use", "cc10", "female")]

# Introduce some missingness into 'female' as well
set.seed(123) # for reproducibility of this step
data_to_impute_orig$female[sample(x = 1:nrow(data_to_impute_orig),
                                  size = round(0.1*nrow(data_to_impute_orig)))] = NA

print("Original data subset with missing values:")
summary(data_to_impute_orig)
print(paste("Number of rows with any NA:", sum(!complete.cases(data_to_impute_orig))))

# --- 3. Manual MICE Loop Demonstration with Parameter Re-estimation ---
cat("\n--- Starting Manual MICE Loop with Parameter Re-estimation ---\n")
cat("NOTE: This shows parameters being re-estimated based on current data,\n")
cat("      but does NOT fully replicate Bayesian posterior draws for parameters.\n")
cat("      Loop may stop if model fitting fails.\n")


# Create a copy to work with
data_imputed_manual <- data_to_impute_orig

# Simple initial imputation (e.g., mean/median/mode) for starting point
mode_female <- as.numeric(names(which.max(table(data_imputed_manual$female))))
data_imputed_manual$female[is.na(data_imputed_manual$female)] <- mode_female
data_imputed_manual$perf[is.na(data_imputed_manual$perf)] <- mean(data_imputed_manual$perf, na.rm = TRUE)
data_imputed_manual$use[is.na(data_imputed_manual$use)] <- mean(data_imputed_manual$use, na.rm = TRUE)
data_imputed_manual$cc10[is.na(data_imputed_manual$cc10)] <- mean(data_imputed_manual$cc10, na.rm = TRUE)

# Define number of iterations for the manual loop
n_iter_manual <- 10
set.seed(456) # Seed for the imputation process itself

# Get indices of missing values for each variable
missing_idx <- lapply(data_to_impute_orig, function(col) which(is.na(col)))
vars_to_impute <- names(missing_idx)[sapply(missing_idx, length) > 0]

# Manual MICE loop
for (iter in 1:n_iter_manual) {
  cat(paste("Manual Iteration:", iter, "\n"))
  
  # Iterate through variables in a specified order
  for (var in vars_to_impute) {
    
    # Define predictors (all other variables in the current imputed dataset)
    predictors <- setdiff(names(data_imputed_manual), var)
    formula_str <- paste(var, "~", paste(predictors, collapse = " + "))
    formula_obj <- as.formula(formula_str)
    
    # Indices for originally missing values in this variable
    idx_miss <- missing_idx[[var]]
    if (length(idx_miss) == 0) next # Skip if no missing values in original
    
    # --- Parameter Re-estimation Step ---
    # Fit model using the *current* state of the data_imputed_manual
    current_data_for_fit <- data_imputed_manual
    current_data_for_fit[[var]][idx_miss] <- NA # Temporarily remove values to be imputed
    current_data_for_fit <- current_data_for_fit[complete.cases(current_data_for_fit[, c(var, predictors)]), ] # Use cases complete for this model
    
    # Data frame with current predictor values for cases needing imputation
    data_to_predict_on <- data_imputed_manual[idx_miss, predictors, drop = FALSE]
    
    # Check if data_to_predict_on has valid rows after potential NA removal above
    if(nrow(data_to_predict_on) == 0 || nrow(current_data_for_fit) == 0) {
      cat(paste("  Skipping imputation for", var, "- insufficient data after handling NAs for model fitting.\n"))
      next
    }
    
    # Fit the appropriate model using the current data
    imputed_values <- NULL # Initialize
    
    if (var == "female") {
      # Fit Logistic Regression
      fit <- glm(formula_obj, data = current_data_for_fit, family = binomial(link = "logit"))
      # Print Coefficients (Demonstration)
      if(iter %% 5 == 0) { # Print every 5 iterations
        cat(paste("  Re-estimated coefs for female (Iter", iter, "):",
                  paste(round(coef(fit), 2), collapse=", "), "\n"))
      }
      # Predict probabilities & draw imputations
      pred_probs <- predict(fit, newdata = data_to_predict_on, type = "response")
      pred_probs <- pmax(pmin(pred_probs, 1 - 1e-8), 1e-8)
      imputed_values <- rbinom(length(idx_miss), 1, pred_probs)
      
    } else {
      # Fit Linear Regression
      fit <- lm(formula_obj, data = current_data_for_fit)
      # Print Coefficients (Demonstration)
      if(iter %% 1 == 0) { # Print every iteration
        cat(paste("  Re-estimated coefs for", var, "(Iter", iter, "):",
                  paste(round(coef(fit), 2), collapse=", "), "\n"))
      }
      # Predict means & draw imputations (add random noise)
      pred_means <- predict(fit, newdata = data_to_predict_on)
      residual_sd <- max(summary(fit)$sigma, 1e-4) # Ensure sd > 0
      imputed_values <- rnorm(length(idx_miss), mean = pred_means, sd = residual_sd)
    }
    
    # Update the imputed dataset ONLY if imputation was successful
    if (!is.null(imputed_values) && length(imputed_values) == length(idx_miss)) {
      data_imputed_manual[[var]][idx_miss] <- imputed_values
    } else {
      # This condition might be less likely without tryCatch, but kept for safety
      cat(paste("  Imputation yielded NULL or wrong length for var:", var, "at iteration", iter, "- keeping previous values.\n"))
    }
    # ***** End of code previously inside tryCatch *****
    
  } # End loop over variables
} # End loop over iterations

cat("--- Manual MICE Loop (with Re-estimation) Finished ---\n")
print("Summary of manually imputed dataset (after last iteration):")
summary(data_imputed_manual)
print(paste("Number of NAs remaining:", sum(is.na(data_imputed_manual)))) # Should be 0


# --- 4. Perform Multiple Imputation using mice package (for comparison) ---
# (Same as before - keeping for reference)
cat("\n--- Starting MICE Package Imputation (for comparison) ---\n")
md.pattern(data_to_impute_orig)
methods <- make.method(data_to_impute_orig)
methods["female"] <- "logreg"
data_to_impute_orig$female = as.factor(data_to_impute_orig$female)

print("Imputation methods to be used by mice package:")
print(methods)
imputed_data_mice <- mice(data_to_impute_orig,
                          m = 5, maxit = 10, method = methods,
                          seed = 500, printFlag = TRUE)
cat("--- MICE Package Imputation Finished ---\n")
summary(imputed_data_mice)

# Compare summary of manually imputed data vs first mice imputed dataset
cat("\nSummary of Manually Imputed Data (Iteration 10, Re-estimation):\n")
print(summary(data_imputed_manual))
cat("\nSummary of MICE Package Imputed Data (Dataset 1):\n")
print(summary(complete(imputed_data_mice, 1)))

# comparing original data with imputed data
hist(data_to_impute_orig$use)
hist(imputed_data_mice$imp$use[,1])

# Plotting the complete data
plot(y = data_to_impute_orig$use, x = data_to_impute_orig$cc10, main = "Original Data")
abline(lm(data_to_impute_orig$use ~ data_to_impute_orig$cc10), col = "red", lwd = 2)

plot(y = complete(imputed_data_mice, 1)$use, x = complete(imputed_data_mice, 1)$cc10, main = "MICE Imputed Data (one iteration)")
abline(lm(complete(imputed_data_mice, 1)$use ~ complete(imputed_data_mice, 1)$cc10), col = "blue", lwd = 2)
abline(lm(data_to_impute_orig$use ~ data_to_impute_orig$cc10), col = "red", lwd = 2)

## Overall example -- predicting PERF and USE from CC10 and Female with HSL, MSC, MAS, and MSE as auxiliary variables using MICE

data_analysis = as.data.frame(allx)
data_analysis$cc10 = data_analysis$cc - 10 # Create centered cc10

# Select ALL variables needed for imputation AND analysis
vars_for_mice <- c("perf", "use", "female", "cc10", "hsl", "msc", "mas", "mse")
data_to_impute <- data_analysis[, vars_for_mice]

# Optionally introduce some missingness into 'female' and auxiliaries for demo
set.seed(123)
data_to_impute$female[sample(1:nrow(data_to_impute), 30)] <- NA
data_to_impute$hsl[sample(1:nrow(data_to_impute), 20)] <- NA
data_to_impute$msc[sample(1:nrow(data_to_impute), 25)] <- NA

cat("\n--- Missing Data Pattern Before Imputation ---\n")
md.pattern(data_to_impute, plot=FALSE)

# --- 2. Imputation Step ---
cat("\n--- Starting Imputation using mice ---\n")

# Set imputation methods
methods <- make.method(data_to_impute)
methods["female"] <- "logreg"
data_to_impute$female = as.factor(data_to_impute$female)

print("Imputation methods:")
print(methods)

# Perform imputation
n_imputations <- 100
mice_output <- mice(data_to_impute,
                    m = n_imputations,
                    maxit = 10, # Increase maxit for real analysis
                    method = methods,
                    seed = 500,
                    printFlag = TRUE)

cat("\n--- Imputation Finished ---\n")

# --- 3. Estimation Step (using lavaanList) ---
cat("\n--- Starting Estimation using lavaanList ---\n")

# Define the lavaan analysis model syntax
analysis_model_syntax <- "
  # Regressions:
  perf ~ 1 + female + cc10
  use  ~ 1 + female + cc10

  # Residual Variances:
  perf ~~ perf
  use  ~~ use

  # Residual Covariance:
  perf ~~ use

"

# Fit the model to the imputed datasets using lavaanList
# This function iterates through the imputed datasets stored within the mice object
# Note: lavaanList requires the mice object directly, not a list of data frames
fit_lavaanList <- lavaan.mi(model = analysis_model_syntax, 
                            data = mice_output
                            ) # Pass the mids object
        
summary(fit_lavaanList, fit.measures = TRUE, standardize=TRUE)
parameterEstimates.mi(fit_lavaanList) 

# compare pooled estimates to data with missing values
missingData.fit = lavaan(analysis_model_syntax, data=data_to_impute, estimator = "MLR", mimic = "MPLUS")
summary(missingData.fit, fit.measures = TRUE, standardize=TRUE)
parameterEstimates(missingData.fit) 

# Define the lavaan analysis model syntax -- adding predictors to the likelihood function
analysis_model_syntax_ML <- "
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

# remove factor from female variable
data_to_impute$female = as.numeric(data_to_impute$female) - 1


ml.fit = lavaan(analysis_model_syntax_ML, data=data_to_impute, estimator = "MLR", mimic = "MPLUS")
summary(missingData.fit, fit.measures = TRUE, standardize=TRUE)
parameterEstimates(missingData.fit) 



