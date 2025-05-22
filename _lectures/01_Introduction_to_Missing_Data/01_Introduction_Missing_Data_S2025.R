# clear environment
rm(list = ls())

# load pain.dat file
painData = read.table("data/pain.dat", header = FALSE, na.strings = "999")

# label variables
names(painData) = c("id", "txgrp", "male", "age", "edugroup", "workhrs", "exercise", "paingrps", 
                "pain", "anxiety", "stress", "control", "depress", "interfere", "disability",
                paste0("dep", seq(1:7)), paste0("int", seq(1:6)), paste0("dis", seq(1:6)))

# create plot of depression vs pain with points shown
plot(painData$pain, painData$depress, pch = 19, xlab = "Pain", ylab = "Depression", main = "Depression vs Pain")

# auxiliary variables and semipartial correlation

# generate data
library(ppcor)
set.seed(42)
n = 100
y = rnorm(n, mean=50, sd=10)   # Outcome variable
x1 = rnorm(n, mean=10, sd=5)    # Predictor 1
x2 = rnorm(n, mean=20, sd=5)    # Predictor 2
aux = .05*y + 0.5*x1 + 0.3*x2 + rnorm(n, mean=0, sd=2)  # Auxiliary variable

data = data.frame(y, x1, x2, aux)

# Step 1: Regress the auxiliary variable on other predictors and obtain residuals
aux_resid = resid(lm(aux ~ x1 + x2, data=data))

# Step 2: Compute the correlation between residuals and the outcome
spr_lm = cor(aux_resid, data$y)
print(spr_lm)

# Compute semi-partial correlation using spcorr package
spr_spcorr = spcor(data)
?spcor.test
print(spr_spcorr)

