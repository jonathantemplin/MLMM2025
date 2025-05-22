#SET WORKING DIRECTORY TO LOCATION OF DATA FILE
setwd(myDataLocation)

#set scientific notation to off
options(scipen=999)

#AUTOMATING PACKAGES NEEDED FOR ANALYSES--------------------------------------------------------------------
needed_packages = c("lavaan","semPlot", "psych", "semTools")

for (i in 1:length(needed_packages)){
  haspackage = require(needed_packages[i], character.only = TRUE)
  if (haspackage==FALSE){
    install.packages(needed_packages[i])
    library(needed_packages[i], character.only = TRUE)
  }
  haspackage = require(needed_packages[i], character.only = TRUE)
}

#Showing Distribution of Sum Scores------------------------------------------------------------------------
data01 = read.csv(file="gamblingdataR.csv", na.strings = "99")

data01$GRIsum = data01$GRI1+data01$GRI3+data01$GRI5
hist(data01$GRIsum, main="Three-item GRI Sum Score", xlab="Sum Score")

#Factor score distribution plot------------------------------------------------------------------------------
par(mfrow = c(1,2))
plot(seq(-4,4,.01), dnorm(seq(-4,4,.01)), type="l", main="Posterior Distribution: Normal", ylab="f(Y|F)", xlab="Factor Score F", lwd=2)
lines(c(0,0), c(0,dnorm(0)), lty = 2, col=2, lwd=3)
lines(c(-5,0), c(dnorm(0),dnorm(0)), lty = 2, col=4, lwd=3)

plot(seq(0,8,.01)-4, dchisq(seq(0,8,.01), df=4), type="l", main="Posterior Distribution: Not Quite Normal", ylab="f(Y|F)", xlab="Factor Score F", lwd=2)
lines(c(0,0), c(0,dchisq(4, df=4)), lty = 2, col=2, lwd=3)
lines(c(0,5), c(dchisq(4,df=4),dchisq(4,df=4)), lty = 2, col=2, lwd=3)
lines(c(-2,-2), c(-2,dchisq(2,df=4)), lty = 2, col=4, lwd=3)
lines(c(-2,5), c(dchisq(2,df=4),dchisq(2,df=4)), lty = 2, col=4, lwd=3)


#FIRST EXAMPLE: COMPARING SUM SCORE WITH FACTOR SCORE--------------------------------------------------------

model01.lavaan = "
  GAMBLING =~ (LOADING)*GRI1+(LOADING)*GRI3+(LOADING)*GRI5
  
  GRI1 ~~ (UVAR)*GRI1
  GRI3 ~~ (UVAR)*GRI3
  GRI5 ~~ (UVAR)*GRI5

  GAMBLING ~~ GAMBLING

  rho := ( (3*LOADING)^2 )/(((3*LOADING)^2)+3*UVAR)
"
model01.fit = sem(model01.lavaan, data=data01, estimator = "MLR", mimic="Mplus", fixed.x=FALSE, std.lv=TRUE)
summary(model01.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

#get factor score estimates from the predict function:
model01.factorscores = predict(model01.fit)

#compare both on plot:
par(mfrow = c(1,1))
plot(model01.factorscores, data01$GRIsum, type ="o", lwd=3)

#compare with correlation
cor(model01.factorscores, data01$GRIsum)

#Comparing scores----------------------------------------------------------------------------------------------
score_mat = NULL

#Factor score of a person with GRI1==1, GRI3==1, & GRI5==1 ---- Sum Score = 3
person111_id = data01[data01$GRI1==1 &data01$GRI3==1 & data01$GRI5==1,]$ID[1]
model01.factorscores[person111_id]
score_mat = cbind(1, 1, 1, 3, model01.factorscores[person111_id])

#Factor score of a person with GRI1==1, GRI3==1, & GRI5==2 ---- Sum Score = 4
person112_id = data01[data01$GRI1==1 &data01$GRI3==1 & data01$GRI5==2,]$ID[1]
model01.factorscores[person112_id]
score_mat = rbind(score_mat,cbind(1, 1, 2, 4, model01.factorscores[person112_id]))

#Factor score of a person with GRI1==1, GRI3==2, & GRI5==1 ---- Sum Score = 4
person121_id = data01[data01$GRI1==1 &data01$GRI3==2 & data01$GRI5==1,]$ID[1]
model01.factorscores[person121_id]
score_mat = rbind(score_mat,cbind(1, 2, 1, 4, model01.factorscores[person121_id]))

#Factor score of a person with GRI1==2, GRI3==1, & GRI5==1 ---- Sum Score = 4
person211_id = data01[data01$GRI1==2 &data01$GRI3==1 & data01$GRI5==1,]$ID[1]
model01.factorscores[person211_id]
score_mat = rbind(score_mat,cbind(2, 1, 1, 4, model01.factorscores[person211_id]))

#Difference between factor scores for sum score 3 vs. sum score 4:
model01.factorscores[person112_id] - model01.factorscores[person111_id]

#Factor score of a person with GRI1==1, GRI3==1, & GRI5==3 ---- Sum Score = 5
person113_id = data01[data01$GRI1==1 &data01$GRI3==1 & data01$GRI5==3,]$ID[1]
model01.factorscores[person113_id]
score_mat = rbind(score_mat,cbind(1, 1, 3, 5, model01.factorscores[person113_id]))

#Difference between factor scores for sum score 3 vs. sum score 4:
model01.factorscores[person113_id] - model01.factorscores[person112_id]

#loading matrices for factor scores----------------------------------------------------------------------------

#getting more decimal places from model estimates estimates
parameterEstimates(model01.fit)$est 

#saving into matrices:
lambda = matrix(.5759246, nrow=3, ncol=1)
psi = diag(rep(.5860709,times=3))
mu = matrix(c(1.8226048, 1.5479042, 1.5928144),nrow=3, ncol = 1)
phi = matrix(1, nrow=1, ncol=1)
mu_f = matrix(0, nrow=1, ncol=1)

sigma = lambda%*%t(lambda) + psi

x = matrix(cbind(data01$GRI1, data01$GRI3, data01$GRI5), ncol=3)

#getting mean and variance of factor scores from slide 35:
scores = t(phi %*% t(lambda) %*% solve(sigma)%*%(t(x) - mu%*%matrix(1,nrow=1, ncol=dim(x)[1])))
varscores = phi - phi %*% t(lambda) %*% solve(sigma) %*% lambda %*% phi

#the standard error of the factor score:
sqrt(varscores)

#showing they match with lavaan's estimated scores
plot(scores, model01.factorscores)

#factor score reliability
print("#factor score reliability")
1/(1+varscores)

#Tau Equivalent Model-----------------------------------------------------------------------
model02.lavaan = "
  GAMBLING =~ (LOADING)*GRI1+(LOADING)*GRI3+(LOADING)*GRI5
  
  GRI1 ~~ (UVAR1)*GRI1
  GRI3 ~~ (UVAR3)*GRI3
  GRI5 ~~ (UVAR5)*GRI5

  GAMBLING ~~ GAMBLING

  rho := ( (3*LOADING)^2 )/(((3*LOADING)^2)+UVAR1+UVAR3+UVAR5)
"
model02.fit = sem(model02.lavaan, data=data01, estimator = "MLR", mimic="Mplus", fixed.x=FALSE, std.lv=TRUE)
summary(model02.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

#get factor score estimates from the predict function:
model02.factorscores = predict(model02.fit)

#compare both on plot:
plot(model02.factorscores, data01$GRIsum)

#compare with correlation
cor(model02.factorscores, data01$GRIsum)

score_mat2 = NULL

#Factor score of a person with GRI1==1, GRI3==1, & GRI5==1 ---- Sum Score = 3
person111_id = data01[data01$GRI1==1 &data01$GRI3==1 & data01$GRI5==1,]$ID[1]
model02.factorscores[person111_id]
score_mat2 = model02.factorscores[person111_id]

#Factor score of a person with GRI1==1, GRI3==1, & GRI5==2 ---- Sum Score = 4
person112_id = data01[data01$GRI1==1 &data01$GRI3==1 & data01$GRI5==2,]$ID[1]
model02.factorscores[person112_id]
score_mat2 = rbind(score_mat2,model02.factorscores[person112_id])

#Factor score of a person with GRI1==1, GRI3==2, & GRI5==1 ---- Sum Score = 4
person121_id = data01[data01$GRI1==1 &data01$GRI3==2 & data01$GRI5==1,]$ID[1]
model02.factorscores[person121_id]
score_mat2 = rbind(score_mat2,model02.factorscores[person121_id])

#Factor score of a person with GRI1==2, GRI3==1, & GRI5==1 ---- Sum Score = 4
person211_id = data01[data01$GRI1==2 &data01$GRI3==1 & data01$GRI5==1,]$ID[1]
model02.factorscores[person211_id]
score_mat2 = rbind(score_mat2,model02.factorscores[person211_id])

#Factor score of a person with GRI1==1, GRI3==1, & GRI5==3 ---- Sum Score = 5
person113_id = data01[data01$GRI1==1 &data01$GRI3==1 & data01$GRI5==3,]$ID[1]
model02.factorscores[person113_id]
score_mat2 = rbind(score_mat2,model02.factorscores[person113_id])

scoremat = cbind(score_mat,score_mat2)
rownames(scoremat) = rep("",5)
colnames(scoremat) = c("GRI1", "GRI3", "GRI5", "SUMSC", "FS-PI", "FS-TE")
scoremat
#Compare reliabilty to that of coefficient alpha
alpha(data.frame(x))

#getting more decimal places from model estimates estimates
parameterEstimates(model02.fit)$est 

#saving into matrices:
lambda = matrix(.5666655, nrow=3, ncol=1)
psi = diag(c(0.7025372, 0.4679966, 0.6034985))
mu = matrix(c(1.8226048, 1.5479042, 1.5928144),nrow=3, ncol = 1)
phi = matrix(1, nrow=1, ncol=1)
mu_f = matrix(0, nrow=1, ncol=1)

sigma = lambda%*%t(lambda) + psi

x = matrix(cbind(data01$GRI1, data01$GRI3, data01$GRI5), ncol=3)

#getting mean and variance of factor scores from slide 35:
scores = t(phi %*% t(lambda) %*% solve(sigma)%*%(t(x) - mu%*%matrix(1,nrow=1, ncol=dim(x)[1])))
varscores = phi - phi %*% t(lambda) %*% solve(sigma) %*% lambda %*% phi

#the standard error of the factor score:
sqrt(varscores)

#showing they match with lavaan's estimated scores
plot(scores, model02.factorscores)

#factor score reliability
print("#factor score reliability")
1/(1+varscores)




#Unrestricted CFA Model-----------------------------------------------------------------------
model03.lavaan = "
GAMBLING =~ (LOADING1)*GRI1+(LOADING3)*GRI3+(LOADING5)*GRI5

GRI1 ~~ (UVAR1)*GRI1
GRI3 ~~ (UVAR3)*GRI3
GRI5 ~~ (UVAR5)*GRI5

GAMBLING ~~ GAMBLING

rho := ( (LOADING1+LOADING3+LOADING5)^2 )/(((LOADING1+LOADING3+LOADING5)^2)+UVAR1+UVAR3+UVAR5)
info1 := LOADING1^2/UVAR1
info3 := LOADING3^2/UVAR3
info5 := LOADING5^2/UVAR5
"

model03.fit = sem(model03.lavaan, data=data01, estimator = "MLR", mimic="Mplus", fixed.x=FALSE, std.lv=TRUE)
summary(model03.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

#compare with models 1 and 2
anova(model01.fit, model02.fit, model03.fit)

#get factor score estimates from the predict function:
model03.factorscores = predict(model03.fit)

#compare both on plot:
plot(model03.factorscores, data01$GRIsum)

#compare with correlation
cor(model03.factorscores, data01$GRIsum)

score_mat2 = NULL

#Factor score of a person with GRI1==1, GRI3==1, & GRI5==1 ---- Sum Score = 3
person111_id = data01[data01$GRI1==1 &data01$GRI3==1 & data01$GRI5==1,]$ID[1]
model03.factorscores[person111_id]
score_mat2 = model03.factorscores[person111_id]

#Factor score of a person with GRI1==1, GRI3==1, & GRI5==2 ---- Sum Score = 4
person112_id = data01[data01$GRI1==1 &data01$GRI3==1 & data01$GRI5==2,]$ID[1]
model03.factorscores[person112_id]
score_mat2 = rbind(score_mat2,model03.factorscores[person112_id])

#Factor score of a person with GRI1==1, GRI3==2, & GRI5==1 ---- Sum Score = 4
person121_id = data01[data01$GRI1==1 &data01$GRI3==2 & data01$GRI5==1,]$ID[1]
model03.factorscores[person121_id]
score_mat2 = rbind(score_mat2,model03.factorscores[person121_id])

#Factor score of a person with GRI1==2, GRI3==1, & GRI5==1 ---- Sum Score = 4
person211_id = data01[data01$GRI1==2 &data01$GRI3==1 & data01$GRI5==1,]$ID[1]
model03.factorscores[person211_id]
score_mat2 = rbind(score_mat2,model03.factorscores[person211_id])

#Factor score of a person with GRI1==1, GRI3==1, & GRI5==3 ---- Sum Score = 5
person113_id = data01[data01$GRI1==1 &data01$GRI3==1 & data01$GRI5==3,]$ID[1]
model03.factorscores[person113_id]
score_mat2 = rbind(score_mat2,model03.factorscores[person113_id])

scoremat = cbind(scoremat,score_mat2)
rownames(scoremat) = rep("",5)
colnames(scoremat) = c("GRI1", "GRI3", "GRI5", "SUMSC", "FS-PI", "FS-TE", "FS-CFA")
scoremat

#getting more decimal places from model estimates estimates
parameterEstimates(model03.fit)$est 

#saving into matrices:
lambda = matrix(c(.6380921, .4631759, .6351938), nrow=3, ncol=1)
psi = diag(c(.6474473, 0.5346703, 0.5459981))
mu = matrix(c(1.8226048, 1.5479042, 1.5928144),nrow=3, ncol = 1)
phi = matrix(1, nrow=1, ncol=1)
mu_f = matrix(0, nrow=1, ncol=1)

sigma = lambda%*%t(lambda) + psi

x = matrix(cbind(data01$GRI1, data01$GRI3, data01$GRI5), ncol=3)

#getting mean and variance of factor scores from slide 35:
scores = t(phi %*% t(lambda) %*% solve(sigma)%*%(t(x) - mu%*%matrix(1,nrow=1, ncol=dim(x)[1])))
varscores = phi - phi %*% t(lambda) %*% solve(sigma) %*% lambda %*% phi

#the standard error of the factor score:
sqrt(varscores)

#showing they match with lavaan's estimated scores
plot(scores, model03.factorscores)

#factor score reliability
print("#factor score reliability")
1/(1+varscores)



#EXAMPLE OF SECONDARY ANALYSIS INCORPORATING SCORES WITH MULTIPLE IMPUTATION------------------------
data02 = read.csv(file="happinessdataR.csv")

#Simultaneous Analysis with SEM------------------------------------------------
example_sem_analysis.syntax = "
happiness =~ X1 + X2 + X3 + X4 + X5
happiness ~ married
"

example_sem_analysis.fit = sem(example_sem_analysis.syntax, data=data02, estimator = "MLR", mimic="Mplus", fixed.x=FALSE, std.lv=TRUE)
summary(example_sem_analysis.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)
standardizedSolution(example_sem_analysis.fit, type="std.nox")


#Sum Score Secondary Analysis  ------------------------------------------------
example_sumscore_analysis.syntax = "
happiness_sumscore ~ married
"

example_sumscore_analysis.fit = sem(example_sumscore_analysis.syntax, data=data02, 
                                    estimator = "MLR", mimic="Mplus", fixed.x=FALSE, std.lv=TRUE)
summary(example_sumscore_analysis.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)
standardizedSolution(example_sumscore_analysis.fit, type="std.nox")

#Factor Score Secondary Analysis  ------------------------------------------------

#part 1: obtaining the factor score with only the measurement model
example_factorscore_pre_analysis.syntax = "
happiness =~ X1 + X2 + X3 + X4 + X5
"


example_factorscore_pre_analysis.fit = sem(example_factorscore_pre_analysis.syntax, data=data02, 
                                           estimator = "MLR", mimic="Mplus", fixed.x=FALSE, std.lv=TRUE)
#check model fit
summary(example_factorscore_pre_analysis.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

semPaths(example_factorscore_pre_analysis.fit,intercepts = FALSE, residuals = TRUE, style="mx", layout="tree", rotation=1,
         optimizeLatRes = TRUE, whatLabels="std")

#model fit is good--get factor scores
data02$happiness_factorscore = predict(example_factorscore_pre_analysis.fit)

#use factor scores as observed variables
example_factorscore_analysis.syntax = "
happiness_factorscore ~ married
"

example_factorscore_analysis.fit = sem(example_factorscore_analysis.syntax, data=data02, 
                                       estimator = "MLR", mimic="Mplus", fixed.x=FALSE, std.lv=TRUE)

summary(example_factorscore_analysis.fit, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)

standardizedSolution(example_factorscore_analysis.fit, type="std.nox")

#Multiple Imputation with Factor Scores------------------------------------------

#preliminary, we need factor score variance for each person (should be the same by CFA...but if missing data is present won't be)

#saving into matrices:
lambda = inspect(example_factorscore_pre_analysis.fit, what="coef")$lambda
psi = inspect(example_factorscore_pre_analysis.fit, what="coef")$theta
mu = inspect(example_factorscore_pre_analysis.fit, what="coef")$nu
phi = inspect(example_factorscore_pre_analysis.fit, what="coef")$psi
mu_f = inspect(example_factorscore_pre_analysis.fit, what="coef")$alpha

sigma = lambda%*%t(lambda) + psi

x = matrix(cbind(data02$X1, data02$X2, data02$X3, data02$X4, data02$X5), ncol=5)

#getting mean and variance of factor scores from slide 35:
scores = t(phi %*% t(lambda) %*% solve(sigma)%*%(t(x) - mu%*%matrix(1,nrow=1, ncol=dim(x)[1])))
varscores = phi - phi %*% t(lambda) %*% solve(sigma) %*% lambda %*% phi

#checking accuracy
plot(scores, data02$happiness_factorscore)

#Step #1: impute multiple factor scores per person into new data frames
n_imputations = 1000

nobs = dim(data02)[1]
imputed_data_list = list()
imputed_data = list()

example_observation = matrix(NA, nrow=n_imputations, ncol=1)
for (i in 1:n_imputations){
  #copy data from data02
  newimputeddata = data02
  
  #draw random normal variables:
  imputed_factor_score = rnorm(n = nobs, mean=0, sd=1)
  
  #transform factor score onto distribution for each person:
  imputed_factor_score = scores + imputed_factor_score*sqrt(varscores)
  
  #change original factor score to imputed value
  newimputeddata$happiness_factorscore = imputed_factor_score 
  
  #add data to list
  imputed_data[[i]] = newimputeddata
  
  #add data to example observation:
  example_observation[i,1] = newimputeddata$happiness_factorscore[1]
}

par(mfrow=c(2,1))
#Sequence of Imputed Factor Scores for Observation #1
plot(1:1000,example_observation[,1], type="l", xlab="Imputation Number", ylab="Imputed Value", main="Factor Score Imputation: Obs #1")

plot(density(example_observation[,1]), main="Density of Imputed Factor Scores: Obs #1", lwd=3)
lines(seq(-4,4,.01),dnorm(seq(-4,4,.01), mean=scores[1], sd=sqrt(varscores)), lwd=3, col=2, lty=2)
legend(x=-1.75,y=.55, c("Density of Imputed", "Density of Theoretical"), lwd=c(3,3), col=c(1,2), lty=c(1,2))

#Step #2 and 3: analyze data and pool results--from semTools package
mi_analysis = runMI(data=imputed_data, model=example_factorscore_analysis.syntax, fun="sem")

#finally: report results
summary(mi_analysis, fit.measures=TRUE, rsquare=TRUE, standardized=TRUE)
standardizedSolution(mi_analysis, type="std.nox")

