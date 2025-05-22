if (!require(psych)) install.packages("psych")
library(psych)

#LOAD INTERNAL R DATA SET
load("dataSexHeightWeight.rda")

#ASSIGN A DIFFERENT NAME FOR THE DATA SET TO MAKE CHANGES
data01 = dataSexHeightWeight

# Creating New variables for weight and height
data01$weightKG = 0.453 * data01$weightLB

# show descriptives
describe(data01)

#show covariance matrix
cov(data01[c("heightIN", "weightKG", "weightLB")])

#show correlation matrix
cor(data01[c("heightIN", "weightKG", "weightLB")])

# centering data for model 05:
#CREATE NEW VARIABLE WHERE FEMALE=1 (DUMMY CODED)
data01$female[data01$sex=="F"]=1
data01$female[data01$sex=="M"]=0

#center height at mean height
data01$heightIN_MC = data01$heightIN - mean(data01$heightIN)

#model #5: adding interaction to model #4:
model05 = lm(weightLB~heightIN_MC+female+female*heightIN_MC,data=data01)

summary(model05)
anova(model05)

#plotting distributions
xbar1 = 127.5466
sigma = 4.73
xbar2 = 248.021

weight1 = seq(xbar1-15,xbar1+15,.1)
weight2 = seq(xbar2-15,xbar2+15,.1)
fy1 = dnorm(weight1,xbar1,sqrt(sigma))
fy2 = dnorm(weight2,xbar2,sqrt(sigma))

par(mfrow=c(1,2))
plot(weight1,fy1,type="l",xlab="Weight",ylab="f(Y|X)",main="Height=62 Female=1")
text(x=120,y=.15,expression(paste(hat(Y)[p],"=127.5")))
plot(weight2,fy2,type="l",xlab="Weight",ylab="f(Y|X)",main="Height=76 Female=0")
text(x=240,y=.15,expression(paste(hat(Y)[p],"=248.0")))

# plotting residuals of model 05:
plot(model05)

# testing normality with Shapiro test
shapiro.test(model05$residuals)

