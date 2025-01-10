rm(list=ls()) 

source("gamlss_BP.R")
source("Residual_H_Log_Like_BP.R") 

# Packages
library(gamlss)       # Generalized additive models for location, scale, and shap
library(extraDistr)   # Additional Univariate and Multivariate Distributions

############### black cherry tree dataset ##########################

#dataset available at Weisberg [28] on the average rent paid (in US dollars) per acre planted to alfalfa in 67 
#counties in Minnesota in 1977
#Weisberg, S. (2014). Applied Linear Regression, 4th edition. Hoboken NJ: Wiley. 

library("alr4")

#X1: average rent paid for land for other agricultural purposes
#X2: density of dairy cows in number per square mile
#X3: proportion of agricultural farmland used for pasture
#X4: 1 if liming is necessary for alfalfa cultivation; 0 otherwise
#Y: average rent paid per acre planted to alfalfa 

data(landrent)
head(landrent)
landrent <- subset(landrent, X4 != 0)
head(landrent)
attach(landrent)

par(mar=c(3.5, 3.2, 1, 1.5)) 
par(mgp=c(2.2, 0.8, 0))
library(AdequacyModel)
TTT(Y, grid = F)

media <- mean(Y)
mediana <- median(Y)
maximo <- max(Y)
minimo <- min(Y)
q1 <- quantile(Y, probs = 0.25)
q3 <- quantile(Y, probs = 0.75)

par(mar=c(2.5, 2.2, 1, 1.5)) 
par(mgp=c(2.2, 0.8, 0))

boxplot(Y, ylim = c(minimo, maximo),col = "white")  
points(0.74, media, pch = 19, col = "black",cex = 0.7)
text(1, media, substitute(bar(Y) == mu, list(mu = format(media, digits = 4))), pos = 2, col = "black",offset = 5.3,cex = 0.8)
text(1, mediana, sprintf("Md = %.2f", mediana), pos = 4, col = "black",offset = 5.2,cex = 0.8)
text(1, maximo, sprintf("Max = %.2f", maximo), pos = 4, col = "black",offset = 5.2,cex = 0.8)
text(1, minimo, sprintf("Min = %.2f", minimo), pos = 4, col = "black",offset = 5.2,cex = 0.8)
text(1, q1, sprintf("Q1 = %.2f", q1), pos = 4, col = "black",offset = 5.2,cex = 0.8)
text(1, q3, sprintf("Q3 = %.2f", q3), pos = 4, col = "black",offset = 5.2,cex = 0.8)

# candidate models
fit0 <- gamlss(Y~I(log(X1)),family = BP(mu.link = "log"),trace=FALSE)
summary(fit0)
envelope.BP(fit0, k=100, type = "quantile", link=c("log","log"))

fit1 <- gamlss(Y~I(log(X1))+I(log(X2)),family = BP(mu.link = "sqrt"),trace=FALSE)
summary(fit1)
envelope.BP(fit1, k=100, type = "quantile", link=c("sqrt","log"))

fit2 <- gamlss(Y~I(log(X1))+I(log(X2)),family = BP(mu.link = "log"),trace=FALSE)
summary(fit2)
envelope.BP(fit2, k=100, type = "quantile", link=c("log","log"))

fit3 <- gamlss(Y~I(log(X1))+I(log(X2)),~X2,family = BP(mu.link = "log", sigma.link = "log"),trace=FALSE) 
summary(fit3)

envelope.BP(fit3, k=100, type = "quantile", link=c("log","log"))
plot.BP(fit3, which = 1, type = "quantile", -3,3, pos1 = c(2,3), pos2 = c(3), cex.lab=1.5)
plot.BP(fit3, which = 2, type = "quantile", -3,3, pos1 = c(2,3), pos2 = c(3), cex.lab=1.5)
plot.BP(fit3, which = 3, type = "quantile", -3,3, pos2 = c(3), cex.lab=1.5)
plot.BP(fit3, which = 4, type = "quantile", -2,2, pos1 = c(2,3), pos2 = c(1), cex.lab=1.5)

shapiro.test(residuals.BP(fit3, type="quantile"))
summary(fit3$mu.fv)
summary(fit2$sigma.fv)
