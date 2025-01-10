rm(list=ls()) 

source("gamlss_BP.R")
source("Residual_H_Log_Like_BP.R") 

# Packages
library(gamlss)       # Generalized additive models for location, scale, and shap
library(extraDistr)   # Additional Univariate and Multivariate Distributions

############### black cherry tree dataset ##########################

#The data provide measurements of diameter (in inches at 4.5 feet above ground), height (feet),
#and volume (cubic feet) of wood in 31 black cherry trees felled in the Allegheny National Forest in Pennsylvania. 
#Ryan et al. (1976, p. 329).

cherry <- read.table("cherry.dat",header = TRUE)
attach(cherry)

par(mar=c(3.5, 3.2, 1, 1.5)) 
par(mgp=c(2.2, 0.8, 0))
library(AdequacyModel)
TTT(Volume, grid = F)

media <- mean(Volume)
mediana <- median(Volume)
maximo <- max(Volume)
minimo <- min(Volume)
q1 <- quantile(Volume, probs = 0.25)
q3 <- quantile(Volume, probs = 0.75)

par(mar=c(2.5, 2.2, 1, 1.5)) 
par(mgp=c(2.2, 0.8, 0))

boxplot(Volume, ylim = c(minimo, maximo),col = "white")  
points(0.74, media, pch = 19, col = "black",cex = 0.7)
text(1, media, substitute(bar(Y) == mu, list(mu = format(media, digits = 4))), pos = 2, col = "black",offset = 5.3,cex = 0.8)
text(1, mediana, sprintf("Md = %.2f", mediana), pos = 4, col = "black",offset = 5.2,cex = 0.8)
text(1, maximo, sprintf("Max = %.2f", maximo), pos = 4, col = "black",offset = 5.2,cex = 0.8)
text(1, minimo, sprintf("Min = %.2f", minimo), pos = 4, col = "black",offset = 5.2,cex = 0.8)
text(1, q1, sprintf("Q1 = %.2f", q1), pos = 4, col = "black",offset = 5.2,cex = 0.8)
text(1, q3, sprintf("Q3 = %.2f", q3), pos = 4, col = "black",offset = 5.2,cex = 0.8)

# candidate models
fit0 <- gamlss(Volume~Girth+Height,family = BP(mu.link = "log"),trace=FALSE)
summary(fit0)
envelope.BP(fit0, k=100, type = "quantile", link=c("log","log"))

fit1 <- gamlss(Volume~Girth+Height,family = BP(mu.link = "sqrt"),trace=FALSE)
summary(fit1)
envelope.BP(fit1, k=100, type = "quantile", link=c("sqrt","log"))

fit2 <- gamlss(Volume~Girth+Height+I(Girth^2),family = BP(mu.link = "log"),trace=FALSE)
summary(fit2)
envelope.BP(fit2, k=100, type = "quantile", link=c("log","log"))

fit3 <- gamlss(Volume~I(log(Girth))+I(log(Height)),family = BP(mu.link = "sqrt"),trace=FALSE)
summary(fit3)
envelope.BP(fit3, k=100, type = "quantile", link=c("sqrt","log"))

fit4 <- gamlss(Volume~I(log(Girth))+I(log(Height)),family = BP(mu.link = "log"),trace=FALSE)
summary(fit4)
envelope.BP(fit4, k=100, type = "quantile", link=c("log","log"))
plot.BP(fit4, which = 1, type = "quantile", -3,3, pos1 = c(2,3), pos2 = c(3), cex.lab=1.5)
plot.BP(fit4, which = 2, type = "quantile", -2,2, pos1 = c(2,3), pos2 = c(3), cex.lab=1.5)
plot.BP(fit4, which = 3, type = "quantile", -2,2, pos1 = c(2,3), pos2 = c(3), cex.lab=1.5)
plot.BP(fit4, which = 4, type = "quantile", -2,2, pos1 = c(2,3), pos2 = c(1), cex.lab=1.5)

shapiro.test(residuals.BP(fit4,type = "quantile"))
summary(fit4$mu.fv)
summary(fit4$sigma.fv)
