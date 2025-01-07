
# Clear all existing objects from the workspace
rm(list=ls()) 

setwd("C:/Users/euedu/Downloads/git")

source("gamlss_BP.R")
source("Residual_H_Log_Like_BP.R") 

# Packages
library(gamlss)       # Generalized additive models for location, scale, and shap
library(extraDistr)   # Additional Univariate and Multivariate Distributions

NREP <- 10000     # Monte Carlo replicates
n <- 200          # Sample size

# True values of the parameters
beta <- c(2.9, -2.5)         
varphi <- c(1.5,1.8)      

set.seed(2022)  # Seed 
# Generate model regressors
x0 <- rep(1,n)
x1 <- runif(n,-1,1)
x2 <- x1^2
z0 <- rep(1,n)
z1 <- runif(n)  
X <- cbind(x0,x2)     # Mean regressor matrix
Z <- cbind(z0,z1)     # Precision regressor matrix

# Linear predictor vectors 
eta1 <- X%*%beta
eta2 <-  Z%*%varphi

# link fuction - log  
mu <- as.vector(exp(eta1))
phi <- as.vector(exp(eta2))

# Variance of response variable
Vmu <-  mu*(1+mu)   # Variance function
vary <-  Vmu/phi

# Skewness of response variable
ske <- (2*(1+phi)*(1+2*mu) / (phi-1)) * sqrt(phi/(Vmu*((1+phi)^2)))

# Kurtosis of response variable
kurt <- 6*((5*phi -1)/((phi-2)*(phi-1)) + (phi/(Vmu*(phi-2)*(phi-1))))

cat("Summary mu\n"); print(summary(mu))
cat("Summary phi\n"); print(summary(phi))
cat("Summary Variance of y\n"); print(summary(vary))
cat("Summary Skewness of y\n"); print(summary(ske))
cat("Summary kurtosis of y\n"); print(summary(kurt))

# Generate the response variable  
y <- rBP(n,mu,phi)

# Estimate model - correct specification
fit1 <- gamlss(y~x2,~z1,family = BP(mu.link = "log", sigma.link = "log"), trace=FALSE)

#Estimate model - incorrect specification  
fit2 <- gamlss(y~x1,~z1,family = BP(mu.link = "log", sigma.link = "log"), trace=FALSE)

# Simulated envelope plots - correct model
envelope.BP(fit1, k=100, type = "quantile")
envelope.BP(fit1, k=100, type = "pearson")
envelope.BP(fit1, k=100, type = "pearson P")
envelope.BP(fit1, k=100, type = "sweighted1")
envelope.BP(fit1, k=100, type = "sweighted2")
envelope.BP(fit1, k=100, type = "variance")
envelope.BP(fit1, k=100, type = "combined")
envelope.BP(fit1, k=100, type = "deviance")
envelope.BP(fit1, k=100, type = "anscombe")
envelope.BP(fit1, k=100, type = "williams")

# index plots - correct model
plot.BP(fit1, which = 5, type = "quantile") 
plot.BP(fit1, which = 5, type = "pearson")
plot.BP(fit1, which = 5, type = "pearson P")
plot.BP(fit1, which = 5, type = "sweighted1")
plot.BP(fit1, which = 5, type = "sweighted2")
plot.BP(fit1, which = 5, type = "variance")
plot.BP(fit1, which = 5, type = "combined")
plot.BP(fit1, which = 5, type = "deviance")
plot.BP(fit1, which = 5, type = "anscombe")
plot.BP(fit1, which = 5, type = "williams")

# Simulated envelope plots - incorrect model
envelope.BP(fit2, k=100, type = "quantile")
envelope.BP(fit2, k=100, type = "pearson")
envelope.BP(fit2, k=100, type = "pearson P")
envelope.BP(fit2, k=100, type = "sweighted1")
envelope.BP(fit2, k=100, type = "sweighted2")
envelope.BP(fit2, k=100, type = "variance")
envelope.BP(fit2, k=100, type = "combined")
envelope.BP(fit2, k=100, type = "deviance")
envelope.BP(fit2, k=100, type = "anscombe")
envelope.BP(fit2, k=100, type = "williams")

# index plots - incorrect model
plot.BP(fit2, which = 5, type = "quantile") 
plot.BP(fit2, which = 5, type = "pearson")
plot.BP(fit2, which = 5, type = "pearson P")
plot.BP(fit2, which = 5, type = "sweighted1")
plot.BP(fit2, which = 5, type = "sweighted2")
plot.BP(fit2, which = 5, type = "variance")
plot.BP(fit2, which = 5, type = "combined")
plot.BP(fit2, which = 5, type = "deviance")
plot.BP(fit2, which = 5, type = "anscombe")
plot.BP(fit2, which = 5, type = "williams")

# Simulation
cont <- 0       # Convergence failure counts

# Matrices to store the p-values of the simulation - correct model
contp <- matrix(0,NREP,10)
pvalue <- matrix(NA,NREP,10)

# Matrices to store the p-values of the simulation - correct model
contp1 <- matrix(0,NREP,10)
pvalue1 <- matrix(NA,NREP,10)

# Matrices to store the residuals of the simulation - correct model
quantile <- matrix(NA,NREP,n)
pearson <- matrix(NA,NREP,n)
pearsonp <- matrix(NA,NREP,n)
sweighted1 <- matrix(NA,NREP,n)
sweighted2 <- matrix(NA,NREP,n)
variance <- matrix(NA,NREP,n)
combined <- matrix(NA,NREP,n)
deviance <- matrix(NA,NREP,n)
anscombe <- matrix(NA,NREP,n)
williams <- matrix(NA,NREP,n)

# Matrices to store the residuals of the simulation - incorrect model
quantile1 <- matrix(NA,NREP,n)
pearson1 <- matrix(NA,NREP,n)
pearsonp1 <- matrix(NA,NREP,n)
sweighted11 <- matrix(NA,NREP,n)
sweighted21 <- matrix(NA,NREP,n)
variance1 <- matrix(NA,NREP,n)
combined1 <- matrix(NA,NREP,n)
deviance1 <- matrix(NA,NREP,n)
anscombe1 <- matrix(NA,NREP,n)
williams1 <- matrix(NA,NREP,n)

# Monte Carlo Simulation
i <- 1
while(i <= NREP)
{ 
  # Print replica ratio
  prop = (i/NREP)*100 
  if(prop==0 || prop==10 || prop==20 || prop==30 || prop==40 || prop==50 || prop==60 || prop==70 
     || prop==80|| prop== 90 || prop==100) cat(paste(prop,"%"),"\n")
  
  # Generate the response variable 
  y <- rBP(n,mu,phi)
  
  # tryCatch to catch errors
  tryCatch({
    
    # Estimate model - correct specification 
    fit1 <- gamlss(y~x2,~z1,family = BP(mu.link = "log", sigma.link = "log"), trace=FALSE)
    
    # Estimate model - incorrect specification 
    fit2 <- gamlss(y~x1,~z1,family = BP(mu.link = "log", sigma.link = "log"), trace=FALSE)
    
    # if converge    
    if(fit1$converged == TRUE && fit2$converged == TRUE) 
    {
      # Residuals of the correct model
      quantile[i,] <- residuals.BP(fit1, type = "quantile")       # Quantile residual
      pearson[i,] <- residuals.BP(fit1, type = "pearson")         # Pearson residual
      pearsonp[i,] <- residuals.BP(fit1, type = "pearson P")      # Standardized Pearson residual
      sweighted1[i,] <- residuals.BP(fit1, type = "sweighted1")   # Weighted residual
      sweighted2[i,] <- residuals.BP(fit1, type = "sweighted2")   # Standardized weighted residual
      variance[i,] <- residuals.BP(fit1, type = "variance")       # Variance residual
      combined[i,] <- residuals.BP(fit1, type = "combined")       # Combined residual
      deviance[i,] <- residuals.BP(fit1, type = "deviance")       # Deviance residual
      anscombe[i,] <- residuals.BP(fit1, type = "anscombe")       # Anscombe residual
      williams[i,] <- residuals.BP(fit1, type = "williams")       # Williams residual
      
      # Residuals of the incorrect model
      quantile1[i,] <- residuals.BP(fit2, type = "quantile")       # Quantile residual
      pearson1[i,] <- residuals.BP(fit2, type = "pearson")         # Pearson residual
      pearsonp1[i,] <- residuals.BP(fit2, type = "pearson P")      # Standardized Pearson residual
      sweighted11[i,] <- residuals.BP(fit2, type = "sweighted1")   # Weighted residual
      sweighted21[i,] <- residuals.BP(fit2, type = "sweighted2")   # Standardized weighted residual
      variance1[i,] <- residuals.BP(fit2, type = "variance")       # Variance residual
      combined1[i,] <- residuals.BP(fit2, type = "combined")       # Combined residual
      deviance1[i,] <- residuals.BP(fit2, type = "deviance")       # Deviance residual
      anscombe1[i,] <- residuals.BP(fit2, type = "anscombe")       # Anscombe residual
      williams1[i,] <- residuals.BP(fit2, type = "williams")       # Williams residual
      
      # Shapiro-Wilk tests - Residuals of the correct model
      pvalue[i,1] <- shapiro.test(quantile[i,])$p.value
      pvalue[i,2] <- shapiro.test(pearson[i,])$p.value
      pvalue[i,3] <- shapiro.test(pearsonp[i,])$p.value
      pvalue[i,4] <- shapiro.test(sweighted1[i,])$p.value
      pvalue[i,5] <- shapiro.test(sweighted2[i,])$p.value
      pvalue[i,6] <- shapiro.test(variance[i,])$p.value
      pvalue[i,7] <- shapiro.test(combined[i,])$p.value
      pvalue[i,8] <- shapiro.test(deviance[i,])$p.value
      pvalue[i,9] <- shapiro.test(anscombe[i,])$p.value
      pvalue[i,10] <- shapiro.test(williams[i,])$p.value
      
      # Shapiro-Wilk tests - Residuals of the incorrect model
      pvalue1[i,1] <- shapiro.test(quantile1[i,])$p.value
      pvalue1[i,2] <- shapiro.test(pearson1[i,])$p.value
      pvalue1[i,3] <- shapiro.test(pearsonp1[i,])$p.value
      pvalue1[i,4] <- shapiro.test(sweighted11[i,])$p.value
      pvalue1[i,5] <- shapiro.test(sweighted21[i,])$p.value
      pvalue1[i,6] <- shapiro.test(variance1[i,])$p.value
      pvalue1[i,7] <- shapiro.test(combined1[i,])$p.value
      pvalue1[i,8] <- shapiro.test(deviance1[i,])$p.value
      pvalue1[i,9] <- shapiro.test(anscombe1[i,])$p.value
      pvalue1[i,10] <- shapiro.test(williams1[i,])$p.value
      
      # Count of p-values (significance levels of 5%) - Residuals of the correct model
      if(pvalue[i,1]>0.05) contp[i,1] <- 1
      if(pvalue[i,2]>0.05) contp[i,2] <- 1
      if(pvalue[i,3]>0.05) contp[i,3] <- 1
      if(pvalue[i,4]>0.05) contp[i,4] <- 1
      if(pvalue[i,5]>0.05) contp[i,5] <- 1
      if(pvalue[i,6]>0.05) contp[i,6] <- 1
      if(pvalue[i,7]>0.05) contp[i,7] <- 1
      if(pvalue[i,8]>0.05) contp[i,8] <- 1
      if(pvalue[i,9]>0.05) contp[i,9] <- 1
      if(pvalue[i,10]>0.05) contp[i,10] <- 1
      
      # Count of p-values (significance levels of 5%) - Residuals of the incorrect model
      if(pvalue1[i,1]>0.05) contp1[i,1] <- 1
      if(pvalue1[i,2]>0.05) contp1[i,2] <- 1
      if(pvalue1[i,3]>0.05) contp1[i,3] <- 1
      if(pvalue1[i,4]>0.05) contp1[i,4] <- 1
      if(pvalue1[i,5]>0.05) contp1[i,5] <- 1
      if(pvalue1[i,6]>0.05) contp1[i,6] <- 1
      if(pvalue1[i,7]>0.05) contp1[i,7] <- 1
      if(pvalue1[i,8]>0.05) contp1[i,8] <- 1
      if(pvalue1[i,9]>0.05) contp1[i,9] <- 1
      if(pvalue1[i,10]>0.05) contp1[i,10] <- 1
      
      i <- i + 1
      
    }else{ # if it does not converge
      cont <- cont + 1
      print(c("Non-convergence",i,cont))
      i <- i - 1 
    }
  }, error = function(e) {
    # If an error occurs, print the message and repeat the loop
    print(paste("Error in the replica", i, ":", e$message))
    cont <- cont + 1
    i <- i - 1 
  })
}

# percentages of the p-values - correct model
colSums(contp)/NREP*100
# percentages of the p-values - incorrect model 
colSums(contp1)/NREP*100

# Histograms of the p-value - correct model
hist(pvalue[,1], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue[,2], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue[,3], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue[,4], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue[,5], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue[,6], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue[,7], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue[,8], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue[,9], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue[,10], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')

# Histograms of the p-value - incorrect model
hist(pvalue1[,1], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue1[,2], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue1[,3], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue1[,4], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue1[,5], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue1[,6], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue1[,7], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue1[,8], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue1[,9], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
hist(pvalue1[,10], xlab = "p-value", ylab = "Frequency", main = "", col = 'white')
abline(v = 0.05, col = 'red')
