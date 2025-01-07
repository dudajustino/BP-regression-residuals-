
# Clear all existing objects from the workspace
rm(list=ls()) 

setwd("C:/Users/euedu/Downloads/git")

source("gamlss_BP.R")
source("Residual_H_Log_Like_BP.R") 

# Packages
library(gamlss)       # Generalized additive models for location, scale, and shap
library(extraDistr)   # Additional Univariate and Multivariate Distributions

NREP <- 10000     # Monte Carlo replicates
n <- 30        # Sample size

set.seed(2022)  # Seed 
# Generate model regressors
x0 <- rep(1,n)  
x1 <- runif(n)
z0 <- rep(1,n)  
z1 <- runif(n)
X <- cbind(x0,x1)  # Mean regressor matrix
Z <- cbind(z0,z1)  # Precision regressor matrix

# Matrices to store the p-values of the simulation
pvalue <- array(NA,dim=c(NREP,10,5))
contp5 <- array(0,dim=c(NREP,10,5))
contp10 <- array(0,dim=c(NREP,10,5))

# Matrices to store the residuals of the simulation
quantile <- array(NA,dim=c(NREP,n,5))
pearson <- array(NA,dim=c(NREP,n,5))
pearsonp <- array(NA,dim=c(NREP,n,5))
sweighted1 <- array(NA,dim=c(NREP,n,5))
sweighted2 <- array(NA,dim=c(NREP,n,5))
variance <- array(NA,dim=c(NREP,n,5))
combined <- array(NA,dim=c(NREP,n,5))
deviance <- array(NA,dim=c(NREP,n,5))
anscombe <- array(NA,dim=c(NREP,n,5))
williams <- array(NA,dim=c(NREP,n,5))

# Matrix to store the true values of the parameters
par <- matrix(0, 5, 4)    # par = (beta0, beta1, nu0, nu1)

# Five scenarios 
par[1,] <- c(3.2, -5.3, 1, 2)
par[2,] <- c(3.2, -5.3, 3, 1.3)
par[3,] <- c(3, 2, 1, 2)
par[4,] <- c(3, 2, 3, 1.3)
par[5,] <- c(0.3, -3.9, 0.7, 0.2)

for(k in 1:5){
  cont <- 0 # convergence failure counts
  
  # Linear predictor vectors  
  eta1 <- X%*%par[k,1:2]
  eta2 <-  Z%*%par[k,3:4]
  
  # link function - log 
  mu <- as.vector(exp(eta1))
  phi <- as.vector(exp(eta2))
  
  # Variance of response variable
  Vmu <-  mu*(1+mu)   # Variance function
  vary <-  Vmu/phi
  
  # Skewness of response variable
  ske <- (2*(1+phi)*(1+2*mu) / (phi-1)) * sqrt(phi/(Vmu*((1+phi)^2)))
  
  # Kurtosis of response variable
  kurt <- 6*((5*phi -1)/((phi-2)*(phi-1)) + (phi/(Vmu*(phi-2)*(phi-1))))
  
  cat("Scenario:", k, "\n")
  cat("Summary mu\n"); print(summary(mu))
  cat("Summary phi\n"); print(summary(phi))
  cat("Summary Variance of y\n"); print(summary(vary))
  cat("Summary Skewness of y\n"); print(summary(ske))
  cat("Summary kurtosis of y\n"); print(summary(kurt))
  
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

    # Estimate the model
    fit <- gamlss(y~x1, sigma.formula=~z1, family = BP(mu.link="log",sigma.link="log"), trace=F)
    
    # if converge    
    if(fit$converged == TRUE)
    {
      # Residuals
      quantile[i,,k] <- residuals.BP(fit, type = "quantile")       # Quantile residual
      pearson[i,,k] <- residuals.BP(fit, type = "pearson")         # Pearson residual
      pearsonp[i,,k] <- residuals.BP(fit, type = "pearson P")      # Standardized Pearson residual
      sweighted1[i,,k] <- residuals.BP(fit, type = "sweighted1")   # Weighted residual
      sweighted2[i,,k] <- residuals.BP(fit, type = "sweighted2")   # Standardized weighted residual
      variance[i,,k] <- residuals.BP(fit, type = "variance")       # Variance residual
      combined[i,,k] <- residuals.BP(fit, type = "combined")       # Combined residual
      deviance[i,,k] <- residuals.BP(fit, type = "deviance")       # Deviance residual
      anscombe[i,,k] <- residuals.BP(fit, type = "anscombe")       # Anscombe residual
      williams[i,,k] <- residuals.BP(fit, type = "williams")       # Williams residual
      
      # Shapiro-Wilk tests
      pvalue[i,1,k] <- shapiro.test(quantile[i,,k])$p.value
      pvalue[i,2,k] <- shapiro.test(pearson[i,,k])$p.value
      pvalue[i,3,k] <- shapiro.test(pearsonp[i,,k])$p.value
      pvalue[i,4,k] <- shapiro.test(sweighted1[i,,k])$p.value
      pvalue[i,5,k] <- shapiro.test(sweighted2[i,,k])$p.value
      pvalue[i,6,k] <- shapiro.test(variance[i,,k])$p.value
      pvalue[i,7,k] <- shapiro.test(combined[i,,k])$p.value
      pvalue[i,8,k] <- shapiro.test(deviance[i,,k])$p.value
      pvalue[i,9,k] <- shapiro.test(anscombe[i,,k])$p.value
      pvalue[i,10,k] <- shapiro.test(williams[i,,k])$p.value

      # Count of p-values (significance levels of 5%)
      if(pvalue[i,1,k]>0.05) contp5[i,1,k] <- 1
      if(pvalue[i,2,k]>0.05) contp5[i,2,k] <- 1
      if(pvalue[i,3,k]>0.05) contp5[i,3,k] <- 1
      if(pvalue[i,4,k]>0.05) contp5[i,4,k] <- 1
      if(pvalue[i,5,k]>0.05) contp5[i,5,k] <- 1
      if(pvalue[i,6,k]>0.05) contp5[i,6,k] <- 1
      if(pvalue[i,7,k]>0.05) contp5[i,7,k] <- 1
      if(pvalue[i,8,k]>0.05) contp5[i,8,k] <- 1
      if(pvalue[i,9,k]>0.05) contp5[i,9,k] <- 1
      if(pvalue[i,10,k]>0.05) contp5[i,10,k] <- 1
      
      # Count of p-values (significance levels of 10%)
      if(pvalue[i,1,k]>0.10) contp10[i,1,k] <- 1
      if(pvalue[i,2,k]>0.10) contp10[i,2,k] <- 1
      if(pvalue[i,3,k]>0.10) contp10[i,3,k] <- 1
      if(pvalue[i,4,k]>0.10) contp10[i,4,k] <- 1
      if(pvalue[i,5,k]>0.10) contp10[i,5,k] <- 1
      if(pvalue[i,6,k]>0.10) contp10[i,6,k] <- 1
      if(pvalue[i,7,k]>0.10) contp10[i,7,k] <- 1
      if(pvalue[i,8,k]>0.10) contp10[i,8,k] <- 1
      if(pvalue[i,9,k]>0.10) contp10[i,9,k] <- 1
      if(pvalue[i,10,k]>0.10) contp10[i,10,k] <- 1
      
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
}

################ residual analysis ##################

# Descriptive statistics
desc <- function(residuos){
  measures = apply(residuos, 2, function(x) c(mean(x), var(x), skewness(x), kurtosis(x)))
  result = c(mean(measures[1,]), mean(measures[2,]), mean(measures[3,]), mean(measures[4,]))
  result
}

measures <- array(NA, dim=c(10,6,5))

for (l in 1:5) {
  measures[,,l] <- cbind(rbind(desc(quantile[,,l]), desc(pearson[,,l]), desc(pearsonp[,,l]), desc(sweighted1[,,l]), 
                              desc(sweighted2[,,l]), desc(variance[,,l]), desc(combined[,,l]), desc(deviance[,,l]), 
                              desc(anscombe[,,l]), desc(williams[,,l])),
                        colSums(contp5[,,l])/NREP*100, colSums(contp10[,,l])/NREP*100)
}

colnames(measures)<- c("Mean", "Variance", "Skewness", "Kurtosis", "SW 5%", "SW 10%")
rownames(measures)<- c("Quantile", "Pearson", "Standardized Pearson", "Weighted", "Standardized Weighted", 
                      "Variance", "Combined", "Deviance", "Anscombe", "Williams")

cat("Descriptive statistics (on average) of residuals for scenario I \n"); round(measures[,,1],4)
cat("Descriptive statistics (on average) of residuals for scenario II \n"); round(measures[,,2],4)
cat("Descriptive statistics (on average) of residuals for scenario III \n"); round(measures[,,3],4)
cat("Descriptive statistics (on average) of residuals for scenario IV \n"); round(measures[,,4],4)
cat("Descriptive statistics (on average) of residuals for scenario V \n"); round(measures[,,5],4)
