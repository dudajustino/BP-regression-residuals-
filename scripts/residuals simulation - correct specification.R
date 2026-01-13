
# Clear all existing objects from the workspace
rm(list=ls()) 

source("gamlss_BP.R")
source("Residual_H_Log_Like_BP.R") 

# Packages
library(gamlss)       # Generalized additive models for location, scale, and shap
library(extraDistr)   # Additional Univariate and Multivariate Distributions

NREP <- 10000     # Monte Carlo replicates
n <- 20        # Sample size

set.seed(2022)  # Seed 
# Generate model regressors
x0 <- rep(1,n)  
x1 <- runif(n)
z0 <- rep(1,n)  
z1 <- runif(n)
X <- cbind(x0,x1)  # Mean regressor matrix
Z <- cbind(z0,z1)  # Precision regressor matrix

# Matrices to store the residuals of the simulation
quantile <- array(NA,dim=c(NREP,n,5))
pearson <- array(NA,dim=c(NREP,n,5))
pearsonp <- array(NA,dim=c(NREP,n,5))
sweighted1 <- array(NA,dim=c(NREP,n,5))
sweighted2 <- array(NA,dim=c(NREP,n,5))
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
      deviance[i,,k] <- residuals.BP(fit, type = "deviance")       # Deviance residual
      anscombe[i,,k] <- residuals.BP(fit, type = "anscombe")       # Anscombe residual
      williams[i,,k] <- residuals.BP(fit, type = "williams")       # Williams residual
      
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
  result = apply(residuos, 2, 
                  function(x) c(mean(x), var(x), skewness(x), kurtosis(x)))
  result
}

estat_obs <- function(array_res, k, estat = c("mean","var","skew","kurt")){
  estat <- match.arg(estat)
  desc_k <- desc(array_res[,,k])  
  pos <- switch(estat,
                mean = 1,
                var  = 2,
                skew = 3,
                kurt = 4)
  desc_k[pos, ]   
}

table_stat <- function(k, estat){
  tab <- cbind(
    Quantile = estat_obs(quantile,  k, estat),
    Pearson  = estat_obs(pearson,   k, estat),
    PearsonP = estat_obs(pearsonp,  k, estat),
    SW1      = estat_obs(sweighted1,k, estat),
    SW2      = estat_obs(sweighted2,k, estat),
    Deviance = estat_obs(deviance,  k, estat),
    Anscombe = estat_obs(anscombe,  k, estat),
    Williams = estat_obs(williams,  k, estat)
  )

    rownames(tab) <- 1:n
  tab
}

k <- 1
                 
tab_mean <- table_stat(k, "mean")
tab_var  <- table_stat(k, "var")
tab_skew <- table_stat(k, "skew")
tab_kurt <- table_stat(k, "kurt")

xtable(tab_mean, digits = 3)
xtable(tab_var, digits = 3)
xtable(tab_skew, digits = 3)
xtable(tab_kurt, digits = 3)

#  Q-Q plot                 
qqplot.mc <- function(res, main = "", ylab = "") {
  n <- ncol(res)
  NREP <- nrow(res)
  
  w <- ppoints(n)          
  quantis_norm <- qnorm(w) 
  
  quantis_replicas <- matrix(NA, nrow = NREP, ncol = n)
  
  for(i in 1:NREP) {
    res_i <- na.omit(res[i, ])
    if(length(res_i) > 0) {
      quantis_replicas[i, ] <- quantile(res_i, probs = w)
    }
  }
  
  quantis_media <- colMeans(quantis_replicas, na.rm = TRUE)
  
  plot(quantis_norm, quantis_media,
       pch = 16, cex = 1, cex.lab = 1.3, cex.axis = 1, cex.main = 2,
       xlab = "N(0,1) quantiles",
       ylab = ylab,
       main = main,
       xlim = c(-2.3, 2.3), ylim = c(-2.3, 2.3))
  
  abline(0, 1, col = "black", lwd = 2)
  
}

par(mfrow = c(1,1), mar = c(3.2, 3.2, 2, 3.2), oma = c(0.5, 0.5, 0.5, 0.5), mgp = c(2.2, 0.6, 0))

# Resíduo quntilico
qqplot.mc(quantile[,,1], main = "Scenario I", ylab = "Quantile residual")
qqplot.mc(quantile[,,2], main = "Scenario II", ylab = "Quantile residual")
qqplot.mc(quantile[,,3], main = "Scenario III", ylab = "Quantile residual")
qqplot.mc(quantile[,,4], main = "Scenario IV", ylab = "Quantile residual")
qqplot.mc(quantile[,,5], main = "Scenario V", ylab = "Quantile residual")

# Resíduo Pearson
qqplot.mc(pearson[,,1], main = "", ylab = "Pearson residual")
qqplot.mc(pearson[,,2], main = "", ylab = "Pearson residual")
qqplot.mc(pearson[,,3], main = "", ylab = "Pearson residual")
qqplot.mc(pearson[,,4], main = "", ylab = "Pearson residual")
qqplot.mc(pearson[,,5], main = "", ylab = "Pearson residual")

# Resíduo Pearson Padronizado
qqplot.mc(pearsonp[,,1], main = "", ylab = "Standardized Pearson residual")
qqplot.mc(pearsonp[,,2], main = "", ylab = "Standardized Pearson residual")
qqplot.mc(pearsonp[,,3], main = "", ylab = "Standardized Pearson residual")
qqplot.mc(pearsonp[,,4], main = "", ylab = "Standardized Pearson residual")
qqplot.mc(pearsonp[,,5], main = "", ylab = "Standardized Pearson residual")

# Resíduo Ponderado
qqplot.mc(sweighted1[,,1], main = "", ylab = "Weighted residual")
qqplot.mc(sweighted1[,,2], main = "", ylab = "Weighted residual")
qqplot.mc(sweighted1[,,3], main = "", ylab = "Weighted residual")
qqplot.mc(sweighted1[,,4], main = "", ylab = "Weighted residual")
qqplot.mc(sweighted1[,,5], main = "", ylab = "Weighted residual")

# Resíduo Ponderado Padronizado
qqplot.mc(sweighted2[,,1], main = "", ylab = "Standardized weighted residual")
qqplot.mc(sweighted2[,,2], main = "", ylab = "Standardized weighted residual")
qqplot.mc(sweighted2[,,3], main = "", ylab = "Standardized weighted residual")
qqplot.mc(sweighted2[,,4], main = "", ylab = "Standardized weighted residual")
qqplot.mc(sweighted2[,,5], main = "", ylab = "Standardized weighted residual")

# Resíduo Deviance
qqplot.mc(deviance[,,1], main = "", ylab = "Deviance residual")
qqplot.mc(deviance[,,2], main = "", ylab = "Deviance residual")
qqplot.mc(deviance[,,3], main = "", ylab = "Deviance residual")
qqplot.mc(deviance[,,4], main = "", ylab = "Deviance residual")
qqplot.mc(deviance[,,5], main = "", ylab = "Deviance residual")

# Resíduo Anscombe
qqplot.mc(anscombe[,,1], main = "", ylab = "Anscombe residual")
qqplot.mc(anscombe[,,2], main = "", ylab = "Anscombe residual")
qqplot.mc(anscombe[,,3], main = "", ylab = "Anscombe residual")
qqplot.mc(anscombe[,,4], main = "", ylab = "Anscombe residual")
qqplot.mc(anscombe[,,5], main = "", ylab = "Anscombe residual")

# Resíduo Williams
qqplot.mc(williams[,,1], main = "", ylab = "Williams residual")
qqplot.mc(williams[,,2], main = "", ylab = "Williams residual")
qqplot.mc(williams[,,3], main = "", ylab = "Williams residual")
qqplot.mc(williams[,,4], main = "", ylab = "Williams residual")
qqplot.mc(williams[,,5], main = "", ylab = "Williams residual")
