library(MASS)
library(fBasics)
library(moments)
library(pracma)

######### funcao log-verossimilhanca #################

log.like = function(y,mu,sigma)
{
  a <- mu*(1+sigma)
  b <- 2 + sigma 
  lmd.i = (a-1)*log(y) - (a+b)*log(1+y) - lbeta(a,b)
  
  log.theta = sum(lmd.i)
  return(log.theta)
}


#####################################
###  Matriz H - Alavancagem      ####
#####################################
###  Para função de ligação log   ###

influence.BP <- function(model){
  n <- model$N
  y <- model$y
  X <- model$mu.x 
  beta <- model$mu.coefficients
  mu <- model$mu.fv  
  phi <- model$sigma.fv 
  linkmu <- model$mu.link 
  
  onephi <- (1+phi) 
  muonephi <- mu*onephi
  muonephitwo <- muonephi+phi+2
  yast <- log(y/(1+y))
  muast <- digamma(muonephi) - digamma(muonephitwo)
  
  #derivative of the linkmu
  if(linkmu=="log"){
    explink<-expression(log(mu))
    DlinkDmu<-D(explink,"mu")
  } else if(linkmu=="sqrt"){
    explink<-expression(sqrt(mu))
    DlinkDmu<-D(explink,"mu")
  }else{warning("Função link diferente")
  }
  
  Dlink <- (1/eval(DlinkDmu)) 
  DLINK <- diag(Dlink,ncol=n,nrow = n)
  Vi <- trigamma(muonephi) - trigamma(muonephitwo)
  Wi <- onephi*(Vi)*Dlink^2 # esse calculo resulta na diagonal da matrix W
  WI <- diag(Wi,ncol=n,nrow = n)

  PHI <- diag(onephi,ncol=n,nrow = n)
  Inv <- solve(t(X)%*%PHI%*%WI%*%X)
  
  z <- X%*%beta + solve(WI)%*%DLINK%*%(yast-muast)
  y1 <- (PHI^(1/2))%*%(WI^(1/2))%*%z
  
  Hat <- (PHI^(1/2))%*%(WI^(1/2))%*%X%*%Inv%*%t(X)%*%(PHI^(1/2))%*%(WI^(1/2))
  hat <- diag(Hat)
  Result <- list(hat=hat,H=Hat, ychap = y1, zi = z)
  Result
}


residuals.BP <- function(object,
                         type = c("quantile", "pearson", "pearson P", "weighted", "sweighted1", "sweighted2", 
                                  "combined", "variance", "deviance", "deviance P", "response", "anscombe", "williams"), ...)
{
  y <- object$y
  mu <- object$mu.fv
  phi <- object$sigma.fv
  res <- y-mu
  hat <-influence.BP(object)$hat
  quantile <- object$residuals 
  
  onemu <- 1+mu
  onephi <- 1+phi
  muonemu <- mu*onemu
  muonephi <- mu*onephi
  muonephitwo <- muonephi+phi+2
  yast <- log(y/(1+y))
  muast <- digamma(muonephi) - digamma(muonephitwo)
  v <- trigamma(muonephi) - trigamma(muonephitwo)
  
  ystar <- mu*log(y) - onemu*log(1+y)
  mustar <- mu*muast - digamma(muonephitwo) + digamma(phi+2)
  e <- mu^2 * trigamma(muonephi) + trigamma(phi+2) - (onemu^2) * trigamma(muonephitwo)
  z <- (onemu^2) * trigamma(muonephi) - ((2+mu)^2) * trigamma(muonephitwo) + trigamma(phi+2)
  
  type <- match.arg(type)
  
  if(type == "response") return(res)
  
  wts <- weights(object)
  if(is.null(wts)) wts <- 1
  
  res <- switch(type,
                
                "quantile" = {
                  sqrt(wts) * quantile #qnorm(pBP(y,mu,phi))
                },
                
                "pearson" = {
                  sqrt(wts) * (sqrt(phi)*res) / sqrt(muonemu)
                },
                
                "pearson P" = {
                  sqrt(wts) * (sqrt(phi)*res) / (sqrt(muonemu*(1-hat)))
                },
                
                "weighted" = {
                  sqrt(wts) * (yast-muast) / sqrt(phi*v)
                },
                
                "sweighted1" = {
                  sqrt(wts) * (yast-muast) / sqrt(v)
                },
                
                "variance" = {
                  sqrt(wts) * (ystar-mustar) / sqrt(e)
                },
                
                "sweighted2" = {
                  sqrt(wts) * (yast-muast) / sqrt(v*(1-hat))
                },
                
                "combined" = {
                  sqrt(wts) * ((yast-muast) + (ystar-mustar)) / sqrt(z)
                },
                
                "deviance" = {
                  ll <- function(mu, phi)
                    (muonephi-1)*log(y) - muonephitwo*log(1+y) - lbeta((mu*(1+phi)),(2+phi))
                  sqrt(wts) * sign(res) * sqrt(2 * abs(ll(y, phi) - ll(mu, phi)))
                },
                
                "deviance P" = {
                  ll <- function(mu, phi)
                    (muonephi-1)*log(y) - muonephitwo*log(1+y) - lbeta((mu*(1+phi)),(2+phi))
                  sqrt(wts) * (sign(res) * sqrt(2 * abs(ll(y, phi) - ll(mu, phi)))) / sqrt(1-hat)
                },
                
                "anscombe" = {
                  func <- function(x) (1/((x*(1+x))^(1/3)))
                  integralVy <- NULL
                  integralVmu <- NULL
                  for (i in 1:length(mu)) {
                    integralVy[i] <- integral(func,0,y[i])
                    integralVmu[i] <- integral(func,0,mu[i])
                  }
                  Vmu <- mu*(1+mu)
                  sqrt(wts) *((integralVy - integralVmu)*sqrt(phi))/((Vmu)^(1/6))
                  
                },
                
                "williams" = {
                  resP <- residuals.BP(object,type="pearson P")
                  resD <- residuals.BP(object,type="deviance P")
                  sqrt(wts) * sign(res) * sqrt((1-hat)*(resD^2) + hat*(resP^2))
                })
  
  return(res)
}

envelope.BP <- function(model, k=100, link=c("log","log"), type=c("quantile", "pearson", "pearson P", "sweighted1", "sweighted2", 
                                                                  "variance", "combined", "deviance", "anscombe", "williams"),  
                        color="grey50", xlabel="N(0,1) quantiles", ylabel="Empirical quantiles", font="Times")
{
  set.seed(2022) 
  n=model$N
  td <- residuals.BP(model, type)
  sigma = model$sigma.fv
  mu = model$mu.fv
  re <- matrix(0,n,k)
  mu.link <- model$mu.link
  sigma.link <- model$sigma.link
  
  if(length(model$mu.coefficients)==1) x <- model$mu.x else x <- model$mu.x[,-1];
  if(length(model$sigma.coefficients)==1) z <- model$sigma.x else z <- model$sigma.x[,-1];
  
  for(i in 1:k)
  { 
    y1 <- rBP(n,mu,sigma)
    
    if(length(model$mu.coefficients)==1){
      form.x <- y1 ~ x-1 
    }else{
      form.x <- y1 ~ x
    }  
    if(length(model$sigma.coefficients)==1){
      form.z <- y1 ~ z-1 
    }else{
      form.z <- y1 ~ z
    }
    
    if(mu.link=="log" & sigma.link == "identity"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="log",sigma.link="identity"),method=RS(), trace=FALSE)
    } else if(mu.link=="log" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="log",sigma.link="log"),method=RS(), trace=FALSE)
    } else if(mu.link=="log" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="log",sigma.link="sqrt"),method=RS(), trace=FALSE)
    } else if(mu.link=="identity" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="identity",sigma.link="sqrt"),method=RS(), trace=FALSE)
    } else if(mu.link=="identity" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="identity",sigma.link="log"),method=RS(), trace=FALSE)
    } else if(mu.link=="identity" & sigma.link == "identity"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="identity",sigma.link="identity"),method=RS(), trace=FALSE)
    } else if(mu.link=="sqrt" & sigma.link == "sqrt"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="sqrt",sigma.link="identity"),method=RS(), trace=FALSE)
    } else if(mu.link=="sqrt" & sigma.link == "log"){
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="sqrt",sigma.link="log"),method=RS(), trace=FALSE)
    } else{
      model1 <- gamlss(formula=form.x ,sigma.formula=form.z, family=BP(mu.link="sqrt",sigma.link="identity"),method=RS(), trace=FALSE)
    }
    
    rdf <- residuals.BP(model1, type) 
    re[,i] <- sort(rdf)
  }
  
  e1 <- numeric(n)
  e2 <- numeric(n)
  
  for(i in 1:n){
    eo <- sort(re[i,])
    e1[i] <- (eo[2]+eo[3])/2 
    e2[i] <- (eo[97]+eo[98])/2
  }
  
  med <- apply(re,1,mean)
  faixa <- range(td,e1,e2)
  
  q1 <- quantile(med, 0.025)
  q2 <- quantile(med, 0.975)
  
  # points outside the envelopes
  c95 <- 0
  posi <- NULL
  
  r_ord <- sort(td) 
  r_ord_p <- NULL
  
  for(q in 1:n){
    for(w in 1:n){
      if(r_ord[q] == td[w]){ 
        r_ord_p[q] <- w  
        if(td[r_ord_p[q]] < e1[q] || td[r_ord_p[q]] > e2[q]){ 
          c95 <- c95 + 1
          posi[q] <- r_ord_p[q] 
        }
      }
    }
  }
  
  prop95 <- round(c95/n*100, 2) # points outside the envelopes (in %)
  
  ###### plots #####
  par(mar=c(4.8, 4.8, 1, 4.8)) 
  par(mgp=c(3, 1, 0))
  
  r <- qqnorm(td,xlab=xlabel, ylab=ylabel, ylim=faixa, pch=1, main="", cex = 1, cex.axis = 1.4, cex.lab = 1.8)$x
  par(new=T)
  qqnorm(e1,axes=F,xlab="",ylab="",type="l",ylim=faixa,lty=1,main="")
  par(new=T)
  qqnorm(e2,axes=F,xlab="",ylab="", type="l",ylim=faixa,lty=1,main="")
  par(new=T)
  qqnorm(med,axes=F,xlab="", ylab="", type="l",ylim=faixa,lty=2,main="")
  legend("topleft", legend=c(paste("Total of points:", n), paste("Points outside the envelope:", c95, "(", prop95, "%)")), bty="n", cex = 1)
  #topleft bottomright
  
  # name the points outside the envelopes (plot the index on the graph)
  #if(length(posi)!=0) text(r[posi], td[posi], posi, pos=1, cex = 0.75) 
  assign("q1",q1[[1]], envir = globalenv())
  assign("q2",q2[[1]], envir = globalenv())
}


plot.BP <- function(model, which = 1:5, q1, q2, pos1, pos2, xlabcov,
                    caption = c("Residuals vs. Index", "Residuals vs. Linear predictor", "Residuals vs. Adjusted values",
                                "Adjusted values vs Observed values"), main = "", 
                    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ylab,
                    ..., type = c("quantile", "pearson", "pearson P", "sweighted1", "sweighted2", "variance", "combined", "deviance", "anscombe", "williams"), 
                    nsim = 100, level = 0.95)
{
  if(!is.numeric(which) || any(which < 1) || any(which > 5)) 
    stop("`which' must be in 1:5")
  
  res <- residuals.BP(model, type = type)
  n <- length(res)
  show <- rep(FALSE, 5)
  show[which] <- TRUE
  one.fig <- prod(par("mfcol")) == 1
  
  if(type=="quantile") ylab = "Quantile"
  if(type=="pearson") ylab = "Pearson"
  if(type=="pearson P") ylab = "Standardized Pearson"
  if(type=="sweighted1") ylab = "Weighted"
  if(type=="sweighted2") ylab = "Standardized weighted"
  if(type=="variance") ylab = "Variance"
  if(type=="combined") ylab = "Combined"
  if(type=="deviance") ylab = "Deviance"
  if(type=="anscombe") ylab = "Anscombe"
  if(type=="williams") ylab = "Williams"
  
  par(mar=c(4.8, 4.8, 1, 4.8)) 
  par(mgp=c(3, 1, 0))
  
  if(ask) {
    op <- par(ask = TRUE)
    on.exit(par(op))
  }
  
  if(show[1]) {
    plot(1:n, res, xlab="Index", ylab=ylab, pch=1, cex=1, cex.lab=1.8, cex.axis=1.3, ylim=c(min(res)-0.5,max(res)+0.5),...)
    #abline(h = q1, lty=2); abline(h = q2, lty=2)
    #I <- q1; S <- q2
    #idI <- which(res < I)
    #idS <- which(res > S)
    #if(length(idS)!= 0) text(idS,res[idS],seq(1:n)[idS],pos=pos1,cex = 1.2)
    #if(length(idI)!= 0) text(idI,res[idI],seq(1:n)[idI],pos=pos2,cex = 1.2)
    if(one.fig) 
      abline(h = 0, lty = 2, col = "gray30")
    
  }
  
  if(show[2]) {
    plot(predict(model, type ="link"), res, xlab="Linear predictor", ylab=ylab, pch=1, cex=1, cex.lab=1.8, cex.axis=1.3, ylim=c(min(res)-0.5,max(res)+0.5),...)
    #abline(h = q1, lty=2); abline(h = q2, lty=2)
    #I <- q1; S <- q2
    #idI <- which(res < I)
    #idS <- which(res > S)
    #if(length(idS)!= 0) text(predict(model)[idS],res[idS],idS,pos=pos1,cex = 1.2)
    #if(length(idI)!= 0) text(predict(model)[idI],res[idI],idI,pos=pos2,cex = 1.2)
    if(one.fig)
      abline(h = 1, lty = 2, col = "gray30")
  }
  
  if(show[3]) {
    plot(fitted(model), res, xlab="Adjusted values", ylab=ylab, pch=1, cex=1, cex.lab=1.8, cex.axis=1.3, ylim=c(min(res)-0.5,max(res)+0.5),...)
    #abline(h = q1, lty=2); abline(h = q2, lty=2)
    #I <- q1; S <- q2
    #idI <- which(res < I)
    #idS <- which(res > S)
    #if(length(idS) != 0) text(fitted(model)[idS],res[idS],idS,pos=pos1,cex = 1.2)
    #if(length(idI)!= 0) text(fitted(model)[idI],res[idI],idI,pos=pos2,cex = 1.2)
    if(one.fig) 
      abline(h = 0, lty = 2, col = "gray30")
  }
  
  if(show[4]) {
    y <- if(is.null(model$y)) model.response(model.frame(model)) else model$y
    plot(y, fitted(model), xlab = "Observed values", ylab = "Adjusted values", pch = "+", ...)
    if(one.fig) 
      abline(0, 1, lty = 2, col = "gray30")
  }
  
  if(show[5]) {
    plot(x1, res, xlab="x", ylab=ylab, pch=1, cex=1, cex.lab=1.8, cex.axis=1.3, ylim=c(min(res)-0.5,max(res)+0.5),...)
    #abline(h = q1, lty=2); abline(h = q2, lty=2)
    #I <- q1; S <- q2
    #idI <- which(res < I)
    #idS <- which(res > S)
    #if(length(idS) != 0) text(x1[idS],res[idS],idS,pos=pos1,cex = 1.2)
    #if(length(idI) != 0) text(x1[idS],res[idI],idI,pos=pos2,cex = 1.2)
    if(one.fig) 
      abline(h = 0, lty = 2, col = "gray30")
  }
  
  if(!one.fig && par("oma")[3] >= 1) mtext(sub.caption, outer = TRUE, cex = 1.25)
  invisible()
}
