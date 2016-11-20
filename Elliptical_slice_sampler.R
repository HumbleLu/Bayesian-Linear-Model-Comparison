library(car)
library(mvtnorm)
detach("package:datasets", unload=TRUE)
attach(Prestige)
source(file = "Linear_model.R")

## Initilize
alpha<- 1
beta<- c(0, 0, 0, 0)
tau<- 1

## Sampling posterior distribution
iter = 20000

## Store trace
alpha.ess.trace<- c(alpha)
beta.ess.trace<- beta
tau.ess.trace<- c(tau)

## Step setting
##radius<- 25
radius.alpha<- 1
radius.tau<- 1

for (i in 1:iter){
  ## update alpha (MH)
  ## draw proposal
  alpha.prop<- alpha + rnorm(1, 0, radius.alpha)
  while(alpha.prop <= 0){
    alpha.prop<- alpha + rnorm(1, 0, radius.alpha)
  }
  
  ## acceptance
  ratio.alpha<- dmvnorm(beta, rep(0, 4), 1/alpha.prop * diag(4), log = T) -
    dmvnorm(beta, rep(0, 4), 1/alpha * diag(4), log = T) +
    dgamma(alpha.prop, 10, 1/10) - dgamma(alpha, 10, 1/10)
  
  if(ratio.alpha > log(runif(1))){
    alpha<- alpha.prop
  }
  
  ## trace
  alpha.ess.trace<- append(alpha.ess.trace, alpha)
  
  ## update beta (Elliptical slice sampler)
  ## chooseellipse
  v<- rmvnorm(1, rep(0, 4), 1/alpha * diag(4))
  
  ## threshold
  logy<- logp(beta, tau) + dmvnorm(beta, rep(0, 4), 1/alpha * diag(4), log = T) + log(runif(1))
  
  ## draw proposal
  theta<- runif(1, 0, 2*pi)
  theta.min<- theta - 2 * pi
  theta.max<- theta
  beta.prop<- beta * cos(theta) + v * sin(theta)
  logbeta<- logp(beta.prop, tau) + dmvnorm(beta.prop, rep(0, 4), 1/alpha * diag(4), log = T)
  
  ## shrink bracket
  while(logbeta <= logy){
    if(theta < 0){
      theta.min<- theta
    }else{
      theta.max<- theta
    }
    theta<- runif(1, theta.min, theta.max)
    
    beta.prop<- beta * cos(theta) + v * sin(theta)
    logbeta<- logp(beta.prop, tau) + dmvnorm(beta.prop, rep(0, 4), 1/alpha * diag(4), log = T)
  }
  beta<- beta.prop
  
  ## trace
  beta.ess.trace<- rbind(beta.ess.trace, beta)
  
  ## update tau (MH)
  ## draw proposal
  tau.prop<- tau + rnorm(1, 0, radius.tau)
  while(tau.prop <= 0){
    tau.prop<- tau + rnorm(1, 0, radius.tau)
  }
  
  ## acceptance
  ratio.tau<- logp(beta, tau.prop) - logp(beta, tau) +
    dgamma(tau.prop, 50, 1/10) - dgamma(tau, 50, 1/10)
  
  if(ratio.tau > log(runif(1))){
    tau<- tau.prop
  }
  
  ## trace
  tau.ess.trace<- append(tau.ess.trace, tau)
}

