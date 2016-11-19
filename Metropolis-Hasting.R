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
iter = 10000

## Store trace (re-run FROM HERE)
alpha.mh.trace<- c(alpha)
beta.mh.trace<- beta
tau.mh.trace<- c(tau)

## Step setting
##radius<- 25
radius.beta<- 10
radius.alpha<- 10
radius.tau<- 10

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
  alpha.mh.trace<- append(alpha.mh.trace, alpha)
  
  ## update beta (MH)
  ## threshold
  logy<- logp(beta, tau) + dmvnorm(beta, rep(0, 4), 1/alpha * diag(4), log = T) + log(runif(1))
  
  ## draw proposal
  beta.prop<- beta + rmvnorm(1, rep(0, 4), radius.beta * diag(4))
  logbeta<- logp(beta.prop, tau) + dmvnorm(beta.prop, rep(0, 4), 1/alpha * diag(4), log = T)
  
  ## shrink bracket
  if(logbeta > logy){
    beta<- beta.prop
  }
  
  ## trace
  beta.mh.trace<- rbind(beta.mh.trace, beta)
  
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
  tau.mh.trace<- append(tau.mh.trace, tau)
}