## MCMC with MH samplers
########################################
##                                    ##
##  Prior:                            ##
##    alpha ~ Gamma(1, 10)            ##
##    beta_(0-3) ~ MN(0, 1/alpha * I) ##
##    tau ~ Gamma(1, 10)              ##
##                                    ##
##  Model                             ##
##    Y ~ N(beta*X, 1/tau)            ##
##                                    ##
########################################

logp<- function(beta, tau){
  mu<- beta%*%t(cbind(1, education, women, prestige))
  sum(dnorm(income, mu, sqrt(1/tau), log = TRUE))
}