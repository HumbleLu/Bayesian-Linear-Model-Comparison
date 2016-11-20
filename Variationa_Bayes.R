## Variational Bayes
iter = 10000

################################
##                            ##
## alpha ~ gamma(a, b)        ##
##                            ##
## beta ~ mvnorm(mu, sigma)   ##
##                            ##
## tau ~ gamma(c, d)          ##
##                            ##
##                            ##
################################

## initilize
a<- 1 + 2
b<- 10
c<- 1 + nrow(Prestige)/2
d<- 10
mu<- c(0, 0, 0, 0)
sigma<- diag(4)

## trace
b.trace<- c(b)
d.trace<- c(d)
mu.trace<- mu
sigma.trace<- list(sigma)

X<- as.matrix(cbind(1, Prestige[, c("education", "women", "prestige")]))
Y<- Prestige[, "income"]

iter = 1000

for(i in 1:iter){
  
  ## update alpha
  b<- as.numeric(b + (mu^2 + diag(sigma))/2)
  b.trace<- append(b.trace, b)
  
  ## update beta
  sigma<- solve(a/b * diag(4) + c * (t(X) %*% X)/d)
  mu<- as.vector(c/d * sigma %*% t(X) %*% Y)
  mu.trace<- rbind(mu.trace, mu)
  sigma.trace[[length(sigma.trace) + 1]]<- sigma
    
  ## update tau
  d<- as.numeric(10 + (t(Y) %*% Y)/2 - (t(mu) %*% t(X) %*% Y) + (t(mu) %*% t(X)%*%X %*% mu)/2)
  d.trace<- append(d.trace, d)
}

