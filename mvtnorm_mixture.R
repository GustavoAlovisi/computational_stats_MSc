######mvt normal copula mixtures


library(mvnfast)


mu <- matrix(c(0,0,0,
               3,3,3), ncol = 3, nrow =2, byrow =T)

covmat <- diag(3)

weights <- c(0.7, 0.3)

normalmvtmix <- rmixn(n=40, mu = mu, sigma = list(covmat, covmat), w = weights, retInd = T)

normalmvtmix1 <- as.data.frame(normalmvtmix)

ggplot(normalmvtmix1, aes(x=V1)) + geom_density()
#########

library(mvtnorm)
p1 <- 0.7
p2 <- 0.3

mvtmixLL_optim <- function(m,x){
  m1 = m[1:3]
  m2 = m[4:6]
  #p1 = m[7] ##estimating weights too
  #p2 = m[8] ##w2 
  -sum(log(p1*mvtnorm::dmvnorm(x = x, mean = m1)+ p2*mvtnorm::dmvnorm(x = x, mean = m2)))
}

##################maximizing log-likelihoods using optim() and nelder mead
#we use nelder-mead because it is a simple unconstrained optmization. If we had equality constraints on weights, 
#we could use L-BFGS-B with two weights (p1, 1-p1) and bound p1 from 0 to 1
optim(par = c(-1,-1,-1,4,4,4), fn = mvtmixLL_optim, x = normalmvtmix, method = 'Nelder-Mead')
#sensivel ao ponto inicial! 

#################using Simulated Annealing

mvtmixLL <- function(m,x){
  m1 = m[1:3]
  m2 = m[4:6]
  #p1 = m[7] ##estimating weights too
  #p2 = m[8] ##w2 
  sum(log(p1*mvtnorm::dmvnorm(x = x, mean = m1)+ p2*mvtnorm::dmvnorm(x = x, mean = m2)))
}


n.iter <- 500
teta <- matrix(0,nrow=n.iter,ncol=6) ##6 parametros, 3 de cada mvt normal
teta[1,] <- c(0,-1,-1, 4,-4,4) #initial guesses para a primeira iteração
f.eval <- rep(0,n.iter)

f.eval[1] <- mvtmixLL(teta[1,],x = normalmvtmix) #x = normalmvtmix


for(i in 2:n.iter){
  temp <- 1/log(i+1)
  #temp <- 1/log(1+i) 
  d <- 0.1*sqrt(temp)
  e <-rnorm(6)*d #runif(6, -d, d)
  te1 <- teta[(i-1),] + e
  u <- runif(1)
  prob <- min(exp((mvtmixLL(m = te1, x = normalmvtmix)-mvtmixLL(m = teta[(i-1),], x = normalmvtmix))/temp), 1)
  print(prob)
  teta[i,] <- (u<=prob)*te1+(u>prob)*teta[(i-1),]
  f.eval[i] <- mvtmixLL(m = teta[i,], x = normalmvtmix)
}

#rnorm(6)*0.6
pos <- which.max(f.eval)
teta.ot <- teta[pos,]
teta.ot

#alternativamente, poderiamos ter feito
optim(par = c(-1,-1,-1,4,4,4), fn = mvtmixLL_optim, x = normalmvtmix, method = 'SANN', control = list(temp = 15, tmax = 30))
#note que os resultados são bem parecidos: nosso d parece ser o adequado. 




