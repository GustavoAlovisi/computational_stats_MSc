######mvt normal copula mixtures


library(mvnfast)
#mvtnorm::dmvnorm(x = normalmvtmix, mean = c(0,3,2))

mu <- matrix(c(0,0,0,
               3,3,3), ncol = 3, nrow =2, byrow =T)

covmat <- diag(3)

weights <- c(0.7, 0.3)
w1 <- 0.7
w2 <- 0.3
N <- 4000
normalmvtmix <- rmixn(n=4000, mu = mu, sigma = list(covmat, covmat), w = weights, retInd = T)

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
#note que os resultados são bem parecidos: nosso 'd' parece ser o adequado. 



##############usando EM - vamos tambem classificar se cada ponto (vetor) pertence a N1 ou N2 

##função do passo E. Poderiamos generalizar para k misturas.
expect_step <- function(o, y){
  o1 <- o[1:ncol(y)]
  o2 <- o[(ncol(y)+1):(2*ncol(y))]
  r1 <- w1 * dmvnorm(x = y, mean = o1) / (w1 * dmvnorm(x = y, mean = o1) + w2 * dmvnorm(x = y, mean = o2))
  r2 <- w2 * dmvnorm(x = y, mean = o2) / (w1 * dmvnorm(x = y, mean = o1) + w2 * dmvnorm(x = y, mean = o2))
  return(cbind(r1, r2))
}


##função do passo M. Poderiamos generalizar para k misturas. 
max_step <- function(pi, y){
  teta1.n <- sum(pi[,1])
  teta2.n <- sum(pi[,2])
  tetak1 <- tetak2 <- c(rep(0, ncol(y)))
  
  tetak1 <- sapply(1:ncol(y), function(i){
    1/teta1.n * sum(pi[,1]*y[,i])
  })
  
  tetak2 <- sapply(1:ncol(y), function(i){
    1/teta2.n * sum(pi[,2]*y[,i])
  })
  return(c(tetak1, tetak2))
}

##função que realiza a estimação via algoritmo EM. 
GMM_EM <- function(eps, nrep, Y, t0){
  cc <- 1
  rep <- 1
  t <- matrix(ncol = 1000, nrow = 2*ncol(Y))
  t[,1] <- t0
  while(cc > eps & rep < nrep){
    #etapa E
    pi <- expect_step(o = t0, y = Y)
    
    #etapa M
    t1 <- max_step(pi = pi, y = Y)
    
    rep <- rep + 1
    cc <- sum((t0 - t1)^2)
    t0 <- t1
    t[,rep] <- t0
  }
  t <- t[, 1:rep] #salvando os valores de todos os passos realizados para analise do algoritmo excluindo NA
  return(list(t, pi))
}

t0 <- c(-1,-2,-1,6,6,4) ##chute inicial para os parâmetros da mistura. 

EM_result <- GMM_EM(eps = 1e-8, nrep = 100, Y = normalmvtmix, t0 = t0) ##EM

theta_EM <- as.matrix(EM_result[[1]]) ##salvando os valores dos parâmetros (mu's)
#podemos ver que para n=40 os resultados estimados são próximos dos parâmetros reais.
#Com aumento do tamanho da amostra, porém, a estimação fica cada vez mais próxima dos parâmetros reais. 


pi_EM <- EM_result[[2]] ##obtendo os Pi's para cada y_i. Como temos dois elementos da mistura, podemos usar uma regra de decisão 
#em que se pi_k=1 >= 0.5, o vetor de obs. y_i pertence a mistura 1. se pi_k=1 <=0, o vetor pertece a mistura k=2. 

##vamos agora analisar a classificação feita pelo algoritmo EM em relação a classificação real que mvnfast::rmixn nos provê
yi_classificado <- cbind(normalmvtmix, ifelse(pi_EM[,1] >= 0.5, 1, 2)) ##classificando y_i conforme pertencer a mistura 1 ou 2
yi_classificado <- cbind(yi_classificado, attributes(normalmvtmix)$index)
colnames(yi_classificado) <- c("Y1", "Y2", "Y3", "Mist_Est", "Mist_Real")
head(yi_classificado)
#analisando os dados visualmente, podemos ver que o algoritmo faz um bom trabalho em classificar os dados. 

##podemos também calcular a porcentagem de classificações erradas:
class_errada <- abs(yi_classificado[,"Mist_Est"] - yi_classificado[,"Mist_Real"])/N
sum(class_errada)
#o calculo mostra uma porcentagem baixa de classificações erradas. Usando n=4000, cerca de 0,2% dos dados foram classificados de forma errada apenas






####Para visualizarmos melhor a classifcação, vamos considerar o caso de uma mistura de normais bivariadas: 
library(ggplot2)
library(dplyr)

mu <- matrix(c(0,0,
               3,3), ncol = 2, nrow =2, byrow =T)

covmat <- diag(2)

weights <- c(0.7, 0.3)
w1 <- 0.7
w2 <- 0.3
normalbvtmix <- mvnfast::rmixn(n=100, mu = mu, sigma = list(covmat, covmat), w = weights, retInd = T)

t0 <- c(0, -1, 2, 4) ##chute inicial

EM_result <- GMM_EM(eps = 1e-8, nrep = 100, Y = normalbvtmix, t0 = t0) ##EM

theta_EM <- as.matrix(EM_result[[1]]) 
pi_EM <- EM_result[[2]]

yi_classificado <- 0
yi_classificado <- cbind(normalbvtmix, ifelse(pi_EM[,1] >= 0.5, 1, 2)) ##classificando y_i conforme pertencer a mistura 1 ou 2
yi_classificado <- cbind(yi_classificado, attributes(normalbvtmix)$index)
colnames(yi_classificado) <- c("Y1", "Y2", "Mist_Est", "Mist_Real")

##visualizando nossa classificacao: 
yi_classificado <- dplyr::as.tbl(yi_classificado)
ggplot(yi_classificado, aes(x = Y1, y = Y2, color = Mist_Est)) + geom_point()

##por ultimo, vamos visualizar quais pontos foram erroneamente classificados:
yi_classificado <- yi_classificado %>% mutate(Class_Errado = abs(Mist_Est - Mist_Real))

ggplot(yi_classificado, aes(x = Y1, y = Y2, color = as.factor(Class_Errado))) + 
  geom_point() 






