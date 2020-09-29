run_MCMC <- function(n_MCMC, ponto_inicial, proposta, prob_aceita, funcao_p){
  it=1
  nro_algarismos <- length(ponto_inicial)
  MCMC_d <- matrix(ncol = n_MCMC, nrow = nro_algarismos)
  MCMC_d[,1] <- ponto_inicial
  while(it < (n_MCMC)){    ### inicia o algoritmo
      x_old <- MCMC_d[,it]
      x_new <- x_old
      x_new <- proposta(x_old)                ## cria uma proposta
      
      if(runif(1) < prob_aceita(x_old,x_new, funcao_p)){
        MCMC_d[,it+1] <- x_new                         #com probabilidade prob_aceita registra o valor da proposta no MCMV 
      }else{
        MCMC_d[,it+1] <- x_old                         #caso contrario repete o valor antigo 
      }
      it <- it+1
      if (it == n_MCMC){
        break }
  }
  #return(MCMC_d[,-n_MCMC])
  return(MCMC_d)
}

runifdisc<-function(n=2, min=0, max=1) sample(min:max, n, replace=T) ##função para amostrar de uma uniforme discreta

#######Problema Original da aula: P(x) ~ (x[1]+x[2]+x[3])^5 para x = (000, 001...999)
#vamos definir nossas funções para o problema:

funcao_p <- function(x){      ## função proporcional a densidade que queremos amostrar (aqui X é o vetor de dados)
  prop <- exp((-x[1]^2)/2 - (x[2]^2)/2)
  return(prop)
}

proposta <- function(xi){  #neste caso, geramos uma proposta bivariada, onde (p1,p2) ~ (p1_anterior,p2_anterior) + (runifdiscreta[-5,5])
  a <- c(0,0)
  while(var(a) == 0){ #nao queremos o caso que a amostragem resulta em proposta igual a anterior
    a <- runifdisc(2,-5,5)
  }
  p1 <- xi[1] + a[1]
  p2 <- xi[2] + a[2]
  return(c(p1,p2))
}

prob_aceita <- function(x_old,x_new, funcao_p){
  min(funcao_p(x_new)/funcao_p(x_old), 1)
}

#rodando o MH
MCMC_d <- run_MCMC(n_MCMC = 50000, c(0,0), proposta = proposta, prob_aceita = prob_aceita, funcao_p = funcao_p)

#################  Estimar E(Xi) e Var(Xi)
mean(MCMC_d[1,]) 
#[1] 0.0038

mean(MCMC_d[2,])
#[1] -0.01218

var(MCMC_d[1,])
#[1] 0.9984855
var(MCMC_d[2,])
#[1] 1.036452

cov(MCMC_d[1,],MCMC_d[2,])
#[1] -0.002533767

#Ou seja, podemos notar que (X1,X2) ~ N((0,0),I)


##Vamos agora plotar as distribuições marginais;
hist(MCMC_d[1,]) #par1 
hist(MCMC_d[2,]) #par2

##distribuição conjunta bivariada discreta: 
library(tidyverse)

MCMC_d <- as_tibble(as.data.frame(t(MCMC_d)))
glimpse(MCMC_d)

ggplot(MCMC_d, aes(x = V1, y = V2)) + geom_density2d()


##Analise do trace plot:
library(coda)
MCMC_coda <- rbind(MCMC_d)
x <- coda::mcmc(MCMC_coda)  
summary(x, na.rm = TRUE)
plot(x) ##podemos ver que o Trace Plot está satisfatoriamente estacionário para ambas as variáveis. 

accept_rate <- 1 - coda::rejectionRate(x)
accept_rate ##a taxa de aceitação ficou em ~7% 

coda::effectiveSize(x) ##Para RE = 50000, o tamanho amostral efetivo ficou em media ~3000 para as variáveis. 
coda::autocorr.plot(x) 
