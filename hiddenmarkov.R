#Yt|Xt = i ~ N(0, sigma²i) ou seja Retorno em t dado estado da cadeia em t segue N(0, sigma_i) 
#onde sigma_i é dado pelo estado da cadeia

#Assim, precisamos simular uma normal com runif onde sigma é dado pelo estado da cadeia no tempo t
#p_i = 1/4 para todas as cadeias, ou seja, elas tem a mesma prob de acontecer

##boxmuller para gerar normal(mu, sigma) a partir de runif
lsboxmuller <- function(mu, sigma, N=2){ #box-muller para normal com media mu e desvio sigma
  u1 = runif(N/2) 
  u2 = runif(N/2)
  norm = sqrt(-2*log(u1))*cos(2*pi*u2)
  return(mu + sigma*norm)
}
#ks.test(a, 'pnorm', 0, 3)


##a) gerando os 100 retornos da cadeia oculta

sim_ret <- function(N){
  yret <- matrix(NA, N ,2)
  for (t in 1:N){
    u <- runif(1)
    if(u <= 0.25){  #se unif <0.25, cadeia 1 e geramos retorno y_i ~ N(0, 1)
      yret[t, 1] <- lsboxmuller(mu=0, sigma=1) #rnorm(1, 0, 1)
      yret[t, 2] <- 1 #retorno gerado sob a hidden markov chain 1
    }
    else if(u <= 0.5){  #se unif entre 0.25 e 0.5, geramos retorno y_i ~ N(0, 3)
      yret[t, 1] <- lsboxmuller(mu=0, sigma=3) #rnorm(1, 0, 3)
      yret[t, 2] <- 2 #retorno gerado sob a hidden markov chain 2
    }
    else if(u <= 0.75){ #...
      yret[t, 1] <- lsboxmuller(mu=0, sigma=6) #rnorm(1, 0, 6)
      yret[t, 2] <- 3 #...
    }
    else{
      yret[t, 1] <- lsboxmuller(mu=0, sigma=9) #rnorm(1, 0, 9)
      yret[t, 2] <- 4
    }
  }
  yret <- as.data.frame(yret)
  colnames(yret) <- c('Retornos', 'EstadoReal')
  return(yret)
}

yret <- sim_ret(100)

library(ggplot2)
ggplot(yret, aes(x = EstadoReal, y = Retornos)) + geom_point() #investigando o que foi gerado
plot(yret$Retornos,type = 'l') #grafico dos retornos que seguem uma cadeia de markov oculta



##b) vamos estimar com o código da aula a cadeia em que y_t se encontra, bem como as probabilidades dele trocar de cadeia em t+1
library(depmixS4)
mod <- depmix(Retornos ~ 1, family = gaussian(), nstates = 4, data = yret)
set.seed(1)
fm2 <- fit(mod, verbose = FALSE)
#
summary(fm2)
#Com 100 retornos simulados, a estimativa das normais não é muito boa. Porém, ao aumentarmos o tamanho da simulação (numero de retornos da serie),
#nossa estimativa das normais dos estados melhora consideravelmente.

print(fm2) #info sobre convergencia/critérios de info

probs <- round(posterior(fm2),3)    #probabilidade de y_t estar em cada estado no tempo t 
head(probs)

yret <- cbind(yret, probs)
head(yret) #analisando as trocas de estado dado o estado real simulado ao longo do tempo 


##c) Vamos aproximar a esperança dos estimadores de EM das médias e desvios das normais utilizando simulação de monte carlo: 
#Para cada iteração vamos samplear de nossas distribuições normais propostas e ajustar o modelo de hmm e depois tomar a media
 
RE <- 100  #100 simulações
par_aux <- matrix(NA, 4,2)
sigma <- matrix(NA, RE, 4)
mu <- matrix(NA, RE, 4)

for(i in 1:RE){
  yret <- sim_ret(100) #simulando retornos
  mod <- depmix(Retornos ~ 1, family = gaussian(), nstates = 4, data = yret) #ajustando o modelo para os retornos
  fm2 <- fit(mod, verbose = F)
  for(j in 1:4){
    par_aux[j,]<-unlist(fm2@response[[j]][[1]]@parameters) #obtendo os parametros de mu e sigma
  }
  sigma[i,] <- par_aux[order(par_aux[,2]),2] #ordenando por ordem de desvio padrao
  mu[i,] <- par_aux[order(par_aux[,2]),1]
}

apply(sigma, 2, mean) 
#1] 1.047266 2.155644 3.904259 7.693259

apply(mu, 2, mean)
#1] -0.4653428  0.9744585  0.6638414 -1.1275777

#percebemos que com uma simulação de monte carlo os valores estimados melhoram bastante 
#a estimação também melhora conforme aumentamos o numero de retornos simulados (padrão = 100) e o numero de repetições da sim. de monte carlo (RE).

