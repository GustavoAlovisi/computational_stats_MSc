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


yret <- matrix(NA, 100 ,2)
for (t in 1:100){
  u <- runif(1)
  if(u <= 0.25){  #se unif <0.25, cadeia 1 e geramos retorno y_i ~ N(0, 1)
    yret[t, 1] <- lsboxmuller(mu=0, sigma=1) #rnorm(1, 0, 1)
    yret[t, 2] <- 1 #retorno gerado sob a hidden markov chain 1
  }
  else if(u <= 0.5){  #se unif entre 0.25 e 0.5, geramos retorno y_i ~ N(0, 3)
    yret[t, 1] <- lsboxmuller(mu=0, sigma=5) #rnorm(1, 0, 3)
    yret[t, 2] <- 2 #retorno gerado sob a hidden markov chain 2
  }
  else if(u <= 0.75){ #...
    yret[t, 1] <- lsboxmuller(mu=0, sigma=10) #rnorm(1, 0, 6)
    yret[t, 2] <- 3 #...
  }
  else{
    yret[t, 1] <- lsboxmuller(mu=0, sigma=20) #rnorm(1, 0, 9)
    yret[t, 2] <- 4
  }
}
yret <- as.data.frame(yret)
colnames(yret) <- c('Retornos', 'EstadoReal')

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
print(fm2)

# Classification (inference task)
probs <- round(posterior(fm2),3)             # Compute probability of being in each state
head(probs)

yret <- cbind(yret, probs)


