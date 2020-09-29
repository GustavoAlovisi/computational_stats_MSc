

####Função para rodar o MCMC dado nro de iterações, ponto inicial, proposta, prob de aceitar e função prop. a densidade
run_MCMC <- function(n_MCMC, ponto_inicial, proposta, prob_aceita, funcao_p){
  it=1
  nro_algarismos <- length(ponto_inicial)
  MCMC_d <- matrix(ncol = n_MCMC, nrow = nro_algarismos)
  MCMC_d[,1] <- ponto_inicial
  while(it < (n_MCMC)){    ### inicia o algoritmo
    for (i in 1:nro_algarismos){         ### para cada um dos dígitos do número
      x_old <- MCMC_d[,it]
      x_new <- x_old
      x_new[i] <- proposta(x_old[i])                ## cria uma proposta
      
      if(runif(1) < prob_aceita(x_old,x_new, funcao_p)){
        MCMC_d[,it+1] <- x_new                         #com probabilidade prob_aceita registra o valor da proposta no MCMV 
      }else{
        MCMC_d[,it+1] <- x_old                         #caso contrario repete o valor antigo 
      }
      it <- it+1
      if (it == n_MCMC){
        break }
    }
  }
  #return(MCMC_d[,-n_MCMC])
  return(MCMC_d)
}



#######Problema Original da aula: P(x) ~ (x[1]+x[2]+x[3])^5 para x = (000, 001...999)
#vamos definir nossas funções para o problema:

funcao_p <- function(x){      ## função proporcional a densidade que queremos amostrar (aqui X é o vetor de dados)
  (x[1]+x[2]+x[3])^(5)
}

proposta <- function(xi){    # aqui xi eh o elemento que queremos amostrar 
  candidatos=c(0:9)
  candidatos=candidatos[-which(candidatos==xi)]#lista posíveis mudanças
  sample(candidatos,1)         # escolhe um elemento da lista de candidatos com dist uniforme
}

#   p_proposta=1/9   - todas as peopostas tem a mesma prob. 

prob_aceita <- function(x_old,x_new, funcao_p){
  min(funcao_p(x_new)/funcao_p(x_old), 1)
}



p=vector()                 #declara vetor que vai guardar resultado da função densidade para cada observação

MCMC_d <- run_MCMC(n_MCMC = 1700, c(0,0,0), proposta = proposta, prob_aceita = prob_aceita, funcao_p = funcao_p)

#################  Estimar E(X)
X <- MCMC_d[1,]*100 + MCMC_d[2,]*10 + MCMC_d[3,]
hist(X,nclass=100)                         ## Olhar a distribuição de X
hist(X,nclass=500)
mean(X,na.rm=T)
#[1] 709.9659

sd(X,na.rm=T)
#[1] 234.4805

library(coda)
MCMC_coda <- rbind(MCMC_d, MCMC_d[1,]*100+MCMC_d[2,]*10+MCMC_d[3,] ,(MCMC_d[1,]+MCMC_d[2,]+MCMC_d[3,])^5)
row.names(MCMC_coda) <- c("centena","dezena","unidade","valor x","f(x)")

#a <- is.na(MCMC_d)
x <- coda::mcmc(t(MCMC_coda))   #cria elemento coda.mcmc
summary(x, na.rm = TRUE)
plot(x) 

accept_rate <- 1 - coda::rejectionRate(x)
accept_rate

coda::effectiveSize(x) ##cerca de 1700 iterações produziram um valor amostral efetivo de ~200 
coda::autocorr.plot(x) #para a centena, dezena, unidade, x e f(x)



########Problema da lista: 
########P(x) ~ (x[1]+x[2]+x[3]+x[4]+x[5])^5 para x = (00000, 00001, ..., 99999)

#vamos definir a função proporcional a densidade que queremos amostrar
funcao_p <- function(x){      ## função proporcional a densidade que queremos amostrar (aqui X é o vetor de dados)
  (x[1]+x[2]+x[3]+x[4]+x[5])^(5)
}

##rodando o MH 
MCMC_d <- run_MCMC(n_MCMC = 2500, c(0,0,0,0,0), proposta = proposta, prob_aceita = prob_aceita, funcao_p = funcao_p)

X <- MCMC_d[1,]*10000 + MCMC_d[2,]*1000 + 100*MCMC_d[3,] + 10*MCMC_d[4,] + MCMC_d[5,]
hist(X,nclass=100)                         ## Olhar a distribuição de X
hist(X,nclass=500)
mean(X,na.rm=T) #media de X
#[1] 66581.37

sd(X,na.rm=T) #desvio de X
#[1] 25677.71


###Analise do tamanho amostral efetivo:
MCMC_coda <- NULL
MCMC_coda <- rbind(MCMC_d,MCMC_d[1,]*10000 + MCMC_d[2,]*1000 + 100*MCMC_d[3,] + 10*MCMC_d[4,] + MCMC_d[5,],
                   (MCMC_d[1,]+MCMC_d[2,]+MCMC_d[3,]+MCMC_d[4,]+MCMC_d[5,])^5)

row.names(MCMC_coda) <- c("10k","1k","centena","dezena","unidade",'x','f(x)')

#a <- is.na(MCMC_d)
x <- coda::mcmc(t(MCMC_coda))   #cria elemento coda.mcmc
summary(x, na.rm = TRUE)

accept_rate <- 1 - coda::rejectionRate(x)
accept_rate

coda::effectiveSize(x) #entre 2000-3000 iterações são necessárias para gerar um tamanho efetivo de ~200 para 
coda::autocorr.plot(x) #10k, 1k, centena, dezena e unidade.
#neste caso, foram necessárias mais iterações do que no problema original. 
 
#10k       1k  centena   dezena  unidade        x     f(x) 
#236.5231 294.1065 234.6794 235.9561 290.4770 252.4362 397.9397 




