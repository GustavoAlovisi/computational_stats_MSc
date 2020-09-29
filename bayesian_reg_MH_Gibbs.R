run_MCMC <- function(n_MCMC, ponto_inicial, proposta_betas, proposta_phi, prob_aceita, funcao_p, sd){
  it=1
  MCMC_d <- matrix(ncol = n_MCMC, nrow = 3)
  MCMC_d[,1] <- ponto_inicial
  while(it < (n_MCMC)){    ### Vamos iniciar nosso MH com uma proposta para Beta0, Beta1 e phi 
    x_old <- MCMC_d[,it]
    x_new <- x_old
    x_new[1] <- proposta_betas(x_old[1], sd) #Beta0              
    x_new[2] <- proposta_betas(x_old[2], sd) #Beta1
    x_new[3] <- proposta_phi(x_old[3], sd) #phi
    
        
    if(runif(1) < prob_aceita(x_old,x_new, funcao_p)){
      MCMC_d[,it+1] <- x_new   #com probabilidade prob_aceita registra o valor da proposta
    }else{
      MCMC_d[,it+1] <- x_old   #caso contrario repete o valor antigo 
    }
    it <- it+1
    if (it == n_MCMC){
      break }
  }
  return(MCMC_d)
}


#vamos definir nossas funções para o problema:

funcao_p <- function(x){  ## função proporcional a densidade que queremos amostrar (encontrada em a)) 
  beta0 <- x[1]
  beta1 <- x[2]
  phi <- x[3]
  prop <- sum(dnorm(Yi, beta0+beta1*Xi, sqrt(1/phi), log = TRUE) + dnorm(beta0, 0, sqrt(1/Tao0), log = TRUE) + 
              dnorm(beta1, 0, sqrt(1/Tao1), log = TRUE) + dgamma(phi, Alpha, Beta))
  return(prop)
}

proposta_betas <- function(par, sd = 1){  ## proposta para os betas: N(i, s²)
  proposta <- rnorm(1, par, sd = sd)
  return(proposta)
}

proposta_phi <- function(par, sd = 1){ ##proposta para phi: |N(i,s²)|, pois a variância deve ser positiva 
  proposta <- abs(rnorm(1, par, sd = sd))
  return(proposta)
}

prob_aceita <- function(x_old,x_new, funcao_p){
  min(exp(funcao_p(x_new))/exp(funcao_p(x_old)), 1)
}


##########simulando os dados
Tao0 <- 1/200
Tao1 <- 1/200
Alpha <- 1
Beta <- 0.2

set.seed(13)
Xi=runif(50,1,9)
Yi=5-0.5*Xi +rnorm(50)


###########b) Rodando Metropolis-Hastings para Reg. Bayesiana:

MCMC_d <- run_MCMC(50000, c(0.5,0.5,0.5), proposta_betas = proposta_betas, 
                   proposta_phi = proposta_phi, prob_aceita = prob_aceita, funcao_p = funcao_p, sd = 0.5)

mean(MCMC_d[1,])
mean(MCMC_d[2,])
mean(MCMC_d[3,])


library(coda)
MCMC_d <- rbind(MCMC_d[,200:50000]) ##burn-in 

x <- coda::mcmc(t(MCMC_d))   #cria elemento coda.mcmc
summary(x, na.rm = TRUE)
plot(x) #Apesar de parecer estacionário, podemos notar certa autocorrelação no Trace Plot dos parâmetros. 

accept_rate <- 1 - coda::rejectionRate(x)
accept_rate #A taxa de aceitação foi baixa para esta proposta. 

coda::effectiveSize(x) ## 50000 replicações resultaram em um tamanho efetivo de +- 500 observações para cada parâmetro
coda::autocorr.plot(x) 



#########c) vamos comparar a taxa de aceitação para s² entre 0.01 e 1:

s2 <- seq(0.01, by = 0.05)

tam_ef <- matrix(NA, 20, 3)
for(i in 1:20){
  MCMC_d <- run_MCMC(50000, c(0.5,0.5,0.5), proposta_betas = proposta_betas, 
                     proposta_phi = proposta_phi, prob_aceita = prob_aceita, funcao_p = funcao_p, sd = s2[i])
  x <- coda::mcmc(t(MCMC_d))
  tam_ef[i,] <- coda::effectiveSize(x)
}

#vamos considerar a mediana dos tamanhos efetivos:
apply(tam_ef, 1, median)

#Analisando a lista, a mediana mais alta para s² foi encontrada em tam_ef[9,]. 
s2[9]
# [1] 0.41
#equivalente a um s² de 0.41. 

#Podemos assim concluir que valores muito baixos ou muito altos para s² diminuem a taxa de aceitação do MH.  


########d) a estimativa para a média a posteriori dos paramêtros resultou em 
mean(MCMC_d[1,])
mean(MCMC_d[2,])
mean(MCMC_d[3,])

#Beta0 [1] 4.622311
#Beta1 [1] -0.464727
#Phi [1] 1.136956  e 1/Phi = 0.884



#########e) vamos implementar o Gibbs-Sampler

N <- length(Xi)
X2 <- sum(Xi^2)

###distribuições condicionais:
b0_prop <- function(par_old_vec){
  s <- 1/(Tao0 + par_old_vec[3]*N)
  m <- par_old_vec[3]*sum(Yi-par_old_vec[2]*Xi)*s
  par_new_vec <- par_old_vec
  par_new_vec[1] <- rnorm(1, m, sqrt(s))
  return(par_new_vec[1])
}

b1_prop <- function(par_old_vec){
  s <- 1/(Tao1 + par_old_vec[3]*X2)
  m <- par_old_vec[3]*sum((Yi-par_old_vec[1])*Xi)*s
  par_new_vec <- par_old_vec
  par_new_vec[2] <- rnorm(1, m, sqrt(s))
  return(par_new_vec[2])
}

sig_prop <- function(par_old_vec){
  beta <- Beta + sum((Yi - par_old_vec[1] - par_old_vec[2]*Xi)^2)/2
  par_new_vec <- par_old_vec
  par_new_vec[3] <- rgamma(1, shape = Alpha + N/2, rate = beta)
  return(par_new_vec[3])
}


##gibbs-sampler: 
#Para o parâmetro beta0_t, utilizamos beta1_t-1 e phi_t-1 
#Para o parâmetro beta1_t, utilizamos beta0_t e phi_t-1 
#Para o parâmetro phi_t, utilizamos beta0_t e beta1_t 

gibbs_sampler = function(par0, n_MCMC){
  
  MCMC_par = matrix(ncol = n_MCMC, nrow = 3)     
  MCMC_par[,1] = par0
  
  #parametro pra ser atualizado pela proposta
  for(it in 1:(n_MCMC-1)){    ### inicia o algoritmo
    if(it == 1){
      MCMC_par[1,1] <- b0_prop(MCMC_par[,1])

      MCMC_par[2,1] <- b1_prop(MCMC_par[,1])

      MCMC_par[3,1] <- sig_prop(MCMC_par[,1])

    } else{
    MCMC_par[1,it] <- b0_prop(MCMC_par[,it-1])
    MCMC_par[2,it] <- b1_prop(c(MCMC_par[1,it],MCMC_par[2,it-1],MCMC_par[3,it-1]))
    MCMC_par[3,it] <- sig_prop(MCMC_par[,it])}
  }
  return(MCMC_par)
}

n_MCMC = 50000
par0 = c(0.5,0.5,0.5)
set.seed(13)
MCMC_gibbs<- gibbs_sampler(par0, n_MCMC)

mean(MCMC_gibbs[1,], na.rm = T)
#[1] 4.762613

mean(MCMC_gibbs[2,], na.rm = T)
#[1] -0.4875055

mean(MCMC_gibbs[3,], na.rm = T)
#[1] 1.22853 


###########f) 

MCMC_d <- rbind(MCMC_gibbs[,200:49999]) ##burn-in 
x <- coda::mcmc(t(MCMC_d))   #cria elemento coda.mcmc
summary(x, na.rm = TRUE)
plot(x) #Notamos que o Trace Plot do Gibbs-Sampler parece totalmente estacionário. 

accept_rate <- 1 - coda::rejectionRate(x)
accept_rate #A taxa de aceitação ficou em 1, resultado esperado teoricamente. 

coda::effectiveSize(x) ## 50000 replicações resultaram em um tamanho efetivo significativamente maior do que o Metropolis. 
coda::autocorr.plot(x) ## Notamos também a ausência de uma autocorrelação persistente. 

#Assim, podemos concluir que o algoritmo de Gibbs-Sampling nos fornece estimativas melhores que MH neste caso.
#Para o Metropolis, devemos nos preocupar com a proposta e eventuais parâmetros, como foi o caso do s². 
#Ainda, notamos a baixa taxa de aceitação e tamanho amostral efetivo.
#Vale notar porém que precisamos obter as distribuições condicionais dos parâmetros para utilizar o método de Gibbs, o que nem sempre
#pode ser trivial. 



