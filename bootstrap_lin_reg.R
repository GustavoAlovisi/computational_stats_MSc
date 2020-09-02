
############linear regression bootstrap
library(tidyverse)


Y <- c(6.9, 5.1, 7.5, 11.1, 10.9, 4.2, 10.5, 6.8, 12.3, 14.3)
X <- c(6.2, 5.1, 7.6, 2.5, 3.5, 9.4, 4.1, 6.3, 3.0, 0.8)
data <- as_tibble(cbind(Y,X))

#a) construa um IC para b0 e b1 com 95% onde Y = b0 + b1*X + e

B <- 10000 #10000 bootstraps 

beta0b <- numeric(10000)
beta1b <- numeric(10000)

for (i in 1:B){
  sample <- dplyr::sample_n(data, 10, replace = T)
  coefs <- coef(lm(sample$Y ~ sample$X))
  beta0b[i] <- coefs[1]
  beta1b[i] <- coefs[2]
}

hist(beta0b)
hist(beta1b)

quantile(beta0b, c(0.025, 0.975))
quantile(beta1b, c(0.025, 0.975))



#b) Assumindo X dado e erro ~ N(0, sigma²) utilize bootstrap parametrico p/ obter as mesmas estimativas

#Queremos estimar sigma do erro ~ N(0, sigma²) 

model <- lm(Y ~ X)
coefs <- model$coefficients #betas iniciais estimados
sigma <- sigma(model) #estimando sigma dos residuos (erro)

coefs_boot <- c(0,0)

for (i in 1:B){
  y_n <- coefs[1] + coefs[2]*X + rnorm(10, mean = 0, sd = sigma) #geramos Y_boot com X de nossos dados, betas estimados e N(0, sigma^)
  #sample <- dplyr::sample_n(data, 10, replace = T)
  coefs_boot <- coef(lm(y_n ~ X)) #com Y_boot, estimamos beta_boot 
  beta0b[i] <- coefs_boot[1]
  beta1b[i] <- coefs_boot[2]
}

hist(beta0b)
hist(beta1b)

quantile(beta0b, c(0.025, 0.975))
quantile(beta1b, c(0.025, 0.975))


###c) Utilize um teste de permutação para verificar H0: B1 = 0 
#Note que sob H0, Y = Beta0 + e, ou seja, X nao tem efeito algum sob Y
#Assim, podemos permutar os valores Yobs|Xobs e X nao deveria ter efeito algum em Y. 

#Vamos estimar beta1^ 
beta1hat <- coef(lm(Y ~ X))[2]

#Agora, vamos permutar os dados e calcular a dist de Beta1^ sob H0 
est_test <- numeric(5000) #B = 5000 permutações

for (i in 1:5000){
  Y_perm <- sample(Y)
  est_test[i] <- coef(lm(Y_perm ~ X))[2]
}

#histograma da dist de beta1_perm em H0 
hist(est_test,50)

#vamos agora calcular o p-valor de H0 ter gerado o valor beta1^ = -1.15
round(mean(est_test<=beta1hat),2)
#ou seja, aproximadamente 0% das permutações geraram um valor tão ou mais extremo que beta1^ sob H0. 
#assim, podemos rejeitar a 5% a hipótese que B1 = 0. 

#quantile(est_test, c(0.025, 0.975))

###d) Este teste pode ser usado para testar B0 = 0 ?
#Não, ele faz sentido para B1 pois ao setar B1 = 0 temos um efeito na relação entre X e Y: X nao afeta Y e assim podemos testar se isto é verdade.
#Porém ao setar B0 = 0, Y = Beta1*X + e
#Ou seja, ainda temos o efeito entre as variáveis presente. 
#Assim, gerar valores permutados nesta distribuição não iria ser significativo para Beta0. 

