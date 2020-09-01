
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
