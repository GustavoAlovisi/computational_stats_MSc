##############k means clustering k=2: otimização 

library(tidyverse)

#f_obj <- function(xi, mu, zi){ #mediax c2 mediax c1 mediay c2 mediay c1
#  mu1 <- as.numeric(mu[seq(from = 2, to = 2*(ncol(data)-1), by = 2 )])
#  mu2 <- as.numeric(mu[seq(from = 1, to = 2*(ncol(data)-1), by = 2 )])
#  f <- sum(zi * (xi-mu1)^2) + sum((1-zi)*(xi - mu2)^2) #queremos minimizar esta função 
#  return(f) #note que em cada chamada da função, f corresponderá apenas no primeiro ou segundo termo devido a zi ser 0 ou 1 (indicadora)
#}

f_obj <- function(xi, mu, zi){ #mediax c2 mediax c1 mediay c2 mediay c1
  mu1 <- as.numeric(mu[2, ])
  mu2 <- as.numeric(mu[1, ])
  f <- sum(zi * (xi-mu1)^2) + sum((1-zi)*(xi - mu2)^2) #queremos minimizar esta função 
  return(f) #note que em cada chamada da função, f corresponderá apenas no primeiro ou segundo termo devido a zi ser 0 ou 1 (indicadora)
}


opt_k2means <- function(data, mu){
  nrep <- 1
  conv <- FALSE
  while(conv == FALSE & nrep < 1000){
    data_temp <- data
    for(i in 1:nrow(data)){
      if(f_obj(xi = data[i,-which(colnames(data)=='which_cl')], mu = mu, zi = 1) < f_obj(xi = data[i,-which(colnames(data)=='which_cl')], mu = mu, zi = 0)){
        data[i,]$which_cl <- 1
      } else{
        data[i,]$which_cl <- 0 ##0 =  cluster 2 = m2 
      }
    }
    print('a')
    mu <- data %>% group_by(which_cl) %>% summarise_all(mean) %>% select(-which_cl) # %>% pivot_wider(names_from = which_cl, values_from = -which_cl) #computando as medias de cada cluster
    #note que ao derivarmos a f.obj em relacao aos vetores de medias m1 e m2, chegamos a conclusao que m1^ é a soma das observações pertencentes ao cluster 1 dividido por a quantidade 
    #de observações presente no cluster 1. Para o m2^, vale a mesma ideia.
    
    if(dplyr::all_equal(data, data_temp) == T){ ##criterio de convergencia: a atualização dos clusters não resultou em nenhuma mudança nova que minimizasse a f. objetivo
      conv <- TRUE
    }
    nrep <- nrep + 1  
  }
  return(list(data, mu))
}


###vamos criar 2 clusters que diferem bastante entre si para checarmos se nosso algoritmo funcionou:
x <- 0
x[1:50] <- rnorm(50, 7, 2)
x[51:100] <- rnorm(50, 2, 1)

y <- 0 
y[1:50] <- rnorm(50, -1, 1)
y[51:100] <- rnorm(50, 3, 1)

z <- 0
z[1:50] <- rnorm(50, 10, 1)
z[51:100] <- rnorm(50, -3, 1)
data <- as.tbl(as.data.frame(cbind(x,y)))
 
###vamos usar um chute inicial para os centroides, que constituem em valores aleatorios amostrados de nossos dados
m1 <- sample_n(data, 1) ###cluster 1 inicial
m2 <- sample_n(data, 1) ###cluster 2 inicial
mu <- rbind(m1, m2)

data$which_cl <- rep(NA, nrow(data))  ##criar coluna para atribuirmos o cluster a cada vetor de obs

##visualizando nossos dados simulados
ggplot(data, aes(x = x, y = y)) + geom_point()

##otimizando com kmeans
opt <- opt_k2means(data, mu = mu) ###k2 means clustering
a <- opt[[1]] ###salvando os dados e seu cluster
b <- opt[[2]] ###salvando o vetor de medias (centroides)

##visualizando a classificação
ggplot(a, aes(x = x, y= y, color = which_cl)) + geom_point()

##vamos agora comparar nossa classificação com a classificação do kmeans():
model <- kmeans(data[,1:2], centers = 2)
data[,3] <- model$cluster
ggplot(data, aes(x=x, y=y, color = which_cl)) + geom_point()
#podemos ver que ambas classifcações são muito parecidas
