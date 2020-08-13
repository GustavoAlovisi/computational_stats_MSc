######Mistura de Copulas Normais Multivariadas: otimização por EM, Sim. Annealing e Nelder-Mead######

set.seed(182)
library(mvnfast)
library(ggplot2)
library(mvtnorm)

######a) Vamos simular um conjunto de dados (V1,V2,V3) provenientes de duas misturas de normais multivariadas
#(V1,V2,V3) = 0.7*mvtnorm(mu1, I) + 0.3*mvtnorm(mu2, I)
#Onde mu1 = t(0,0,0) e mu2 = t(3,3,3), I = matriz identidade 3x3

mu <- matrix(c(0,0,0,
               3,3,3), ncol = 3, nrow =2, byrow =T) #matriz de medias (para esta função eh transposto)

covmat <- diag(3) #a covariancia eh a matriz identidade (var 1 e cov 0)

#pesos da mistura 
w1 <- 0.7
w2 <- 0.3
N <- 40  #tamanho da amostra a ser simulada

#simulando uma mistura de normais multivariadas. retInd = T nos da um atributo que indica se o elemento (V1i, V2i, V3i)
#veio da mistura 1 ou 2 que iremos utilizar para análises posteriores:
normalmvtmix <- mvnfast::rmixn(n = N, mu = mu, sigma = list(covmat, covmat), w = c(w1,w2), retInd = T)

#plotando a marginal V1, podemos ver que de fato ela parece ser uma mistura de duas distribuicoes
normalmvtmix1 <- as.data.frame(normalmvtmix)
ggplot(normalmvtmix1, aes(x=V1)) + geom_density()
normalmvtmix1 <- NULL


#####b) Vamos obter os parâmetros das misturas de normais utilizando o metodo de Simulated Annealing
#problema: devemos obter um 'd' adequado para nosso algoritmo

#defininindo a função de log-verossimilhança de nosso problema
mvtmixLL <- function(m,x){
  m1 = m[1:3]
  m2 = m[4:6]
  #p1 = m[7] ##estimating weights too
  #p2 = m[8] ##w2 
  sum(log(w1*mvtnorm::dmvnorm(x = x, mean = m1)+ w2*mvtnorm::dmvnorm(x = x, mean = m2)))
} #soma está positiva, devemos utilizar fnscale = -1 no optim() 

##definindo a função de Simulated Anealing
#através de testes realizados, d <- 0.1sqrt(temp) parece ser um bom valor para d neste problema. 
#valores mais altos retornam piores resultados na otimização
GMM_SA <- function(n.iter, y, x0){
  teta <- matrix(NA, nrow = n.iter, ncol = length(x0))
  teta[1,] <- x0
  f.eval <- rep(0, n.iter)
  f.eval[1] <- mvtmixLL(teta[1,], x = y)
  for(i in 2:n.iter){
    temp <- 1/log(i+1)
    #temp <- 1/log(1+i) 
    d <- 0.1*sqrt(temp)
    e <-rnorm(6)*d #runif(6, -d, d)
    te1 <- teta[(i-1),] + e
    u <- runif(1)
    prob <- min(exp((mvtmixLL(m = te1, x = y)-mvtmixLL(m = teta[(i-1),], x = y))/temp), 1)
    #print(prob) #printar a probabilidade de aceitação para calibragem do algoritmo
    teta[i,] <- (u<=prob)*te1+(u>prob)*teta[(i-1),]
    f.eval[i] <- mvtmixLL(m = teta[i,], x = y)
  }
  #rnorm(6)*0.6
  pos <- which.max(f.eval)
  teta.ot <- teta[pos,]
  return(teta.ot)
}


x0 <- c(-1, 1, 0, 4, 2, 8) #initial guesses para a primeira iteração
x1 <- c(8, 9, -1, 4, 3, 6) #initial guess mais distante

teta_SA <- GMM_SA(n.iter = 500, y = normalmvtmix, x0 = x0) #otimizando via SA
teta_SA
#[1]  0.1239047  0.2901662 -0.1826373  3.0737463  2.8966653  3.2500398


teta_SA <- GMM_SA(n.iter = 500, y = normalmvtmix, x0 = x1) #otimizando via SA
teta_SA 
#[1]  7.3154158  8.6761766 -2.3489667  0.8210402  0.9318625  0.687138
##analisando chutes inicias, podemos ver que os resultados para teta são próximos dos tetas reais apenas se o chute inicial for um bom chute
##para o Simulated Annealing neste caso. 


#alternativamente, poderiamos ter feito Simulated Annealing utilizando optim:
optim(par = x0, fn = mvtmixLL, x = normalmvtmix, method = 'SANN', control = list(temp = 15, tmax = 30, fnscale = -1))$par
#[1]  0.09850878  0.37002510 -0.11494675  2.96382203  2.83863939  3.19747390

optim(par = x1, fn = mvtmixLL, x = normalmvtmix, method = 'SANN', control = list(temp = 15, tmax = 30, fnscale = -1))$par
#[1] 9.2777542 9.0986352 0.2013623 0.8537164 0.9482063 0.7030627
#com a funcao do optim, a otimizacao parece ser igualmente sensivel aos parametros iniciais. 
#A simulação estocastica em algumas vezes rodando a otimização produziu valores bem proximos aos reais, e as vezes valores distantes.


#####c) Vamos agora derivar o algoritmo EM para o problema de k=2 misturas de normais multivariadas 
#Vamos definir os passos E,M, o algoritmo e fazer algumas análises baseadas na clusterização que o algoritmo resulta
#Note que no algoritmo, pi_i = E(Z_i | X_i, mu1, mu2) é a probabilidade de cada vetor X_i pertencer a mistura 1
#Analogamente, 1 - pi_i é a probabilidade de cada X_i pertencer a mistura 2. 
#Assim, ao realizarmos o cálculo de pi_i, i = 1..n, podemos classificar cada vetor de dados em relação a qual mistura ele pertence. 
#Como ao gerar as 'n' observacoes do vetor da mistura, obtivemos tambem a qual mistura o dado real pertence, podemos determinar quais elementos
#foram classificados de forma correta ou errada. 



##função do passo E. Poderiamos generalizar para k misturas desconhecidas.
#Note que o EM neste caso também está estimando os pesos das misturas que será utilizado posteriormente na g). 
expect_step <- function(o, y, weights){
  o1 <- o[1:ncol(y)] #salva parâmetros mu da primeira mistura 
  o2 <- o[(ncol(y)+1):(2*ncol(y))] #salva parâmetros mu da segunda mistura
  w1 <- weights[1] #salvando pesos 1 e 2
  w2 <- weights[2]
  r1 <- w1 * mvtnorm::dmvnorm(x = y, mean = o1) / (w1 * mvtnorm::dmvnorm(x = y, mean = o1) + w2 * mvtnorm::dmvnorm(x = y, mean = o2))
  r2 <- w2 * mvtnorm::dmvnorm(x = y, mean = o2) / (w1 * mvtnorm::dmvnorm(x = y, mean = o1) + w2 * mvtnorm::dmvnorm(x = y, mean = o2))
  return(cbind(r1, r2))
}


##função do passo M. Poderiamos generalizar para k misturas desconhecidas. 
max_step <- function(pi, y){
  teta1.n <- sum(pi[,1])  #total alocado para a mistura 1
  teta2.n <- sum(pi[,2])  #total alocado para a mistura 2
  tetak1 <- tetak2 <- c(rep(0, ncol(y)))
  
  tetak1 <- sapply(1:ncol(y), function(i){ #estimando vetor de médias da mistura 1
    1/teta1.n * sum(pi[,1]*y[,i]) 
  })
  
  tetak2 <- sapply(1:ncol(y), function(i){ #estimando vetor de médias da mistura 2
    1/teta2.n * sum(pi[,2]*y[,i])
  })
  
  w1 <- teta1.n/(teta1.n + teta2.n) #estimando o peso da mistura 1. Note que o peso da mistura 1 eh apenas o que foi alocado para a mistura 1/soma de tudo alocado
  w2 <- teta2.n/(teta1.n + teta2.n) #estimando o peso da mistura 2
  
  return(list(c(tetak1, tetak2), c(w1,w2)))
}

##função que realiza a estimação via algoritmo EM. 
GMM_EM <- function(eps, nrep, Y, t0, weights0){
  cc <- 1
  rep <- 1
  t <- matrix(ncol = 1000, nrow = 2*ncol(Y))
  t[,1] <- t0
  while(cc > eps & rep < nrep){
    #etapa E
    pi <- expect_step(o = t0, y = Y, weights = weights0)
    
    #etapa M
    max_res <- max_step(pi = pi, y = Y)
    t1 <- max_res[[1]]
    weights0 <- max_res[[2]]
    rep <- rep + 1
    cc <- sum((t0 - t1)^2)
    t0 <- t1
    t[,rep] <- t0
  }
  t <- t[, 1:rep] #salvando os valores de todos os passos realizados para analise do algoritmo excluindo NA
  return(list(t, pi, weights0))
}

w0 <- c(0.4, 0.6) #chute inicial para os pesos das misturas

EM_result <- GMM_EM(eps = 1e-8, nrep = 100, Y = normalmvtmix, t0 = x0, weights0 = w0) ##EM

theta_EM <- as.matrix(EM_result[[1]]) ##salvando os valores dos parâmetros (mu's) em relação a todos os passos
round(theta_EM,3)
#podemos ver que para n=40 os resultados estimados são relativamente próximos dos parâmetros reais.
#Com aumento do tamanho da amostra, porém, a estimação fica cada vez mais próxima dos parâmetros reais. 
#Tambem, o algoritmo eh robusto em relacao aos pontos iniciais, ao contrario de nossa implementacao via Sim. Annealing. 

pi_EM <- EM_result[[2]] ##obtendo os Pi's para cada y_i. Como temos dois elementos da mistura, podemos usar uma regra de decisão 
#em que se pi_k=1 >= 0.5, o vetor de obs. y_i pertence a mistura 1. se pi_k=1 <0.5, o vetor pertece a mistura 2. 
round(head(pi_EM),3)

pesos_EM <- EM_result[[3]] ##obtendo os pesos estimados das misturas
pesos_EM #para N=40, o peso estimado é proximo, mas a estimação melhora bastante conforme o tamanho da amostra aumenta 
#em relação aos pesos, o algoritmo parece estimar bem independente do chute inicial para os mesmos.
#[1] 0.6499474 0.3500526


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
#atenção: existem problemas de métodos locais, e dependendo do chute inicial, existem casos em que a maximização resulta em 
#um vetor estimado ~(3, 3, 3, 0, 0, 0) ao invés do contrário (0,0,0,3,3,3). Isto pode levar a classificação a produzir resultados errôneos. 

####Para visualizarmos melhor a classifcação, vamos considerar o caso de uma mistura de normais bivariadas: 
library(dplyr)

w1 <- 0.7
w2 <- 0.3
mu_bvt <- matrix(c(0,0,
                   3,3), ncol = 2, nrow =2, byrow =T) #gerando uma mistura k=2 de normais bivariadas

covmat_bvt <- diag(2)

normalbvtmix <- mvnfast::rmixn(n=100, mu = mu_bvt, sigma = list(covmat_bvt, covmat_bvt), w = c(w1,w2), retInd = T)

t0 <- c(0, -1, 2, 4) ##chute inicial
w0 <- c(0.5, 0.5)
EM_result <- GMM_EM(eps = 1e-8, nrep = 100, Y = normalbvtmix, t0 = t0, weights0 = w0) ##EM

theta_EM <- as.matrix(EM_result[[1]]) 
pi_EM <- EM_result[[2]] #salvando os pi's para realizar algumas análises visuais 
peso_EM <- EM_result[[3]] #salvando os pesos

yi_classificado <- 0
yi_classificado <- cbind(normalbvtmix, ifelse(pi_EM[,1] >= 0.5, 1, 2)) ##classificando y_i conforme pertencer a mistura 1 ou 2
yi_classificado <- cbind(yi_classificado, attributes(normalbvtmix)$index)
colnames(yi_classificado) <- c("Y1", "Y2", "Mist_Est", "Mist_Real")

##visualizando nossa classificação de dois clusters: 
yi_classificado <- dplyr::as.tbl(as.data.frame(yi_classificado))
ggplot(yi_classificado, aes(x = Y1, y = Y2, color = Mist_Est)) + geom_point()

######Vamos comparar a classificação de nosso algoritmo EM com a classificação de um pacote de misturas gaussianas (flexmix)
#para checar se estamos de fato classificando corretamente:

#install.packages('flexmix')
#install.packages('ellipse')
library(flexmix)
fit_with_covariance <- flexmix(normalbvtmix ~ 1,
                               k = 2, 
                               model = FLXMCmvnorm(diag = T),
                               control = list(tolerance = 1e-08, iter.max = 1000))

plotEll(fit_with_covariance, data = normalbvtmix)
###podemos ver que no caso da normal bivariada, a classificação parece ter sido muito parecida em ambos os casos.



##por ultimo, vamos visualizar quais pontos foram erroneamente classificados por nosso algoritmo:
yi_classificado <- yi_classificado %>% mutate(Class_Errado = abs(Mist_Est - Mist_Real))

ggplot(yi_classificado, aes(x = Y1, y = Y2, color = as.factor(Class_Errado))) + 
  geom_point() 





#####d) Vamos utilizar optim() para resolver o mesmo problema.
#Podemos utilizar o metodo de Simulated Annealing ('SANN') do optim como fizemos na letra b) 
#Ou ainda Nelder-Mead uma vez que temos apenas uma otimização sem restrição da log-verossimilhança. 
#Nota: Ao utilizarmos pesos variáveis, eh necessario a restrição em que w1+w2 = 1. 
#Desta forma, seria interessante usar 'L-BFGS-B' e limitar w1 entre 0 e 1. 
#Com mais de 2 misturas, porém, devemos pensar em como restringir a soma dos pesos, penalizando a função de verossimilhança
#ou utilizando outro algoritmo como Augmented-Lagrange. 

#chute inicial melhor:
optim(par = x0, fn = mvtmixLL, x = normalmvtmix, method = 'Nelder-Mead', control = list(fnscale = -1))$par
#sensível ao ponto inicial! 
#[1] 0.53246960 0.22688230 0.06694153 3.62802198 3.33742254 3.20681230

#chute inicial pior:
optim(par = x1, fn = mvtmixLL, x = normalmvtmix, method = 'Nelder-Mead', control = list(fnscale = -1))$par
#sensível ao ponto inicial! 
#[1]  8.7926401 12.8008444  4.0747302  0.8426125  0.9472383  0.6966989


#####e) Com a implementação e teste dos três métodos de otimização (Sim. Ann, EM, Nelder-Mead), 
#podemos observar que para nossos dados simulados, o método EM é o mais robusto em relação aos parâmetros iniciais,
#ainda, com aumento da amostra, os parâmetros estimados se aproximam bastante dos parametros reais. 
#O metodo de EM ainda nos fornece uma ferramenta capaz de classificar cada vetor de observações (i = 1..n) em relação a probabilidade
#de pertencer em cada mistura. 
#Os algoritmos de SA e Nelder-Mead oferecem bons resultados apenas se dermos bons chutes iniciais. 
#Existe ainda a questão de ótimos locais em nosso problema. 
#Para Nelder-Mead, com um determinado chute inicial [c(3, 4, 5,-1,-2,-6)], os parâmetros encontrados foram ~c(3,3,3,0,0,0) 
#ao invés de c(0,0,0,3,3,3). Isto sugere que ele acabou convergindo para um otimo local. 
#O mesmo acontece para Simulated Annealing e também para o algoritmo EM. 
#Analisando a convergência, o método EM é drasticamente mais rápido. Em uma estimação com chute inicial ruim, o algoritmo convergiu em ~7 iterações. 
#Com chutes iniciais mais próximos, ele chegou a convergir em 3 iterações. 
#O método de Simulated Annealing acaba demorando mais que os outros para rodar todas as suas iterações.
#Com Nelder-Mead, a convergência ocorreu entre 200 e 600 iterações, dependendo da precisão do chute inicial. 

#####f) Vamos agora realizar a mesma análise de e) porém com vetor de médias da segunda mistura mu2 = t(1,1,0)
#gerando duas misturas com médias mais parecidas do que da primeira vez:

mu <- matrix(c(0,0,0,
               1,1,0), ncol = 3, nrow =2, byrow =T) #matriz de medias

covmat <- diag(3) #a covariancia eh a matriz identidade (var 1 e cov 0)

normalmvtmix <- mvnfast::rmixn(n = N, mu = mu, sigma = list(covmat, covmat), w = c(w1,w2), retInd = T)

x1 <- c(-1, 1, 0, 4 ,2 ,8)#mesmos chutes iniciais de antes, para compararmos

#Simulated Annealing
teta_SA <- GMM_SA(n.iter = 1000, y = normalmvtmix, x0 = x1) #otimizando via SA
teta_SA
#[1]  7.78488771  8.79780235 -0.89209191  0.26988808  0.28823402  0.01122817

#Optim Simulated Annealing
optim(par = x1, fn = mvtmixLL, x = normalmvtmix, method = 'SANN', control = list(temp = 15, tmax = 30, fnscale = -1))$par
#[1] 14.140905244  3.708064862 -4.750182133  0.268029406  0.281554807  0.002570182

#EM
w0 <- c(0.5, 0.5) #chute inicial dos pesos
EM_result <- GMM_EM(eps = 1e-8, nrep = 100, Y = normalmvtmix, t0 = x1, weights0 = w0) ##EM
theta_EM <- as.matrix(EM_result[[1]]) ##salvando os valores dos parâmetros (mu's)
round(theta_EM,3)
#[1,] 0.522 0.566 -0.004 -0.311 -0.357 0.033

#Note que neste caso o algoritmo EM rodou até a última repetição, apesar de fornecer uma estimação relativamente próxima dos parâmetros.
#Assim, podemos aumentar o nrep para obter convergência. Aumentando para 300, o algoritmo convergiu em 180 passos.
#Podemos ver que misturas com médias parecidas são mais difíceis de serem classificadas, ainda mais com poucas observações. 

#Nelder-Mead
optim(par = x1, fn = mvtmixLL, x = normalmvtmix, method = 'Nelder-Mead', control = list(fnscale = -1))$par
#[1] 11.733734894 11.775296159  0.021060133  0.265933067  0.282902967  0.007179512
#sensivel ao ponto inicial


#Podemos ver que quando as misturas possuem parâmetros parecidos, o algoritmo EM resulta em uma otimização melhor que os outros métodos, apeasr de convergir com bem mais passos.
#Ao utilizarmos chutes iniciais distantes, a estimação por EM foi consideravelmente melhor para nossos dados especificos simulados. 
#Conforme o tamanho da amostra aumenta, a precisão da estimação do método também aumenta, o que não foi observado com clareza para Sim. Annealing e Nelder Mead neste caso. 
#O problema de máximos locais também continua neste exemplo. Ao mudarmos os chutes iniciais, podemos ver que os algoritmos convergem para diferentes valores. 
#Em relação ao tempo de convergência, o algoritmo EM passou a levar mais iterações para convergir, (~20) o que ainda é bastante baixo e consideravelmente mais rapido que os outros métodos. 
#Com um bom chute, o método de Nelder-Mead convergiu em cerca de 500 iterações. 



####g) vamos usar nosso EM para classificar os dados genéticos em duas misturas normais multivariadas 
dados_gen <- read.table('Dadosgeneticos.txt', header = T) #lendo os dados
dados_gen <- as.tbl(as.data.frame(t(dados_gen))) #transformando em tbl 

##vamos remover colunas que sao todas 0 
ind <- sapply(dados_gen, function(x) sum(x == 0)) != nrow(dados_gen)
dados_gen <- dados_gen[,ind]


w0 <- c(0.5, 0.5) #chute inicial para a mistura
x1 <- c(rep(-0.1, 78), rep(0.6, 78)) #chute inicial para as medias. Removendo as colunas de 0, temos 78 variaveis. 
EM_result_gen <- GMM_EM(eps = 1e-8, nrep = 100, Y = dados_gen, t0 = x1, weights0 = w0)

theta_EM_gen <- EM_result_gen[[1]] #estimação dos parâmetros 
head(theta_EM_gen[,3])

pi_EM_gen <- EM_result_gen[[2]] #estimação dos pi's (a qual mistura pertence)

pesos_EM_gen <- EM_result_gen[[3]] #estimação dos pesos
pesos_EM_gen
#[1] 0.9993578202 0.0006421798
#No caso deste exemplo, o algoritmo EM atribui peso ~1 para a primeira mistura. Ao checar a estimação dos parâmetros para cada mistura, podemos ver que eles são os mesmos.
#Assim, o algoritmo EM rejeitou a hipótese de existir duas misturas (clusters) para os dados. 

#O algoritmo K-means ao utilizar dois clusters foi capaz de diferenciar as populações europeias das americanas:
kmeans(dados_gen, centers = 2)[1]


