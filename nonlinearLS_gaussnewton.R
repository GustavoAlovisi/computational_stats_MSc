#############questão 2
######Minimos Quadrados Não Lineares / Newton-Raphson e Gauss-Newton 

#vamos simular a regressão 
#Y_i = b0 * exp(b1 * x_i) + e_i 

n <- 200
beta <- c(60)
beta[2] <- -0.05
xi <- runif(n, 2, 40)
ei <- rnorm(n, 0, 2)
yi <- beta[1] * exp(beta[2] * xi) + ei


#queremos utilizar os algoritmos de otimização para minimzar a soma dos quadrados dos residuos da regressão 
#onde fun <- sum(yi - b0 * exp(b1 * x_i))^2

########resolvendo através do método de Newton-Rhapson: 
gradf <- function(beta, xi, yi){ #grad_beta
  g <- rbind(0,0)
  g[1] <- sum(2*(yi - beta[1]*exp(beta[2]*xi))*(-exp(beta[2]*xi)))
  g[2] <- sum(2*(yi - beta[1]*exp(beta[2]*xi))*(-beta[1]*exp(beta[2]*xi)*xi))
  return(g)
}

hessf <- function(beta, xi, yi){ #hessian matrix of f wrt beta[1], beta[2]
  b0 <- beta[1]
  b1 <- beta[2]
  A <- matrix(0, 2,2)
  A[1,1] <- sum(-2*exp(b1*xi)*(-exp(b1*xi)))
  A[2,2] <- sum(2*(-b0*xi*exp(b1*xi))*(-b0*xi*exp(b1*xi))) + sum(2*(yi-b0*exp(b1*xi))*(-b0*xi^2*exp(b1*xi)))
  A[2,1] <- A[1,2] <- sum(2*(-b0*xi*exp(b1*xi))*(-exp(b1*xi))) + sum(2*(yi-b0*exp(b1*xi))*(-xi*exp(b1*xi)))
  return(A)
}

f <- function(beta, xi ,yi){
  res <- sum((yi - beta[1] * exp(beta[2] * xi))^2)
  attr(res, "gradient") <- gradf(beta, xi, yi)
  #attr(res, 'hessian') <- hessf(beta, xi, yi)
  return(res)
}

opt_NR <- nlm(f, c(40,-2), xi, yi) ##otimizando por Newton Rhapson
opt_NR$estimate
#[1] 60.29774970 -0.05050371
#podemos ver que as estimativas para os betas foram bem próximas de seus valores reais 



############resolvendo com gauss-newton: 

#residuos da regresão:
r <- function(beta, xi, yi){
  res <- yi - beta[1]*exp(beta[2]*xi)
  return(res)
}

#matriz do jacobiano:
J_mat <- function(b, x, y){
  Jac <- cbind(-exp(beta[2]*xi), -beta[1]*xi*exp(beta[2]*xi))
  return(Jac)
}


newton_rhapson <- function(set_tol, x0, xi, yi){
  tol <- 1
  i <- 0
  z <- x0
  while(tol > set_tol & i < 1000){
    x1 = x0 - solve(t(J_mat(x0,xi,yi))%*%J_mat(x0,xi,yi)) %*% t(J_mat(x0,xi,yi)) %*% r(x0,xi,yi)
    tol <- sum((x1-x0)^2)                                                         
    x0 <- x1
    i <- i+1
    z <- cbind(z,x1)
  }
  return(z)
}

opt <- newton_rhapson(set_tol = 1e-8, x0 = rbind(40,-2), xi = xi, yi = yi)
opt[ ,ncol(opt)]
#[1] 60.2961746 -0.0505015
#Os valores estimados por Newton-Raphson e Gauss-Newton foram muito próximos. 
#Isto sugere que a aproximação da matriz hessiana pelo método de GN com derivadas primeiras foi uma boa aproximação. 


#Esta aproximação porém é sensível a problemas computacionais
#Ao utilizarmos o método de GN com um chute inicial de c(70, -2), o algoritmo não converge 
opt <- newton_rhapson(set_tol = 1e-8, x0 = rbind(70,-2), xi = xi, yi = yi)
opt[ ,ncol(opt)]

#Para este problema, utilizando Newton-Raphson o algoritmo converge: 
opt_NR <- nlm(f, c(70,-2), xi, yi) ##otimizando por Newton Rhapson
opt_NR$estimate
#Porém, este problema é um problema de baixa dimensão. Para diversas dimensões torna-se dificil o calculo da derivada analítica das funções e da inversa da matriz hessiana/derivadas segundas. 



