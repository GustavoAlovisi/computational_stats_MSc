#############questão 2
######Minimos quadrados nao lineares / Newton Rapson e Gauss-Newton 

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


########resolvendo através do método de Newton-Rhapson: 
f <- function(beta, xi ,yi){
  res <- sum((yi - beta[1] * exp(beta[2] * xi))^2)
  attr(res, "gradient") <- gradf(beta, xi, yi)
  
  attr(res, 'hessian') <- hessf(beta, xi, yi)
  return(res)
}
nlm(f, c(40,-2), xi, yi)
#podemos ver que as estimativas para os betas foram bem próximas de seus valores reais 




############resolvendo com gauss-newton: nao funcionando por enquanto

f <- function(beta, xi ,yi){
  sum((yi - beta[1] * exp(beta[2] * xi))^2)
}

hessapprox <- function(beta, xi){
  grad1 <- rbind(0,0)
  grad1[1] <- sum(exp(beta[2]*xi))
  grad1[2] <- sum(-beta[1]*xi*exp(beta[2]*xi))
  return(grad1 %*% t(grad1))
}

newton_rhapson <- function(fun, set_tol, gradf, hessapprox, x0, xi, yi){
  tol <- 1
  i <- 0
  z <- x0
  while(tol > set_tol){
    x1 <- x0 - as.numeric(solve(2*t(gradf(x0, xi,yi)) %*% gradf(x0,xi,yi))) * gradf(x0, xi, yi)
    tol <- sum((x1-x0)^2)                                                         
    x0 <- x1
    i <- i+1
    z <- cbind(z,x1)
   # print(z)
  }
  return(z)
}

opt <- newton_rhapson(fun = f, set_tol = 1e-15, gradf = gradf, 
                      hessapprox = hessapprox, x0 = c(30,-5), xi = xi, yi = yi)
