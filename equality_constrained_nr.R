

################solving non-linear equality constrained min:
#min f(x1, y) = x^4 + y^2 + 4xy 
#s.t g(x, y) = x^2 + y^2 = 1 (unit circle)
################


#####penalty method for equality: we set a very high positive value (penalty) for deviating from g(x,y)
f_penalty <- function(x){
  z <- x[1]^4 + x[2]^2 + 4*x[1]*x[2] + 1e20*abs(x[1]^2 +x[2]^2 - 1)
  return(z) 
}

##using built-in base R function to solve with Newton-Rapson method
results_f_penalty <- nlm(f_penalty, c(1000,1000),steptol  = 1e-10, gradtol = 1e-10)
results_f_penalty$estimate
#[1] -0.7071068 -0.7071068 = c(x1, x2)



####using classes constructed algorithm, with numerical hessian/gradient calculation:
library(pracma) 

fx <- function(x){
  z <- x[1]^4 + x[2]^2 + 4*x[1]*x[2] + (x[1]^2 +x[2]^2 - 1)
  return(z) 
}

newton_rhapson <- function(fun, set_tol, x0){
  tol <- 1
  i <- 0
  z <- x0
  while(tol > set_tol){
    x1 <- x0 - solve(as.matrix(pracma::hessian(f=fun, x0=x0))) %*% (as.matrix(pracma::grad(f=fun, x0=x0)))
    tol <- sum((x1-x0)^2)                                                         
    x0 <- x1
    i <- i+1
    z <- cbind(z,x1)
  }
  return(z)
}

results_f_penalty2 <- newton_rhapson(fun = fx, set_tol = 1e-15, x0 = c(10,10))
round(results_f_penalty2[,ncol(results_f_penalty2)],6)
#[1]  0.707107 -0.707107 = c(x1,x2)



####bonus: Using Rsolnp (Augmented Lagrange + SQP) with equality constraints that we define
#(this method has given me good results in different occasions)
library(Rsolnp)

f <- function(x){
  X <- x[1]
  Y <- x[2]
  z <- X^4 + Y^2 + 4*X*Y + Y^2 + X^2 -1 
  return(z)
}

equality_const <- function(x){
  z <- x[1]^2 + x[2]^2
  return(z)
}

results_auglag <- Rsolnp::solnp(c(1000,1000), fun = f, eqfun = equality_const, eqB = c(1))
results_auglag$pars
#[1] -0.7071068  0.7071068 = c(x1,x2)



##Given the equality constrained problem on a closed and bounded set (unit circle), we know 
##that the obtained minimum points must be global minimum points because of Extreme Value Theorem.
##Hence, c(x1, y1) = (0.707107, -0.707107) and c(x2, y2) = (-0.707107, 0.707107) are both 
##Global Min. points. This suggests NR algorithm in this case is sensitive to starting points x0,
##since different starting points can lead to one of the two global min. points. 
