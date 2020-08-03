########################################Poisson Sports Ranking################################
#This script fits a Poisson Ranking model to Brazilian soccer data. 
#We define X_ij ~ poisson(lambda) where lambda = o_i - d_j = expected goals from i against j on a match
#o_i is the ofensive power of team i. This is usually computed as average(goals from i)\average(goals from league)
#d_i is the defensive power of team j. Usually computed as average(goals against j)\average(goals against league)

#Due to nature of this problem, log-likelihood function uses information from d_j to estimate o_i 
#And information from o_i to estimate d_j
#A typical approach to this kind of optimization is Block Relaxation, where first we optimize o_i, then d_j.
#We also consider L-BFGS-B/optim and NR with inequality constraint. 


##Entering real data from Campeonato Brasileiro 2018
#The i_j entry indicates goals from i against j 
X <- matrix(c(0,1,3,1,1,0,0,0,1,2,0,1,2,0,0,2,1,3,2,2,
              0,0,3,1,1,2,3,1,1,0,5,0,0,1,2,3,1,5,0,2,
              4,1,0,2,2,2,5,1,2,3,3,2,2,1,3,2,0,4,1,4,
              1,2,0,0,3,2,1,1,0,0,2,0,0,1,2,1,2,2,3,4,
              1,0,2,0,0,0,1,1,1,2,2,2,1,1,2,0,2,2,1,1,
              2,2,0,0,0,0,3,2,0,0,1,0,1,2,1,1,0,1,0,2,
              1,1,2,1,0,2,0,2,2,3,1,1,2,1,1,0,1,2,1,0,
              1,1,0,2,2,1,0,0,2,0,2,0,1,1,1,1,1,2,1,0,
              3,0,2,1,1,0,3,1,0,0,2,0,0,1,3,2,0,2,1,3,
              2,2,1,2,2,0,2,1,1,0,3,2,2,1,2,1,0,4,1,1,
              1,1,2,1,1,0,3,1,1,0,0,0,0,1,4,0,1,0,0,0,
              1,2,0,2,4,3,2,1,1,2,0,0,0,0,2,5,2,3,2,4,
              2,1,2,2,3,1,3,2,0,2,2,1,0,0,1,2,3,0,3,2,
              4,3,2,3,2,2,0,1,3,1,3,2,1,0,3,3,3,2,1,3,
              1,0,0,1,1,0,1,0,1,0,2,0,1,1,0,0,1,1,1,1,
              0,3,1,2,1,2,0,1,0,1,3,0,1,1,3,0,0,3,1,5,
              1,2,0,1,3,1,2,3,1,2,1,1,0,0,1,1,0,0,2,3,
              0,3,1,2,1,1,1,1,0,0,1,0,2,0,1,2,1,0,2,0,
              4,2,1,2,1,1,3,1,2,1,1,1,1,0,1,0,2,3,0,2,
              1,1,1,2,3,2,1,2,1,2,1,0,2,0,1,0,0,1,1,0),ncol=20,nrow=20,byrow=T)




#Block Relaxation/Cyclic Coordinate Descent Implementation
##########################################################

n=20; #20 teams

d = rep(1,20) #we set ofensive and defensive initial guess values as 1 
o = rep(1,20)

O = matrix(ncol=20,nrow=n) #to visualize algorithm, we construct a matrix were column = optization step
D = matrix(ncol=20,nrow=n) #and row = team parameter
O[,1] = o
D[,1] = d


epslon = 1e-15 #we set a convergence criteria

i = 1  

#####Cyclick Coordinate Descent: 
#block updating: we estimate oMLE1 using dMLE1 and initial guesses, and then dMLE2 using oMLE1. 

poisson_block_opt <- function(epslon, X, o, d, O, D){ #block optimization of X log-lik function
  i <- 1
  d_old <- rep(0, n)
  o_old <- rep(0, n)
  while(sum((d-d_old)^2) + sum((o-o_old)^2) > epslon){   #convergence criteria
    o_old <- o                          
    o <- log(rowSums(X)/sum(exp(-d)))  
    d_old <- d
    d <- -log(colSums(X)/sum(exp(o)))  
    i <- i+1
    O[,i] <- o
    D[,i] <- d
  }
  return(list(O,D))
}



results <- poisson_block_opt(epslon = 1e-15, X = X, o=o, d=d, O=O, D=D) #poisson sports ranking com block opt 
#para jogos em Casa dos times
O_est_BO <- results[[1]] ##ofensive estimate
D_est_BO <- results[[2]] ##defensive estimate

O_est_BO <- O_est_BO[,3]
D_est_BO <- D_est_BO[,3]

names(O_est_BO) <- names(D_est_BO) <- c("AMM","ATM","ATP","BAH","BOT","CEA","CHA","COR","CRU","FLA",
                                        "FLU","GRE","INT","PAL","PAR","SAN","SPA","SPT","VAS","VIT")


#Direct -loglikelihood minimization using optin and L-BFGS-B
#We use L-BFGS-B to add lower>0 bounds on our parameters
##########################################################

### simple loglikelihood function
logL = function(t,n){
  temp = 0
  for(i in 1:n){
    for(j in 1:n){
      temp = temp - exp(t[i]-t[n+j]) + X[i,j] * (t[i]-t[n+j]) #elements to be added from each pair of observed teams
    }
  }
  -temp
}


n <- 20 #20 teams
o <- numeric(n)  
d <- numeric(n)
t = c(o, d) #reunites ofensive and defensive parameters in a single vector

# creating initial guesses from each parameter
t0 = c(rep(2,20), rep(1,20)) 

control <- list(factr = 1e3) #setting L-BFGS-B convergence with smaller eps

est_pars <- optim(par = t0, fn = logL, method = 'L-BFGS-B', 
                  lower =  c(rep(0.1,40), upper = c(rep(5, 40))),
                  control = control, n=20)$par  #optmization

O_est_LBFGSB <- est_pars[1:20]
D_est_LBFGSB <- est_pars[21:40]
names(O_est_LBFGSB) <- names(D_est_LBFGSB) <- c("AMM","ATM","ATP","BAH","BOT","CEA","CHA","COR","CRU","FLA",
                                                "FLU","GRE","INT","PAL","PAR","SAN","SPA","SPT","VAS","VIT")


#Direct -loglikelihood minimization using Newton-Rapson 
#There is no lambda*g(x) but we need to positive-bound our parameters
##########################################################

#problem: we need to set inequality bounds on our parameters (>0)
#lets update our log-likelihood function to greatly penalize less/equal to 0 parameters
mod_logL <- function(t,n){
  temp = 0
  for(i in 1:n){
    for(j in 1:n){
      if((t[i] <= 0) | (t[n+j] <= 0)){ #if parameter is <= 0, returns a big positive value to minimization of neg. log.lik
        temp <- -1e10
        return(-temp)
      }
      else{
        temp <- temp - exp(t[i]-t[n+j]) + X[i,j]*(t[i]-t[n+j]) #else, continue optimization
        #due to nature of the problem, this will never be bigger than 1e10 for bad <= 0 parameters
      }# elementos a serem somados para cada par de times registrado
    }
  }
  -temp
}

est_pars_NR <- nlm(mod_logL, p = t0, n = 20)$estimate ##optimization
O_est_NR <- est_pars_NR[1:20]
D_est_NR <- est_pars_NR[21:40]
names(O_est_NR) <- names(D_est_NR) <- c("AMM","ATM","ATP","BAH","BOT","CEA","CHA","COR","CRU","FLA",
                                        "FLU","GRE","INT","PAL","PAR","SAN","SPA","SPT","VAS","VIT")


###########################################
#now, lets compare the three methods
library(tidyverse)

##First, lets compare each team's 'absolute' score, being ofensive + defensive score for each method
OD_est3 <- matrix(0, 20, 3)
OD_est3[,1] <- O_est_BO + D_est_BO
OD_est3[,2] <- O_est_LBFGSB + D_est_LBFGSB
OD_est3[,3] <- O_est_NR + D_est_NR
OD_est3 <- as.data.frame(OD_est3)
colnames(OD_est3) <- c("BO", 'LBFGSB','NR')
OD_est3 <- OD_est3 %>% gather('Method', 'parvalue')
median <- OD_est3 %>% group_by(Method) %>% summarize(mediana = median(parvalue)) 

ggplot(OD_est3, aes(x = parvalue, color = Method, fill= Method)) + geom_density() + facet_wrap(~Method)+
  geom_vline(data = median, aes(xintercept=mediana),
             linetype="dashed")

#We can see the 3 methods produced different optimizations, where Block-Opt produced Ofensive+Defensive parameter
#distribution with lowest median, and L-BFGS-B produced O+D param. distr. with highest median. 
#However, the three distributions seem to have a similar shape.


###Lets now check for the 5 better teams in each method:

tail(sort(O_est_BO + D_est_BO),5) ##sorting for defensive + ofensive power
#     FLA      INT      ATP      GRE      PAL 
#2.628609 2.793231 2.920629 3.290503 3.444654 

tail(sort(O_est_LBFGSB + D_est_LBFGSB),5)
#FLA      INT      ATP      GRE      PAL 
#4.107937 4.272560 4.399958 4.769832 4.923983 

tail(sort(O_est_NR + D_est_NR),5)
#FLA      INT      ATP      GRE      PAL 
#3.355717 3.520335 3.647738 4.017618 4.171771 

#We can see the rank is the same for each method, what seems to change is the parameters absolute value. 


###Now, lets compute lambda = expected goals from i against j 

##expected value of goals when ATP faces SPT:
exp(O_est_BO["ATP"] - D_est_BO["SPT"]) #lambda
#3.352381 goals expected from ATP against SPT

exp(O_est_LBFGSB["ATP"] - D_est_LBFGSB["SPT"]) #lambda
#3.35238

exp(O_est_NR["ATP"] - D_est_NR["SPT"]) #lambda
#3.352392 

#Here, we can see the expected goals are the same in every method, even though the different parameter values. 

#This suggests the three algorithms can be used to estimate this sports rank example. 
#However, L-BFGS-B and Block Optimization can be thought of more natural approaches, since the derivation of 
#log-likelihood function is straightforward. 
#In Newton-Rapson log-lik function, we had to add a conditional statement to manually bound the parameters. 
#In this case, we could see adding this condition wouldn't give us bad funtion values, but in more complicated functions
#this could cause problems. Also, this makes analytical gradient/hessian calculation less intuitive. 


