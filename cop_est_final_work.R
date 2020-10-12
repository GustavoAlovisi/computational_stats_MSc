

############################################################################
#A Sequential Quadratic Programming Estimation of Mixture Copula Parameters#
############################################################################
#Code by Gustavo Alovisi - PPGEst UFRGS


#####################################################
#Simulating a 8-d Mixture of Clayton-t-Gumbel Copula:

library(copula)

#Clayton par = pi 
#Gumbel par = 2.5
#tCop par = 0.7
#weights = 1/3 each

mC <- mixCopula(list(claytonCopula(pi, dim=8),  #defining the mixture object
                     gumbelCopula(2.5, dim=8),
                     tCopula(0.7, dim=8, df = 5)),
                c(1/3, 1/3, 1/3))
mC

cop_mix_sim <- copula::rCopula(1200, mC)  #generating 1200 8d observations from the mixture
plot(cop_mix_sim)

U <- copula::pobs(cop_mix_sim)  #transforming the observations to pseudo-observations (this is needed for semiparametric copula estimation)
plot(U)


#Sequential Quadratic Programming for maximization of mixture copula pseudo-likelihood:
library(Rsolnp) #Ye's algorithm implemented by Alexios Ghalanos 

#Defining the pseudo log-likelihood function we want to maximize:
LLCG <- function(params,U, copC, copG, copt){ #negative log likelihood function for estimating copulae weights and parameters
  slot(copC, "parameters") <- params[1]    #Initial Clayton parameter provided
  slot(copG, "parameters") <- params[2]    #Initial Gumbel parameter provided 
  slot(copt, "parameters") <- params[3:4]  #Initial t parameters provided (correlation and Degrees of Freedom)
  pi1 <- params[5] #weight of Clayton copula
  pi2 <- params[6] #weight of Gumbel copula
  pi3 <- params[7] #weight of t copula
  
  opt <- log(pi1 * dCopula(U, copC) + pi2 * dCopula(U, copG)
             + pi3 * dCopula(U, copt))  ##log-likelihood function to be optimized
  
  if(any(is.infinite(opt))){    #this is used to avoid computational issues
    opt[which(is.infinite(opt))]<-0
  }
  -sum(opt) #RSOLNP minimizes by default so we minimized negative loglik
}


#Defining the sum(weights) = 1 constraint: 
eqfun <- function(params,U,copC,copG,copt){
  z <- params[5]+params[6]+params[7]
  return(z)
}

#initializing 8d t-copula
copt <- copula::tCopula(param = 0.5, dim = 8) #initial parameters here aren't useful: we will later estimate them for initial guess 
# Initializing archimedian copula objects
copC <- copula::claytonCopula(2, dim = 8) # delta= 2
copG <- copula::gumbelCopula(2, dim = 8)  # theta= 2

##Creating elliptical copula objects and estimating "initial guesses" for each copula parameter. 
#Each copula initial guess is just the max. of log-likelihood imagining the mixture has only one component (the copula itself).
par1 <- copula::fitCopula(copC, U, "itau", estimate.variance = T)@estimate #inversion of Kendall's tau for Clayton 
par2 <- copula::fitCopula(copG, U,"itau", estimate.variance = T)@estimate #inversion of Kendall's tau for Gumbel 
par3 <- copula::fitCopula(copt, U,"mpl", estimate.variance = F)@estimate ###mpl para poder estimar tambem DF. Na documentacao diz que nao pode usar 'itau' pois ele n estima DF.
par4 <- 1/3 #initial guesses for weights = 1/3 each
par5 <- 1/3
par6 <- 1/3

##Setting parameter bounds for the copulas to be used in optimization. 
lower <- c(0.1, 1, -0.9,(2+.Machine$double.eps), 0,0,0)     #lower and upper bounds of the parameters and weights for bounded non linear opt.
upper <- c(copC@param.upbnd, copG@param.upbnd, 1,100, 1,1,1) #2+eps so that variance of t copula is defined


##non linear constrained optimization of the mixture following Augmented Lagrange SQP method 
opt <- Rsolnp::solnp(pars = c(par1,par2,par3,par4,par5,par6), fun = LLCG, LB = lower, 
                     UB = upper, copt=copt,copC = copC, copG = copG, U=U,eqfun = eqfun, 
                     eqB=c(1)) ####RSOLNP
opt



###################################################################################################
#running a Monte Carlo opt. to compute mixture copula parameters expected value and standard error#

MC_ests <- matrix(NA, 100, 7) #we ran 100 iterations due to computation intensity of the method

for(i in 1:100){
  cop_mix_sim <- copula::rCopula(1200, mC)
  U <- copula::pobs(cop_mix_sim)
  
  par1 <- copula::fitCopula(copC, U, "itau", estimate.variance = T)@estimate #inversion of Kendall's tau for Clayton 
  par2 <- copula::fitCopula(copG, U,"itau", estimate.variance = T)@estimate #inversion of Kendall's tau for Gumbel 
  par3 <- copula::fitCopula(copt, U,"mpl", estimate.variance = F)@estimate ###mpl para poder estimar tambem DF. Na documentacao diz que nao pode usar 'itau' pois ele n estima DF.
  par4 <- par5 <- par6 <- 1/3 #initial guesses for weights = 1/3 each
  
  MC_ests[i,] <- Rsolnp::solnp(pars = c(par1,par2,par3,par4,par5,par6), fun = LLCG, LB = lower, 
                               UB = upper, copt=copt,copC = copC, copG = copG, U=U,eqfun = eqfun, 
                               eqB=c(1))$pars
  
}

apply(MC_ests, 2, mean)
#[1] 3.0965449 2.4577772 0.7132326 4.7238911 0.3228383 0.3302576 0.3469041

apply(MC_ests, 2, sd)
#[1] 0.14430112 0.08292706 0.02255114 0.76873906 0.02709860 0.03742296 0.03527428




###########################################################################################################
#Real world application: We use the above Clayton-t-Gumbel copula mixture to model 8 U.S. assets dependence 
#while marginals are modeled using an ARMA-GARCH(1,1). 
#Then, with copula dependence and univariate modeling, we can approximate the assets mvt probability density
#Required on CVaR optimization approach for the portfolio of stocks.




######################Now lets define the numerous functions for application: 

###############Modeling each assets marginals following ARMA(p,q)-GARCH(1,1)
#(p,q) are selected from the combination that maximizes BIC for each optimization

arima_est <- function(data){
  ##ARMA estimation
  xx <- forecast::auto.arima(y=data, max.p=2, max.d = 0,    #auto.arima fit, where we select (p,q) with better BIC criteria
                             max.q=2, seasonal=FALSE, stationary = TRUE,  #max (p,q) = 2, assumes series are already stationary (wont dif)
                             ic = c('bic'), allowmean =FALSE, allowdrift = FALSE)  #wont allow mean or drift since we expect mean(returns) = 0
  ordem <- forecast::arimaorder(xx) #getting arma orders 
  ordem <- c(ordem[1], ordem[3])  #getting [p,q]
  return(ordem)
}

garch_est <- function(data, ordem){
  ##GARCH esitmation 
  armagarchspec <- rugarch::ugarchspec(list(model="sGARCH", garchOrder = c(1,1), variance.targeting = TRUE),
                                       mean.model=list(armaOrder=ordem,include.mean = FALSE), 
                                       distribution.model = 'sstd') #creating garch(1,1) with skewed t specification
  garch_fit <- rugarch::ugarchfit(armagarchspec, data = data, solver = "hybrid") #fitting garch for the j asset
  #we use variance.targeting for faster estimation
  return(garch_fit)
}

coefs_vector <- function(ordem, garch_fit){ ###creating an easy-to-deal matrix of coefs [this is ugly, i know :( ]
  aux_list <- vector('list', 9)
  if(ordem[1] == 0 && ordem[2] == 0){           
    aux_list[1:4] = 0                           ##from the form (AR1, AR2, MA1, MA2,  GARCH coefs)
    aux_list[5:9] = garch_fit@fit$coef[1:5]  
  }else
    if(ordem[1] == 1 && ordem[2] == 0){    #if p=1 and q=0, (ar1 coef, 0, 0, 0, garch coefs)
      aux_list[2:4]=0
      aux_list[1]=garch_fit@fit$coef[1]
      aux_list[5:9]= garch_fit@fit$coef[2:6]
    } else
      if(ordem[1] == 2 && ordem[2] == 0){  #if p=2, q=0, (ar1 coef, ar2 coef, 0, 0 garch coefs)
        aux_list[3:4]=0
        aux_list[1:2]=garch_fit@fit$coef[1:2]
        aux_list[5:9]= garch_fit@fit$coef[3:7]
      }
  else
    if(ordem[1]==2 && ordem[2]==1){  #if p=2, q=1, (ar1 coef, ar2 coef, ma1 coef, 0 garch coefs)
      aux_list[4]=0
      aux_list[1:3]=garch_fit@fit$coef[1:3]
      aux_list[5:9]= garch_fit@fit$coef[4:8]
    }
  else
    if(ordem[1]==2 && ordem[2] ==2){ #if p=2, q=2, (ar1 coef, ar2 coef, ma1 coef, ma2 coef,  garch coefs)
      aux_list[1:4]=garch_fit@fit$coef[1:4]
      aux_list[5:9]= garch_fit@fit$coef[5:9]
    }
  else
    if(ordem[1]==0 && ordem[2]==1){  #...
      aux_list[1:2]=0
      aux_list[4]=0
      aux_list[3]=garch_fit@fit$coef[1]
      aux_list[5:9]= garch_fit@fit$coef[2:6]
    }
  else
    if(ordem[1]==0 && ordem[2]==2){
      aux_list[1:2]=0
      aux_list[3:4]=garch_fit@fit$coef[1:2]
      aux_list[5:9]= garch_fit@fit$coef[3:7]
    }
  else
    if(ordem[1]==1 && ordem[2] ==2){
      aux_list[2]=0
      aux_list[1]=garch_fit@fit$coef[1]
      aux_list[3:4]=garch_fit@fit$coef[2:3]
      aux_list[5:9]= garch_fit@fit$coef[4:8]
    } else
      if(ordem[1]==1 && ordem[2]==1){
        aux_list[2]=aux_list[4]=0
        aux_list[1]=garch_fit@fit$coef[1]
        aux_list[3]=garch_fit@fit$coef[2]
        aux_list[5:9]=garch_fit@fit$coef[3:7]
      } 
  return(aux_list)
}


#######ARMA GARCH ESTIMATION FUNCTION
arma_garch_fit_store <- function(returns, per_est){ #per_est = 1260
  residuos<- vector('list', 100)     #a list of list to store residuals, 
  coeficientes <- vector('list', 100) #arma-garch coeficients and sigma from each asset estimation
  sigma <- vector('list', 100)
  for (i in 1:100){  ##for period 1:100, we estimate ARMA-GARCH for each asset
    ret_data <- returns[i:per_est+(i-1),] #goes from window_sup-per_est to window_sup  (per_est=1260)
    for (j in 1:ncol(returns))   #for period i, fit the model for every j>1 asset, where j=1 is Date column. 
    {
      ordem <- arima_est(ret_data[,j]) ##ARMA estimation
     # print(ordem)
      garch_fit <- garch_est(ret_data[,j], ordem) ##Garch estimation
     # print(garch_fit@fit$coef)
      ###storing data from the estimations
      aux_list <- coefs_vector(ordem, garch_fit)
    # print(aux_list)
      coeficientes[[i]][[j]] <- aux_list  ##storing coeficients vector resulting for the j asset
      residuos[[i]][[j]] <- garch_fit@fit$residuals ##storing residuals vector for the j asset
      sigma[[i]][[j]] <- garch_fit@fit$sigma ##storing sigma vector for the j asset
      cat("\nFit Col/Per: ",j,i)
    }
  }
  
  ############storing resulting data in R data files
  saveRDS(coeficientes, file = "coef.Rds")
  saveRDS(residuos, file = "residuos.Rds")
  saveRDS(sigma, file = "sigma.Rds")
}


#function for transforming observations to pseudo-uniform [0,1] data needed to copula estimation
gen_emp_dist <- function(){ 
  unif <- vector('list', 100 ) ##initializing an empty list to save data
  residuos <- readRDS("residuos.Rds") ##reading our data
  sigma <- readRDS("sigma.Rds")
  garch_coef <- readRDS("coef.Rds") 
  
  for(i in 1:length(residuos)){          
    for(j in 1:8){                       
      unif[[i]][[j]] = copula::pobs(residuos[[i]][[j]]/sigma[[i]][[j]])
    }
  }
  saveRDS(unif, file = "Unif_ParDist.Rds") ##saving cdf data in R data file
}
  

#Function that reads our pseudo-uniform transformed data and estimate mixture copula parameters using our SQP approach
cop_mix_fit <- function(){
  par_cdf <- readRDS("unif_ParDist.RDS") #reading our parametric [0,1] cdf
  weight_t <- vector('list', length(par_cdf))  #creating a list to store each period's parameters

  #initializing 8d t copula
  copt <- copula::tCopula(param = 0.5, dim = 8) # 
  
  ## Initializing archimedian copula objects
  copC <- copula::claytonCopula(2, dim = 8) # delta= 2
  copG <- copula::gumbelCopula(2, dim = 8)  # theta= 2
  
  lower <- c(0.1, 1, -0.9,(2+.Machine$double.eps), 0,0,0)     #lower and upper bounds of the parameters and weights for bounded non linear opt.
  upper <- c(copC@param.upbnd, copG@param.upbnd, 1,100, 1,1,1) #2+eps so that variance of t copula is defined
  
  
  for(i in 1:length(par_cdf)){  ####for each (1:100) optimization, we estimate copula parameters for data
    v<-as.matrix(do.call(cbind, par_cdf[[i]])) ##transforming in Matrix 
    U<-v[,1:8]  ##pseudo-uniform [0,1] observations for each asset
    
    ##Creating elliptical copula objects and estimating "initial guesses" for each copula parameter. 
    #Then, we maximize loglikelihood of the linear combination of the three copulas
    par1 <- copula::fitCopula(copC, U, "itau", estimate.variance = TRUE)@estimate #inversion of Kendall's tau for Clayton 
    par2 <- copula::fitCopula(copG, U,"itau", estimate.variance = TRUE)@estimate #inversion of Kendall's tau for Gumbel 
    par3 <- copula::fitCopula(copt, U,"mpl", estimate.variance = FALSE)@estimate ###mpl para poder estimar tambem DF. Na documentacao diz que nao pode usar 'itau' pois ele n estima DF.
    par4 <- 1/3 #initial guesses for weights = 1/3 each
    par5 <- 1/3
    par6 <- 1/3
    
    ##non linear constrained optimization 
    opt <- Rsolnp::solnp(pars = c(par1,par2,par3,par4,par5,par6), fun = LLCG, LB = lower, 
                         UB = upper, copt=copt,copC = copC, copG = copG, U=U,eqfun = eqfun, 
                         eqB=c(1)) ####RSOLNP
    
    weight_t[[i]]<-opt$pars  ##saving optimization parameters in a list
    print(i)
  }
  saveRDS(weight_t, file = "copulaParams.Rds")  #saving resulting the parameters in a file
}


#From the estimated copula parameters, we generate 'nsim' 8-d mixture copula variates as a Monte Carlo way
#to approximate parameters.
mix_cop_sim <- function(cop_pars, i, nsim){
  Cc <- Cg <- Ct <- matrix(0,nrow = nsim, ncol = 8)   #clayton, t, gumbel variates matrix
  ctg <- matrix(0, nrow = nsim, ncol = 8)              #copula mixture variates matrix
  ##generating copula variates
  Cc[,]<- cop_pars[[i]][[5]]*copula::rCopula(n = nsim, copula = claytonCopula(param = cop_pars[[i]][[1]], dim = 8))   
  Cg[,]<- cop_pars[[i]][[6]]*copula::rCopula(n = nsim, copula = gumbelCopula(param = cop_pars[[i]][[2]], dim = 8))
  Ct[,]<- cop_pars[[i]][[7]]*copula::rCopula(n = nsim, copula = tCopula(param = cop_pars[[i]][[3]], 
                                                                        df = cop_pars[[i]][[4]], dim = 8))
  ctg <- Cc + Ct + Cg #linear combination of them 
  return(ctg)
}



z_cop_sim <- function(garch_coefs, ctg, i, nsim){
  zsim <- matrix(0, nrow = nsim, ncol = 8) #'z' copula matrix, resulting from quantile of mixture 
  for(j in 1:8){  #for each asset, generate copula 'z' dependence structure, applying the Quantile function 
    #in the mixture of copulas
    # if(j==8 && anyNA(fGarch::qsstd(ctg[,j], nu = garch_coefs[[i]][[j]][[8]], xi =garch_coefs[[i]][[j]][[7]]))){
    #   garch_coefs[[i]][[j]][[7]] <- 1
    # }
    nu <- garch_coefs[[i]][[j]][[8]]
    xi <- garch_coefs[[i]][[j]][[7]]
    
    zsim[,j] <- fGarch::qsstd(ctg[,j], nu = nu, xi = xi) / # garch_coefs[[i]][[j]][[7]]
      sd(fGarch::qsstd(ctg[,j], nu = nu, xi = xi)) #nu = t's DF, xi = t's skew
  }
  return(zsim)
}

##simulating returns (monte carlo) using copula dependence and ARMA-GARCH marginals
ret_cop_sim <- function(garch_coefs, zsim, sigma_per, resid_per, returns, i, nsim){
  ret_sim <- matrix(0, nrow = nsim, ncol = 8)    #simulated returns matrix
  RZ <- as.matrix(returns[i:(1259+i),1:8]) #getting real returns we'll use in one-step forward armaGarch's AR forecasting 
  for(j in 1:8){   #K scenarios generation for each asset
    sigma_f_t1 <- tail(sigma_per[[j]],1) ##(t-1) sigma for GARCH forecasting
    e_f_t2_t1 <- tail(resid_per[[j]],2) ##(t-2, t-1) residuals for MA forecasting
    for(z in 1:nsim){
      ret_sim[z,j] <- ((garch_coefs[[i]][[j]][[1]] * RZ[1260,j]) + (garch_coefs[[i]][[j]][[2]] * RZ[1259,j]) +  #AR1 * R_t-1, AR2 * R_t-2
                        (garch_coefs[[i]][[j]][[3]] * e_f_t2_t1[1]) + (garch_coefs[[i]][[j]][[4]] * e_f_t2_t1[2]) +  #MA1*e_t-1, MA2*e_t-2
                        (zsim[z,j] * (sqrt(garch_coefs[[i]][[j]][[9]]) +  ##alfa0 (omega)
                                        sqrt(garch_coefs[[i]][[j]][[5]]) * e_f_t2_t1[2] + ##alfa1 * e_t-1
                                        sqrt(garch_coefs[[i]][[j]][[6]]) * sigma_f_t1))) ##beta1  * s_t-1
    }
  }
  return(ret_sim)
}


###Portfolio Optimization Function: for each optimization step we minimize portfolio CVaR subject to a portfolio expected return objective
##The returns used in the optimization are the Monte Carlo estimated returns for one-day ahead forecast (see Pfaff 2012 or Rockafellar and Uryasev 2002)
cop_portf_opt <- function(targetReturn, filename, nsim, type, returns){
  returns <- returns
  if(type == 'mixture'){
    cop_pars <- readRDS("copulaParams.Rds")
  }
  else if(type == 'gaussian'){
    cop_pars <- readRDS("GausscopulaParams.Rds")
  }
   #reading copula parameters
  garch_coefs <- readRDS("coef.Rds") #reading ArmaGarch parameters
  sigma_fit <- readRDS("sigma.Rds") #reading armaGarch sigma. We need this to estimate J one-steap ahead returns 
  residual_fit <- readRDS("residuos.Rds")  #reading armaGarch fitted residuals 
  
  #####setting up a fGarch min CVaR optimization
  frontierSpec <- portfolioSpec() 
  setType(frontierSpec) <- "CVaR"  
  setSolver(frontierSpec) <- "solveRglpk.CVAR"  #solving as a linear program as Rockafellar & Uryasev (2000)
  setAlpha(frontierSpec) <- 0.05   #CVaR_0.95  
  setTargetReturn(frontierSpec) <- targetReturn  #daily target return constrain
  
  ####initializing returns and weights matrixes
  ret_sim <- matrix(0, nrow = nsim, ncol = 8) #initializaing sim.ret matrix
  cvar_opt <- matrix(0, nrow = nsim, ncol = 8) #initializing optimized portfolio weights matrix
  
  for(i in 1:length(sigma_fit)){  #we do everything for each optimization
    ##generating copula variates 
    if(type == 'mixture'){         #if we want mixture of copula, type should be 'mixture'
      cop_sim <- mix_cop_sim(cop_pars, i, nsim)
    #  print(head(cop_sim))
    } else if(type == 'gaussian'){ #else, 'gaussian'
      cop_sim <- gauss_cop_sim(cop_pars, i, nsim)
    }
    ##generating zsim 
    zsim <- z_cop_sim(garch_coefs, cop_sim, i, nsim)
   # print(head(zsim))
    #sigma_per <- sigma_fit[[i]]  #getting fitted sigma and residuals we'll use in MA and Garch's forecasting
    #resid_per <- residual_fit[[i]] 
    
    ##K return scenarios generation, using z and armaGarch
    ret_sim <- ret_cop_sim(garch_coefs, zsim, sigma_fit[[i]], residual_fit[[i]], returns, i, nsim)
   # print(head(ret_sim))
    ##optimizing portfolio using K simulated return for each assets, for optimization period i 
    retornofPort <- as.timeSeries(ret_sim[, 1:8])
    frontier1g <- fPortfolio::efficientPortfolio(data = retornofPort, spec = frontierSpec, constraints = "LongOnly")
    cvar_opt[i,1:8] <- fPortfolio::getWeights(frontier1g)   #storing resulting weights   
  #  print(cvar_opt)
  }
  saveRDS(cvar_opt, file = filename) ##saving weights data, we repeat this for 0,3,6 and 9%
}




#################################################################
#Now that our functions are defined, we will run the optimization
library(xts)
library(quantmod)
library(forecast)
library(rugarch)
library(Rsolnp)
library(fPortfolio)
library(fGarch)
library(copula)

#Downloading the 8 U.S. stocks 
getSymbols(Symbols = c("AAPL",'TSLA','BAC', 'MSFT', 'WFC', "T", "XOM", 'INTC'), 
           env = parent.frame(),
           reload.Symbols = FALSE,
           verbose = FALSE,
           warnings = TRUE,
           src = "yahoo",
           symbol.lookup = TRUE,
           auto.assign = getOption('getSymbols.auto.assign',T),
           from = '2015-01-01', to = '2020-06-02') 

stockprice <- Cl(merge.xts(AAPL, BAC, INTC, MSFT, T, TSLA, WFC, XOM))

#calculating LOG returns from prices
returns <- Return.calculate(stockprice, method = 'log') #1362 no total
returns <- returns[-1,] #dropping first line 

AAPL <- BAC <- INTC <- MSFT <- TSLA <- WFC <- XOM <- NULL

per_est <- 1260 #using 1260 days for model estimation 

##fitting arma-garch models
arma_garch_fit_store(returns[1:1360,], per_est) #ok

##generating empirical cdf
gen_emp_dist() #ok

##estimating copula parameters
cop_mix_fit() #ok

##generating MC simulated returns using copula-arma-garch and optimizing portfolios weights 
cop_portf_opt(targetReturn = 0.00072, filename = "mix_cop_wgt_18pct.Rds",
              nsim= 10000, type = 'mixture', returns = returns)
#we set targetReturn to 0.00072 daily and run 10k MC estimations



###############################################################################
#Now, lets run a benchmark comparing our portfolio to an Equal Weight portfolio

#loading our estimated weights
estw <- readRDS('mix_cop_wgt_12pct.Rds')


library(PerformanceAnalytics) 
estw <- estw[1:100,]




#transforming back from log-returns to normal returns
bench_ret <- exp(returns[1261:1360,]) - 1

#calculating forecasted returns from estimated weights and real returns
out_sam_ret <- estw * bench_ret

#calculating the returns of the portfolio for mixture and for 1/n
out_sam_p_ret <- PerformanceAnalytics::Return.portfolio(out_sam_ret, geometric = TRUE, wealth.index = FALSE)
charts.PerformanceSummary(out_sam_p_ret) ##wealth + dd plot

out_sam_1n_ret <- PerformanceAnalytics::Return.portfolio(bench_ret, geometric = TRUE, wealth.index = FALSE)
charts.PerformanceSummary(out_sam_1n_ret)

retursnmatrix <- NULL
returnsmatrix <- cbind(out_sam_p_ret, out_sam_1n_ret)
charts.PerformanceSummary(returnsmatrix) #####plotting performance summary


#####calculating risk statistics (we didnt use this)

bench <- rbind(Return.annualized(returnsmatrix), sd.annualized(returnsmatrix),
               VaR(returnsmatrix, method = "historical"), CVaR(returnsmatrix, method = "historical"), 
               SemiDeviation(returnsmatrix), CDD(returnsmatrix, method = "historical"), 
               maxDrawdown(returnsmatrix, method = "historical"),
               AverageDrawdown(returnsmatrix, method = "historical"), AverageLength(returnsmatrix, method = "historical"),
               SharpeRatio.annualized(returnsmatrix, Rf = 0), 
               BurkeRatio(returnsmatrix, Rf = 0),
               SortinoRatio(returnsmatrix, MAR = 0), UpsidePotentialRatio(returnsmatrix, MAR = 0), 
               DownsideFrequency(returnsmatrix, MAR = 0),
               CalmarRatio(returnsmatrix),
               DrawdownDeviation(returnsmatrix),
               Kappa(returnsmatrix,MAR=0, 1), 
               OmegaSharpeRatio(returnsmatrix, MAR=0))
