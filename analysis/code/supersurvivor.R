# Meta =====================
# Title: supersurvivor
# Author: Sanghee
# Description: haha

# Last Edit/Editor: Feb-15-2024 / Wonjun
# Description: Re-write the script as a function.

# Load packages
my_SuperSurvivor <- function(data,
                             unit_variable,
                             time_variable,
                             duration_variable,
                             censored_indicator,
                             logit_regressors,
                             survival_regressors,
                             survival_function_type,
                             optimization_method='Nelder-Mead',
                             return_params=FALSE){
  # unit_variable: str; variable name for unit (in ours, i)
  # time_variable: str; variable name for time (in ours, t)
  # duration_variable: str; variable name for duration (in ours, G)
  # logit_regressors: list(str); list of variable names for logit (x)
  # survival_regressors: list(str); list of variable names for survival (z)
  # survival_function_type: 'Weibull', 'LogNormal'
  
  # t_max <- data %>% pull({{time_variable}}) %>% max()
  ones <- rep(1, nrow(data))
  X <- data[,logit_regressors] %>% as.matrix()
  Z <- data[,survival_regressors] %>% as.matrix()
  G <- data %>% pull({{duration_variable}}) %>% as.numeric()
  C <- data %>% pull({{censored_indicator}}) %>% as.numeric()
  
  if (survival_function_type == 'Weibull') {
    # hahaha
  } else if (survival_function_type == 'LogNormal') {
    # log-likelihood function for cured model
    # L(\theta;X,Z) = \Prod_{i=1}^n [p_uncured * f(dur|Z)]^{uncensored} * \Prod_{i=1}^n [p_cured + p_uncured*S(dur|Z)]^{censored}
    # Thus, l(\tehta;X,Z) = \Sum_{i=1}^n uncensored * [log(p_uncured) + log(f)] + \Sum_{i=1}^n censored[log(p_cured + p_uncured*S)]
    loglik <- function(theta) {
      # ones <- rep(1, nrow(data))
      
      gamma <- theta[1:(length(logit_regressors)+1)]
      # X <- data[,logit_regressors] %>% as.matrix()
      Xg <- cbind(ones,X) %*% t(t(gamma))
      p_cured <- exp(Xg)/(1+exp(Xg))
      p_uncured <- 1 - p_cured
      
      lambda <- theta[(length(logit_regressors)+2):(length(theta)-1)]
      s <- theta[length(theta)]
      # Z <- data[,survival_regressors] %>% as.matrix()
      Zl <- cbind(ones,Z) %*% t(t(lambda))
      f <- dnorm((log(G)-Zl)/s)/(G*s) # density
      S <- 1 - pnorm((log(G)-Zl)/s)
      
      # C <- data %>% pull(censored_indicator) %>% as.numeric()
      ret <- (1-C) * (log(p_uncured) + log(f)) + C * (log(p_cured + p_uncured*S))
      return(-sum(ret)) # for minimization
      }
    } else if (survival_function_type == 'LogNormal_discrete') {
    # In discrete duration case, the density can be replaced by a probability mass
    
    # log-likelihood function
    # L(\theta;X,Z) = \Prod_{i=1}^n [p_uncured * f(dur|Z)]^{uncensored} * \Prod_{i=1}^n [p_cured + p_uncured*S(dur|Z)]^{censored}
    # Thus, l(\tehta;X,Z) = \Sum_{i=1}^n uncensored * [log(p_uncured) + log(f)] + \Sum_{i=1}^n censored[log(p_cured + p_uncured*S)]
    loglik <- function(theta) {
      # ones <- rep(1, nrow(data))
      
      gamma <- theta[1:(length(logit_regressors)+1)]
      # X <- data[,logit_regressors] %>% as.matrix()
      Xg <- cbind(ones,X) %*% t(t(gamma))
      p_cured <- exp(Xg)/(1+exp(Xg))
      p_uncured <- 1 - p_cured
      
      lambda <- theta[(length(logit_regressors)+2):(length(theta)-1)]
      s <- theta[length(theta)]
      # Z <- data[,survival_regressors] %>% as.matrix()
      Zl <- cbind(ones,Z) %*% t(t(lambda))
      f <- pnorm((log(G+1)-Zl)/s) - pnorm((log(G)-Zl)/s) # point mass instead of density
      S <- 1 - pnorm((log(G)-Zl)/s)
      
      # C <- data %>% pull(censored_indicator) %>% as.numeric()
      ret <- (1-C) * (log(p_uncured) + log(f)) + C * (log(p_cured + p_uncured*S))
      return(-sum(ret)) # for minimization
    }
  } else {
    stop('survival_function_type should be one of "Weibull", "LogNormal", "LogNormal_discrete"')
  }
  
  # optimization
  initial_value = c(-1,2,0,0,0,-1,2,0,0,2)
  # initial_value = c(-1,2,-2,2,3, -1,2,-3,4, 1)  # sigma_V
  res <- optim(par = initial_value,
                 fn = loglik,
                 method = optimization_method)
  gamma_hat <- res$par[1:(length(logit_regressors)+1)]
  lambda_hat <- res$par[(length(logit_regressors)+2):(length(res$par)-1)]
  # print('gamma_hat: ')
  print(gamma_hat)
  # print('lambda_hat: ')
  print(lambda_hat)
  num <- exp(cbind(ones,X) %*% t(t(gamma_hat)))
  denom <- 1+exp(cbind(ones,X) %*% t(t(gamma_hat)))
  phat <- num/denom
  
  if (return_params == TRUE) {
    return(c(phat, gamma_hat, lambda_hat))
  }
  return(phat)
}

if (sys.nframe() == 0) {
  res = my_SuperSurvivor(data = df,
                         unit_variable = 'i',
                         time_variable = 't',
                         duration_variable = 'G',
                         censored_indicator = 'C',
                         logit_regressors = c('x1','x2','x3','x4'),
                         survival_regressors = c('z1','z2','z3'),
                         survival_function_type = 'LogNormal_discrete',
                         optimization_method='Nelder-Mead')
}

# 
# df$phat <- res
# df %>% select(phat,p_cured) %>% View()

# input
# data <- df
# unit_variable = 'i'
# time_variable = 't'
# duration_variable = 'G'
# censored_indicator = 'C'
# logit_regressors = c('x1','x2','x3','x4')
# survival_regressors = c('z1','z2','z3')
# survival_function_type = 'LogNormal_discrete'
# theta = c(1,0,0,0,0,1,0,0,0,3) # gamma (logit), lambda (survival), s_v (survival)

# # loglikelihood function
# ones <- rep(1, nrow(data))
# 
# gamma <- theta[1:(length(logit_regressors)+1)]
# X <- data[,logit_regressors] %>% as.matrix()
# Xg <- cbind(ones,X) %*% t(t(gamma))
# p_cured <- exp(Xg)/(1+exp(Xg))
# p_uncured <- 1 - p_cured
# 
# lambda <- theta[(length(logit_regressors)+2):(length(theta)-1)]
# s <- theta[length(theta)]
# Z <- data[,survival_regressors] %>% as.matrix()
# Zl <- cbind(ones,Z) %*% t(t(lambda))
# f <- pnorm((log(G+1)-Zl)/s) - pnorm((log(G)-Zl)/s) # point mass instead of density
# S <- 1 - pnorm((log(G)-Zl)/s)
# 
# C <- data %>% pull(censored_indicator) %>% as.numeric()
# 
# ret <- (1-C) * (log(p_uncured) + log(f)) + C * (log(p_cured + p_uncured*S))
# sum(ret, rm.na=T)


######## old code ########
# df$time_to_treat<-df$G
# df$time_to_treat[df$time_to_treat > t_max] <- t_max # this is unnecessary
# 
# #generate status
# df$status <- ifelse(df$C == 0, 1, ifelse(df$C == 1, 0, NA))
# 
# super_survivor_logit <- fit.cure.model(Surv(time_to_treat,status) ~ z1+z2, formula.surv=list(~x1+x2+1), data=df, type='mixture',
#                                        bhazard=NULL,dist='weibull', link='logit')
# 
# summary(super_survivor_logit)  # pi: logit, theta: weibull
# gammahat <- super_survivor_logit$coefs[['1']]
# 
# ###### pb of being supersurvivor
# df$phat<-exp(gammahat[1]+gammahat[2]*df$x1+gammahat[3]*df$x2)/(1+exp(gammahat[1]+gammahat[2]*df$x1+gammahat[3]*df$x2))
# 
# rm(list = c("gammahat", "super_survivor_logit", "t_max"))
# 
