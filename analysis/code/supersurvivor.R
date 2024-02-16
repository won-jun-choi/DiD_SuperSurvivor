# Meta =====================
# Title: supersurvivor
# Author: Sanghee
# Description: haha

# Last Edit/Editor: Feb-15-2024 / Wonjun
# Description: Re-write the script as a function.

# Load packages
library(cuRe)

my_SuperSurvivor <- function(data,
                             unit_variable,
                             time_variable,
                             duration_variable,
                             uncensored_indicator,
                             logit_regressors,
                             survival_regressors,
                             survival_function_type){
  # unit_variable: str; variable name for unit (in ours, i)
  # time_variable: str; variable name for time (in ours, t)
  # duration_variable: str; variable name for duration (in ours, G)
  # logit_regressors: list(str); list of variable names for logit (x)
  # survival_regressors: list(str); list of variable names for survival (z)
  # survival_function_type: 'Weibull', 'LogNormal'
  t_max <- data %>% pull({{time_variable}}) %>% max()
  
  if (survival_function_type == 'Weibull') {
    # generate time_to_treat and status for cuRe package
    data['time_to_treat'] <- data %>% pull({{duration_variable}})
    data['status'] <- ifelse(data %>% pull({{uncensored_indicator}}) == 0,
                             1,
                             ifelse(data %>% 
                                      pull({{uncensored_indicator}}) == 1, 
                                    0, NA))
    # From c('x1','x2') how can I make x1+x2 and define it as formula
    logit_formula <- as.formula(paste("Surv(time_to_treat,status) ~", paste(survival_regressors,collapse = '+')))
    surv_formula <- as.formula(paste(paste(logit_regressors,collapse = '+'),"+1)"))
    super_survivor_logit <- fit.cure.model(logit_formula,
                                           formula.surv = list(~surv_formula),
                                           data = data,
                                           type = 'mixture',
                                           bhazard = NULL,
                                           dist = 'weibull',
                                           link = 'logit')
    gammahat <- super_survivor_logit$coefs[['1']]
    # generate X which is a matrix of 1 and logit_regressors
    X <- cbind(1, data %>% pull({{logit_regressors}}))
    phat <- exp(gammahat %*% X) / (1 + exp(gammahat %*% X))
  } else if (survival_function_type == 'LogNormal') {
    # Do MLE using the code in log_linear_dgp
  }
  return(phat)
}

my_SuperSurvivor(data = df,
                 unit_variable = i,
                 time_variable = t,
                 duration_variable = G,
                 uncensored_indicator = C,
                 logit_regressors = c('x1','x2','x3','x4'),
                 survival_regressors = c('z1','z2','z3'),
                 survival_function_type = 'Weibull')

as.formula(paste(c('z1','z2','z3'),collapse = '+')) # z1+z2+z3

df$time_to_treat<-df$G
df$time_to_treat[df$time_to_treat > t_max] <- t_max # this is unnecessary

#generate status
df$status <- ifelse(df$C == 0, 1, ifelse(df$C == 1, 0, NA))

super_survivor_logit <- fit.cure.model(Surv(time_to_treat,status) ~ z1+z2, formula.surv=list(~x1+x2+1), data=df, type='mixture',
                                       bhazard=NULL,dist='weibull', link='logit')

summary(super_survivor_logit)  # pi: logit, theta: weibull
gammahat <- super_survivor_logit$coefs[['1']]

###### pb of being supersurvivor
df$phat<-exp(gammahat[1]+gammahat[2]*df$x1+gammahat[3]*df$x2)/(1+exp(gammahat[1]+gammahat[2]*df$x1+gammahat[3]*df$x2))

rm(list = c("gammahat", "super_survivor_logit", "t_max"))

