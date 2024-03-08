if (sys.nframe() == 0){
  rm(list = ls())
}
#'
#'
#'
#'
DiD_supersurvivor <- function(data,
                              unit_variable,
                              time_variable,
                              group_variable,
                              outcome_variable,
                              outcome_regressors,
                              censored_indicator,
                              logit_regressors,
                              survival_regressors,
                              survival_function_type,
                              survival_MLE_optimizer,
                              MLE_init_values=NULL
) {
  # generate dataframe for DD =====================
  df_DD <- data
  df_DD['t_'] <- df_DD[time_variable]
  df_DD['i_'] <- df_DD[unit_variable]
  df_DD['G_'] <- df_DD[group_variable]
  df_DD['C_'] <- df_DD[censored_indicator]
  df_DD['Y_'] <- df_DD[outcome_variable]
  
  # Run MLE for the cured model to get phat =====================
  source('analysis/code/supersurvivor.R')
  df_DD$phat <- my_SuperSurvivor(data = df_DD,
                                 unit_variable = unit_variable,
                                 time_variable = time_variable,
                                 duration_variable = group_variable,
                                 censored_indicator = censored_indicator,
                                 logit_regressors = logit_regressors,
                                 survival_regressors = survival_regressors,
                                 survival_function_type = survival_function_type,
                                 optimization_method = survival_MLE_optimizer,
                                 MLE_init_values = MLE_init_values,
                                 return_params = FALSE)
  # Doesn't support control variables for DD yet.
  
  # Estimate DiD using reweighting =====================
  t_min <- df_DD['t_'] %>% min()
  t_max <- df_DD['t_'] %>% max()  # we can get ATT from (t_min+1)~t_max
  G_max <- df_DD[df_DD['G_']<9000, 'G_'] %>% max()
  results <- data.frame(G = rep(c(1:G_max), each=t_max),
                        t = rep(c(1:t_max), times=G_max))
  results <- results %>% filter(G > t_min, t >= G)
  results$ATT_SS <- NA
  
  # # DiD estimator
  for (i in 1:nrow(results)) {
    # ATT(g__,t__)
    g__ <- results[i, group_variable]
    t__ <- results[i, time_variable]
    # Select T group and generate dY
    T_group <- df_DD %>% filter(G_==g__)
    T_group <- T_group %>%
      mutate(Yt = ifelse(t_==t__,Y_,NA),
             Ygm1 = ifelse(t_==g__-1,Y_,NA)) %>%
      group_by(i_) %>%
      mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
      ungroup() %>%
      distinct(i_, .keep_all=TRUE)
    # Select C group and generate dY
    C_group <- df_DD %>% filter(C_==1)
    C_group <- C_group %>%
      mutate(Yt = ifelse(t_==t__,Y_,NA),
             Ygm1 = ifelse(t_==g__-1,Y_,NA)) %>%  # gm1: g minus 1
      group_by(i_) %>%
      mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>% # just to remove NA
      ungroup() %>%
      distinct(i_, .keep_all=TRUE)
    
    # re-weighted DiD
    Em1 <- T_group$dY %>% mean()
    formula <- paste0('dY ~ ', paste0(outcome_regressors, collapse='+'))
    m0_reweight <- lm(formula, data=C_group, weights=phat)
    # Em0_reweight <- predict(m0_reweight, newdata=C_group) %>% mean()
    Em0_reweight <- predict(m0_reweight, newdata=T_group) %>% mean()
    results[i, 'ATT_SS'] <- Em1 - Em0_reweight
  }  # end of for loop of DiD estimators
  return(results)
}

if (sys.nframe() == 0) {
  # test DGP to check code validity
  source('analysis/code/testDGP.R')
  df <- read_csv('analysis/temp/testDGP.csv')
  res <- DiD_supersurvivor(data = df,
                           unit_variable = 'unit',
                           time_variable = 't',
                           group_variable = 'G',
                           censored_indicator = 'C',
                           outcome_variable= 'Y',
                           outcome_regressors = c('x1','x2','x4'),
                           logit_regressors = c('x1', 'x2', 'x3'),
                           survival_regressors = c('z1', 'z2', 'z3'),
                           survival_function_type = 'LogNormal_discrete',
                           survival_MLE_optimizer = 'Nelder-Mead',
                           MLE_init_values = c(-1,2,0,0,-1,2,0,0,2))
  res
}

################################################################################
# input
# rm(list =ls())
# data <- read_csv('analysis/temp/testDGP.csv')
# unit_variable <- 'unit'
# time_variable <- 't'
# group_variable <- 'G'
# outcome_variable <- 'Y'
# outcome_regressors <- c('x1', 'x2', 'x4')
# censored_indicator <- 'C'
# logit_regressors <- c('z1', 'z2', 'z3')
# survival_regressors <- c('x1', 'x2', 'x3')
# survival_function_type <- 'LogNormal_discrete'
# survival_MLE_optimizer <- 'Nelder-Mead'
# 
# # code
# df_DD <- data
# df_DD['i_'] <- df_DD[unit_variable]
# df_DD['t_'] <- df_DD[time_variable]
# df_DD['G_'] <- df_DD[group_variable]
# df_DD['C_'] <- df_DD[censored_indicator]
# df_DD['Y_'] <- df_DD[outcome_variable]
# 
# source('analysis/code/supersurvivor.R')
# df_DD$phat <- my_SuperSurvivor(data = df_DD,
#                                unit_variable = unit_variable,
#                                time_variable = time_variable,
#                                duration_variable = group_variable,
#                                censored_indicator = censored_indicator,
#                                logit_regressors = logit_regressors,
#                                survival_regressors = survival_regressors,
#                                survival_function_type = survival_function_type,
#                                optimization_method = survival_MLE_optimizer,
#                                return_params = FALSE)
# 
# t_min <- df_DD['t_'] %>% min()
# t_max <- df_DD['t_'] %>% max()  # we can get ATT from (t_min+1)~t_max
# G_max <- df_DD[df_DD['G_']<9000, 'G_'] %>% max()
# results <- data.frame(G = rep(c(1:G_max), each=t_max),
#                       t = rep(c(1:t_max), times=G_max))
# results <- results %>% filter(G > t_min, t >= G)
# results$ATT_SS_noX <- NA
# 
# # DiD
# g__ <- 2
# t__ <- 2
# 
# T_group <- df_DD %>% filter(G_==g__)
# T_group <- T_group %>%
#   mutate(Yt = ifelse(t_==t__,Y_,NA),
#          Ygm1 = ifelse(t_==g__-1,Y_,NA)) %>%
#   group_by(i_) %>%
#   mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
#   ungroup() %>%
#   distinct(i_, .keep_all=TRUE)
# 
# C_group <- df_DD %>% filter(C_==1)
# C_group <- C_group %>%
#   mutate(Yt = ifelse(t_==t__,Y_,NA),
#          Ygm1 = ifelse(t_==g__-1,Y_,NA)) %>%  # gm1: g minus 1
#   group_by(i_) %>%
#   mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>% # just to remove NA
#   ungroup() %>%
#   distinct(i_, .keep_all=TRUE)
# 
# # reweighted DiD
# Em1 <- T_group$dY %>% mean()
# formula <- paste0('dY ~ ', paste0(outcome_regressors, collapse='+'))
# m0_reweight <- lm(formula, data=C_group, weights=phat)
# Em0_reweight <- predict(m0_reweight, newdata=T_group) %>% mean()
# Em0_reweight1 <- predict(m0_reweight, newdata=C_group) %>% mean()
# Em1 - Em0_reweight
