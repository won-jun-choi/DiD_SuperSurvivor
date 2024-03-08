# Meta =====================
# Title: DiD using supersurvivor
# Author: Wonjun
# Last Edit: Dec-03-2023
# Description: haha
if (sys.nframe() == 0) {
 rm(list = ls())
}
    
DiD_infeasible <- function(data,
                           unit_variable,
                           time_variable,
                           group_variable,
                           supersurvivor_indicator,
                           outcome_variable,
                           outcome_regressors) {
  # generate dataframe for DD
  df_DD <- data
  df_DD['i_'] <- df_DD[unit_variable]
  df_DD['t_'] <- df_DD[time_variable]
  df_DD['G_'] <- df_DD[group_variable]
  df_DD['C_tilde_'] <- df_DD[supersurvivor_indicator]
  
  t_min <- df_DD['t_'] %>% min()
  t_max <- df_DD['t_'] %>% max()  # we can get ATT from (t_min+1)~t_max
  G_max <- df_DD[df_DD['G_']<1000, 'G_'] %>% max()
  
  # dataframe that contains all possible ATTs
  results <- data.frame(G = rep(c(1:G_max), each=t_max),
                        t = rep(c(1:t_max), times=G_max))
  results <- results %>% filter(G > t_min, t >= G)
  results['ATT_infeasible'] <- NA
  
  # DD using C_tilde as control
  for (i in 1:nrow(results)) {
    # ATT(g__,t__)
    g__ <- results[i, group_variable]
    t__ <- results[i, time_variable]
    T_group <- df_DD %>% filter(G_==g__)
    C_group <- df_DD %>% filter(C_tilde_==1)
    
    # T group and C group
    T_group <- T_group %>%
      mutate(Yt = ifelse(t_==t__,Y,NA),
             Ygm1 = ifelse(t_==g__-1,Y,NA)) %>%
      group_by(i_) %>%
      mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
      ungroup() %>%
      distinct(i_, .keep_all=TRUE)
    C_group <- C_group %>%
      mutate(Yt = ifelse(t_==t__,Y,NA),
             Ygm1 = ifelse(t_==g__-1,Y,NA)) %>%  # gm1: g minus 1
      group_by(i_) %>%
      mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
      ungroup() %>%
      distinct(i_, .keep_all=TRUE)
    
    # ATT(g,t)
    Em1 <- mean(T_group$dY)
    formula <- paste0('dY ~ ', paste0(outcome_regressors, collapse = '+'))
    m0 <- lm(formula, data=C_group)
    Em0 <- predict(m0, newdata=T_group) %>% mean()
    results[i,'ATT_infeasible'] <- Em1 - Em0
  }
  return(results)
}

if (sys.nframe() == 0) {
  # test DGP for checking code validity
  source('analysis/code/testDGP.R')
  test_data <- read_csv('analysis/temp/testDGP.csv')
  res <- DiD_infeasible(data=test_data,
                        unit_variable='unit',
                        time_variable='t',
                        group_variable='G',
                        outcome_variable='Y',
                        outcome_regressors = c('x1','x2','x4'),
                        supersurvivor_indicator='C_tilde')
  res
}

################################################################################
# data <- read_csv('analysis/temp/testDGP.csv')
# unit_variable <- 'unit'
# time_variable <- 't'
# group_variable <- 'G'
# supersurvivor_indicator <- 'C_tilde'
# outcome_regressors <- c('x1', 'x2', 'x4')
# outcome_regressors <- '1'
# 
# formula <- paste0('dY ~ ', paste0(outcome_regressors, collapse = '+'))
# formula
