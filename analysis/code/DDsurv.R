#' DiD using supersurvivor
#' gagaga 
#' 
#' 
#'
DiD_supersurvivor <- function(data,
                              unit_variable,
                              time_variable,
                              group_variable,
                              censored_indicator,
                              logit_regressors,
                              survival_regressors,
                              survival_function_type,
                              survival_MLE_optimizer
                              ) {
  # generate dataframe for DD =====================
  df_DD <- df
  df_DD['t_'] <- df_DD[time_variable]
  df_DD['i_'] <- df_DD[unit_variable]
  df_DD['G_'] <- df_DD[group_variable]
  
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
                                  optimization_method = survival_MLE_optimizer)
  # Doesn't support control variables for DD yet.
  
  # Estimate DiD using reweighting =====================
  # Do ATT_RI of Callaway and Sant'Anna (2020) and two reweight version of it.
  # CS(2020) get ATT_gt of all possible g and t separately.
  # The first reweighted version is using phat as a weight in the OLS of m_hat.
  
  # dataframe that will contains all possible ATTs of CS(2020)
  t_min <- df_DD['t_'] %>% min()
  t_max <- df_DD['t_'] %>% max()  # we can get ATT from (t_min+1)~t_max
  G_max <- df_DD[df_DD['G_']<9000, 'G_'] %>% max()
  G_censored <- df_DD['G'] %>% max()
  
  results <- data.frame(G = rep(c(1:G_max), each=t_max),
                        t = rep(c(1:t_max), times=G_max))
  results <- results %>% filter(G > t_min, t >= G)
  results$ATT <- NA
  results$ATT_reweight <- NA
  
  # DiD estimators
  for (i in 1:nrow(results)) {
    g__ <- results[i, group_variable]
    t__ <- results[i, time_variable]
    T_group <- df_DD %>% filter(G_==g__)
    C_group <- df_DD %>% filter(G_==G_censored)
    
    # m(g,t|X)
    C_group <- C_group %>%
      mutate(Yt = ifelse(t_==t__,Y,NA),
             Ygm1 = ifelse(t_==g__-1,Y,NA)) %>%  # gm1: g minus 1
      group_by(i_) %>%
      mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
      ungroup() %>%
      distinct(i_, .keep_all=TRUE)
    # weighted OLS of dY on x1, x2 using phat as a weight
    formula <- paste('dY ~', paste(logit_regressors, collapse = '+'))
    m0 <- lm(formula, data=C_group)
    m0_reweight <- lm(formula, data=C_group, weights=phat)
    
    # regression imputation DD
    T_group <- T_group %>%
      mutate(Yt = ifelse(t_==t__,Y,NA),
             Ygm1 = ifelse(t_==g__-1,Y,NA)) %>%
      group_by(i_) %>%
      mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
      ungroup() %>%
      distinct(i_, .keep_all=TRUE)
    
    Em1 <- mean(T_group$dY)
    T_group$Em0 <- predict(m0, newdata=T_group)
    Em0 <- predict(m0, newdata=T_group) %>% mean()
    Em0_reweight <- predict(m0_reweight, newdata=T_group) %>% mean()
    
    # ATT(g,t)
    results[i,'ATT'] <- Em1 - Em0
    results[i,'ATT_reweight'] <- Em1 - weighted.mean(T_group$Em0, T_group$phat)
    results[i,'ATT_reweight2'] <- Em1 - Em0_reweight
  }  # end of for loop of DiD estimators
  return(results)
}

if (sys.nframe() == 0) {
  res <- DiD_supersurvivor(data = df,
                           unit_variable = 'unit',
                           time_variable = 't',
                           group_variable = 'G',
                           censored_indicator = 'C',
                           logit_regressors = c('x1', 'x2', 'x3', 'x4'),
                           survival_regressors = c('z1', 'z2', 'z3'),
                           survival_function_type = 'LogNormal_discrete',
                           survival_MLE_optimizer = 'Nelder-Mead')
}

# # potential inputs
# time_variable <- 't'
# unit_variable <- 'unit'
# group_variable <- 'G'
# 
# # generate dataframe for DD =====================
# df_DD <- df
# df_DD['t_'] <- df_DD[time_variable]
# df_DD['i_'] <- df_DD[unit_variable]
# df_DD['G_'] <- df_DD[group_variable]
# 
# t_min <- df_DD['t_'] %>% min()
# t_max <- df_DD['t_'] %>% max()  # we can get ATT from (t_min+1)~t_max
# G_max <- df_DD[df_DD['G_']<100, 'G_'] %>% max()
# G_censored <- df_DD['G'] %>% max()
# 
# # dataframe that contains all possible ATTs
# results <- data.frame(G = rep(c(1:G_max), each=t_max),
#                       t = rep(c(1:t_max), times=G_max))
# results <- results %>% filter(G > t_min, t >= G)
# results$ATT <- NA
# results$ATT_reweight <- NA
# 
# # DD =====================
# for (i in 1:nrow(results)) {
#   g__ <- results[i, group_variable]
#   t__ <- results[i, time_variable]
#   T_group <- df_DD %>% filter(G_==g__)
#   C_group <- df_DD %>% filter(G_==G_censored)
#   
#   # m(g,t|X)
#   C_group <- C_group %>%
#     mutate(Yt = ifelse(t_==t__,Y,NA),
#            Ygm1 = ifelse(t_==g__-1,Y,NA)) %>%  # gm1: g minus 1
#     group_by(i_) %>%
#     mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
#     ungroup() %>%
#     distinct(i_, .keep_all=TRUE)
#   # weighted OLS of dY on x1, x2 using phat as a weight
#   m0 <- lm('dY ~ x1+x2', data=C_group)
#   m0_reweight <- lm('dY ~ x1+x2', data=C_group, weights=phat)
#   
#   
#   # regression imputation DD
#   T_group <- T_group %>%
#     mutate(Yt = ifelse(t_==t__,Y,NA),
#            Ygm1 = ifelse(t_==g__-1,Y,NA)) %>%
#     group_by(i_) %>%
#     mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
#     ungroup() %>%
#     distinct(i_, .keep_all=TRUE)
#   
#   Em1 <- mean(T_group$dY)
#   T_group$Em0 <- predict(m0, newdata=T_group)
#   Em0 <- predict(m0, newdata=T_group) %>% mean()
#   Em0_reweight <- predict(m0_reweight, newdata=T_group) %>% mean()
#   
#   # ATT(g,t)
#   results[i,'ATT'] <- Em1 - Em0
#   results[i,'ATT_reweight'] <- Em1 - weighted.mean(T_group$Em0, T_group$phat)
#   results[i,'ATT_reweight2'] <- Em1 - Em0_reweight
# }
# 
# DDsurv <- results
# 
# rm(list=c('time_variable','unit_variable','group_variable',
#           't_min','t_max','G_max','G_censored','i','g__','t__',
#           'T_group','C_group','m0', 'm0_reweight', 'Em1','Em0', 'Em0_reweight',
#           'results','df_DD'))
