#'
#'
#'
#'
DiD_CS <- function(data,
                   unit_variable,
                   time_variable,
                   group_variable,
                   censored_indicator,
                   outcome_variable,
                   outcome_regressors
) {
  # generate dataframe for DD =====================
  df_DD <- data
  df_DD['t_'] <- df_DD[time_variable]
  df_DD['i_'] <- df_DD[unit_variable]
  df_DD['G_'] <- df_DD[group_variable]
  df_DD['C_'] <- df_DD[censored_indicator]
  df_DD['Y_'] <- df_DD[outcome_variable]
  
  ## ATT_RI of Callaway and Sant'Anna (2020)
  # dataframe that will contains all possible ATTs of CS(2020)
  t_min <- df_DD['t_'] %>% min()
  t_max <- df_DD['t_'] %>% max()  # we can get ATT from (t_min+1)~t_max
  G_max <- df_DD[df_DD['G_']<9000, 'G_'] %>% max()
  
  results <- data.frame(G = rep(c(1:G_max), each=t_max),
                        t = rep(c(1:t_max), times=G_max))
  results <- results %>% filter(G > t_min, t >= G)
  results$ATT_CS <- NA
  
  # ATT_RI
  for (i in 1:nrow(results)) {
    # ATT(g__,t__)
    g__ <- results[i, group_variable]
    t__ <- results[i, time_variable]
    # Select T group and generate dY
    T_group <- df_DD %>% filter(G_==g__)
    T_group <- T_group %>% # generate dY
      mutate(Yt = ifelse(t_==t__,Y_,NA),  # select Yt
             Ygm1 = ifelse(t_==g__-1,Y_,NA)) %>%  # select Y_{g-1}
      group_by(i_) %>%
      mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%  # generate dY
      ungroup() %>%
      distinct(i_, .keep_all=TRUE)  # keep only one dY for each i_
    
    # Select C group and generate dY
    C_group <- df_DD %>% filter(C_==1)
    C_group <- C_group %>%
      mutate(Yt = ifelse(t_==t__,Y_,NA),
             Ygm1 = ifelse(t_==g__-1,Y_,NA)) %>%
      group_by(i_) %>%
      mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
      ungroup() %>%
      distinct(i_, .keep_all=TRUE)
    
    # CS(2020)
    formula <- paste('dY ~', paste(outcome_regressors, collapse = '+'))
    m0 <- lm(formula, data=C_group)
    Em1 <- mean(T_group$dY)
    Em0 <- predict(m0, newdata=T_group) %>% mean()
    results[i,'ATT_CS'] <- Em1 - Em0
  }  # end of for loop of DiD estimators
  return(results)
}

if (sys.nframe() == 0) {
  # test DGP to check code validity
  source('analysis/code/testDGP.R')
  df <- read_csv('analysis/temp/testDGP.csv')
  # DiD_CS
  res <- DiD_CS(data = df,
                unit_variable = 'unit',
                time_variable = 't',
                group_variable = 'G',
                censored_indicator = 'C',
                outcome_variable = 'Y',
                outcome_regressors = c('x1','x2','x4')
  )
}

################################################################################
# inputs

# Select T group and generate dY
# g__ = 2
# t__ = 2
# df_DD <- df
# df_DD['G_'] <- df_DD['G']
# df_DD['t_'] <- df_DD['t']
# df_DD['Y_'] <- df_DD['Y']
# df_DD['C_'] <- df_DD['C']
# df_DD['i_'] <- df_DD['unit']
# 
# T_group <- df_DD %>% filter(G_==g__)
# T_group <- T_group %>% # generate dY
#   mutate(Yt = ifelse(t_==t__,Y_,NA),  # select Yt
#          Ygm1 = ifelse(t_==g__-1,Y_,NA)) %>%  # select Y_{g-1}
#   group_by(i_) %>%
#   mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%  # generate dY
#   ungroup() %>%
#   distinct(i_, .keep_all=TRUE)  # keep only one dY for each i_
# 
# C_group <- df_DD %>% filter(C_==1)
# C_group <- C_group %>%
#   mutate(Yt = ifelse(t_==t__,Y_,NA),
#          Ygm1 = ifelse(t_==g__-1,Y_,NA)) %>%
#   group_by(i_) %>%
#   mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
#   ungroup() %>%
#   distinct(i_, .keep_all=TRUE)
# 
# outcome_regressors <- c('x1','x2','x4')
# formula <- paste('dY ~', paste(outcome_regressors, collapse = '+'))
# m0 <- lm(formula, data=C_group)
# Em1 <- mean(T_group$dY)
# Em0 <- predict(m0, newdata=T_group) %>% mean()
# Em1 - Em0
