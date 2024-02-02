# Meta =====================
# Title: DiD using supersurvivor
# Author: Wonjun
# Last Edit: Dec-03-2023
# Description: haha

# potential input =====================
time_variable <- 't'
unit_variable <- 'unit'
group_variable <- 'G'

# generate dataframe for DD =====================
df_DD <- df
df_DD['t_'] <- df_DD[time_variable]
df_DD['i_'] <- df_DD[unit_variable]
df_DD['G_'] <- df_DD[group_variable]

t_min <- df_DD['t_'] %>% min()
t_max <- df_DD['t_'] %>% max()  # we can get ATT from (t_min+1)~t_max
G_max <- df_DD[df_DD['G_']<100, 'G_'] %>% max()
G_censored <- df_DD['G'] %>% max()  # 10000

# dataframe that contains all possible ATTs
results_inf <- data.frame(G = rep(c(1:G_max), each=t_max),
                          t = rep(c(1:t_max), times=G_max))
results_inf <- results_inf %>% filter(G > t_min, t >= G)
results_inf$ATT <- NA

# DD =====================
for (i in 1:nrow(results_inf)) {
  g__ <- results_inf[i, group_variable]
  t__ <- results_inf[i, time_variable]
  T_group <- df_DD %>% filter(G_==g__)
  C_group <- df_DD %>% filter(C_tilde==1)
  
  # m(g,t|X)
  C_group <- C_group %>%
    mutate(Yt = ifelse(t_==t__,Y,NA),
           Ygm1 = ifelse(t_==g__-1,Y,NA)) %>%  # gm1: g minus 1
    group_by(i_) %>%
    mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
    ungroup() %>%
    distinct(i_, .keep_all=TRUE)
  m0 <- lm('dY ~ x1+x2', data=C_group)
  
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
  
  # ATT(g,t)
  results_inf[i,'ATT'] <- Em1 - Em0
}
results_inf <- results_inf %>% rename(ATT_inf = ATT)

rm(list=c('time_variable','unit_variable','group_variable',
          't_min','t_max','G_max','G_censored','i','g__','t__',
          'T_group','C_group','m0', 'm0_reweight', 'Em1','Em0', 'Em0_reweight',
          'results','df_DD'))