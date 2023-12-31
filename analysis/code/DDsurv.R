# Meta =====================
# Title: DiD using supersurvivor
# Author: Wonjun
# Last Edit: Dec-03-2023
# Description: haha

time_variable <- 't'
unit_variable <- 'unit'
group_variable <- 'G'

df_DD <- df
df_DD['t_'] <- df_DD[time_variable]
df_DD['i_'] <- df_DD[unit_variable]
df_DD['G_'] <- df_DD[group_variable]

t_min <- df_DD['t_'] %>% min()
t_max <- df_DD['t_'] %>% max()  # we can get ATT from (t_min+1)~t_max
G_max <- df_DD[df_DD['G_']<100, 'G_'] %>% max()
G_censored <- df_DD['G'] %>% max()

results <- data.frame(G = rep(c(1:G_max), each=t_max),
                      t = rep(c(1:t_max), times=G_max))
results <- results %>% filter(G > t_min, t >= G)
results$ATT <- NA
results$ATT_reweight <- NA

# ATT(g,t)
for (i in 1:nrow(results)) {
  g__ <- results[i, group_variable]
  t__ <- results[i, time_variable]
  T_group <- df_DD %>% filter(G_==g__)
  C_group <- df_DD %>% filter(G_==G_censored)
  
  # m(g,t|X)
  C_group <- C_group %>%
    mutate(Yt = ifelse(t_==t__,Y,NA),
           Ygm1 = ifelse(t_==g__-1,Y,NA)) %>%
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
  results[i,'ATT'] <- Em1 - Em0
  results[i,'ATT_reweight'] <- Em1 - weighted.mean(T_group$Em0, T_group$phat)
}

# loopppp
# t__ <- 2
# g__ <- 2
# T_group <- df_DD %>% filter(G_==g__)
# C_group <- df_DD %>% filter(G_==G_censored)  # why different?
# C_group <- df_DD %>% filter(C==1)  # why different?
# 
# # m(g,t|X)
# C_group <- C_group %>%
#   mutate(Yt = ifelse(t_==t__,Y,NA),
#          Ygm1 = ifelse(t_==g-1,Y,NA)) %>%
#   group_by(i_) %>%
#   mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
#   ungroup() %>%
#   distinct(i_, .keep_all=TRUE)
# m0 <- lm('dY ~ x1+x2', data=C_group)
# 
# T_group <- T_group %>%
#   mutate(Yt = ifelse(t_==t__,Y,NA),
#          Ygm1 = ifelse(t_==g-1,Y,NA)) %>%
#   group_by(i_) %>%
#   mutate(dY = sum(Yt, na.rm=T) - sum(Ygm1, na.rm=T)) %>%
#   ungroup() %>%
#   distinct(i_, .keep_all=TRUE)
# 
# Em1 <- mean(T_group$dY)
# T_group$Em0 <- predict(m0, newdata=T_group)
# Em0 <- predict(m0, newdata=T_group) %>% mean()
# 
# TE <- Em1 - Em0
# print(TE)
# 
# # Reweight
# Em0_reweright <- weighted.mean(T_group$Em0, T_group$phat)
# print(Em1 - Em0_reweright)
