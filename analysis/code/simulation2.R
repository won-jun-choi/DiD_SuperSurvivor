# Meta =====================
# Title: Simulation2
# Author: Sanghee
# Last Edit/Editor: Jan-19-2024 / Wonjun
# Description: Main file to run Monte Carlo simulations.

rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here)

###### DGP ######
source('analysis/code/simDGP2.R')  # generate simDGP2.csv in temp folder
df <- read_csv('analysis/temp/simDGP2.csv')

# histrogram of G other than 10000
hist(df$G[df$G != 10000])

df %>% group_by(G) %>% summarise(n()) # the number of units in each group
df %>% filter(C_tilde==0,C==1) %>% nrow()  # the number of censored units
df %>% filter(C_tilde==1,C==1) %>% nrow()  # the number of supersurvivors

# Plot average Y grouped by G
ggplot(df, aes(x = t, y = Y, group = factor(G), color = factor(G))) +
  stat_summary(fun = mean, geom = 'line') +
  labs(title = "Time Trend of Y by Group G",
       x = "Time (t)",
       y = "Y",
       color = "Group (G)")
ggsave('analysis/output/fig_DGP2.png', width = 6, height = 4)

###### Supersurvivor regression ######
source('analysis/code/supersurvivor.R')  # generate super_survivor in dataframe

###### DiD using reweighting #######
source('analysis/code/DDsurv.R')  # generate DDsurv in dataframe

###### print the result ###### 
DDsurv <- DDsurv %>%
  mutate(TE = ifelse(G==t, 3, 0),
         TE = ifelse(G==t-1, 2, TE),
         TE = ifelse(G==t-2, 1, TE),
         TE = ifelse(G==t-3, 0.8, TE),
         TE = ifelse(G==t-4, 0.6, TE),
         TE = ifelse(G==t-5, 0.3, TE)) %>%
  select(G,t,TE,ATT,ATT_reweight)

# add DiD using C_tilde (infeasible)
source('analysis/code/InfeasibleDD.R')
# merge DDsurv and results_inf using G and t as keys
DDsurv <- DDsurv %>%
  left_join(results_inf, by=c('G','t')) %>%
  select(G,t,TE,ATT,ATT_reweight,ATT_inf)

# print the xtable result without row number
print(xtable(DDsurv), include.rownames = FALSE)

df_censored <- df %>% filter(C==1) %>% select(unit,t,C,C_tilde,p_cured,phat)
