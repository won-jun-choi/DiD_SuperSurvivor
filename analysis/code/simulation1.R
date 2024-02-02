# Meta =====================
# Title: Simulation1
# Author: Wonjun
# Last Edit/Editor: Jan-19-2024 / Wonjun
# Description: MC Simulation 1

rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, xtable)

###### DGP ######
source('analysis/code/simDGP1.R')  # generate simDGP.csv
df <- read_csv('analysis/temp/simDGP1.csv')

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

# save the plot
ggsave('analysis/output/fig_DGP1.png', width = 6, height = 4)

###### Supersurvivor regression ######
source('analysis/code/supersurvivor.R')

###### DiD using reweighting #######
source('analysis/code/DDsurv.R')

###### print the results ###### 
# add true TE to the results
DDsurv <- DDsurv %>%
  mutate(TE = ifelse(G==t, 1, 0),
         TE = ifelse(G==t-1, 0.8, TE),
         TE = ifelse(G==t-2, 0.6, TE)) %>%
  select(G,t,TE,ATT,ATT_reweight,ATT_reweight2)

# add DiD using C_tilde (infeasible)
source('analysis/code/InfeasibleDD.R')
# merge DDsurv and results_inf using G and t as keys
DDsurv <- DDsurv %>%
  left_join(results_inf, by=c('G','t')) %>%
  select(G,t,TE,ATT,ATT_reweight,ATT_reweight2,ATT_inf)

# print the xtable result without row number
print(xtable(DDsurv), include.rownames = FALSE)

df$p_cured %>% mean()
df$phat %>% mean()

mean((df$p_cured - df$phat)^2)*10000

df %>% group_by(C,C_tilde) %>% summarise(mean(phat))
df_censored <- df %>% filter(G==10000) %>% select(unit, t, p_cured, phat)
