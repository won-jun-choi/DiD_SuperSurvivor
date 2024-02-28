# Meta =====================
# Title: Simulation1
# Author: Wonjun
# Last Edit/Editor: Feb-15-2024 / Wonjun
# Description: MC Simulation 1

rm(list=ls())
library(tidyverse)
library(here)
library(xtable)

###### DGP ######
source('analysis/code/simDGP1.R')  # generate simDGP.csv in temp folder.
df <- read_csv('analysis/temp/simDGP1.csv')

# histrogram of G other than 10000
hist(df$G[df$G != 10000])

df %>% group_by(G) %>% summarise(n()) # the number of units in each group
df %>% filter(C_tilde==0,C==1) %>% nrow()  # the number of censored units
df %>% filter(C_tilde==1,C==1) %>% nrow()  # the number of supersurvivors

# Histogram of G_star other than 10000.
# The observations with G_star > 5 are considered as C=1 in naive DD.
hist(df$G_star[(df$G_star != 10000) & (df$G_star <=10)])

# Plot average Y grouped by G -- do it in a separate code.
ggplot(df, aes(x = t, y = Y, group = factor(G), color = factor(G))) +
  stat_summary(fun = mean, geom = 'line') +
  labs(title = "Time Trend of Y by Group G",
       x = "Time (t)",
       y = "Y",
       color = "Group (G)")

# ggsave('analysis/output/fig_DGP1.png', width = 6, height = 4)

###### Supersurvivor regression ######
# no need to run unless want to diagnose the supersurvivor estimation
# source('analysis/code/supersurvivor.R')
# res <- my_SuperSurvivor(data = df,
#                         unit_variable = 'unit',
#                         time_variable = 't',
#                         duration_variable = 'G',
#                         censored_indicator = 'C',
#                         logit_regressors = c('x1','x2','x3','x4'),
#                         survival_regressors = c('z1','z2','z3'),
#                         survival_function_type = 'LogNormal',
#                         optimization_method = 'Nelder-Mead',
#                         return_params = TRUE)
# df$phat <- res[1]
# gamma_hat <- res[2]
# lambda_hat <- res[3]
# df %>% select(p_cured,phat, C_tilde) %>% View()

###### DiD using reweighting #######
source('analysis/code/DDsurv.R')
res = DiD_supersurvivor(data = df,
                       unit_variable = 'unit',
                       time_variable = 't',
                       group_variable = 'G',
                       outcome_variable = 'Y',
                       outcome_regressors = c('x1','x2','x3','x4'),
                       censored_indicator = 'C',
                       logit_regressors = c('x1','x2','x3','x4'),
                       survival_regressors = c('z1','z2','z3'),
                       survival_function_type = 'LogNormal_discrete',
                       survival_MLE_optimizer='Nelder-Mead')
res$TE <- 1

# add DiD using C_tilde (infeasible)
source('analysis/code/InfeasibleDD.R')
# merge DDsurv and results_inf using G and t as keys
res_infeasible <- DD_infeasible(data=df,
                                unit_variable = 'unit',
                                time_variable = 't',
                                group_variable = 'G',
                                supersurvivor_indicator = 'C_tilde')

res <- res %>%
  left_join(res_infeasible, by=c('G','t')) %>%
  select(G,t,TE,ATT,ATT_reweight,ATT_infeasible)

# print the xtable result without row number
print(xtable(DDsurv), include.rownames = FALSE)

df$p_cured %>% mean()
df$phat %>% mean()

mean((df$p_cured - df$phat)^2)*10000

df %>% group_by(C,C_tilde) %>% summarise(mean(phat))
df_censored <- df %>% filter(G==10000) %>% select(unit, t, p_cured, phat)
