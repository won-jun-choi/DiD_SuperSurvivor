# This is a sample code to run a simulation

rm(list=ls())
library(tidyverse)
library(here)
library(xtable)
library(dplyr)

###### DGP ######
source('analysis/code/testDGP.R')  # generate testDGP.csv in temp folder.
df <- read_csv('analysis/temp/testDGP.csv')

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

# DiD estimators
source('analysis/code/DiD_infeasible.R')
res_inf <- DiD_infeasible(data=df,
                          unit_variable = 'unit',
                          time_variable = 't',
                          group_variable = 'G',
                          supersurvivor_indicator = 'C_tilde',
                          outcome_variable = 'Y',
                          outcome_regressors = c('x1','x2','x4'))

source('analysis/code/DiD_CS.R')
res_CS <- DiD_CS(data=df,
                 unit_variable = 'unit',
                 time_variable = 't',
                 group_variable = 'G',
                 censored_indicator = 'C',
                 outcome_variable = 'Y',
                 outcome_regressors = c('x1','x2','x4'))

source('analysis/code/DiD_supersurvivor_NoControl.R')
res_SS_noX <- DiD_supersurvivor_NoControl(data=df,
                                          unit_variable = 'unit',
                                          time_variable = 't',
                                          group_variable = 'G',
                                          outcome_variable = 'Y',
                                          censored_indicator = 'C',
                                          logit_regressors = c('z1','z2','z3'),
                                          survival_regressors = c('x1','x2','x3'),
                                          survival_function_type = 'LogNormal_discrete',
                                          survival_MLE_optimizer = 'Nelder-Mead',
                                          MLE_init_values = c(-1,2,0,0,-1,2,0,0,2))

source('analysis/code/DiD_supersurvivor.R')
res_SS <- DiD_supersurvivor(data=df,
                            unit_variable = 'unit',
                            time_variable = 't',
                            group_variable = 'G',
                            outcome_variable = 'Y',
                            outcome_regressors = c('x1','x2','x4'),
                            censored_indicator = 'C',
                            logit_regressors = c('z1','z2','z3'),
                            survival_regressors = c('x1','x2','x3'),
                            survival_function_type = 'LogNormal_discrete',
                            survival_MLE_optimizer = 'Nelder-Mead',
                            MLE_init_values = c(-1,2,0,0,-1,2,0,0,2))                                          

# Put everything in one table
res <- res_inf %>% 
  left_join(res_CS, by = c('t','G')) %>%
  left_join(res_SS_noX, by = c('t','G')) %>%
  left_join(res_SS, by = c('t','G')) %>% 
  mutate(TE=1)
rm(list = c('res_CS','res_inf','res_SS','res_SS_noX',
            'my_SuperSurvivor',
            'DiD_CS','DiD_infeasible','DiD_supersurvivor','DiD_supersurvivor_NoControl'))
  