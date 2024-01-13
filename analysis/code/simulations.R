# Meta =====================
# Title: Simulations
# Author: Wonjun, Juhyun, Sanghee
# Last Edit/Editor: Dec-16-2023 / Wonjun
# Description: Main file to run Monte Carlo simulations.

rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here)

###### DGP ######
source('analysis/code/simDGP.R')  # generate simDGP.csv
df <- read_csv('analysis/temp/simDGP.csv')

# histrogram of G other than 10000
hist(df$G[df$G != 10000])

df %>% group_by(G) %>% summarise(n()) # the number of units in each group
df %>% filter(C_tilde==0,C==1) %>% nrow()  # the number of censored units
df %>% filter(C_tilde==1,C==1) %>% nrow()  # the number of supersurvivors

cor(df$V, df$e)  # corr between surv. ftn. err. and observed Y

# Plot average Y grouped by G
ggplot(df, aes(x = t, y = Y, group = factor(G), color = factor(G))) +
  stat_summary(fun = mean, geom = 'line') +
  labs(title = "Time Trend of Y by Group G",
       x = "Time (t)",
       y = "Y",
       color = "Group (G)")

###### Supersurvivor regression ######
source('analysis/code/supersurvivor.R')

# show df of G=10000 (supersurvivors)
View(df %>% filter(G==10000) %>% select(i,G_star,phat))

df %>% filter(G==10000) %>% group_by(C_tilde) %>% summarise(mean(phat))

###### DiD using reweighting #######
source('analysis/code/DDsurv.R')
library(xtable)
results <- results %>%
  mutate(TE = ifelse(G==t, 1, 0),
         TE = ifelse(G==t-1, 0.8, TE),
         TE = ifelse(G==t-2, 0.6, TE)) %>%
  select(G,t,TE,ATT,ATT_reweight)

# print the xtable result without row number
print(xtable(results), include.rownames = FALSE)
