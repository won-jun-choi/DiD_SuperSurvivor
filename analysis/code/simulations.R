# Meta =====================
# Title: Simulations
# Author: Wonjun, Juhyun, Sanghee
# Last Edit/Editor: Sep-11-2023 / Wonjun
# Description: haha
rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here)

# here::here()
# here::i_am('analysis/code/simulations.R')
here::here()

###### DGP ######
source('analysis/code/simDGP.R')

# hist(df$G[df$G != 10000])
# df %>% filter(C==1) %>% nrow()  # the number of control units (ss+censored)
# df %>% filter(C_tilde==0,C==1) %>% nrow()  # the number of censored units
# 
# df %>% group_by(G) %>% summarise(n()) # the number of units in each group
# 
# cor(df$V, df$Y)  # corr between surv. ftn. err. and observed Y

# Create a plot of Y grouped by G
ggplot(df, aes(x = t, y = Y, group = factor(G), color = factor(G))) +
  stat_summary(fun = mean, geom = "line") +
  labs(title = "Mean Time Trend of Y by Group G",
       x = "Time (t)",
       y = "Mean Y",
       color = "Group (G)")

###### Supersurvivor regression ######
source('analysis/code/supersurvivor.R')
View(df %>% filter(G==10000) %>% select(i,G_star,phat))
df %>% filter(G==10000) %>% group_by(C_tilde) %>% summarise(mean(phat))

###### DiD using reweighting #######
source('analysis/code/DDsurv.R')
