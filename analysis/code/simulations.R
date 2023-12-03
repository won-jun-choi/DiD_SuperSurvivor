# Meta =====================
# Title: Simulations
# Author: Wonjun
# Last Edit: Sep-11-2023
# Description: haha
rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here)

# DGP
n <- 1000
t_max <- 5
gamma <- c(1, 2, 0.5)  # supersurvivor
lambda <- c(0.1, 3, 5)  # survival function
TE <- c(1, 0.8, 0.6)  # time varying treatment effect
df <- tibble(unit = rep(c(1:n), each = t_max),
             t = rep(c(1:t_max), times = n),
             ones = 1,
             x1 = rep(runif(n)*(-2) + 1, each = t_max),  # X: SS regressors
             x2 = rep(runif(n)*(-2) + 1, each = t_max),
             z1 = x1,  # Z: survival function regressors
             z2 = x2) %>%
  group_by(unit) %>%
  mutate(p_uncured = exp(gamma[1]+gamma[2]*x1+gamma[3]*x2)/(1+exp(gamma[1]+gamma[2]*x1+gamma[3]*x2)),
         C_tilde = ifelse(runif(1)<p_uncured,0,1),  # SS indicator
         G = ifelse(C_tilde==1,10000,-1),
         u = runif(1),
         V = rnorm(1), # survival function error
         G_star = exp(-lambda[1] - lambda[2]*z1 - lambda[3]*z2 - V) * (-log(u)),
         G_star = ifelse(C_tilde==1,10000,G_star),
         G_star = floor(G_star+2),
         G = ifelse(G_star>t_max,10000,G_star),  # censored
         d = ifelse(G_star>t_max,1,0),  # censored indicator
         C = ifelse(G_star>t_max,1,0)) %>%
  ungroup() %>%
  mutate(tau = t-G) %>%
  group_by(unit) %>%
  mutate(unitFE = rnorm(1),
         Y_base = (V*runif(G))+rnorm(1) + unitFE) %>%
  ungroup() %>%
  mutate(Y0 = Y_base + 1*t + 0.5*rnorm(1),
         Y1 = Y0 + TE[1]*(tau==0) + TE[2]*(tau==1) + TE[3]*(tau==2),
         Y = Y0*(G>t) + Y1*(G<=t))

hist(df$G[df$G != 10000])
df %>% filter(C==1) %>% nrow()
df %>% filter(C_tilde==0,C==1) %>% nrow()

df %>% group_by(G) %>% summarise(n())

cor(df$V, df$Y)

# Create a boxplot of Y grouped by G
ggplot(df, aes(x = t, y = Y, group = factor(G), color = factor(G))) +
  stat_summary(fun = mean, geom = "line") +
  labs(title = "Mean Time Trend of Y by Group G",
       x = "Time (t)",
       y = "Mean Y",
       color = "Group (G)")

library(cuRe)

super_survivor_logit <- fit.cure.model(Surv(G,C) ~ x1+x2, data=df, type='mixture',
                                       bhazard='bhaz',dist='weibull', link='logit')
