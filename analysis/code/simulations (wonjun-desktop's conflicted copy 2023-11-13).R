# Meta =====================
# Title: Simulations
# Author: Wonjun
# Last Edit: Sep-11-2023
# Description: haha

library(tidyverse)
library(here)
# DGP
n <- 1000
t_max <- 5
df <- tibble(unit = rep(c(1:n), each = t_max),
             t = rep(c(1:t_max), times = n),
             ones = 1,
             x1 = rep(runif(n)*(-2) + 1, each = t_max),  # X: SS regressors
             x2 = rep(runif(n)*(-2) + 1, each = t_max),
             z1 = x1,  # Z: survival function regressors
             z2 = x2,
             p_uncured = exp(1 + 2*x1 + 0.5*x2)/(1+exp(1 + 2*x1 + 0.5*x2))
             ) %>%
  group_by(unit) %>%
  mutate(C = ifelse(runif(1)<p_uncured,0,1),
         U = runif(1),  # for survival function
         G_star = exp(-0.5 - 0.5*z1 - 3*z2)*(-log(U)),  # survival function
         G_star = floor(G_star),
         G = ifelse(C==1,10000,G_star),
         G = ifelse(G_star>=t_max, t_max,G_star)) %>%
  ungroup()
df <- df %>% mutate(Y_base = U+rnorm(1),
         Y_t1 = Y_base + 0.5 + 1*(t-G==1) + 0.8*(t-G==2) +
           0.6*(t-G==3) + 0.4*(t-G==4) + rnorm(1),
         Y_t2 = Y_t1 - 0.5 + rnorm(1,0,1)+U-0.5)
df

# potential outcomes and treatment effect
e = rnorm(n)
Y_base = U + e - 0.5
