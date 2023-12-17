# Meta =====================
# Title: simulation DGP
# Author: Wonjun
# Last Edit/Editor: Dec-16-2023 / Wonjun
# Description: DGP for Monte Carlo simulations

n <- 1000
t_max <- 5
gamma <- c(1, 2, 0.5)  # params for supersurvivor logit
lambda <- c(1, 2, 3)  # parmas for survival function
TE <- c(1, 0.8, 0.6)  # time varying treatment effect
df <- tibble(unit = rep(c(1:n), each = t_max),
             t = rep(c(1:t_max), times = n),
             ones = 1,
             x1 = rep(runif(n)*(-2) + 1, each = t_max),  # X: SS regressors
             x2 = rep(runif(n)*(-2) , each = t_max),
             z1 = x1,  # Z: survival function regressors
             z2 = x2) %>%
  group_by(unit) %>%
  mutate(p_uncured = exp(gamma[1]+gamma[2]*x1+gamma[3]*x2)/(1+exp(gamma[1]+gamma[2]*x1+gamma[3]*x2)),
         p_cured = exp(gamma[1]+gamma[2]*x1+gamma[3]*x2)/(1+exp(gamma[1]+gamma[2]*x1+gamma[3]*x2)),
         C_tilde = ifelse(runif(1)<p_uncured,0,1),  # SS indicator
         C_tilde = ifelse(runif(1)<p_cured,1,0),  # SS indicator 바꿈 check
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
  group_by(G_star) %>%
  mutate(groupFE = rnorm(1)) %>%
  ungroup() %>%
  group_by(unit) %>%
  mutate(unitFE = rnorm(1),
         e = V + rnorm(1),  # Y0 error correlated with V
         Y_base = e + unitFE + groupFE) %>%
  ungroup() %>%
  mutate(Y0 = Y_base + 1*t + 0.5*rnorm(1) + 0.2*t*(G==10000)*(G_star<10000),
         Y1 = Y0 + TE[1]*(tau==0) + TE[2]*(tau==1) + TE[3]*(tau==2),
         Y = Y0*(G>t) + Y1*(G<=t))

# Save data as csv in temp folder
write_csv(df, "analysis/temp/simDGP.csv")

# clear local variables
rm(list=c('n','t_max','gamma','lambda','TE','df'))


