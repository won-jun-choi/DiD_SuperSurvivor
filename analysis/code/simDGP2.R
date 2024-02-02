# Meta =====================
# Title: simulation DGP
# Author: Wonjun
# Last Edit/Editor: Dec-16-2023 / Wonjun
# Description: DGP for Monte Carlo simulations

n <- 1000
t_max <- 10
gamma <- c(1, 2, 0.5)  # params for supersurvivor logit
lambda <- c(1, 2, 3)  # parmas for survival function
TE <- c(3,2,1, 0.8, 0.6,0.3)  # time varying treatment effect
df <- tibble(unit = rep(c(1:n), each = t_max),
             t = rep(c(1:t_max), times = n),
             ones = 1,
             x1 = rep(runif(n)*(-2) + 1, each = t_max),  # X: SS regressors
             x2 = rep(runif(n)*(-2) , each = t_max),
             x3 = rep(runif(n) , each = t_max),
             x4 = rep(runif(n)*(-1.5) , each = t_max),
             z1 = x1,  # Z: survival function regressors
             z2 = x2) %>%
     mutate(
      x1 = (x1 - mean(x1)) / sd(x1),
      x2 = (x2 - mean(x2)) / sd(x2),
      x3 = (x3 - mean(x3)) / sd(x3),
      x4 = (x4 - mean(x4)) / sd(x4),
      z1 = x1,  # Z: survival function regressors
      z2 = x2
    ) %>%
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
         e_base = 20*(V^2)+10*V,
         Y_base =  unitFE + groupFE) %>%
  ungroup() %>%
  mutate(cov=210+x1+ 30*x2+12*x3+10*x4, 
         e=e_base+rnorm(1), ###
         Y0 = Y_base +cov+ t + 0.5*rnorm(1)-2*log(t)*(G==10000)*(G_star<10000)+0.5*t*(G_star==10000)+e,
         Y1 = Y0 + TE[1]*(tau==0) + TE[2]*(tau==1) + TE[3]*(tau==2)+TE[4]*(tau==3)+TE[5]*(tau==4)+TE[6]*(tau==5),
         Y = Y0*(G>t) + Y1*(G<=t))
  
# Save data as csv in temp folder
write_csv(df, "analysis/temp/simDGP2.csv")

# clear local variables
rm(list=c('n','t_max','gamma','lambda','TE','df'))


