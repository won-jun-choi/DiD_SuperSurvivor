# Meta =====================
# Title: supersurvivor
# Author: Sanghee
# Last Edit/Editor: Dec-03-2023 / Sanghee
# Description: haha

library(cuRe)

###### test DGP ######
# set.seed(1)
# n <- 1000
# t_max <- 5
# gamma <- c(1, 2, 0.5)  # supersurvivor
# lambda <- c(0.1, 3, 5)  # survival function
# TE <- c(1, 0.8, 0.6)  # time varying treatment effect
# df <- tibble(unit = rep(c(1:n), each = t_max),
#              t = rep(c(1:t_max), times = n),
#              ones = 1,
#              x1 = rep(runif(n)*(-2) + 1, each = t_max),  # X: SS regressors
#              x2 = rep(runif(n)*(-2) + 1, each = t_max),
#              z1 = x1,  # Z: survival function regressors
#              z2 = x2) %>%
#   group_by(unit) %>%
#   mutate(p_cured = exp(gamma[1]+gamma[2]*x1+gamma[3]*x2)/(1+exp(gamma[1]+gamma[2]*x1+gamma[3]*x2)), #package seem to use p_cured 
#          C_tilde = ifelse(runif(1)<p_cured,1,0),  # SS indicator 바꿈 check
#          G = ifelse(C_tilde==1,10000,-1),
#          u = runif(1),
#          V = rnorm(1), # survival function error
#          G_star = exp(-lambda[1] - lambda[2]*z1 - lambda[3]*z2 - V) * (-log(u)),
#          G_star = ifelse(C_tilde==1,10000,G_star),
#          G_star = floor(G_star+2),
#          G = ifelse(G_star>t_max,10000,G_star),  # censored
#          d = ifelse(G_star>t_max,1,0),  # censored indicator
#          C = ifelse(G_star>t_max,1,0)) %>%
#   ungroup() %>%
#   mutate(tau = t-G) %>%
#   group_by(unit) %>%
#   mutate(unitFE = rnorm(1),
#          Y_base = (V*runif(G))+rnorm(1) + unitFE) %>%
#   ungroup() %>%
#   mutate(Y0 = Y_base + 1*t + 0.5*rnorm(1),
#          Y1 = Y0 + TE[1]*(tau==0) + TE[2]*(tau==1) + TE[3]*(tau==2),
#          Y = Y0*(G>t) + Y1*(G<=t))

#generate time_to_treat
t_max <- 10

df$time_to_treat<-df$G-1
df$time_to_treat[df$time_to_treat > t_max] <- t_max
#df$time_to_treat<-df$G_star
#df$time_to_treat[df$time_to_treat > 100] <- 100
df<-df%>% mutate(C = ifelse(G_star>100,1,0))

#generate status
df$status <- ifelse(df$C == 0, 1, ifelse(df$C == 1, 0, NA))

super_survivor_logit <- fit.cure.model(Surv(time_to_treat,status) ~ z1+z2, formula.surv=list(~x1+x2+1), data=df, type='mixture',
                                       bhazard=NULL,dist='weibull', link='logit')

summary(super_survivor_logit)  # pi: logit, theta: weibull
gammahat <- super_survivor_logit$coefs[['1']]

#pb of being supersurvivor

df$phat<-exp(gammahat[1]+gammahat[2]*df$x1+gammahat[3]*df$x2)/(1+exp(gammahat[1]+gammahat[2]*df$x1+gammahat[3]*df$x2))

# to compare with p_uncured:
df1 <- df %>% select(C, C_tilde, phat, p_cured)
df1<-df1 %>% mutate(dif=phat-p_cured)
