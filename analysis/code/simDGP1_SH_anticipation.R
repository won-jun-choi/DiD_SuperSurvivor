# Meta =====================
# Title: simulation DGP
# Author: Wonjun
# Last Edit/Editor: Jan-18-2024 / Wonjun
# Description: DGP for Monte Carlo simulations
# load libary
install.packages("truncnorm")
library(truncnorm)
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
n <- 1000
t_max <- 10


# generate a balanced panel with n units and t_max time periods
df <- tibble(unit = rep(c(1:n), each = t_max),
             t = rep(c(1:t_max), times = n),
             ones = 1)

# generate regressors X: supersurvivor logit, Z: survival function
df <- df %>% mutate(x1 = rep(20*(runif(n)*(-2) + 1), each = t_max),
                    x2 = rep(15*(runif(n)*(-2) + 1), each = t_max),
                    z1 = x1,
                    z2 = x2)

# supersurvivor logit
gamma <- c(0.1, 2, 0.5)  # supersurvivor logit parameters
df <- df %>% 
  group_by(unit) %>%
  mutate(p_uncured = exp(gamma[1] + gamma[2]*x1 + gamma[3]*x2) / (1 + exp(gamma[1] + gamma[2]*x1 + gamma[3]*x2))) %>%
  mutate(p_cured = 1 - p_uncured) %>%
  mutate(C_tilde = ifelse(rtruncnorm(1, a = 0.5, mean = 0, sd = 1)<p_uncured,0,1), # supersurvivor indicator
         C_tilde = ifelse(rtruncnorm(1, a = 0.5, mean = 0, sd = 1)<p_cured,1,0)) %>%
  ungroup()

# average p_cured and C_tilde (should be similar)
# df$p_cured %>% mean()
# df$C_tilde %>% mean()

# survival function
lambda <- c(0.05, 0.15, 0.75)  # Adjust these values experimentally
# parameters for the survival function
df <- df %>% 
  group_by(unit) %>%
  mutate(u = runif(1),
         V = rnorm(1), # survival function error
         G_star = exp(-lambda[1] - lambda[2]*z1 - lambda[3]*z2 - V) * (-log(u)),
         G_star = ifelse(C_tilde==1,10000-2,G_star),  # -2 just to make it look better
         G_star = floor(G_star+2),  # shift to right so that the first adoption t==2
         G = ifelse(G_star>t_max,10000,G_star),  # censored units
         C = ifelse(G_star>t_max,1,0), # C group indicator
         d = C) %>% # for cured model estimation
  ungroup()
count_censor<- sum(df$G_star != 10000 & df$C == 1)
count_g2<- sum(df$G_star==2)
count_g10<- sum(df$G_star==10)
count_sup<- sum(df$G_star == 10000 & df$C == 1)
print(count_censor)
print(count_g10)
print(count_sup)
print(count_g2)
# histogram of G
# hist(df$G[df$G != 10000])
df %>% group_by(G) %>% summarise(n=n())

# In C==1, the proportion of supersurvivors is?
num = df[df$C==1, 'C_tilde'] %>% sum()
denom = df[df$C==1, 'C'] %>% sum()
num/denom  # kind of high...

df %>% filter(C==1) %>% summarise(n=sum(C_tilde))  # proportion of supersur

# generate Y0 and Y1
beta <- c(1,2)
TE <- c(5,4,2,0.8,0.6)  # time varying treatment effect
AE <- c(0.7, 2,5) # add anticipated effects
df <- df %>% 
  mutate(tau = t-G) %>% # time varying treatment
  mutate(delta_gt = ifelse(G_star>t_max & G_star <= 10, 3*sin(2*t), 0)) %>% # group specific trend
  mutate(e = rnorm(1)) %>% # Y0 error
  mutate(Y0 = x1*beta[1] + x2*beta[2] + 1*t + delta_gt + e) %>%
  mutate(Y1 = Y0 + AE[1]*(tau==-3)+AE[2]*(tau==-2)+AE[3]*(tau==-1)+TE[1]*(tau==0) + TE[2]*(tau==1) + TE[3]*(tau==2)+TE[4]*(tau==3)+TE[5]*(tau==4)) %>%
  mutate(Y = Y0*(G>t) + Y1*(G<=t))

# time trend by G
ggplot(df, aes(x = t, y = Y, group = factor(G), color = factor(G))) +
  stat_summary(fun = mean, geom = 'line') +
  labs(title = "Time Trend of Y by Group G",
       x = "Time (t)",
       y = "Y",
       color = "Group (G)")

# Save data as csv in temp folder
write_csv(df, "analysis/temp/simDGP1.csv")

# clear local variables
rm(list=c('n','t_max','gamma','lambda','TE','df','beta','num','denom'))


