# Description: Generate test DGP to check the validity of the estimators.
#              The outcome data has an ideal data without confounders or
#              endogeneity or what so ever.

n <- 1000
t_max <- 5
# generate a balanced panel with n units and t_max time periods
df <- tibble(unit = rep(c(1:n), each = t_max),
             t = rep(c(1:t_max), times = n),
             ones = 1)

### Generate regressors.
# For supersurvivor logit: x1, x2, x3
# For duration model: z1, z2, z3 (where z1=x1, z2=x2)
# For outcome model: x1, x2, x4
df <- df %>% mutate(x1 = rep(runif(n)*(-2) + 1, each = t_max),  # x1 ~ U(-1,1)
                    x2 = rep(runif(n), each = t_max),  # x2 ~ U(0,1)
                    x3 = rep(rnorm(n), each = t_max),  # x3 ~ N(0,1)
                    x4 = rep(rnorm(n), each = t_max),  # x4 ~ N(0,1)
                    z1 = x1,
                    z2 = x2,
                    z3 = rep(rnorm(n), each = t_max))  # z3 ~ N(0,1) 

### Supersurvivor logit
gamma <- c(-1, 2, -2, 2)  # supersurvivor logit parameters (g0,g1,g2,g3)
print('gamma: ')
print(gamma)
df <- df %>% 
  group_by(unit) %>%
  mutate(U = rlogis(1,0,1)) %>%
  mutate(Xg = gamma[1] + gamma[2]*x1 + gamma[3]*x2 + gamma[4]*x3) %>%
  mutate(p_cured = exp(Xg) / (1 + exp(Xg)) ) %>%
  mutate(p_uncured = 1 - p_cured) %>%
  mutate(C_tilde = ifelse(Xg+U>0,1,0)) %>% # supersurvivor indicator
  ungroup()

# df$C_tilde %>% mean()



### Duration model
print("Duration model: log-normal")
lambda <- c(4, 2, -3, 1)  # parameters (z0, z1, z2, z3)
print('lambda: ')
print(lambda)
df <- df %>% 
  group_by(unit) %>%
  mutate(V = rnorm(n=1)) %>% # survival function error
  mutate(Zl = lambda[1] + lambda[2]*z1 + lambda[3]*z2 + lambda[4]*z3) %>%
  mutate(G_star = exp(Zl + V)) %>% 
  mutate(G_star = ifelse(C_tilde==1,10000-2,G_star)) %>%
  mutate(G_star = floor(G_star+2)) %>%  # so that first adoption is at t=2
  mutate(G = ifelse(G_star>t_max,10000,G_star)) %>%  # censored and supersurvivors
  mutate(C = ifelse(G==10000,1,0)) %>% # C group indicator
  ungroup()

# Weibull survival function
# lambda <- c(-1, 2, -3, 4)  # parameters (z0~z3)
# print('lambda: ')
# print(lambda)
# df <- df %>% 
#   group_by(unit) %>%
#   mutate(u = runif(1)) %>% # uniform rv for probability inversion
#   mutate(V = rnorm(n=1)) %>% # survival function error
#   mutate(G_star = exp(-lambda[1] - lambda[2]*z1 - lambda[3]*z2 - V) * (-log(u))) %>%
#   mutate(G_star = ifelse(C_tilde==1,10000-2,G_star)) %>% # -2 just to make it look better
#   mutate(G_star = floor(G_star+2)) %>%  # shift to right so that the first adoption t==2
#   mutate(G = ifelse(G_star>t_max,10000,G_star)) %>%  # censored and supersurvivors
#   mutate(C = ifelse(G_star>t_max,1,0)) %>% # C group indicator
#   ungroup()

# histogram of G
if (sys.nframe() == 0) {
  hist(df$G[df$G != 10000])
  df %>% group_by(G) %>% summarise(n=n())
}

# In C==1, the proportion of supersurvivors is?
if (sys.nframe() == 0) {
  num = df[df$C==1, 'C_tilde'] %>% sum()
  denom = df[df$C==1, 'C'] %>% sum()
  print("The proportion of supersurvivor in C==1:")
  print(num/denom)
  rm(list = c('num', 'denom'))
}

## Outcome model
beta <- c(10,3,4,5)  # x0,x1,x2,x4
df$e <- stats::rnorm(n=n*t_max, mean=0, sd=1)  # Y0 error
df <- df %>% 
  mutate(Y0 = beta[1] + x1*beta[2] + x2*beta[3] + x4*beta[4] + 1*t + 0.5*e) %>%
  mutate(Y1 = Y0 + 1) %>%
  mutate(Y = Y0*(G>t) + Y1*(G<=t))

# print true ATT
print("True ATT: ")
print(1)

# time trend by G
if (sys.nframe() == 0) {
  ggplot(df, aes(x = t, y = Y, group = factor(G), color = factor(G))) +
    stat_summary(fun = mean, geom = 'line') +
    labs(title = "Time Trend of Y by Group G",
         x = "Time (t)",
         y = "Y",
         color = "Group (G)")
}

# Save data as csv in temp folder
write_csv(df, "analysis/temp/testDGP.csv")

# clear local variables
rm(list=c('n','t_max','gamma','lambda','df','beta'))


