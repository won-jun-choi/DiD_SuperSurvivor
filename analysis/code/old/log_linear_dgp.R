# Generate super survivor logit
#X: supersurvivor logit covariates
#Z: Duration covariates
# X and Z maybe identical except that x has intercept (Handbook of Survival Analysis. pg 115)
library(dplyr)
library(MASS)

n<-2000
t_max<-20
gamma <- c(1, 1, 0.5,1,2)  # params for supersurvivor logit
lambda <- c( 1, 2,3)  # parmas for duration
df <- tibble(unit = rep(c(1:n), each = t_max),
             t = rep(c(1:t_max), times = n),
             ones = 1,
             x1 = rep(runif(n)*(-199) + 100, each = t_max),  
             x2 = rep(runif(n)*(-200)+120 , each = t_max),
             x3 = rep(runif(n) , each = t_max),
             x4 = rep(runif(n)*(-1.5) , each = t_max),
             z1 = x1,   
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
  mutate(U = rlogis(1,0,1)) %>%
  mutate(p_uncured = exp(gamma[1]+gamma[2]*x1+gamma[3]*x2+gamma[4]*x3+gamma[5]*x4)/(1+exp(gamma[1]+gamma[2]*x1+gamma[3]*x2+gamma[4]*x3+gamma[5]*x4)),
         C_tilde = ifelse(U + gamma[1]+gamma[2]*x1+gamma[3]*x2+gamma[4]*x3+gamma[5]*x4 > 0,0,1)  
         ) 

# number of supersurvivor: 
sum(df$C_tilde == 1)


#log-linear
#u : Duration error term
#ln(y)=z'b+u <=> y=exp(z'b+u)
#y: duration! (not outcome; don't be confused)
#question: do we need constant term?
df <- df %>%
  mutate(u = rnorm(n(), mean = 0, sd = 1),
         y = exp(lambda[1]*z1+ lambda[2]*z2 +lambda[3]+ u))  # Example formula for y, you can replace it with your desired formula

df<-df%>%
  group_by(unit)%>%
  mutate(G_star=floor(y),
         G_star=G_star+2,
         G_star=ifelse(C_tilde==1,10000,G_star),
         G = ifelse(G_star>t_max,t_max,G_star),
        C = ifelse(G_star>t_max,0,1),  #Not censored indicator
         )

# proportion of supersur C not censored indicator!!
cat("Number of occurrences where C = 0 and C_tilde = 0:", sum(df$C == 0 & df$C_tilde == 0), "\n")
cat("Number of occurrences where C = 0 and C_tilde = 1:", sum(df$C == 0 & df$C_tilde == 1), "\n")


#log likelihood
# problem in continuous=>discrete: we observe G_star (integer), but y is continuous
log_likelihood <- function(parameters, c, y, x1, x2, x3, x4, z1, z2) {
  lambda1 <- parameters[1]
  lambda2 <- parameters[2]
  gamma1 <- parameters[3]
  gamma2 <- parameters[4]
  gamma3 <- parameters[5]
  gamma4 <- parameters[6]
  gamma5 <- parameters[7]
  mu <- parameters[8]
  sigma <- parameters[9]
  lambda3 <- parameters[10]
  p_uncured <-exp(gamma1 + gamma2 * x1 + gamma3 * x2 + gamma4 * x3 + gamma5 * x4) / (1 + exp(gamma1 + gamma2 * x1 + gamma3 * x2 + gamma4 * x3 + gamma5 * x4))

  
  S_u <- 1 - pnorm((log(y) - (lambda1 * z1 + lambda2 * z2+lambda3)) / sigma, mean = mu, sd = sigma)
  
  f_u <- dnorm((log(y) - (lambda1 * z1 + lambda2 * z2+lambda3)) / sigma, mean = mu, sd = sigma) / (sigma * y) #G should not be zero!
  
  log_likelihood_values <- c * log(p_uncured) + c * log(f_u) + (1 - c) * log(1 - p_uncured + (p_uncured * S_u))
  
  return(-sum(log_likelihood_values))  # Sum the log-likelihood values
}


# Initial guess for parameters
initial_guess <- c(rep(1, 7), 2, 3,4)  # Initial guess for lambda1, lambda2, ..., sigma
initial_guess <- c(1,2,2,1,0.5,1,2,0,1,3) # True value

# Optimization
result <- optim(par = initial_guess, fn = log_likelihood, c = df$C, y = df$G,
                x1 = df$x1, x2 = df$x2, x3 = df$x3, x4 = df$x4,
                z1 = df$z1, z2 = df$z2, method="Nelder-Mead")

# Estimated parameters
estimated_params <- result$par
print(estimated_params)

# gamma1_ <- estimated_params[3]
# gamma2_ <- estimated_params[4]
# gamma3_ <- estimated_params[5]
# gamma4_ <- estimated_params[6]
# gamma5_ <- estimated_params[7]
# 
# df$phat <- exp(gamma1_ + gamma2_ * df$x1 + gamma3_ * df$x2 + gamma4_ * df$x3 + gamma5_ * df$x4) / (1 + exp(gamma1_ + gamma2_ * df$x1 + gamma3_ * df$x2 + gamma4_ * df$x3 + gamma5_ * df$x4))
# 
# df %>% select(p_uncured, phat) %>% View()
