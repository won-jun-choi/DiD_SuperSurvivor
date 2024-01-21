# Meta =====================
# Title: supersurvivor
# Author: Sanghee
# Last Edit/Editor: Dec-03-2023 / Sanghee
# Description: haha

library(cuRe)

t_max <- df %>% pull(t) %>% max()

df$time_to_treat<-df$G
df$time_to_treat[df$time_to_treat > t_max] <- t_max

#generate status
df$status <- ifelse(df$C == 0, 1, ifelse(df$C == 1, 0, NA))

super_survivor_logit <- fit.cure.model(Surv(time_to_treat,status) ~ z1+z2, formula.surv=list(~x1+x2+1), data=df, type='mixture',
                                       bhazard=NULL,dist='weibull', link='logit')

summary(super_survivor_logit)  # pi: logit, theta: weibull
gammahat <- super_survivor_logit$coefs[['1']]

###### pb of being supersurvivor
df$phat<-exp(gammahat[1]+gammahat[2]*df$x1+gammahat[3]*df$x2)/(1+exp(gammahat[1]+gammahat[2]*df$x1+gammahat[3]*df$x2))

rm(list = c("gammahat", "super_survivor_logit", "t_max"))

