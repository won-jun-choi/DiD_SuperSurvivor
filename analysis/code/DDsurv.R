# Meta =====================
# Title: DiD using supersurvivor
# Author: Wonjun
# Last Edit: Dec-03-2023
# Description: haha

time_variable <- 't'
unit_variable <- 'unit'

t_min <- df %>% select(t) %>% min()
t_max <- df %>% select(t) %>% max()  # we can get ATT from (t_min+1)~t_max

results <- data.frame()

# ATT(g=2,t=2)?
t <- 2
g <- 2
T_group <- df %>% filter(G==g, t==2)
C_group <- df %>% filter(C==1, t==2)

m1 <- lm('Y ~ x1+x2', data=T_group)
m0 <- lm('Y ~ x1+x2', data=C_group)

Em1 <- mean(T_group$Y)
T_group$Em0 <- predict(m0, newdata=T_group)
Em0 <- predict(m0, newdata=T_group) %>% mean()

TE <- Em1 - Em0
print(mean(T_group$Y1) - mean(T_group$Y0))
print(TE)

# Reweight
Em0_reweright <- weighted.mean(T_group$Em0, T_group$phat)
print(Em1 - Em0_reweright)
