# Meta =====================
# Title: Simulations
# Author: Wonjun, Juhyun, Sanghee
# Last Edit/Editor: Dec-16-2023 / Wonjun
# Description: Main file to run Monte Carlo simulations.

rm(list=ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here)

###### DGP ######
source('/Users/sangheemun/Library/Mobile Documents/com~apple~CloudDocs/Sanghee Mun/ Academic/RESEARCH/Staggered_DiD/2024Jan/simDGP.R')  # generate simDGP.csv
df <- read_csv('/Users/sangheemun/Library/Mobile Documents/com~apple~CloudDocs/Sanghee Mun/ Academic/RESEARCH/Staggered_DiD/2024Jan/temp/simDGP.csv')

# histrogram of G other than 10000
hist(df$G[df$G != 10000])

df %>% group_by(G) %>% summarise(n()) # the number of units in each group
df %>% filter(C_tilde==0,C==1) %>% nrow()  # the number of censored units
df %>% filter(C_tilde==1,C==1) %>% nrow()  # the number of supersurvivors

cor(df$V, df$e)  # corr between surv. ftn. err. and observed Y

# Plot average Y grouped by G
ggplot(df, aes(x = t, y = Y, group = factor(G), color = factor(G))) +
  stat_summary(fun = mean, geom = 'line') +
  labs(title = "Time Trend of Y by Group G",
       x = "Time (t)",
       y = "Y",
       color = "Group (G)")

###### Supersurvivor regression ######
source('/Users/sangheemun/Library/Mobile Documents/com~apple~CloudDocs/Sanghee Mun/ Academic/RESEARCH/Staggered_DiD/2024Jan/supersurvivor.R')

# show df of G=10000 (supersurvivors)
View(df %>% filter(G==10000) %>% select(i,G_star,phat))

df %>% filter(G==10000) %>% group_by(C_tilde) %>% summarise(mean(phat))
#TE <- c(3,2,1, 0.8, 0.6,0.3) 
###### DiD using reweighting #######
source('/Users/sangheemun/Library/Mobile Documents/com~apple~CloudDocs/Sanghee Mun/ Academic/RESEARCH/Staggered_DiD/2024Jan/DDsurv.R')
library(xtable)
results <- results %>%
  mutate(TE = ifelse(G==t, 3, 0),
         TE = ifelse(G==t-1, 2, TE),
         TE = ifelse(G==t-2, 1, TE),
         TE = ifelse(G==t-3, 0.8, TE),
         TE = ifelse(G==t-4, 0.6, TE),
         TE = ifelse(G==t-5, 0.3, TE)) %>%
  select(G,t,TE,ATT,ATT_reweight)

# print the xtable result without row number
print(xtable(results), include.rownames = FALSE)
