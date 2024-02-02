# install.packages("devtools")
devtools::install_github("jonathandroth/staggered")
library(staggered)
library(ggplot2)
#Calculate Callaway and Sant'Anna estimator for the simple weighted average
staggered_cs(df = df, 
             i = "unit",
             t = "t",
             g = "G",
             y = "Y", 
             estimand = "eventstudy")
eventPlotResults <- staggered_cs(df = df, 
                                 i = "unit",
                                 t = "t",
                                 g = "G",
                                 y = "Y", 
                                 estimand = "eventstudy",eventTime = 0:8)


list(eventPlotResults)