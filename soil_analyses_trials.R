##### Coskun Kucukkaragoz
#### 25th of April 2024
### R script to analyse soil data

# set-up ####

# packages
rm(list = ls())
library(tidyverse)
library(weights)
library(emmeans)
library(EnvStats)
library(lubridate)
library(NADA)
library(NADA2)

#### analyses ####

soil <- read.csv("./soil/soil_data_2022_censored.csv")
soil <- soil %>% mutate(across(c(site, replicate), factor))
soil <- soil %>% mutate(across(c(Phosphorous_censored, Potassium_censored, Sulfur_censored, NO3_censored), as.logical))
str(soil)

equivalent_n(soil$Phosphorus, soil$Phosphorous_censored)
cenboxplot(soil$Phosphorus, soil$Phosphorous_censored, soil$site)


m0 <- cenfit(Cen(soil$Phosphorus, soil$Phosphorous_censored))
m1 <- cenfit(Cen(soil$Phosphorus, soil$Phosphorous_censored)~soil$site)
censummary(soil$Phosphorus, soil$Phosphorous_censored, soil$site)



plot(m1)
summary(m1)

summary(m0)
