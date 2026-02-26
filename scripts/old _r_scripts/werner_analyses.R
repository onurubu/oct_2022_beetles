##### Coskun Kucukkaragoz
#### 25th of April 2024
### Start of the R script to analyse the data from my mounting, pinning, and sorting, of beetles for the October 2022 season from the Cederburg

# set-up ####

# packages
rm(list = ls())
library(tidyverse)
library(weights)
library(emmeans)
library(EnvStats)
library(sf)
library(hms)
library(vegan)

# reading in raw data with correct data types

{abundance_2022 <- read.csv("werner_ID_with_3_species.csv", stringsAsFactors = TRUE) %>% mutate(across(c(year, site, replicate), factor)) %>% filter(year == "2022")

beet_02_raw <- read.csv("oct_2002_beetles_coskun_edit.csv", stringsAsFactors = TRUE) %>% mutate(across(c(label,site, replicate, year), factor)) %>% filter(year == "2002")

abundance_2002 <- beet_02_raw %>% select(-label) %>% select(year, site, replicate, everything())

rm(beet_02_raw)}

# Analyses start here

# beta diversity


factors_2022 <- abundance_2022[,1:3]
values_2022 <- abundance_2022[,4:45]

dist_2022 <- vegdist(values_2022)
ano_2022 <- with(factors_2022, anosim(dist_2022, site))
summary(ano_2022)

## nmds trials

# switch to presence-absence data and make it a matrix
species_data <- values_2022 %>% decostand(.,"pa") %>% as.matrix(.)
