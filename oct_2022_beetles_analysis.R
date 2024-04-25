#### Coskun Kucukkaragoz
### 25th of April 2024
## Start of the R script to analyse the data from my mounting, pinning, and sorting, of beetles for the October 2022 season from the Cederburg

# Packages
library(tidyverse)

# Data wrangling

{beet_22_raw <- read.csv("ceder_beetle_diversity_october_2022.csv", stringsAsFactors = TRUE)
  beet_22_raw$ID <- as.factor(beet_22_raw$ID)}
summary(beet_22_raw)

beet_22_no_replicates <- beet_22_raw %>% group_by(ID, site, replicate, trap, label, season, year, full_label) %>% summarise(across(c(anthia_decemguttata,stenocara_dentata,starts_with("morphospecies")), sum)) %>% ungroup()

beet_22_condensed <- beet_22_no_replicates %>% select(-c(site, replicate, trap, season, year, full_label))


