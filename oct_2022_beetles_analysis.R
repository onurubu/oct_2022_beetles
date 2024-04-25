#### Coskun Kucukkaragoz
### 25th of April 2024
## Start of the R script to analyse the data from my mounting, pinning, and sorting, of beetles for the October 2022 season from the Cederburg

# packages
library(tidyverse)

# reading in raw data with correct data types

{beet_22_raw <- read.csv("ceder_beetle_diversity_october_2022.csv", stringsAsFactors = TRUE)
  beet_22_raw$ID <- as.factor(beet_22_raw$ID)}
summary(beet_22_raw)

# reducing all unique specimen labels into single row entries

beet_22_no_replicates <- beet_22_raw %>% group_by(ID, site, replicate, trap, label, season, year, full_label) %>% summarise(across(c(anthia_decemguttata,stenocara_dentata,starts_with("morphospecies")), sum)) %>% ungroup()

# creating a data frame excluding currently unnecesary columns

beet_22_condensed <- beet_22_no_replicates %>% select(-c(site, replicate, trap, season, year, full_label))

# calculating the number of individual specimens found for each morphospecies/species

{x <- mapply(sum,beet_22_condensed[,-c(1:2)])
species_total <- cbind(read.table(text = names(x)), x)
rownames(species_total) <- NULL
colnames(species_total) <- c("morphospecies", "individuals")
rm(x)}
