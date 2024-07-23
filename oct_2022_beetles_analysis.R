##### Coskun Kucukkaragoz
#### 25th of April 2024
### Start of the R script to analyse the data from my mounting, pinning, and sorting, of beetles for the October 2022 season from the Cederburg

# set-up ####

# packages
rm(list = ls())
library(tidyverse)

# reading in raw data with correct data types

beet_22_raw <- read.csv("ceder_beetle_diversity_october_2022.csv", stringsAsFactors = TRUE) %>% mutate(across(c(ID, site, replicate, trap, label, season, year, full_label), factor))

# reducing all unique specimen labels into single row entries

beet_22_condensed <- beet_22_raw %>% group_by(ID, site, replicate, trap, label, season, year, full_label) %>% summarise(across(c(anthia_decemguttata, stenocara_dentata, starts_with("morphospecies")), sum)) %>% ungroup()

# calculating the number of individual specimens found for each morphospecies/species

{x <- mapply(sum,beet_22_condensed[,-c(1:8)])
species_specimen_count <- cbind(read.table(text = names(x)), x)
rownames(species_specimen_count) <- NULL
colnames(species_specimen_count) <- c("morphospecies", "individuals")
rm(x)}

species_specimen_count %>% arrange(desc(individuals)) # ordered descendingly
species_specimen_count %>% summarise(individuals = sum(individuals))

# calculating number of individual insects found per altitudinal site

site_count <- beet_22_condensed %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(c(ID, site, label, individuals)) %>% group_by(site) %>% summarise(individuals = sum(individuals))

# calculating species diversity per site

site_diversity <- beet_22_condensed %>% group_by(site) %>% summarise(across(c(anthia_decemguttata, stenocara_dentata, starts_with("morphospecies")), sum)) %>% mutate_if(is.numeric, ~1 * (. > 0)) %>% mutate(diversity = rowSums(across(where(is.numeric)))) %>% select(c(site, diversity))

# graphs ####

# species diversity per altitudinal site

ggplot(site_diversity, aes(x = site, y = diversity, group = 1)) + geom_point() + geom_line()

# individual insects per altitudinal site

ggplot(site_count, aes (x = site, y = individuals, group = 1)) + geom_point() + geom_line()

# working data ####

species_specimen_count %>% arrange(individuals)

# edit the morphospecies numbers and sites for the specimens of interest

beet_22_condensed %>% group_by(site) %>% summarise(individuals = sum(morphospecies_135))

# individual site samples

(species_by_site <- beet_22_raw %>% select(ID, site, label, morphospecies_135) %>% filter(site == 1, morphospecies_135 > 0))

