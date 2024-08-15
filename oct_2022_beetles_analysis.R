##### Coskun Kucukkaragoz
#### 25th of April 2024
### Start of the R script to analyse the data from my mounting, pinning, and sorting, of beetles for the October 2022 season from the Cederburg

# set-up ####

# packages
rm(list = ls())
library(tidyverse)

# reading in raw data with correct data types

beet_22_raw <- read.csv("family_seperated_ceder_beetle_diversity_october_2022.csv", stringsAsFactors = TRUE) %>% mutate(across(c(ID, site, replicate, trap, label, season, year, full_label), factor))

beet_02_raw <- read.csv("oct_2002_beetles_coskun_edit.csv", stringsAsFactors = TRUE) %>% mutate(across(c(label,site, replicate, year), factor))

# reducing all unique specimen labels into single row entries

beet_22_condensed <- beet_22_raw %>% group_by(ID, site, replicate, trap, label, season, year, full_label) %>% summarise(across(c(anthia_decemguttata, stenocara_dentata, starts_with("morphospecies")), sum)) %>% ungroup()

# calculating the number of individual specimens found for each morphospecies/species

{x <- mapply(sum,beet_22_condensed[,-c(1:8)])
species_specimen_count <- cbind(read.table(text = names(x)), x)
rownames(species_specimen_count) <- NULL
colnames(species_specimen_count) <- c("morphospecies", "individuals")
rm(x)}

{y <- mapply(sum,beet_02_raw[,-c(1:4)])
  species_02 <- cbind(read.table(text = names(y)), y)
  rownames(species_02) <- NULL
  colnames(species_02) <- c("morphospecies", "individuals")
  rm(y)}

species_specimen_count %>% arrange(individuals) # ordered descendingly
species_specimen_count %>% summarise(individuals = sum(individuals))

species_02 %>%  arrange(individuals)
species_02 %>% summarise(individuals = sum(individuals))

# calculating number of individual insects found per altitudinal site

site_count <- beet_22_condensed %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(c(ID, site, label, individuals)) %>% group_by(site) %>% summarise(individuals = sum(individuals))

abundance_2002 <- beet_02_raw %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(site, individuals) %>% group_by(site) %>% summarise(individuals = sum(individuals))

# calculating species diversity per site

site_diversity <- beet_22_condensed %>% group_by(site) %>% summarise(across(c(anthia_decemguttata, stenocara_dentata, starts_with("morphospecies")), sum)) %>% mutate_if(is.numeric, ~1 * (. > 0)) %>% mutate(diversity = rowSums(across(where(is.numeric)))) %>% select(c(site, diversity))

(diversity_2002 <- beet_02_raw %>% group_by(site) %>% mutate_if(is.numeric, ~1 * (. > 0)) %>% mutate(diversity = rowSums(across(where(is.numeric)))) %>% select(c(site, diversity)) %>% summarise(diversity = sum(diversity)))

# graphs ####

alt_labels <- c("0m W","200m W","300m W","500m W","700m W","900m W","1100m W","1300m W","1500m W","1700m W","1900m W","1700m E","1500m E","1300m E","1100m E","900m E","500m E")

# species diversity per altitudinal site

ggplot(site_diversity, aes(x = site, y = diversity, group = 1, colour = "black")) + geom_point() + geom_line() + xlab("Altitudinal site") + ylab("Total number of morphospecies") + scale_x_discrete(label = alt_labels) + scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) + theme(axis.text = element_text(size = 8)) + geom_point(data = diversity_2002, aes(colour = "grey")) + geom_line(data = diversity_2002, aes(colour = "grey")) + scale_colour_manual(name = NULL, values =c("black"="black","grey"="grey"), labels = c("2022","2002")) + theme(legend.position = c(.9, .9))

# individual insects per altitudinal site

ggplot(site_count, aes (x = site, y = individuals, group = 1, colour = "black")) + geom_point() + geom_line() + xlab("Altitudinal site") + ylab("Total number of individual beetles") + scale_x_discrete(label = alt_labels) + scale_y_continuous(expand = c(0, 0), limits = c(0, 450)) + theme(axis.text = element_text(size = 8)) + geom_point(data = abundance_2002, aes(colour = "grey")) + geom_line(data = abundance_2002, aes(colour = "grey")) + scale_colour_manual(name = NULL, values =c("grey"="grey", "black"="black"), labels = c("2022","2002")) + theme(legend.position = c(.9, .9))

# ICSZ specific analyses ####

reduced_2022 <- beet_22_raw %>% filter(site %in% c("8","9","10")) %>% group_by(ID, site, replicate, trap, label, season, year, full_label) %>% summarise(across(c(anthia_decemguttata, stenocara_dentata, starts_with("morphospecies")), sum)) %>% ungroup()

reduced_2002 <- read.csv("oct_2002_beetles_coskun_edit.csv", stringsAsFactors = TRUE) %>% mutate(across(c(label,site, replicate, year), factor)) %>% filter(site %in% c("8","9","10"))

{x <- mapply(sum,reduced_2022[,-c(1:8)])
  red_species_specimen_count <- cbind(read.table(text = names(x)), x)
  rownames(red_species_specimen_count) <- NULL
  colnames(red_species_specimen_count) <- c("morphospecies", "individuals")
  rm(x)}

{y <- mapply(sum,reduced_2002[,-c(1:4)])
  red_species_02 <- cbind(read.table(text = names(y)), y)
  rownames(red_species_02) <- NULL
  colnames(red_species_02) <- c("morphospecies", "individuals")
  rm(y)}


red_species_specimen_count %>% filter(individuals > 0) %>% arrange(individuals) # ordered descendingly
red_species_specimen_count %>% summarise(individuals = sum(individuals))

red_species_02 %>% filter(individuals > 0) %>% arrange(individuals)
red_species_02 %>% summarise(individuals = sum(individuals))

# selecting the three primary species to conduct hypothesis tests

(hyp_red_2022 <- reduced_2022 %>% rename(zophosis_gracilicornis = morphospecies_6) %>% select(site,replicate,year,anthia_decemguttata,stenocara_dentata,zophosis_gracilicornis) %>% group_by(site,replicate,year) %>% summarise(across(c(anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis), sum)) %>% ungroup())

(hyp_red_2002 <- reduced_2002 %>% rename(zophosis_gracilicornis = Zophosis.sp.1, anthia_decemguttata = Thermophilum.decemguttatum, stenocara_dentata = Stenocara.dentata) %>% select(site,replicate,year,anthia_decemguttata,stenocara_dentata,zophosis_gracilicornis))

hyp_red_all <- rbind(hyp_red_2002,hyp_red_2022)

t.test(x = hyp_red_2002$anthia_decemguttata, y = hyp_red_2022$anthia_decemguttata, alternative = "two.sided", paired = TRUE)

t.test(x = hyp_red_2002$stenocara_dentata, y = hyp_red_2022$stenocara_dentata, alternative = "two.sided", paired = TRUE)

t.test(x = hyp_red_2002$zophosis_gracilicornis, y = hyp_red_2022$zophosis_gracilicornis, alternative = "two.sided", paired = TRUE)

(1 - mean(hyp_red_all$zophosis_gracilicornis[hyp_red_all$year=="2022"])/mean(hyp_red_all$zophosis_gracilicornis[hyp_red_all$year=="2002"]))*100

# climate analyses begin here

for (k in 1:length(unique(beet_22_raw$site))){
  dat <-  read.csv(file = paste0("./climate_data/climate_site_", k , ".csv"))
  assign(paste0("temp_",k), dat)
  {if (k==length(unique(beet_22_raw$site))){rm(dat, k)}}
}





# end of ICSZ analyses



# sorting working data ####
# used to sort morphospecies into boxes

# species_specimen_count %>% filter(individuals < 2) %>%  arrange(individuals)

# edit the morphospecies numbers and sites for the specimens of interest

# beet_22_condensed %>% filter(trap != 3) %>% group_by(site) %>% 
#   summarise(individuals = sum(morphospecies_12)) %>%
#   filter(individuals > 0)

# individual site samples

# (species_by_site <- beet_22_condensed %>% select(ID, site, label, morphospecies_69) %>%
#   filter(morphospecies_69 > 0))

