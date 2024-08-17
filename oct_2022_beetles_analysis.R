##### Coskun Kucukkaragoz
#### 25th of April 2024
### Start of the R script to analyse the data from my mounting, pinning, and sorting, of beetles for the October 2022 season from the Cederburg

# set-up ####

# packages
rm(list = ls())
library(tidyverse)
library(weights)

# reading in raw data with correct data types

beet_22_raw <- read.csv("family_seperated_ceder_beetle_diversity_october_2022.csv", stringsAsFactors = TRUE) %>% mutate(across(c(ID, site, replicate, trap, label, season, year, full_label), factor))

beet_02_raw <- read.csv("oct_2002_beetles_coskun_edit.csv", stringsAsFactors = TRUE) %>% mutate(across(c(label,site, replicate, year), factor))

# calculating the number of individual specimens found for each morphospecies/species

{x <- mapply(sum, beet_22_raw[,-c(1:8)])
species_22 <- cbind(read.table(text = names(x)), x)
rownames(species_22) <- NULL
colnames(species_22) <- c("morphospecies", "individuals")
rm(x)}

{x <- mapply(sum,beet_02_raw[,-c(1:4)])
  species_02 <- cbind(read.table(text = names(x)), x)
  rownames(species_02) <- NULL
  colnames(species_02) <- c("morphospecies", "individuals")
  rm(x)}

species_22 %>% arrange(individuals) # ordered descendingly
species_22 %>% summarise(individuals = sum(individuals))

species_02 %>%  arrange(individuals)
species_02 %>% summarise(individuals = sum(individuals))

# calculating number of individual insects found per altitudinal site

abundance_2022 <- beet_22_raw %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(c(ID, year, site, label, individuals)) %>% group_by(year, site) %>% summarise(individuals = sum(individuals), .groups = "drop")

abundance_2002 <- beet_02_raw %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(year, site, individuals) %>% group_by(year, site) %>% summarise(individuals = sum(individuals), .groups = "drop")

# calculating species diversity per site

diversity_2022 <- beet_22_raw %>% group_by(site) %>% summarise(across(where(is.numeric), sum)) %>% mutate_if(is.numeric, ~1 * (. > 0)) %>% mutate(diversity = rowSums(across(where(is.numeric)))) %>% select(c(site, diversity))

diversity_2002 <- beet_02_raw %>% group_by(site) %>% summarise(across(where(is.numeric), sum)) %>% mutate_if(is.numeric, ~1 * (. > 0)) %>% mutate(diversity = rowSums(across(where(is.numeric)))) %>% select(c(site, diversity))

# graphs ####

alt_labels <- c("0m W","200m W","300m W","500m W","700m W","900m W","1100m W","1300m W","1500m W","1700m W","1900m W","1700m E","1500m E","1300m E","1100m E","900m E","500m E")

# species diversity per altitudinal site

ggplot(diversity_2002, aes(x = site, y = diversity, group = 1, colour = "grey")) + geom_point() + geom_line() + xlab("Altitudinal site") + ylab("Total number of morphospecies") + scale_x_discrete(label = alt_labels) + scale_y_continuous(expand = c(0, 0), limits = c(0, 40)) + theme(axis.text = element_text(size = 8)) + geom_point(data = diversity_2022, aes(colour = "black")) + geom_line(data = diversity_2022, aes(colour = "black")) + scale_colour_manual(name = NULL, values =c("black"="black", "grey"="grey"), labels = c("2022","2002")) + theme(legend.position = "inside", legend.position.inside = c(.9, .9))

# individual insects per altitudinal site

ggplot(abundance_2002, aes (x = site, y = individuals, group = 1, colour = "grey")) + geom_point() + geom_line() + xlab("Altitudinal site") + ylab("Total number of individual beetles") + scale_x_discrete(label = alt_labels) + scale_y_continuous(expand = c(0, 0), limits = c(0, 450)) + theme(axis.text = element_text(size = 8)) + geom_point(data = abundance_2022, aes(colour = "black")) + geom_line(data = abundance_2022, aes(colour = "black")) + scale_colour_manual(name = NULL, values =c("black"="black", "grey"="grey"), labels = c("2022","2002"))  + theme(legend.position = "inside", legend.position.inside = c(.9, .9))

# ICSZ specific analyses ####

# poisson models for abundance

{pois_2022 <- beet_22_raw %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(c(ID, year, site, replicate, label, individuals)) %>% group_by(year, site, replicate) %>% summarise(individuals = sum(individuals), .groups = "drop")

pois_2002 <- beet_02_raw %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(year, site, replicate, individuals) %>% group_by(year, site, replicate) %>% summarise(individuals = sum(individuals), .groups = "drop")

pois_all <- rbind(pois_2002,pois_2022)
  
  vec <- c("1" = 0, "2" = 200, "3" = 300, "4" = 500, "5" = 700, "6" = 900, "7" = 1100, "8" = 1300, "9" = 1500, "10" = 1700, "11" = 1900, "12" = 1700, "13" = 1500, "14" = 1300, "15" = 1100, "16" = 900, "17" = 500)
  
  pois_all <- left_join(pois_all, enframe(vec), by = c("site" = "name")) %>% rename(altitude = value) %>% mutate(site = as.numeric(site)) %>% mutate(site = as.factor(site)) %>% mutate(year = as.numeric(as.character(year)))
  rm(vec, pois_2002, pois_2022)}

m_0 <- glm(individuals ~ year, data = pois_all, family = poisson)

for (k in 1:length(unique(pois_all$site))){
  m <- glm(individuals[site == k] ~ year[site == k], data = pois_all, family = poisson)
  assign(paste0("m_", k), m)
  rm(m)
}

summary(m_0) # *** --
summary(m_1) # *** --
summary(m_2) # ** ++
summary(m_3) # * ++
summary(m_4) # x
summary(m_5) # ** --
summary(m_6) # x
summary(m_7) # x
summary(m_8) # *** --
summary(m_9) # *** --
summary(m_10) # *** --
summary(m_11) # *** ++
summary(m_12) # x
summary(m_13) # x
summary(m_14) # *** ++
summary(m_15) # *** ++
summary(m_16) # x
summary(m_17) # x

predict(m_0, newdata = data.frame(year = "2022"), type = "response")

# reducing data to most affected sites

reduced_2022 <- beet_22_raw %>% filter(site %in% c("8","9","10"))

reduced_2002 <- beet_02_raw %>% filter(site %in% c("8","9","10"))

{x <- mapply(sum,reduced_2022[,-c(1:8)])
  red_species_22 <- cbind(read.table(text = names(x)), x)
  rownames(red_species_22) <- NULL
  colnames(red_species_22) <- c("morphospecies", "individuals")
  rm(x)}

{x <- mapply(sum,reduced_2002[,-c(1:4)])
  red_species_02 <- cbind(read.table(text = names(x)), x)
  rownames(red_species_02) <- NULL
  colnames(red_species_02) <- c("morphospecies", "individuals")
  rm(x)}


red_species_22 %>% filter(individuals > 0) %>% arrange(individuals) # ordered descendingly
red_species_22 %>% summarise(individuals = sum(individuals))

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


# altitudinal range shifts

{shifts_22 <- beet_22_raw %>% rename(zophosis_gracilicornis = morphospecies_6) %>% select(site,replicate,year,anthia_decemguttata,stenocara_dentata,zophosis_gracilicornis) %>% group_by(site,replicate,year) %>% summarise(across(c(anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis), sum), .groups = "drop")

shifts_02 <- beet_02_raw %>% rename(zophosis_gracilicornis = Zophosis.sp.1, anthia_decemguttata = Thermophilum.decemguttatum, stenocara_dentata = Stenocara.dentata) %>% select(site,replicate,year,anthia_decemguttata,stenocara_dentata,zophosis_gracilicornis)

shifts_all <- rbind(shifts_02,shifts_22)

vec <- c("1" = 0, "2" = 200, "3" = 300, "4" = 500, "5" = 700, "6" = 900, "7" = 1100, "8" = 1300, "9" = 1500, "10" = 1700, "11" = 1900, "12" = 1700, "13" = 1500, "14" = 1300, "15" = 1100, "16" = 900, "17" = 500)

shifts_all <- left_join(shifts_all, enframe(vec), by = c("site" = "name")) %>% rename(altitude = value) %>% mutate(site = as.numeric(site)) %>% mutate(site = as.factor(site))
rm(shifts_02, shifts_22, vec)}

shifts_all %>% group_by(year, site) %>%  summarise(anthia_decemguttata = sum(anthia_decemguttata), .groups = "drop") %>% ggplot(aes (x = site, y = anthia_decemguttata, colour = year, group = year)) + geom_point() + geom_line() + xlab("Altitudinal site") + ylab("Total number of individual beetles") + scale_x_discrete(label = alt_labels) + theme(axis.text = element_text(size = 8))

# shifts_all %>% group_by(year, site) %>%  summarise(stenocara_dentata = sum(stenocara_dentata), .groups = "drop") %>% ggplot(aes (x = site, y = stenocara_dentata, colour = year, group = year)) + geom_point() + geom_line() + xlab("Altitudinal site") + ylab("Total number of individual beetles") + scale_x_discrete(label = alt_labels) + theme(axis.text = element_text(size = 8))
# 
# shifts_all %>% group_by(year, site) %>%  summarise(zophosis_gracilicornis = sum(zophosis_gracilicornis), .groups = "drop") %>% ggplot(aes (x = site, y = zophosis_gracilicornis, colour = year, group = year)) + geom_point() + geom_line() + xlab("Altitudinal site") + ylab("Total number of individual beetles") + scale_x_discrete(label = alt_labels) + theme(axis.text = element_text(size = 8))
# 
# xx <- shifts_all %>% group_by(year, altitude) %>% summarise(across(c(anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis), sum), .groups = "drop") %>% select(year, altitude, anthia_decemguttata) %>% rename(abundance = anthia_decemguttata)
# 
# xx <- shifts_all %>% group_by(year, altitude) %>% summarise(across(c(anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis), sum), .groups = "drop") %>% select(year, altitude, stenocara_dentata) %>% rename(abundance = stenocara_dentata)
# 
# xx <- shifts_all %>% group_by(year, altitude) %>% summarise(across(c(anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis), sum), .groups = "drop") %>% select(year, altitude, zophosis_gracilicornis) %>% rename(abundance = zophosis_gracilicornis)
# 
# wtd.t.test(xx$altitude[xx$year==2002], xx$altitude[xx$year==2002], weight = xx$abundance[xx$year==2002], weighty = xx$abundance[xx$year==2022])

m_a_0 <- glm(anthia_decemguttata[year==2022] ~ altitude[year==2022], data = shifts_all, family = poisson)
summary(m_a_0)

m_s_0 <- glm(stenocara_dentata ~ altitude*year, data = shifts_all, family = poisson)
summary(m_s_0)

m_z_0 <- glm(zophosis_gracilicornis[year==2002] ~ altitude[year==2002], data = shifts_all, family = poisson)
summary(m_z_0)

predict(m_a_0, newdata = data.frame(altitude = 2000), type = "response")

# {shifts_raw <- rbind(abundance_2002,abundance_2022)
# 
# vec <- c("1" = 0, "2" = 200, "3" = 300, "4" = 500, "5" = 700, "6" = 900, "7" = 1100, "8" = 1300, "9" = 1500, "10" = 1700, "11" = 1900, "12" = 1700, "13" = 1500, "14" = 1300, "15" = 1100, "16" = 900, "17" = 500)
# 
# shifts_raw <- left_join(shifts_raw, enframe(vec), by = c("site" = "name")) %>% rename(altitude = value) %>% mutate(site = as.numeric(site)) %>% mutate(site = as.factor(site))
# rm(vec)}
# 
# 
# wtd.t.test(shifts_raw$altitude[shifts_raw$year==2002], shifts_raw$altitude[shifts_raw$year==2022], weight = shifts_raw$individuals[shifts_raw$year==2002], weighty = shifts_raw$individuals[shifts_raw$year==2022])

# more hypothesis tests

fit <- aov()

# climate analyses begin here

for (k in 1:length(unique(beet_22_raw$site))){
  dat <-  read.csv(file = paste0("./climate_data/climate_site_", k , ".csv"))
  dat <- data.frame(site = rep(k, nrow(dat)), dat)
  dat <- dat %>% rename(year = Year, month = Month, day = Day, time = Time, temperature = Temp1) %>%
    mutate(across(c(site, month, day, time), factor))
  {if (k == 1){temp_sites <- dat}}
  {if (k > 1){temp_sites <- (bind_rows(temp_sites, dat))}}
  # assign(paste0("temp_", k), dat)
  # {if (k == 1){common_cols <- names(get(paste0("temp_", k)))}}
  # {if (k > 1){common_cols <- intersect(common_cols, names(get(paste0("temp_", k))))}}
  {if (k == length(unique(beet_22_raw$site))){temp_sites <- temp_sites %>% mutate(across(c(month, day), ordered))}}
  {if (k == length(unique(beet_22_raw$site))){rm(dat, k)}}
}


datetest <- temp_sites %>% unite(datetime, c(year, month,day), sep = "-", remove = FALSE) %>% unite(datetime, c(datetime,time), sep = " ", remove = TRUE) %>% mutate(datetime = as_datetime(datetime, format="%Y-%m-%d %H:%M", tz = "Africa/Johannesburg"))

datetest %>% filter(year < 2022, site == 8) %>% ggplot(aes(x = datetime, y = temperature)) + geom_smooth(aes(datetime, colour = "1"), formula = y ~ x, method = lm, se = T)

datetest_filt <- datetest %>% filter(year < 2022, site == 8)
m1 <- lm(temperature ~ datetime, data = datetest_filt)
summary(m1)


# end of ICSZ analyses








# sorting working data ####
# used to sort morphospecies into boxes

# species_22 %>% filter(individuals < 2) %>%  arrange(individuals)

# edit the morphospecies numbers and sites for the specimens of interest

# beet_22_raw %>% filter(trap != 3) %>% group_by(site) %>% 
#   summarise(individuals = sum(morphospecies_12)) %>%
#   filter(individuals > 0)

# individual site samples

# (species_by_site <- beet_22_raw %>% select(ID, site, label, morphospecies_69) %>%
#   filter(morphospecies_69 > 0))

