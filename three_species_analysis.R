##### Coskun Kucukkaragoz
#### 25th of April 2024
### R script to analyse the data for the beetles from the Cederburg

# set-up ####

# packages
rm(list = ls())
library(tidyverse)
library(weights)
library(emmeans)
library(EnvStats)

# analyses ####

# reading in raw data with correct data types
{
beet_22_fullraw <- read.csv("family_separated_ceder_beetles_2022.csv", stringsAsFactors = TRUE) %>% mutate(across(c(ID, site, replicate, trap, label, season, year, full_label), factor)) %>% filter(year == 2022)

beet_23_fullraw <- read.csv("oct_2023_beetle_abundance.csv", stringsAsFactors = TRUE) %>% mutate(across(c(ID, site, replicate, trap, year, label), factor))

beet_02_raw <- read.csv("oct_2002_beetles_coskun_edit.csv", stringsAsFactors = TRUE) %>% mutate(across(c(label,site, replicate, year), factor)) %>% filter(year == 2002)

beet_03_raw <- read.csv("oct_2002_beetles_coskun_edit.csv", stringsAsFactors = TRUE) %>% mutate(across(c(label,site, replicate, year), factor)) %>% filter(year == 2003)

# reducing all entries of modern samples to 1 row per unique ID and also filling in empty IDs

beet_22_fullraw <- beet_22_fullraw %>% group_by(ID, site, replicate, trap, label, season, year, full_label) %>% summarise(across(where(is.numeric), sum), .groups = "drop")

# section to try and create a function to fill in the gaps where no instances were found

{
  missID <- symdiff(as.numeric(as.character(beet_22_fullraw$ID)), 1:680)
for (k in 1:length(missID)){
  i <- missID[k]
  s <- floor((i-1)/40) + 1
  r <- (ceiling((i-(floor((i-1)/40)*40))/10))
  t <- i-(((s-1)*40)+((r-1)*10))
  {if (k == 1){i_list <- i}}
  {if (k == 1){s_list <- s}}
  {if (k == 1){r_list <- r}}
  {if (k == 1){t_list <- t}}
  {if (k == 1){label_list <- paste0(s, ".", r, ".", t)}}
  {if (k > 1){i_list <- append(i_list, i)}}
  {if (k > 1){s_list <- append(s_list, s)}}
  {if (k > 1){r_list <- append(r_list, r)}}
  {if (k > 1){t_list <- append(t_list, t)}}
  {if (k > 1){label_list <- append(label_list, paste0(s, ".", r, ".", t))}}
  {if (k == length(missID)){miss_labels <- data.frame(ID = as.factor(i_list), site = as.factor(s_list), replicate = as.factor(r_list), trap = as.factor(t_list), label = as.factor(label_list), "year" = as.factor(rep(2022, length(missID))), "anthia_decemguttata" = rep(0, length(missID)), "stenocara_dentata" = rep(0, length(missID)), "morphospecies_6" = rep(0, length(missID)), "morphospecies_7" = rep(0, length(missID)), "morphospecies_12" = rep(0, length(missID)), "morphospecies_13" = rep(0, length(missID)), "morphospecies_18" = rep(0, length(missID)), 
"morphospecies_29" = rep(0, length(missID)), "morphospecies_34" = rep(0, length(missID)), "morphospecies_44" = rep(0, length(missID)), "morphospecies_48" = rep(0, length(missID)), 
"morphospecies_52" = rep(0, length(missID)), "morphospecies_55" = rep(0, length(missID)), "morphospecies_56" = rep(0, length(missID)), "morphospecies_63" = rep(0, length(missID)), 
"morphospecies_69" = rep(0, length(missID)), "morphospecies_70" = rep(0, length(missID)), "morphospecies_72" = rep(0, length(missID)), "morphospecies_73" = rep(0, length(missID)), 
"morphospecies_74" = rep(0, length(missID)), "morphospecies_76" = rep(0, length(missID)), "morphospecies_77" = rep(0, length(missID)), "morphospecies_78" = rep(0, length(missID)), 
"morphospecies_81" = rep(0, length(missID)), "morphospecies_82" = rep(0, length(missID)), "morphospecies_86" = rep(0, length(missID)), "morphospecies_88" = rep(0, length(missID)), 
"morphospecies_104" = rep(0, length(missID)), "morphospecies_108" = rep(0, length(missID)), "morphospecies_114" = rep(0, length(missID)), 
"morphospecies_118" = rep(0, length(missID)), "morphospecies_124" = rep(0, length(missID)), "morphospecies_125" = rep(0, length(missID)), 
"morphospecies_129" = rep(0, length(missID)), "morphospecies_133" = rep(0, length(missID)), "morphospecies_135" = rep(0, length(missID)), 
"morphospecies_136" = rep(0, length(missID)), "morphospecies_140" = rep(0, length(missID)), "morphospecies_141" = rep(0, length(missID)), 
"morphospecies_142" = rep(0, length(missID)), "morphospecies_146" = rep(0, length(missID)), "morphospecies_155" = rep(0, length(missID)), 
"morphospecies_156" = rep(0, length(missID)))}}
  {if (k == length(missID)){rm(i, r, r_list, s, s_list, t, t_list, i_list, label_list)}}
  {if (k == length(missID)){beet_22_raw <- beet_22_fullraw %>% select(-season, -full_label) %>% bind_rows(miss_labels) %>% select(-morphospecies_146, -morphospecies_155, -morphospecies_156) %>% arrange(as.numeric(as.character(ID)))}}
  {if (k == length(missID)){rm(miss_labels, missID, k, beet_22_fullraw)}}
}}
{
  missID <- symdiff(as.numeric(as.character(beet_23_fullraw$ID)), 1:680)
  for (k in 1:length(missID)){
    i <- missID[k]
    s <- floor((i-1)/40) + 1
    r <- (ceiling((i-(floor((i-1)/40)*40))/10))
    t <- i-(((s-1)*40)+((r-1)*10))
    {if (k == 1){i_list <- i}}
    {if (k == 1){s_list <- s}}
    {if (k == 1){r_list <- r}}
    {if (k == 1){t_list <- t}}
    {if (k == 1){label_list <- paste0(s, ".", r, ".", t)}}
    {if (k > 1){i_list <- append(i_list, i)}}
    {if (k > 1){s_list <- append(s_list, s)}}
    {if (k > 1){r_list <- append(r_list, r)}}
    {if (k > 1){t_list <- append(t_list, t)}}
    {if (k > 1){label_list <- append(label_list, paste0(s, ".", r, ".", t))}}
    {if (k == length(missID)){miss_labels <- data.frame(ID = as.factor(i_list), site = as.factor(s_list), replicate = as.factor(r_list), trap = as.factor(t_list), label = as.factor(label_list), "year" = as.factor(rep(2023, length(missID))), "anthia_decemguttata" = rep(0, length(missID)), "stenocara_dentata" = rep(0, length(missID)), "morphospecies_6" = rep(0, length(missID)))}}
    {if (k == length(missID)){rm(i, r, r_list, s, s_list, t, t_list, i_list, label_list)}}
    {if (k == length(missID)){beet_23_raw <- beet_23_fullraw %>% bind_rows(miss_labels) %>% arrange(as.numeric(as.character(ID)))}}
    {if (k == length(missID)){rm(miss_labels, missID, k, beet_23_fullraw)}}
  }}
  beet_02_raw <- beet_02_raw %>% select(1:4, Thermophilum.decemguttatum, Stenocara.dentata, Zophosis.sp.1) %>% rename(anthia_decemguttata = Thermophilum.decemguttatum, stenocara_dentata = Stenocara.dentata, zophosis_gracilicornis = Zophosis.sp.1)
  beet_03_raw <- beet_03_raw %>% select(1:4, Thermophilum.decemguttatum, Stenocara.dentata, Zophosis.sp.1) %>% rename(anthia_decemguttata = Thermophilum.decemguttatum, stenocara_dentata = Stenocara.dentata, zophosis_gracilicornis = Zophosis.sp.1)
  beet_22_raw <- beet_22_raw %>% select(1:6, anthia_decemguttata, stenocara_dentata, morphospecies_6) %>% rename(zophosis_gracilicornis = morphospecies_6)
  beet_23_raw <- beet_23_raw %>% select(1:6, anthia_decemguttata, stenocara_dentata, morphospecies_6) %>% rename(zophosis_gracilicornis = morphospecies_6)
  }

# calculating the number of individual specimens found for each species

{x <- beet_22_raw %>% filter(year == "2022")
  x <- mapply(sum, x[,-c(1:6)])
species_22 <- cbind(read.table(text = names(x)), x)
rownames(species_22) <- NULL
colnames(species_22) <- c("morphospecies", "individuals")
rm(x)}

{x <- beet_23_raw
  x <- mapply(sum, x[,-c(1:6)])
  species_23 <- cbind(read.table(text = names(x)), x)
  rownames(species_23) <- NULL
  colnames(species_23) <- c("morphospecies", "individuals")
  rm(x)}

{x <- beet_02_raw %>% filter(year == "2002")
  x <- mapply(sum, x[,-c(1:4)])
  species_02 <- cbind(read.table(text = names(x)), x)
  rownames(species_02) <- NULL
  colnames(species_02) <- c("morphospecies", "individuals")
  rm(x)}

{x <- beet_03_raw %>% filter(year == "2003")
  x <- mapply(sum, x[,-c(1:4)])
  species_03 <- cbind(read.table(text = names(x)), x)
  rownames(species_03) <- NULL
  colnames(species_03) <- c("morphospecies", "individuals")
  rm(x)}

species_22 %>% arrange(individuals)
species_22 %>% summarise(individuals = sum(individuals))

species_23 %>% arrange(individuals)
species_23 %>% summarise(individuals = sum(individuals))

species_02 %>%  arrange(individuals)
species_02 %>% summarise(individuals = sum(individuals)) # 21% loss in abundance overall

species_03 %>%  arrange(individuals)
species_03 %>% summarise(individuals = sum(individuals))

# calculating number of individual insects found per altitudinal site

{abundance_2022 <- beet_22_raw %>% filter(year == "2022") %>% mutate(individuals = rowSums(across(where(is.numeric)))) %>% select(c(ID, year, site, label, individuals)) %>% group_by(year, site) %>% summarise(individuals = sum(individuals), .groups = "drop")
  
abundance_2023 <- beet_23_raw %>% filter(year == "2023") %>% mutate(individuals = rowSums(across(where(is.numeric)))) %>% select(c(ID, year, site, label, individuals)) %>% group_by(year, site) %>% summarise(individuals = sum(individuals), .groups = "drop")

abundance_2002 <- beet_02_raw %>% filter(year == "2002") %>% mutate(individuals = rowSums(across(where(is.numeric)))) %>% select(year, site, individuals) %>% group_by(year, site) %>% summarise(individuals = sum(individuals), .groups = "drop")

abundance_2003 <- beet_03_raw %>% filter(year == "2003") %>% mutate(individuals = rowSums(across(where(is.numeric)))) %>% select(year, site, individuals) %>% group_by(year, site) %>% summarise(individuals = sum(individuals), .groups = "drop")

abundance_all <- rbind(abundance_2002, abundance_2003, abundance_2022, abundance_2023)}

# calculating mean abundance

m_abund_2022 <- beet_22_raw %>% select(1:9) %>% mutate(individuals = rowSums(across(where(is.numeric)))) %>% group_by(year, site, replicate) %>% summarise(individuals = sum(individuals), .groups = "drop")

m_abund_2023 <- beet_23_raw %>% select(1:9) %>% mutate(individuals = rowSums(across(where(is.numeric)))) %>% group_by(year, site, replicate) %>% summarise(individuals = sum(individuals), .groups = "drop")

m_abund_2002 <- beet_02_raw %>% select(1:4, anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis) %>% mutate(individuals = rowSums(across(where(is.numeric)))) %>% select(year, site, replicate, individuals)

m_abund_2003 <- beet_03_raw %>% select(1:4, anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis) %>% mutate(individuals = rowSums(across(where(is.numeric)))) %>% select(year, site, replicate, individuals)

### Trying models
  
m_abund_all <- rbind(m_abund_2002, m_abund_2022, m_abund_2003, m_abund_2023)

m_ave <- m_abund_all %>% lm(data = ., individuals ~ year + site)
summary(m_ave)

predict(m_ave, newdata = data.frame(year = c("2002")))

m_abund_all %>% plot(data = ., individuals[year == 2022] ~ site[year == 2022])

alt_labels <- c("0m(W)","200m(W)","300m(W)","500m(W)","700m(W)","900m(W)","1100m(W)","1300m(W)","1500m(W)","1700m(W)","1900m(W)","1700m(E)","1500m(E)","1300m(E)","1100m(E)","900m(E)","500m(E)")

m_abund_all %>%  ggplot(aes(x=year, y=individuals, fill=year)) + geom_boxplot() + scale_fill_manual(name = "Year", values = c("grey", "black","red","blue")) + theme_minimal(base_size = 20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "inside", legend.position.inside = c(.9, .8)) + scale_x_discrete(label = alt_labels) + theme(axis.text.x = element_text(size = 10)) + xlab("Altitudinal site") + ylab("Mean abundance per site") + scale_y_continuous(expand = c(0, 0), limits = c(0, 250))

m_abund_all %>%  ggplot(aes(x=site, y=individuals, fill=year)) + geom_boxplot() + scale_fill_manual(name = "Year", values = c("grey", "black","red","blue")) + theme_minimal(base_size = 20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + theme(legend.position = "inside", legend.position.inside = c(.9, .8)) + scale_x_discrete(label = alt_labels) + theme(axis.text.x = element_text(size = 10)) + xlab("Altitudinal site") + ylab("Mean abundance per site") + scale_y_continuous(expand = c(0, 0), limits = c(0, 250))

# poisson models for abundance

{pois_2022 <- beet_22_raw %>% filter(year == "2022") %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(c(ID, year, site, replicate, label, individuals)) %>% group_by(year, site, replicate) %>% summarise(individuals = sum(individuals), .groups = "drop")

  pois_2002 <- beet_02_raw %>% filter(year == "2002") %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(year, site, replicate, individuals) %>% group_by(year, site, replicate) %>% summarise(individuals = sum(individuals), .groups = "drop")
  
  pois_2003 <- beet_03_raw %>% filter(year == "2003") %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(year, site, replicate, individuals) %>% group_by(year, site, replicate) %>% summarise(individuals = sum(individuals), .groups = "drop")
  
  pois_2023 <- beet_23_raw %>% filter(year == "2023") %>% mutate(individuals = rowSums(across(where(is.numeric)), na.rm=TRUE)) %>% select(c(ID, year, site, replicate, label, individuals)) %>% group_by(year, site, replicate) %>% summarise(individuals = sum(individuals), .groups = "drop")

  pois_all <- rbind(pois_2002, pois_2022, pois_2003, pois_2023)
  
  vec <- c("1" = 0, "2" = 200, "3" = 300, "4" = 500, "5" = 700, "6" = 900, "7" = 1100, "8" = 1300, "9" = 1500, "10" = 1700, "11" = 1900, "12" = 1700, "13" = 1500, "14" = 1300, "15" = 1100, "16" = 900, "17" = 500)
  veg_cover <- read.csv("veg_cover.csv", stringsAsFactors = TRUE) %>% mutate(across(c(site, year), as.factor))
  soil_type <- read.csv("soil_type.csv", stringsAsFactors = TRUE) %>% mutate(across(c(site, year), as.factor))
  pois_all <- left_join(pois_all, enframe(vec), by = c("site" = "name")) %>% rename(altitude = value) %>% mutate(site = as.numeric(as.character(site)))
  pois_all <- merge(pois_all, veg_cover) %>% arrange(year, site) %>% mutate(site = as.factor(site)) %>% rename("exposed_rock" = "rock")
  pois_all <- merge(pois_all,soil_type)%>% arrange(year, site) %>% mutate(site = as.factor(site))
  rm(pois_2002, pois_2022, vec, pois_2003, pois_2023)
}



pois_all %>% mutate(site = as.numeric(site)) %>% ggplot(aes(x=site,y=altitude)) + geom_smooth(se=FALSE, linewidth = 3, colour = "brown", method = "loess", formula = "y ~ x") + scale_x_continuous(breaks = c(1:17), labels = alt_labels) + theme(axis.line = element_line(color='black'),plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border = element_blank())

m_0 <- glm(individuals ~ year, data = pois_all, family = poisson)

m_soil <- pois_all %>% filter(site %in% c(5,8,9,10)) %>% glm(data = ., individuals ~ year*clay*sand, family = poisson)
summary(m_soil)

for (k in 1:length(unique(pois_all$site))){
  m <- glm(individuals[site == k] ~ year[site == k], data = pois_all, family = poisson)
  assign(paste0("m_", k), m)
  rm(k, m)
}

summary(m_0) # *** --
summary(m_1) # x
summary(m_2) # *** ++
summary(m_3) # *** ++
summary(m_4) # * ++
summary(m_5) # ** --
summary(m_6) # * --
summary(m_7) # x
summary(m_8) # *** --
summary(m_9) # *** --
summary(m_10) # *** --
summary(m_11) # x
summary(m_12) # x
summary(m_13) # x
summary(m_14) # x
summary(m_15) # * ++
summary(m_16) # x
summary(m_17) # x

# climate analyses begin

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

datetest %>% filter(year > 2002, year < 2020, site %in% c(1,5,8,9,10)) %>% ggplot(aes(x = datetime, y = temperature)) + geom_smooth(aes(datetime, colour = "1"), formula = y ~ x, method = lm, se = T)

datacamp_colors <- list(green = "#74F065", darkblue = "#05192D", blue = "#47BEFB", pink = "#ED6AA8")
datetest %>% filter(year < 2020, year > 2002) %>% group_by(year,site, month) %>% summarise(temperature = mean(temperature)) %>% filter(site == 10) %>% ggplot(aes(x = month, y = temperature, color = year)) +geom_point() + geom_smooth(color = "red", se = FALSE) + scale_color_gradient(name = "ºC", low = datacamp_colors$blue, high = datacamp_colors$pink) + ggtitle("Temperature", subtitle = "Visualization using theme_datacamp()") + xlab("Year") + ylab("Mean Temperature")

datetest %>% group_by(site,year) %>% summarise(temperature = mean(temperature)) %>% filter(site %in% c(15,8,9,10), year < 2020, year > 2002) %>% ggplot(aes(x = year, y = temperature, color = temperature)) + geom_point() + geom_smooth(color = "red", se = FALSE) + scale_color_gradient(name = "ºC", low = datacamp_colors$blue, high = datacamp_colors$pink) + ggtitle("Temperature", subtitle = "Visualization using theme_datacamp()") + xlab("Year") + ylab("Mean Temperature")

for (k in (1:17)){m_t <- datetest %>% filter(year > 2002, year < 2021, site == k) %>% lm(data = ., temperature ~ year)
summary(m_t)
{if (k == 1){x1 <- predict(m_t, newdata = data.frame(year = c(2002,2022)))}}
{if (k == 1){x1 <- diff(x1)}}
{if (k > 1){x1 <- rbind(x1, diff(predict(m_t, newdata = data.frame(year = c(2002,2022)))))}}
{if (k==17){row.names(x1) <- NULL}}
{if (k==17){x1 <- data.frame(x1)}}
{if (k==17){x1$site <- c(1:17)}}
{if (k==17){x1 <- x1 %>% mutate(site = as.factor(site))}}
{if (k==17){rm(k)}}
}

x1 %>% ggplot(aes(x = site, y = X2)) + geom_point(aes(colour = as.factor(sign(X2))), size = 5, show.legend = FALSE, alpha = 0.6) + ylab("Change in MAT (°C)") + xlab("Site") + scale_x_discrete(label = c(1:17)) + geom_line(y=0, group = 1) + scale_colour_manual(name = "Year", values = c("-1" = "blue", "1" = "red")) + theme_minimal(base_size = 20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_discrete(label = alt_labels) + theme(axis.text.x = element_text(size = 10))

x1 %>% ggplot(aes(x = site, y = X2)) + geom_point(aes(colour = as.factor(c(-1,-1,-1,-1,1,-1,-1,1,1,1,-1,-1,-1,-1,-1,-1,-1))), size = 5, show.legend = FALSE, alpha = 0.6) + ylab("Change in MAT (°C)") + xlab("Site") + scale_x_discrete(label = c(1:17)) + geom_line(y= 0, group = 1) + scale_colour_manual(name = "Year", values = c("-1" = "grey", "1" = "red")) + theme_minimal(base_size = 20) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + scale_x_discrete(label = alt_labels) + theme(axis.text.x = element_text(size = 10))

{k <- c(5,8,9,10)
m_t <- datetest %>% filter(year > 2002, year < 2020, site %in% k) %>% lm(data = ., temperature ~ year + site)
rm(k)
summary(m_t)
}
