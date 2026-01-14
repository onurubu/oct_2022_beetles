## Adding collection data to beetle sorting spreadsheet for ID purposes
# Set up
rm(list = ls())
library(tidyverse)

# reading the data

beet_data <- read.csv("ceder_beetles_2022.csv", stringsAsFactors = TRUE) %>% mutate(across(c(ID, site, replicate, trap, label, season, year, full_label), factor)) %>% group_by(ID, site, replicate, trap, label, season, year, full_label) %>% summarise(across(where(is.numeric), sum), .groups = "drop")
collection_data <- read.csv("collection_data.csv") %>% mutate(across(c(site, replicate), factor))

beet_data_full <- left_join(collection_data, beet_data, by = c("site" = "site", "replicate" = "replicate")) %>% select(ID, site, replicate, trap, label, season, year, latitude, longitude, biome, elevation, country, province, locality, collector, method, everything()) %>% select(-full_label) %>% rename("elevation_m" = "elevation")

write.csv(beet_data_full, "./beetle_counts_with_collection_data_CK.csv")
