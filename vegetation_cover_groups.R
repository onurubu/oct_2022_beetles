##### Coskun Kucukkaragoz
#### 6th of August 2025
### R script to analyse the dominant group in the vegetation cover

# set-up ####

# rm(list = ls())

# packages
library(tidyverse)

# vegetation cover composition ####

veg_cover_all <- read.csv("veg_cover_modern.csv", stringsAsFactors = TRUE) %>% mutate(across(c(site, replicate, trap, year), as.factor))

{veg_groups <- veg_cover_all %>% select(c(-bare_ground, -litter, -rock, -tot_veg)) %>% group_by(year, site) %>% summarise(across(c(protea, restio, erica, strandveld, succ_karoo), mean), .groups = "drop") %>% 
  mutate(dom_veg = case_when(
    protea == pmax(protea, restio, erica, strandveld, succ_karoo) ~ "protea",
    restio == pmax(protea, restio, erica, strandveld, succ_karoo) ~ "restio",
    erica == pmax(protea, restio, erica, strandveld, succ_karoo) ~ "erica",
    strandveld == pmax(protea, restio, erica, strandveld, succ_karoo) ~ "strandveld",
    succ_karoo == pmax(protea, restio, erica, strandveld, succ_karoo) ~ "succulent_karoo"))
veg_groups_summary <- veg_groups %>% select(year, site, dom_veg)
write.csv(veg_groups, "veg_groups.csv", row.names = FALSE)}
