##### Coskun Kucukkaragoz
#### 25th of April 2024
### Start of the R script to analyse the data from my mounting, pinning, and sorting, of beetles for the October 2022 season from the Cederburg

# set-up ####

# packages
rm(list = ls())

{

  library(tidyverse)
  library(weights)
  library(emmeans)
  library(EnvStats)
  library(sf)
  library(hms)
  library(vegan)

}

# reading in raw data with correct data types

# species analysis ####
# set-up ###

{

  beet_wern <- read.csv("beet_wern_ID.csv", stringsAsFactors = TRUE) %>% mutate(across(c(ID, site, replicate, trap, label, season, year), factor)) %>% filter(ID != 578) %>% select(-label) %>% rename("cicindela_lurida" = Ropaloteres.lurida)
  
  beet_22 <- read.csv("ceder_beetles_2022.csv", stringsAsFactors = TRUE) %>% mutate(across(c(ID, site, replicate, trap, label, season, year, full_label), factor)) %>% group_by(ID, site, replicate, trap, season, year) %>% summarise(across(where(is.numeric), sum), .groups = "drop") %>% select(1:8, morphospecies_6) %>% rename("zophosis_gracilicornis" = morphospecies_6)
  
  abund_old <- read.csv("oct_2002_beetles_coskun_edit.csv", stringsAsFactors = TRUE) %>% mutate(across(c(label,site, replicate, year), factor)) %>% select(-label) %>% rename("anthia_decemguttata" = Thermophilum.decemguttatum, "stenocara_dentata" = Stenocara.dentata, "zophosis_gracilicornis" = Zophosis.sp.1, "cicindela_lurida" = Cicindela.brevicollis)
  
  # creating new dataframe with ID'd beetles as well as known species
  
  beet_modern <- inner_join(beet_22, beet_wern, by = c("ID", "site", "replicate", "trap", "season", "year"))
  
  abund_modern <- beet_modern %>% group_by(year, site, replicate) %>% summarise(across(where(is.numeric), sum), .groups = "drop") %>% relocate(site, replicate) %>% mutate(site = factor(site, levels = c(1:17)), replicate = factor(replicate, levels = c(1:4)), year = factor(year, levels = c(2022, 2023))) %>% complete(site, replicate, year) %>% arrange(year, site, replicate) %>% replace(is.na(.), 0)
  
  # data.table::fwrite(beet_modern, file = "./raw_identified_beetle_compiled_datasets/modern_families_raw.csv")
  
  # data.table::fwrite(abund_modern, file = "./raw_identified_beetle_compiled_datasets/modern_families_grouped.csv")
  
  # data.table::fwrite(abund_old, file = "./raw_identified_beetle_compiled_datasets/old_families_grouped.csv")
  
  df1 <- abund_modern %>% select(1:3, anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis, cicindela_lurida)
  
  df2 <- abund_old %>% select(1:3, anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis, cicindela_lurida)
  
  abund_all <- rbind(df2, df1) %>%  mutate(period = factor(case_when(year %in% c(2002, 2003) ~ "old", year %in% c(2022, 2023) ~ "modern"))) %>% relocate(site, replicate, year, period)
  
  # data.table::fwrite(abund_all, file = "./raw_identified_beetle_compiled_datasets/all_families_grouped.csv")
  
  rm(abund_modern, abund_old, beet_22, beet_modern, beet_wern, df1, df2)

}

# mds ###

# presence-absence ##
# all sites #

{
  
  x <- abund_all[, 5:8] %>% mutate(decostand(across(where(is.numeric)), "pa")) %>% filter(rowSums(across(where(is.numeric)))!=0)
  y <- abund_all %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(1:4)
  Herb_community.mds <- x %>% metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- y$site
  MDS_xy$replicate <- y$replicate
  MDS_xy$year <- y$year
  MDS_xy$period <- y$period
  Herb_community.mds$stress
  
  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + geom_point() + stat_ellipse() + theme_bw() + ggtitle("Species presence by year") + scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + geom_point() + stat_ellipse() + theme_bw() + ggtitle("Species presence by sampling period") + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"))
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Species presence by year with sites, overlaps jittered") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Species presence by year with sites") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE)
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Species presence by sampling period") + stat_ellipse() + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"), guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE)
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Species presence by site") + scale_color_discrete(name = "Site", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  
}

{

  ggsave("./MDS_plots/species/presence_absence/year_all.png", plot = p1, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/species/presence_absence/period_all.png", plot = p2, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/species/presence_absence/year_all_sites_jitter.png", plot = p3, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/species/presence_absence/year_all_sites.png", plot = p4, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/species/presence_absence/period_all_sites.png", plot = p5, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/species/presence_absence/sites_jitter.png", plot = p6, width = 4000, height = 2160, units = "px", bg = "white")

  rm(x, y, Herb_community.mds, MDS_xy, p1, p2, p3, p4, p5, p6)

}

dist <- vegdist(fam_presence_values, method = "bray")
anosim_type <- anosim(dist,fam_presence_cats$site)
summary(anosim_type)

# site-specific ##

for (k in 2:17){
  
  x <- abund_all %>% filter(site == k) %>% mutate(decostand(across(where(is.numeric)), "pa"))
  x1 <- x %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(5:8)
  x2 <- x %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(1:4)
  Herb_community.mds <- x1 %>% metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- x2$site
  MDS_xy$replicate <- x2$replicate
  MDS_xy$year <- x2$year
  MDS_xy$period <- x2$period
  Herb_community.mds$stress
  
  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + geom_jitter() + theme_bw() + ggtitle(paste0("Site ", k, ", species presence by year")) + scale_color_discrete(name = "Year")
  (pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + geom_jitter(width = 0.02, height = 0.02) + stat_ellipse() + theme_bw() + ggtitle(paste0("Site ", k, ", species presence by sampling period")) + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023")))
  
  ggsave(paste0("./MDS_plots/species/presence_absence/year/site_", k, ".png"), plot = pa, width = 4000, height = 2160, units = "px", bg = "white", create.dir = TRUE)
  ggsave(paste0("./MDS_plots/species/presence_absence/period/site_", k, ".png"), plot = pb, width = 4000, height = 2160, units = "px", bg = "white", create.dir = TRUE)
  
  {if (k == 17){rm(k, x, x1, x2, Herb_community.mds, MDS_xy, pa, pb)}}
  
}

# random

MDS_xy %>% filter(year == "2022") %>% ggplot(., aes(MDS1, MDS2, color = site)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family presence by site") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000) + stat_ellipse()



# abundance ##

{
  
  x <- abund_all[, 5:8] %>% filter(rowSums(across(where(is.numeric)))!=0)
  y <- abund_all %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(1:4)
  Herb_community.mds <- x %>% metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- y$site
  MDS_xy$replicate <- y$replicate
  MDS_xy$year <- y$year
  MDS_xy$period <- y$period
  Herb_community.mds$stress
  
  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + geom_point() + stat_ellipse() + theme_bw() + ggtitle("Species abundance by year") + scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + geom_point() + stat_ellipse() + theme_bw() + ggtitle("Species abundance by sampling period") + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"))
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Species abundance by year with sites, overlaps jittered") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Species abundance by year with sites") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE)
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Species abundance by sampling period") + stat_ellipse() + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"), guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE)
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Species abundance by site") + scale_color_discrete(name = "Site", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  
}

{
  
  ggsave("./MDS_plots/species/abundance/year_all.png", plot = p1, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/species/abundance/period_all.png", plot = p2, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/species/abundance/year_all_sites_jitter.png", plot = p3, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/species/abundance/year_all_sites.png", plot = p4, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/species/abundance/period_all_sites.png", plot = p5, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/species/abundance/sites_jitter.png", plot = p6, width = 4000, height = 2160, units = "px", bg = "white")
  
  rm(x, y, Herb_community.mds, MDS_xy, p1, p2, p3, p4, p5, p6)
  
}

dist <- vegdist(fam_presence_values, method = "bray")
anosim_type <- anosim(dist,fam_presence_cats$site)
summary(anosim_type)

# site-specific ##

for (k in 1:17){
  
  x <- abund_all %>% filter(site == k)
  x1 <- x %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(5:8)
  x2 <- x %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(1:4)
  Herb_community.mds <- x1 %>% metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- x2$site
  MDS_xy$replicate <- x2$replicate
  MDS_xy$year <- x2$year
  MDS_xy$period <- x2$period
  Herb_community.mds$stress
  
  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + geom_jitter() + theme_bw() + ggtitle(paste0("Site ", k, ", species abundance by year")) + scale_color_discrete(name = "Year")
  (pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + geom_jitter(width = 0.02, height = 0.02) + stat_ellipse() + theme_bw() + ggtitle(paste0("Site ", k, ", species abundance by sampling period")) + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023")))
  
  ggsave(paste0("./MDS_plots/species/abundance/year/site_", k, ".png"), plot = pa, width = 4000, height = 2160, units = "px", bg = "white", create.dir = TRUE)
  ggsave(paste0("./MDS_plots/species/abundance/period/site_", k, ".png"), plot = pb, width = 4000, height = 2160, units = "px", bg = "white", create.dir = TRUE)
  
  {if (k == 17){rm(k, x, x1, x2, Herb_community.mds, MDS_xy, pa, pb)}}
  
}

# family analysis ####
# set-up ###

{
  
  fam_carabid_modern <- data.table::fread("./family_groups/carabids_modern.csv") %>% mutate(across(c(site, replicate, year), factor))
  
  fam_tenebrionid_modern <- data.table::fread("./family_groups/tenebrionids_modern.csv") %>% mutate(across(c(site, replicate, year), factor))
  
  fam_cicindelid_modern <- data.table::fread("./family_groups/cicindelids_modern.csv") %>% mutate(across(c(site, replicate, year), factor))
  
  fam_carabid_old <- data.table::fread("./family_groups/carabids_old.csv") %>% mutate(across(c(site, replicate, year), factor))
  
  fam_tenebrionid_old <- data.table::fread("./family_groups/tenebrionids_old.csv") %>% mutate(across(c(site, replicate, year), factor))
  
  fam_cicindelid_old <- data.table::fread("./family_groups/cicindelids_old.csv") %>% mutate(across(c(site, replicate, year), factor))
  
  {

    fam_carabid_modern_sums <- fam_carabid_modern %>% mutate("carabids" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, carabids)
    fam_carabid_old_sums <- fam_carabid_old %>% mutate("carabids" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, carabids)
    fam_carabid_sums <- rbind(fam_carabid_old_sums, fam_carabid_modern_sums)
    rm(fam_carabid_modern_sums, fam_carabid_old_sums)
    
      fam_tenebrionid_modern_sums <- fam_tenebrionid_modern %>% mutate("tenebrionids" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, tenebrionids)
      fam_tenebrionid_old_sums <- fam_tenebrionid_old %>% mutate("tenebrionids" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, tenebrionids)
      fam_tenebrionid_sums <- rbind(fam_tenebrionid_old_sums, fam_tenebrionid_modern_sums)
      rm(fam_tenebrionid_modern_sums, fam_tenebrionid_old_sums)
    
      fam_cicindelid_modern_sums <- fam_cicindelid_modern %>% mutate("cicindelids" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, cicindelids)
      fam_cicindelid_old_sums <- fam_cicindelid_old %>% mutate("cicindelids" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, cicindelids)
      fam_cicindelid_sums <- rbind(fam_cicindelid_old_sums, fam_cicindelid_modern_sums)
      rm(fam_cicindelid_modern_sums, fam_cicindelid_old_sums)
    
    fam_sums <- inner_join(fam_carabid_sums, fam_cicindelid_sums, by = c("site", "replicate", "year"))
    fam_sums <- inner_join(fam_sums, fam_tenebrionid_sums, by = c("site", "replicate", "year"))
    rm(fam_cicindelid_sums, fam_carabid_sums, fam_tenebrionid_sums)

  }

  {
  
    x1 <- fam_carabid_modern %>% mutate(decostand(across(where(is.numeric)), "pa")) %>% mutate("carabid_diversity" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, carabid_diversity)
    x2 <- fam_carabid_old %>% mutate(decostand(across(where(is.numeric)), "pa")) %>% mutate("carabid_diversity" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, carabid_diversity)
    x <- rbind(x2, x1)
    
    y1 <- fam_cicindelid_modern %>% mutate(decostand(across(where(is.numeric)), "pa")) %>% mutate("cicindelid_diversity" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, cicindelid_diversity)
    y2 <- fam_cicindelid_old %>% mutate(decostand(across(where(is.numeric)), "pa")) %>% mutate("cicindelid_diversity" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, cicindelid_diversity)
    y <- rbind(y2, y1)
    
    z1 <- fam_tenebrionid_modern %>% mutate(decostand(across(where(is.numeric)), "pa")) %>% mutate("tenebrionid_diversity" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, tenebrionid_diversity)
    z2 <- fam_tenebrionid_old %>% mutate(decostand(across(where(is.numeric)), "pa")) %>% mutate("tenebrionid_diversity" = rowSums(across(where(is.numeric)))) %>% select(site, replicate, year, tenebrionid_diversity)
    z <- rbind(z2, z1)
    
    fam_presence <- inner_join(x, z, by = c("site", "replicate", "year"))
    fam_presence <- inner_join(fam_presence, y, by = c("site", "replicate", "year"))
    
    rm(x1, x2, x, y1, y2, y, z1, z2, z)
  
  }

  fam_sums <- fam_sums %>%  mutate(period = factor(case_when(year %in% c(2002, 2003) ~ "old", year %in% c(2022, 2023) ~ "modern"))) %>% relocate(site, replicate, year, period)
  
  fam_presence <- fam_presence %>%  mutate(period = factor(case_when(year %in% c(2002, 2003) ~ "old", year %in% c(2022, 2023) ~ "modern"))) %>% relocate(site, replicate, year, period)
  
  rm(fam_carabid_modern, fam_carabid_old, fam_cicindelid_modern, fam_cicindelid_old, fam_tenebrionid_modern, fam_tenebrionid_old)
  
}

# mds trials ###

# presence-absence ##
# all sites #

{
  
  fam_presence_values <- fam_presence[, 5:7] %>% mutate(decostand(across(where(is.numeric)), "pa")) %>% filter(rowSums(across(where(is.numeric)))!=0)
  fam_presence_cats <- fam_presence %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(1:4)
  Herb_community.mds <- fam_presence_values %>% metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress
  
  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + geom_point() + stat_ellipse() + theme_bw() + ggtitle("Family presence by year") + scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + geom_point() + stat_ellipse() + theme_bw() + ggtitle("Family presence by sampling period") + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"))
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family presence by year with sites, overlaps jittered") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family presence by year with sites") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE)
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family presence by sampling period") + stat_ellipse() + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"), guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE)
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family presence by site") + scale_color_discrete(name = "Site", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  
}

{
  ggsave("./MDS_plots/family_grouped/presence_absence/year_all.png", plot = p1, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/family_grouped/presence_absence/period_all.png", plot = p2, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/family_grouped/presence_absence/year_all_sites_jitter.png", plot = p3, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/family_grouped/presence_absence/year_all_sites.png", plot = p4, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/family_grouped/presence_absence/period_all_sites.png", plot = p5, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/family_grouped/presence_absence/sites_jitter.png", plot = p6, width = 4000, height = 2160, units = "px", bg = "white")

  rm(fam_presence_values, fam_presence_cats, Herb_community.mds, MDS_xy, p1, p2, p3, p4, p5, p6)
  
}

dist <- vegdist(fam_presence_values, method = "bray")
anosim_type <- anosim(dist,fam_presence_cats$site)
summary(anosim_type)

# site-specific #

for (k in 1:17){
  
  x <- fam_presence %>% filter(site == k) %>% mutate(decostand(across(where(is.numeric)), "pa"))
  fam_presence_values <- x %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(5:7)
  fam_presence_cats <- x %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(1:4)
  Herb_community.mds <- fam_presence_values %>% metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress
  
  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + geom_jitter() + theme_bw() + ggtitle(paste0("Site ", k, ", family presence by year")) + scale_color_discrete(name = "Year")
  (pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + geom_jitter(width = 0.02, height = 0.02) + stat_ellipse() + theme_bw() + ggtitle(paste0("Site ", k, ", family presence by sampling period")) + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023")))
  
  ggsave(paste0("./MDS_plots/family_grouped/presence_absence/year/site_", k, ".png"), plot = pa, width = 4000, height = 2160, units = "px", bg = "white", create.dir = TRUE)
  ggsave(paste0("./MDS_plots/family_grouped/presence_absence/period/site_", k, ".png"), plot = pb, width = 4000, height = 2160, units = "px", bg = "white", create.dir = TRUE)
  
  {if (k == 17){rm(x, fam_presence_values, fam_presence_cats, Herb_community.mds, MDS_xy, pa, pb)}}
  
  }

# random

MDS_xy %>% filter(year == "2022") %>% ggplot(., aes(MDS1, MDS2, color = site)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family presence by site") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000) + stat_ellipse()



# diversity ##
# all sites #

{
  
  fam_presence_values <- fam_presence[, 5:7] %>% filter(rowSums(across(where(is.numeric)))!=0)
  fam_presence_cats <- fam_presence %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(1:4)
  Herb_community.mds <- fam_presence_values %>% metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress
  
  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + geom_point() + stat_ellipse() + theme_bw() + ggtitle("Family diversity by year") + scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + geom_point() + stat_ellipse() + theme_bw() + ggtitle("Family diversity by sampling period") + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"))
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family diversity by year with sites, overlaps jittered") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family diversity by year with sites") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE)
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family diversity by sampling period") + stat_ellipse() + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"), guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE)
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family diversity by site") + scale_color_discrete(name = "Site", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  
}

{
  
ggsave("./MDS_plots/family_grouped/diversity/year_all.png", plot = p1, width = 4000, height = 2160, units = "px", bg = "white")
ggsave("./MDS_plots/family_grouped/diversity/period_all.png", plot = p2, width = 4000, height = 2160, units = "px", bg = "white")
ggsave("./MDS_plots/family_grouped/diversity/year_all_sites_jitter.png", plot = p3, width = 4000, height = 2160, units = "px", bg = "white")
ggsave("./MDS_plots/family_grouped/diversity/year_all_sites.png", plot = p4, width = 4000, height = 2160, units = "px", bg = "white")
ggsave("./MDS_plots/family_grouped/diversity/period_all_sites.png", plot = p5, width = 4000, height = 2160, units = "px", bg = "white")
ggsave("./MDS_plots/family_grouped/diversity/sites_jitter.png", plot = p6, width = 4000, height = 2160, units = "px", bg = "white")

rm(fam_presence_values, fam_presence_cats, Herb_community.mds, MDS_xy, p1, p2, p3, p4, p5, p6)

}

# dist <- vegdist(fam_presence_values, method = "bray")
# anosim_type <- anosim(dist,fam_presence_cats$site)
# summary(anosim_type)

# site-specific #

for (k in 1:17){
  
  x <- fam_presence %>% filter(site == k)
  fam_presence_values <- x %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(5:7)
  fam_presence_cats <- x %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(1:4)
  Herb_community.mds <- fam_presence_values %>% metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress
  
  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + geom_jitter() + theme_bw() + ggtitle(paste0("Site ", k, ", family diversity by year")) + scale_color_discrete(name = "Year")
  pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + geom_jitter(width = 0.02, height = 0.02) + stat_ellipse() + theme_bw() + ggtitle(paste0("Site ", k, ", family diversity by sampling period")) + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"))
  
  ggsave(paste0("./MDS_plots/family_grouped/diversity/year/site_", k, ".png"), plot = pa, width = 4000, height = 2160, units = "px", bg = "white", create.dir = TRUE)
  ggsave(paste0("./MDS_plots/family_grouped/diversity/period/site_", k, ".png"), plot = pb, width = 4000, height = 2160, units = "px", bg = "white", create.dir = TRUE)
  
  {if (k == 17){rm(k, x, fam_presence_values, fam_presence_cats, Herb_community.mds, MDS_xy, pa, pb)}}
  
}

# abundance ##
# all sites #

{
  
  fam_sums_values <- fam_sums[, 5:7] %>% filter(rowSums(across(where(is.numeric)))!=0)
  fam_sums_cats <- fam_sums %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(1:4)
  Herb_community.mds <- fam_sums_values %>% metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_sums_cats$site
  MDS_xy$replicate <- fam_sums_cats$replicate
  MDS_xy$year <- fam_sums_cats$year
  MDS_xy$period <- fam_sums_cats$period
  Herb_community.mds$stress
  
  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + geom_point() + stat_ellipse() + theme_bw() + ggtitle("Family abundance by year") + scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + geom_point() + stat_ellipse() + theme_bw() + ggtitle("Family abundance by sampling period") + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"))
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family abundance by year with sites, overlaps jittered") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family abundance by year with sites") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE)
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family abundance by sampling period") + stat_ellipse() + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023"), guide = guide_legend(override.aes = list(alpha = 1))) + geom_text(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE)
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Family abundance by site") + scale_color_discrete(name = "Site", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  
}

{
    
  ggsave("./MDS_plots/family_grouped/abundance/year_all.png", plot = p1, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/family_grouped/abundance/period_all.png", plot = p2, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/family_grouped/abundance/year_all_sites_jitter.png", plot = p3, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/family_grouped/abundance/year_all_sites.png", plot = p4, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/family_grouped/abundance/period_all_sites.png", plot = p5, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave("./MDS_plots/family_grouped/abundance/sites_jitter.png", plot = p6, width = 4000, height = 2160, units = "px", bg = "white")
  
  rm(fam_sums_values, fam_sums_cats, Herb_community.mds, MDS_xy, p1, p2, p3, p4, p5, p6)
  
}

dist <- vegdist(fam_sums_values, method = "bray")
anosim_type <- anosim(dist,fam_sums_cats$site)
summary(anosim_type)

# site-specific #

for (k in 1:17){
  
  x <- fam_sums %>% filter(site == k)
  fam_sums_values <- x %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(5:7)
  fam_sums_cats <- x %>% filter(rowSums(across(where(is.numeric)))!=0) %>% select(1:4)
  Herb_community.mds <- fam_sums_values %>% metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_sums_cats$site
  MDS_xy$replicate <- fam_sums_cats$replicate
  MDS_xy$year <- fam_sums_cats$year
  MDS_xy$period <- fam_sums_cats$period
  Herb_community.mds$stress
  
  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + geom_jitter() + theme_bw() + ggtitle(paste0("Site ", k, ", family abundance by year")) + scale_color_discrete(name = "Year")
  (pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) + geom_jitter(width = 0.02, height = 0.02) + stat_ellipse() + theme_bw() + ggtitle(paste0("Site ", k, ", family abundance by sampling period")) + scale_color_discrete(name = "Sampling period", labels = c("old" = "2002/2003", "modern" = "2022/2023")))
  
  ggsave(paste0("./MDS_plots/family_grouped/abundance/year/site_", k, ".png"), plot = pa, width = 4000, height = 2160, units = "px", bg = "white", create.dir = TRUE)
  ggsave(paste0("./MDS_plots/family_grouped/abundance/period/site_", k, ".png"), plot = pb, width = 4000, height = 2160, units = "px", bg = "white", create.dir = TRUE)
  
  {if (k == 17){rm(k, x, fam_sumse_values, fam_sums_cats, Herb_community.mds, MDS_xy, pa, pb)}}
  
}