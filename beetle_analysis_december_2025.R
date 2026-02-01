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
  library(MCMCglmm)
}

# reading in raw data with correct data types

alt_labels <- c(
  "0m(W)",
  "200m(W)",
  "300m(W)",
  "500m(W)",
  "700m(W)",
  "900m(W)",
  "1100m(W)",
  "1300m(W)",
  "1500m(W)",
  "1700m(W)",
  "1900m(W)",
  "1700m(E)",
  "1500m(E)",
  "1300m(E)",
  "1100m(E)",
  "900m(E)",
  "500m(E)"
)

{
  beet_wern <- read.csv(
    "./working_beetle_data/beet_wern_ID.csv",
    stringsAsFactors = TRUE
  ) %>%
    mutate(across(
      c(ID, site, replicate, trap, label, season, year),
      factor
    )) %>%
    filter(ID != 578) %>%
    select(-label) %>%
    rename("cicindela_lurida" = Ropaloteres.lurida) %>%
    arrange(year)

  beet_22 <- read.csv(
    "./working_beetle_data/ceder_beetles_2022.csv",
    stringsAsFactors = TRUE
  ) %>%
    mutate(across(
      c(ID, site, replicate, trap, label, season, year, full_label),
      factor
    )) %>%
    group_by(ID, site, replicate, trap, season, year) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop") %>%
    select(1:8, morphospecies_6) %>%
    rename("zophosis_gracilicornis" = morphospecies_6)

  abund_old <- read.csv(
    "./working_beetle_data/oct_2002_beetles_coskun_edit.csv",
    stringsAsFactors = TRUE
  ) %>%
    mutate(across(c(label, site, replicate, year), factor)) %>%
    select(-label) %>%
    rename(
      "anthia_decemguttata" = Thermophilum.decemguttatum,
      "stenocara_dentata" = Stenocara.dentata,
      "zophosis_gracilicornis" = Zophosis.sp.1,
      "cicindela_lurida" = Cicindela.brevicollis
    ) %>%
    mutate(season = as.factor("October")) %>%
    relocate(season, .before = 4)

  beet_23 <- data.table::fread(
    "./working_beetle_data/oct_2023_beetle_abundance.csv"
  ) %>%
    mutate(across(
      c(ID, site, replicate, trap, label, season, year),
      factor
    )) %>%
    select(-label)

  # creating new dataframe with ID'd beetles as well as known species
  beet_wide <- inner_join(
    beet_22,
    beet_wern,
    by = c("ID", "site", "replicate", "trap", "season", "year")
  ) %>%
    arrange(year)

  beet_modern <- rbind(beet_22, beet_23) %>% arrange(year, season, ID)

  abund_modern <- beet_modern %>%
    group_by(year, season, site, replicate) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop") %>%
    relocate(site, replicate) %>%
    mutate(
      site = factor(site, levels = c(1:17)),
      replicate = factor(replicate, levels = c(1:4)),
      year = factor(year, levels = c(2022, 2023)),
      season = factor(season, levels = c("March", "October"))
    ) %>%
    complete(site, replicate, year, season) %>%
    arrange(year, season, site, replicate) %>%
    replace(is.na(.), 0) %>%
    filter(paste(year, season) != "2022 March")

  # df1 <- abund_modern %>% select(1:3, anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis, cicindela_lurida)
  #
  # df2 <- abund_old %>% select(1:3, anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis, cicindela_lurida)

  df1 <- abund_modern %>%
    select(1:4, anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis)

  df2 <- abund_old %>%
    select(1:4, anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis)

  abund_all <- rbind(df2, df1) %>%
    mutate(
      period = factor(case_when(
        year %in% c(2002, 2003) ~ "old",
        year %in% c(2022, 2023) ~ "modern"
      ))
    ) %>%
    relocate(site, replicate, year, season, period)

  # data.table::fwrite(abund_all, file = "./raw_identified_beetle_compiled_datasets/all_families_grouped.csv")

  rm(abund_modern, abund_old, beet_wern, df1, df2, beet_22, beet_23)
}

abund_all <- abund_all %>% filter(paste(year, season) != "2023 March")

Sys.sleep(0.1)

# species analysis december 2025 ####
# mds ###

# presence-absence ##
# all sites #

abund_all <- abund_all %>% filter(paste(year, season) != "2023 March")


{
  x <- abund_all[, 6:8] %>%
    mutate(decostand(across(where(is.numeric)), "pa")) %>%
    filter(rowSums(across(where(is.numeric))) != 0)
  y <- abund_all %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:5)
  Herb_community.mds <- x %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- y$site
  MDS_xy$replicate <- y$replicate
  MDS_xy$year <- y$year
  MDS_xy$period <- y$period
  Herb_community.mds$stress

  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Species presence by year") +
    scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Species presence by sampling period") +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species presence by year with sites, overlaps jittered") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species presence by year with sites") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species presence by sampling period") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023"),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species presence by site") +
    scale_color_discrete(
      name = "Site",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
}

{
  ggsave(
    "./MDS_plots/species/presence_absence/year_all.png",
    plot = p1,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/species/presence_absence/period_all.png",
    plot = p2,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/species/presence_absence/year_all_sites_jitter.png",
    plot = p3,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/species/presence_absence/year_all_sites.png",
    plot = p4,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/species/presence_absence/period_all_sites.png",
    plot = p5,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/species/presence_absence/sites_jitter.png",
    plot = p6,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )

  rm(x, y, Herb_community.mds, MDS_xy, p1, p2, p3, p4, p5, p6)
}

dist <- vegdist(fam_presence_values, method = "bray")
anosim_type <- anosim(dist, fam_presence_cats$site)
summary(anosim_type)

# site-specific ##

for (k in 1:17) {
  x <- abund_all %>%
    filter(site == k) %>%
    mutate(decostand(across(where(is.numeric)), "pa"))
  x1 <- x %>% filter(rowSums(across(where(is.numeric))) != 0) %>% select(6:8)
  x2 <- x %>% filter(rowSums(across(where(is.numeric))) != 0) %>% select(1:5)
  Herb_community.mds <- x1 %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- x2$site
  MDS_xy$replicate <- x2$replicate
  MDS_xy$year <- x2$year
  MDS_xy$period <- x2$period
  Herb_community.mds$stress

  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_jitter() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", species presence by year")) +
    scale_color_discrete(name = "Year")
  (pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_jitter(width = 0.02, height = 0.02) +
    stat_ellipse() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", species presence by sampling period")) +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    ))

  ggsave(
    paste0("./MDS_plots/species/presence_absence/year/site_", k, ".png"),
    plot = pa,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0("./MDS_plots/species/presence_absence/period/site_", k, ".png"),
    plot = pb,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  {
    if (k == 17) {
      rm(k, x, x1, x2, Herb_community.mds, MDS_xy, pa, pb)
    }
  }
}

# random

MDS_xy %>%
  filter(year == "2022") %>%
  ggplot(., aes(MDS1, MDS2, color = site)) +
  theme_bw() +
  geom_point(alpha = 0) +
  ggtitle("Family presence by site") +
  scale_color_discrete(
    name = "Year",
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  geom_text(
    data = MDS_xy,
    mapping = aes(MDS1, MDS2, label = site),
    show.legend = FALSE,
    max.overlaps = 1000
  ) +
  stat_ellipse()


# abundance ##

{
  x <- abund_all[, 6:8] %>% filter(rowSums(across(where(is.numeric))) != 0)
  y <- abund_all %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:5)
  Herb_community.mds <- x %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- y$site
  MDS_xy$replicate <- y$replicate
  MDS_xy$year <- y$year
  MDS_xy$period <- y$period
  Herb_community.mds$stress

  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Species abundance by year") +
    scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Species abundance by sampling period") +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species abundance by year with sites, overlaps jittered") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species abundance by year with sites") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species abundance by sampling period") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023"),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species abundance by site") +
    scale_color_discrete(
      name = "Site",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
}

{
  ggsave(
    "./MDS_plots/species/abundance/year_all.png",
    plot = p1,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/species/abundance/period_all.png",
    plot = p2,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/species/abundance/year_all_sites_jitter.png",
    plot = p3,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/species/abundance/year_all_sites.png",
    plot = p4,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/species/abundance/period_all_sites.png",
    plot = p5,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/species/abundance/sites_jitter.png",
    plot = p6,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )

  rm(x, y, Herb_community.mds, MDS_xy, p1, p2, p3, p4, p5, p6)
}

dist <- vegdist(fam_presence_values, method = "bray")
anosim_type <- anosim(dist, fam_presence_cats$site)
summary(anosim_type)

# site-specific ##

for (k in 1:17) {
  x <- abund_all %>% filter(site == k)
  x1 <- x %>% filter(rowSums(across(where(is.numeric))) != 0) %>% select(6:8)
  x2 <- x %>% filter(rowSums(across(where(is.numeric))) != 0) %>% select(1:5)
  Herb_community.mds <- x1 %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- x2$site
  MDS_xy$replicate <- x2$replicate
  MDS_xy$year <- x2$year
  MDS_xy$period <- x2$period
  Herb_community.mds$stress

  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_jitter() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", species abundance by year")) +
    scale_color_discrete(name = "Year")
  (pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_jitter(width = 0.02, height = 0.02) +
    stat_ellipse() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", species abundance by sampling period")) +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    ))

  ggsave(
    paste0("./MDS_plots/species/abundance/year/site_", k, ".png"),
    plot = pa,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0("./MDS_plots/species/abundance/period/site_", k, ".png"),
    plot = pb,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  {
    if (k == 17) {
      rm(k, x, x1, x2, Herb_community.mds, MDS_xy, pa, pb)
    }
  }
}

# family analysis december 2025 ####
# set-up ###

{
  fam_carabid_modern <- data.table::fread(
    "./family_groups/carabids_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  fam_tenebrionid_modern <- data.table::fread(
    "./family_groups/tenebrionids_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  fam_cicindelid_modern <- data.table::fread(
    "./family_groups/cicindelids_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  fam_carabid_old <- data.table::fread("./family_groups/carabids_old.csv") %>%
    mutate(across(c(site, replicate, year), factor))

  fam_tenebrionid_old <- data.table::fread(
    "./family_groups/tenebrionids_old.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  fam_cicindelid_old <- data.table::fread(
    "./family_groups/cicindelids_old.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  {
    fam_carabid_modern_sums <- fam_carabid_modern %>%
      mutate("carabids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, carabids)
    fam_carabid_old_sums <- fam_carabid_old %>%
      mutate("carabids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, carabids)
    fam_carabid_sums <- rbind(fam_carabid_old_sums, fam_carabid_modern_sums)
    rm(fam_carabid_modern_sums, fam_carabid_old_sums)

    fam_tenebrionid_modern_sums <- fam_tenebrionid_modern %>%
      mutate("tenebrionids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, tenebrionids)
    fam_tenebrionid_old_sums <- fam_tenebrionid_old %>%
      mutate("tenebrionids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, tenebrionids)
    fam_tenebrionid_sums <- rbind(
      fam_tenebrionid_old_sums,
      fam_tenebrionid_modern_sums
    )
    rm(fam_tenebrionid_modern_sums, fam_tenebrionid_old_sums)

    fam_cicindelid_modern_sums <- fam_cicindelid_modern %>%
      mutate("cicindelids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, cicindelids)
    fam_cicindelid_old_sums <- fam_cicindelid_old %>%
      mutate("cicindelids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, cicindelids)
    fam_cicindelid_sums <- rbind(
      fam_cicindelid_old_sums,
      fam_cicindelid_modern_sums
    )
    rm(fam_cicindelid_modern_sums, fam_cicindelid_old_sums)

    fam_sums <- inner_join(
      fam_carabid_sums,
      fam_cicindelid_sums,
      by = c("site", "replicate", "year")
    )
    fam_sums <- inner_join(
      fam_sums,
      fam_tenebrionid_sums,
      by = c("site", "replicate", "year")
    )
    rm(fam_cicindelid_sums, fam_carabid_sums, fam_tenebrionid_sums)
  }

  {
    x1 <- fam_carabid_modern %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("carabid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, carabid_diversity)
    x2 <- fam_carabid_old %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("carabid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, carabid_diversity)
    x <- rbind(x2, x1)

    y1 <- fam_cicindelid_modern %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("cicindelid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, cicindelid_diversity)
    y2 <- fam_cicindelid_old %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("cicindelid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, cicindelid_diversity)
    y <- rbind(y2, y1)

    z1 <- fam_tenebrionid_modern %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("tenebrionid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, tenebrionid_diversity)
    z2 <- fam_tenebrionid_old %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("tenebrionid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, tenebrionid_diversity)
    z <- rbind(z2, z1)

    fam_presence <- inner_join(x, z, by = c("site", "replicate", "year"))
    fam_presence <- inner_join(
      fam_presence,
      y,
      by = c("site", "replicate", "year")
    )

    rm(x1, x2, x, y1, y2, y, z1, z2, z)
  }

  fam_sums <- fam_sums %>%
    mutate(
      period = factor(case_when(
        year %in% c(2002, 2003) ~ "old",
        year %in% c(2022, 2023) ~ "modern"
      ))
    ) %>%
    relocate(site, replicate, year, period)

  fam_presence <- fam_presence %>%
    mutate(
      period = factor(case_when(
        year %in% c(2002, 2003) ~ "old",
        year %in% c(2022, 2023) ~ "modern"
      ))
    ) %>%
    relocate(site, replicate, year, period)

  rm(
    fam_carabid_modern,
    fam_carabid_old,
    fam_cicindelid_modern,
    fam_cicindelid_old,
    fam_tenebrionid_modern,
    fam_tenebrionid_old
  )
}

# mds trials ###

# presence-absence ##
# all sites #

{
  fam_presence_values <- fam_presence[, 5:7] %>%
    mutate(decostand(across(where(is.numeric)), "pa")) %>%
    filter(rowSums(across(where(is.numeric))) != 0)
  fam_presence_cats <- fam_presence %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_presence_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress

  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family presence by year") +
    scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family presence by sampling period") +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family presence by year with sites, overlaps jittered") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family presence by year with sites") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family presence by sampling period") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023"),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family presence by site") +
    scale_color_discrete(
      name = "Site",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
}

{
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/year_all.png",
    plot = p1,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/period_all.png",
    plot = p2,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/year_all_sites_jitter.png",
    plot = p3,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/year_all_sites.png",
    plot = p4,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/period_all_sites.png",
    plot = p5,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/sites_jitter.png",
    plot = p6,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )

  rm(
    fam_presence_values,
    fam_presence_cats,
    Herb_community.mds,
    MDS_xy,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6
  )
}

dist <- vegdist(fam_presence_values, method = "bray")
anosim_type <- anosim(dist, fam_presence_cats$site)
summary(anosim_type)

# site-specific #

for (k in 1:17) {
  x <- fam_presence %>%
    filter(site == k) %>%
    mutate(decostand(across(where(is.numeric)), "pa"))
  fam_presence_values <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(5:7)
  fam_presence_cats <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_presence_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress

  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_jitter() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family presence by year")) +
    scale_color_discrete(name = "Year")
  (pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_jitter(width = 0.02, height = 0.02) +
    stat_ellipse() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family presence by sampling period")) +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    ))

  ggsave(
    paste0("./MDS_plots/family_grouped/presence_absence/year/site_", k, ".png"),
    plot = pa,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./MDS_plots/family_grouped/presence_absence/period/site_",
      k,
      ".png"
    ),
    plot = pb,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  {
    if (k == 17) {
      rm(
        k,
        x,
        fam_presence_values,
        fam_presence_cats,
        Herb_community.mds,
        MDS_xy,
        pa,
        pb
      )
    }
  }
}

# random

MDS_xy %>%
  filter(year == "2022") %>%
  ggplot(., aes(MDS1, MDS2, color = site)) +
  theme_bw() +
  geom_point(alpha = 0) +
  ggtitle("Family presence by site") +
  scale_color_discrete(
    name = "Year",
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  geom_text(
    data = MDS_xy,
    mapping = aes(MDS1, MDS2, label = site),
    show.legend = FALSE,
    max.overlaps = 1000
  ) +
  stat_ellipse()


# diversity ##
# all sites #

{
  fam_presence_values <- fam_presence[, 5:7] %>%
    filter(rowSums(across(where(is.numeric))) != 0)
  fam_presence_cats <- fam_presence %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_presence_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress

  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family diversity by year") +
    scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family diversity by sampling period") +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family diversity by year with sites, overlaps jittered") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family diversity by year with sites") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family diversity by sampling period") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023"),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family diversity by site") +
    scale_color_discrete(
      name = "Site",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
}

{
  ggsave(
    "./MDS_plots/family_grouped/diversity/year_all.png",
    plot = p1,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/diversity/period_all.png",
    plot = p2,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/diversity/year_all_sites_jitter.png",
    plot = p3,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/diversity/year_all_sites.png",
    plot = p4,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/diversity/period_all_sites.png",
    plot = p5,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/diversity/sites_jitter.png",
    plot = p6,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )

  rm(
    fam_presence_values,
    fam_presence_cats,
    Herb_community.mds,
    MDS_xy,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6
  )
}

# dist <- vegdist(fam_presence_values, method = "bray")
# anosim_type <- anosim(dist,fam_presence_cats$site)
# summary(anosim_type)

# site-specific #

for (k in 1:17) {
  x <- fam_presence %>% filter(site == k)
  fam_presence_values <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(5:7)
  fam_presence_cats <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_presence_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress

  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_jitter() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family diversity by year")) +
    scale_color_discrete(name = "Year")
  pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_jitter(width = 0.02, height = 0.02) +
    stat_ellipse() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family diversity by sampling period")) +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )

  ggsave(
    paste0("./MDS_plots/family_grouped/diversity/year/site_", k, ".png"),
    plot = pa,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0("./MDS_plots/family_grouped/diversity/period/site_", k, ".png"),
    plot = pb,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  {
    if (k == 17) {
      rm(
        k,
        x,
        fam_presence_values,
        fam_presence_cats,
        Herb_community.mds,
        MDS_xy,
        pa,
        pb
      )
    }
  }
}

# abundance ##
# all sites #

{
  fam_sums_values <- fam_sums[, 5:7] %>%
    filter(rowSums(across(where(is.numeric))) != 0)
  fam_sums_cats <- fam_sums %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_sums_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_sums_cats$site
  MDS_xy$replicate <- fam_sums_cats$replicate
  MDS_xy$year <- fam_sums_cats$year
  MDS_xy$period <- fam_sums_cats$period
  Herb_community.mds$stress

  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family abundance by year") +
    scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family abundance by sampling period") +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family abundance by year with sites, overlaps jittered") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family abundance by year with sites") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family abundance by sampling period") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023"),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family abundance by site") +
    scale_color_discrete(
      name = "Site",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
}

{
  ggsave(
    "./MDS_plots/family_grouped/abundance/year_all.png",
    plot = p1,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/abundance/period_all.png",
    plot = p2,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/abundance/year_all_sites_jitter.png",
    plot = p3,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/abundance/year_all_sites.png",
    plot = p4,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/abundance/period_all_sites.png",
    plot = p5,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/abundance/sites_jitter.png",
    plot = p6,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )

  rm(
    fam_sums_values,
    fam_sums_cats,
    Herb_community.mds,
    MDS_xy,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6
  )
}

dist <- vegdist(fam_sums_values, method = "bray")
anosim_type <- anosim(dist, fam_sums_cats$site)
summary(anosim_type)

# site-specific #

for (k in 1:17) {
  x <- fam_sums %>% filter(site == k)
  fam_sums_values <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(5:7)
  fam_sums_cats <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_sums_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_sums_cats$site
  MDS_xy$replicate <- fam_sums_cats$replicate
  MDS_xy$year <- fam_sums_cats$year
  MDS_xy$period <- fam_sums_cats$period
  Herb_community.mds$stress

  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_jitter() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family abundance by year")) +
    scale_color_discrete(name = "Year")
  (pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_jitter(width = 0.02, height = 0.02) +
    stat_ellipse() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family abundance by sampling period")) +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    ))

  ggsave(
    paste0("./MDS_plots/family_grouped/abundance/year/site_", k, ".png"),
    plot = pa,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0("./MDS_plots/family_grouped/abundance/period/site_", k, ".png"),
    plot = pb,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  {
    if (k == 17) {
      rm(
        k,
        x,
        fam_sumse_values,
        fam_sums_cats,
        Herb_community.mds,
        MDS_xy,
        pa,
        pb
      )
    }
  }
}


# January 2026 analyses #####

# werner ID boxplots ####

# all species ###

{
  modern_species_richness <- beet_wide %>%
    group_by(year, season, site, replicate) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop") %>%
    mutate(
      site = factor(site, levels = c(1:17)),
      replicate = factor(replicate, levels = c(1:4)),
      year = factor(year, levels = c(2022, 2023)),
      season = factor(season, levels = c("March", "October"))
    ) %>%
    complete(site, replicate, year, season) %>%
    replace(is.na(.), 0) %>%
    filter(
      paste(year, season) != "2022 March",
      paste(year, season) != "2023 October"
    ) %>%
    relocate(site, replicate) %>%
    arrange(year)

  modern_species_richness <- modern_species_richness %>%
    mutate(richness = specnumber(modern_species_richness[, -(1:4)])) %>%
    select(1:4, richness)

  modern_species_abund <- beet_wide %>%
    select(-ID) %>%
    mutate(
      site = factor(site, levels = c(1:17)),
      replicate = factor(replicate, levels = c(1:4)),
      trap = factor(trap, levels = c(1:10)),
      year = factor(year, levels = c(2022, 2023)),
      season = factor(season, levels = c("March", "October"))
    ) %>%
    complete(site, replicate, trap, year, season) %>%
    replace(is.na(.), 0) %>%
    filter(
      paste(year, season) != "2022 March",
      paste(year, season) != "2023 October"
    ) %>%
    relocate(site, replicate) %>%
    arrange(year) %>%
    mutate(abundance = rowSums(across(where(is.numeric)))) %>%
    group_by(year, season, site, replicate) %>%
    summarise(abundance = mean(abundance), .groups = "drop")

  modern_species_div <- modern_species_richness %>%
    mutate(abundance = modern_species_abund$abundance) %>%
    mutate(
      veg_type = factor(case_when(
        site %in% c(1) ~ "strandveld",
        site %in% c(2, 6, 16) ~ "restioid",
        site %in% c(3, 4, 5) ~ "proteoid",
        site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
        site %in% c(11) ~ "alpine",
        site %in% c(17) ~ "succulent_karoo"
      ))
    ) %>%
    relocate(site, replicate, year, veg_type)
}

# combined years ##

{
  p1 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle("Species richness per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness)))

  p2 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle("Mean abundance per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance)))

  # separate years ##

  p3 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle("Species richness per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p4 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle("Mean abundance per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  # veg types combined year ##

  p5 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p6 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  # veg types separate year ##

  p7 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p8 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))
}


{
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/all_species/combined_years_richness.png"
    ),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/all_species/combined_years_abundance.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/all_species/seasonal_richness.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/all_species/seasonal_abundance.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/all_species/veg_type/combined_years_richness.png"
    ),
    plot = p5,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/all_species/veg_type/combined_years_abundance.png"
    ),
    plot = p6,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/all_species/veg_type/seasonal_richness.png"
    ),
    plot = p7,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/all_species/veg_type/seasonal_abundance.png"
    ),
    plot = p8,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(
    modern_species_abund,
    modern_species_richness,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    p7,
    p8
  )
}

# families ###

# carabids ###

{
  carabid_modern <- data.table::fread(
    "./working_beetle_data/family_groups/carabid_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  modern_species_richness <- carabid_modern %>%
    mutate(richness = specnumber(carabid_modern[, -(1:3)])) %>%
    select(1:3, richness)

  modern_species_abund <- carabid_modern %>%
    mutate(abundance = (rowSums(across(where(is.numeric))) / 10)) %>%
    select(1:3, abundance)

  modern_species_div <- modern_species_richness %>%
    mutate(abundance = modern_species_abund$abundance) %>%
    mutate(
      veg_type = factor(case_when(
        site %in% c(1) ~ "strandveld",
        site %in% c(2, 6, 16) ~ "restioid",
        site %in% c(3, 4, 5) ~ "proteoid",
        site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
        site %in% c(11) ~ "alpine",
        site %in% c(17) ~ "succulent_karoo"
      ))
    ) %>%
    relocate(site, replicate, year, veg_type)

  # combined years ##

  p1 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle("Carabid species richness per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness)))

  p2 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle("Carabid mean abundance per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance)))

  # separate years ##

  p3 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle("Carabid species richness per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p4 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle("Carabid mean abundance per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  # veg types combined year ##

  p5 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Carabid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p6 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Carabid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  # veg types separate year ##

  p7 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Carabid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p8 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Carabid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))
}

{
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/carabid/combined_years_richness.png"
    ),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/carabid/combined_years_abundance.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/carabid/seasonal_richness.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/carabid/seasonal_abundance.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/carabid/veg_type/combined_years_richness.png"
    ),
    plot = p5,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/carabid/veg_type/combined_years_abundance.png"
    ),
    plot = p6,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/carabid/veg_type/seasonal_richness.png"
    ),
    plot = p7,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/carabid/veg_type/seasonal_abundance.png"
    ),
    plot = p8,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(
    modern_species_abund,
    modern_species_div,
    modern_species_richness,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    p7,
    p8
  )
}

# tenebrionids ###

{
  tenebrionid_modern <- data.table::fread(
    "./working_beetle_data/family_groups/tenebrionid_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  modern_species_richness <- tenebrionid_modern %>%
    mutate(richness = specnumber(tenebrionid_modern[, -(1:3)])) %>%
    select(1:3, richness)

  modern_species_abund <- tenebrionid_modern %>%
    mutate(abundance = (rowSums(across(where(is.numeric))) / 10)) %>%
    select(1:3, abundance)

  modern_species_div <- modern_species_richness %>%
    mutate(abundance = modern_species_abund$abundance) %>%
    mutate(
      veg_type = factor(case_when(
        site %in% c(1) ~ "strandveld",
        site %in% c(2, 6, 16) ~ "restioid",
        site %in% c(3, 4, 5) ~ "proteoid",
        site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
        site %in% c(11) ~ "alpine",
        site %in% c(17) ~ "succulent_karoo"
      ))
    ) %>%
    relocate(site, replicate, year, veg_type)

  # combined years ##

  p1 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle(
      "Tenebrionid species richness per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness)))

  p2 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle(
      "Tenebrionid mean abundance per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance)))

  # separate years ##

  p3 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle(
      "Tenebrionid species richness per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p4 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle(
      "Tenebrionid mean abundance per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  # veg types combined year ##

  p5 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Tenebrionid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p6 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Tenebrionid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  # veg types separate year ##

  p7 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Tenebrionid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p8 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Tenebrionid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))
}

{
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/tenebrionid/combined_years_richness.png"
    ),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/tenebrionid/combined_years_abundance.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/tenebrionid/seasonal_richness.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/tenebrionid/seasonal_abundance.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/tenebrionid/veg_type/combined_years_richness.png"
    ),
    plot = p5,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/tenebrionid/veg_type/combined_years_abundance.png"
    ),
    plot = p6,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/tenebrionid/veg_type/seasonal_richness.png"
    ),
    plot = p7,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/tenebrionid/veg_type/seasonal_abundance.png"
    ),
    plot = p8,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(
    modern_species_abund,
    modern_species_div,
    modern_species_richness,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    p7,
    p8
  )
}

# cicindelids ###

{
  cicindelid_modern <- data.table::fread(
    "./working_beetle_data/family_groups/cicindelid_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  modern_species_richness <- cicindelid_modern %>%
    mutate(richness = specnumber(cicindelid_modern[, -(1:3)])) %>%
    select(1:3, richness)

  modern_species_abund <- cicindelid_modern %>%
    mutate(abundance = (rowSums(across(where(is.numeric))) / 10)) %>%
    select(1:3, abundance)

  modern_species_div <- modern_species_richness %>%
    mutate(abundance = modern_species_abund$abundance) %>%
    mutate(
      veg_type = factor(case_when(
        site %in% c(1) ~ "strandveld",
        site %in% c(2, 6, 16) ~ "restioid",
        site %in% c(3, 4, 5) ~ "proteoid",
        site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
        site %in% c(11) ~ "alpine",
        site %in% c(17) ~ "succulent_karoo"
      ))
    ) %>%
    relocate(site, replicate, year, veg_type)

  # combined years ##

  p1 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle(
      "Cicindelid species richness per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness)))

  p2 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle(
      "Cicindelid mean abundance per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance)))

  # separate years ##

  p3 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle(
      "Cicindelid species richness per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p4 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle(
      "Cicindelid mean abundance per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  # veg types combined year ##

  p5 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Cicindelid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p6 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Cicindelid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  # veg types separate year ##

  p7 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Cicindelid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p8 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Cicindelid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$abundance))) +
    scale_fill_manual(values = c("white", "darkgrey"))
}

{
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/cicindelid/combined_years_richness.png"
    ),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/cicindelid/combined_years_abundance.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/cicindelid/seasonal_richness.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/cicindelid/seasonal_abundance.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/cicindelid/veg_type/combined_years_richness.png"
    ),
    plot = p5,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/cicindelid/veg_type/combined_years_abundance.png"
    ),
    plot = p6,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/cicindelid/veg_type/seasonal_richness.png"
    ),
    plot = p7,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/seasonal_comparison/cicindelid/veg_type/seasonal_abundance.png"
    ),
    plot = p8,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(
    modern_species_abund,
    modern_species_div,
    modern_species_richness,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    p7,
    p8
  )
}

# werner ID barcharts ####

# all species ###

{
  modern_species_richness <- beet_wide %>%
    group_by(year, season, site, replicate) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop") %>%
    mutate(
      site = factor(site, levels = c(1:17)),
      replicate = factor(replicate, levels = c(1:4)),
      year = factor(year, levels = c(2022, 2023)),
      season = factor(season, levels = c("March", "October"))
    ) %>%
    complete(site, replicate, year, season) %>%
    replace(is.na(.), 0) %>%
    filter(
      paste(year, season) != "2022 March",
      paste(year, season) != "2023 October"
    ) %>%
    relocate(site, replicate) %>%
    arrange(year)

  modern_species_richness <- modern_species_richness %>%
    mutate(richness = specnumber(modern_species_richness[, -(1:4)])) %>%
    select(1:4, richness)

  modern_species_abund <- beet_wide %>%
    select(-ID) %>%
    mutate(
      site = factor(site, levels = c(1:17)),
      replicate = factor(replicate, levels = c(1:4)),
      trap = factor(trap, levels = c(1:10)),
      year = factor(year, levels = c(2022, 2023)),
      season = factor(season, levels = c("March", "October"))
    ) %>%
    complete(site, replicate, trap, year, season) %>%
    replace(is.na(.), 0) %>%
    filter(
      paste(year, season) != "2022 March",
      paste(year, season) != "2023 October"
    ) %>%
    relocate(site, replicate) %>%
    arrange(year) %>%
    mutate(abundance = rowSums(across(where(is.numeric)))) %>%
    group_by(year, season, site, replicate) %>%
    summarise(abundance = mean(abundance), .groups = "drop")

  modern_species_div <- modern_species_richness %>%
    mutate(abundance = modern_species_abund$abundance) %>%
    mutate(
      veg_type = factor(case_when(
        site %in% c(1) ~ "strandveld",
        site %in% c(2, 6, 16) ~ "restioid",
        site %in% c(3, 4, 5) ~ "proteoid",
        site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
        site %in% c(11) ~ "alpine",
        site %in% c(17) ~ "succulent_karoo"
      ))
    ) %>%
    relocate(site, replicate, year, veg_type)

  # combined years ##

  p1 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle("Species richness per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p2 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle("Mean abundance per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  # separate years ##

  p3 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness, fill = year)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge") +
    scale_fill_discrete(
      name = "Season",
      labels = c("2022" = "October 2022", "2023" = "March 2023")
    )

  p4 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance, fill = year)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle("") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge") +
    scale_fill_discrete(
      name = "Season",
      labels = c("2022" = "October 2022", "2023" = "March 2023")
    )

  # veg types combined year ##

  p5 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p6 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  # veg types separate year ##

  p7 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness, fill = year)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      ""
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge") +
    scale_fill_discrete(
      name = "Season",
      labels = c("2022" = "October 2022", "2023" = "March 2023")
    )

  p8 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance, fill = year)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      ""
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge") +
    scale_fill_discrete(
      name = "Season",
      labels = c("2022" = "October 2022", "2023" = "March 2023")
    )
}

gridExtra::grid.arrange(p3, p4)
gridExtra::grid.arrange(p7, p8)

{
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/all_species/combined_years_richness.png"
    ),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/all_species/combined_years_abundance.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/all_species/seasonal_richness.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/all_species/seasonal_abundance.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/all_species/veg_type/combined_years_richness.png"
    ),
    plot = p5,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/all_species/veg_type/combined_years_abundance.png"
    ),
    plot = p6,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/all_species/veg_type/seasonal_richness.png"
    ),
    plot = p7,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/all_species/veg_type/seasonal_abundance.png"
    ),
    plot = p8,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(
    modern_species_abund,
    modern_species_div,
    modern_species_richness,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    p7,
    p8
  )
}

# families ###

# carabids ###

{
  carabid_modern <- data.table::fread(
    "./working_beetle_data/family_groups/carabid_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  modern_species_richness <- carabid_modern %>%
    mutate(richness = specnumber(carabid_modern[, -(1:3)])) %>%
    select(1:3, richness)

  modern_species_abund <- carabid_modern %>%
    mutate(abundance = (rowSums(across(where(is.numeric))) / 10)) %>%
    select(1:3, abundance)

  modern_species_div <- modern_species_richness %>%
    mutate(abundance = modern_species_abund$abundance) %>%
    mutate(
      veg_type = factor(case_when(
        site %in% c(1) ~ "strandveld",
        site %in% c(2, 6, 16) ~ "restioid",
        site %in% c(3, 4, 5) ~ "proteoid",
        site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
        site %in% c(11) ~ "alpine",
        site %in% c(17) ~ "succulent_karoo"
      ))
    ) %>%
    relocate(site, replicate, year, veg_type)

  # combined years ##

  p1 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle("Carabid species richness per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p2 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle("Carabid mean abundance per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  # separate years ##

  p3 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle("Carabid species richness per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p4 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle("Carabid mean abundance per site (October 2022 and March 2023)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  # veg types combined year ##

  p5 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Carabid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p6 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Carabid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")
  # veg types separate year ##

  p7 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Carabid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p8 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Carabid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")
}

{
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/carabid/combined_years_richness.png"
    ),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/carabid/combined_years_abundance.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/carabid/seasonal_richness.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/carabid/seasonal_abundance.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/carabid/veg_type/combined_years_richness.png"
    ),
    plot = p5,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/carabid/veg_type/combined_years_abundance.png"
    ),
    plot = p6,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/carabid/veg_type/seasonal_richness.png"
    ),
    plot = p7,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/carabid/veg_type/seasonal_abundance.png"
    ),
    plot = p8,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(
    modern_species_abund,
    modern_species_div,
    modern_species_richness,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    p7,
    p8,
    carabid_modern
  )
}

# tenebrionids ###

{
  tenebrionid_modern <- data.table::fread(
    "./working_beetle_data/family_groups/tenebrionid_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  modern_species_richness <- tenebrionid_modern %>%
    mutate(richness = specnumber(tenebrionid_modern[, -(1:3)])) %>%
    select(1:3, richness)

  modern_species_abund <- tenebrionid_modern %>%
    mutate(abundance = (rowSums(across(where(is.numeric))) / 10)) %>%
    select(1:3, abundance)

  modern_species_div <- modern_species_richness %>%
    mutate(abundance = modern_species_abund$abundance) %>%
    mutate(
      veg_type = factor(case_when(
        site %in% c(1) ~ "strandveld",
        site %in% c(2, 6, 16) ~ "restioid",
        site %in% c(3, 4, 5) ~ "proteoid",
        site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
        site %in% c(11) ~ "alpine",
        site %in% c(17) ~ "succulent_karoo"
      ))
    ) %>%
    relocate(site, replicate, year, veg_type)

  # combined years ##

  p1 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle(
      "Tenebrionid species richness per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p2 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle(
      "Tenebrionid mean abundance per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  # separate years ##

  p3 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle(
      "Tenebrionid species richness per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p4 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle(
      "Tenebrionid mean abundance per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  # veg types combined year ##

  p5 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Tenebrionid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p6 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Tenebrionid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")
  # veg types separate year ##

  p7 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Tenebrionid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p8 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Tenebrionid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")
}

{
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/tenebrionid/combined_years_richness.png"
    ),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/tenebrionid/combined_years_abundance.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/tenebrionid/seasonal_richness.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/tenebrionid/seasonal_abundance.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/tenebrionid/veg_type/combined_years_richness.png"
    ),
    plot = p5,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/tenebrionid/veg_type/combined_years_abundance.png"
    ),
    plot = p6,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/tenebrionid/veg_type/seasonal_richness.png"
    ),
    plot = p7,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/tenebrionid/veg_type/seasonal_abundance.png"
    ),
    plot = p8,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(
    modern_species_abund,
    modern_species_div,
    modern_species_richness,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    p7,
    p8,
    tenebrionid_modern
  )
}

# cicindelids ###

{
  cicindelid_modern <- data.table::fread(
    "./working_beetle_data/family_groups/cicindelid_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  modern_species_richness <- cicindelid_modern %>%
    mutate(richness = specnumber(cicindelid_modern[, -(1:3)])) %>%
    select(1:3, richness)

  modern_species_abund <- cicindelid_modern %>%
    mutate(abundance = (rowSums(across(where(is.numeric))) / 10)) %>%
    select(1:3, abundance)

  modern_species_div <- modern_species_richness %>%
    mutate(abundance = modern_species_abund$abundance) %>%
    mutate(
      veg_type = factor(case_when(
        site %in% c(1) ~ "strandveld",
        site %in% c(2, 6, 16) ~ "restioid",
        site %in% c(3, 4, 5) ~ "proteoid",
        site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
        site %in% c(11) ~ "alpine",
        site %in% c(17) ~ "succulent_karoo"
      ))
    ) %>%
    relocate(site, replicate, year, veg_type)

  # combined years ##

  p1 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle(
      "Cicindelid species richness per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p2 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle(
      "Cicindelid mean abundance per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  # separate years ##

  p3 <- modern_species_div %>%
    ggplot(aes(x = site, y = richness, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Species richness") +
    ggtitle(
      "Cicindelid species richness per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p4 <- modern_species_div %>%
    ggplot(aes(x = site, y = abundance, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Altitude and Aspect", y = "Mean Abundance") +
    ggtitle(
      "Cicindelid mean abundance per site (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_discrete(labels = alt_labels) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  # veg types combined year ##

  p5 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness)) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Cicindelid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(modern_species_div$richness))) +
    scale_fill_manual(values = c("white", "darkgrey"))

  p6 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Cicindelid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  # veg types separate year ##

  p7 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = richness, fill = year)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Species richness") +
    ggtitle(
      "Cicindelid species richness per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")

  p8 <- modern_species_div %>%
    ggplot(aes(x = veg_type, y = abundance, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Vegetation Type", y = "Mean Abundance") +
    ggtitle(
      "Cicindelid mean abundance per vegetation type (October 2022 and March 2023)"
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_summary(geom = "col", position = "dodge", fun.data = "mean_se") +
    stat_summary(geom = "errorbar", fun.data = "mean_se", position = "dodge")
}

{
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/cicindelid/combined_years_richness.png"
    ),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/cicindelid/combined_years_abundance.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/cicindelid/seasonal_richness.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/cicindelid/seasonal_abundance.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/cicindelid/veg_type/combined_years_richness.png"
    ),
    plot = p5,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/cicindelid/veg_type/combined_years_abundance.png"
    ),
    plot = p6,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/cicindelid/veg_type/seasonal_richness.png"
    ),
    plot = p7,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_barcharts/seasonal_comparison/cicindelid/veg_type/seasonal_abundance.png"
    ),
    plot = p8,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(
    modern_species_abund,
    modern_species_div,
    modern_species_richness,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6,
    p7,
    p8,
    cicindelid_modern
  )
}

# werner ID MDS ####

# set-up ###

# mds ###

# presence-absence ##
# all sites #

{
  dat <- beet_wide %>%
    select(-ID) %>%
    mutate(
      veg_type = factor(case_when(
        site %in% c(1) ~ "strandveld",
        site %in% c(2, 6, 16) ~ "restioid",
        site %in% c(3, 4, 5) ~ "proteoid",
        site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
        site %in% c(11) ~ "alpine",
        site %in% c(17) ~ "succulent_karoo"
      ))
    ) %>%
    relocate(site, replicate, year, veg_type) %>%
    group_by(year, season, site, replicate, veg_type) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop")
  x <- dat[, 6:47] %>%
    mutate(decostand(across(where(is.numeric)), "pa")) %>%
    filter(rowSums(across(where(is.numeric))) != 0)
  y <- dat %>% filter(rowSums(across(where(is.numeric))) != 0) %>% select(1:5)
  Herb_community.mds <- x %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- y$site
  MDS_xy$replicate <- y$replicate
  MDS_xy$year <- y$year
  MDS_xy$veg_type <- y$veg_type
  Herb_community.mds$stress

  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Species presence by year") +
    scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = veg_type)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Species presence by vegetation type") +
    scale_color_discrete(name = "Vegetation Type")
  # p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Species presence by year with sites, overlaps jittered") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species presence by year with sites") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = veg_type)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species presence by vegetation type") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Vegetation Type",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Species presence by site") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Site",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
}

{
  ggsave(
    "./MDS_plots/modern_only/species/presence_absence/year_all.png",
    plot = p1,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/modern_only/species/presence_absence/vegtype_all.png",
    plot = p2,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/modern_only/species/presence_absence/year_all_sites.png",
    plot = p4,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/modern_only/species/presence_absence/vegtypes_all_sites.png",
    plot = p5,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/modern_only/species/presence_absence/sites.png",
    plot = p6,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )

  rm(x, y, Herb_community.mds, MDS_xy, p1, p2, p4, p5, p6)
}

dist <- vegdist(x, method = "bray")
anosim_type <- anosim(dist, y$veg_type)
summary(anosim_type)

# site-specific ##

for (k in 1:17) {
  x <- dat %>% filter(site == k)
  x1 <- x %>% filter(rowSums(across(where(is.numeric))) != 0) %>% select(6:47)
  x2 <- x %>% filter(rowSums(across(where(is.numeric))) != 0) %>% select(1:5)
  Herb_community.mds <- x1 %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- x2$site
  MDS_xy$replicate <- x2$replicate
  MDS_xy$year <- x2$year
  MDS_xy$veg_type <- x2$veg_type
  Herb_community.mds$stress

  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_jitter() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", species abundance by year")) +
    scale_color_discrete(name = "Year")

  ggsave(
    paste0(
      "./MDS_plots/modern_only/species/presence_absence/sites/site_",
      k,
      ".png"
    ),
    plot = pa,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  {
    if (k == 17) {
      rm(k, x, x1, x2, Herb_community.mds, MDS_xy, pa)
    }
  }
}

# abundance ##

{
  dat <- beet_wide %>%
    select(-ID) %>%
    mutate(
      veg_type = factor(case_when(
        site %in% c(1) ~ "strandveld",
        site %in% c(2, 6, 16) ~ "restioid",
        site %in% c(3, 4, 5) ~ "proteoid",
        site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
        site %in% c(11) ~ "alpine",
        site %in% c(17) ~ "succulent_karoo"
      ))
    ) %>%
    relocate(site, replicate, year, veg_type) %>%
    group_by(year, season, site, replicate, veg_type) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  x <- dat[, 6:47] %>% filter(rowSums(across(where(is.numeric))) != 0)
  y <- dat %>% filter(rowSums(across(where(is.numeric))) != 0) %>% select(1:5)
  Herb_community.mds <- x %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- y$site
  MDS_xy$replicate <- y$replicate
  MDS_xy$year <- y$year
  MDS_xy$veg_type <- y$veg_type
  Herb_community.mds$stress

  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Mean species abundance by year") +
    scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = veg_type)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Mean species abundance by vegetation type") +
    scale_color_discrete(name = "Vegetation type")
  # p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) + theme_bw() + geom_point(alpha = 0) + ggtitle("Species abundance by year with sites, overlaps jittered") + scale_color_discrete(name = "Year", guide = guide_legend(override.aes = list(alpha = 1))) + ggrepel::geom_text_repel(data = MDS_xy, mapping = aes(MDS1, MDS2, label = site), show.legend = FALSE, max.overlaps = 1000)
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Mean species abundance by year with sites") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = veg_type)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Mean species abundance by vegetation type") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Vegetation Type",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Mean species abundance by site") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Site",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
}

{
  ggsave(
    "./MDS_plots/modern_only/species/abundance/year_all.png",
    plot = p1,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/modern_only/species/abundance/period_all.png",
    plot = p2,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  # ggsave("./MDS_plots/species/abundance/year_all_sites_jitter.png", plot = p3, width = 4000, height = 2160, units = "px", bg = "white")
  ggsave(
    "./MDS_plots/modern_only/species/abundance/year_all_sites.png",
    plot = p4,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/modern_only/species/abundance/period_all_sites.png",
    plot = p5,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/modern_only/species/abundance/sites_jitter.png",
    plot = p6,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )

  rm(x, y, Herb_community.mds, MDS_xy, p1, p2, p4, p5, p6)
}

dist <- vegdist(fam_presence_values, method = "bray")
anosim_type <- anosim(dist, fam_presence_cats$site)
summary(anosim_type)

# site-specific ##

for (k in 1:17) {
  x <- dat %>% filter(site == k)
  x1 <- x %>% filter(rowSums(across(where(is.numeric))) != 0) %>% select(6:47)
  x2 <- x %>% filter(rowSums(across(where(is.numeric))) != 0) %>% select(1:5)
  Herb_community.mds <- x1 %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- x2$site
  MDS_xy$replicate <- x2$replicate
  MDS_xy$year <- x2$year
  MDS_xy$veg_type <- x2$veg_type
  Herb_community.mds$stress

  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_jitter() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", mean species abundance by year")) +
    scale_color_discrete(name = "Year")

  ggsave(
    paste0("./MDS_plots/modern_only/species/abundance/sites/site_", k, ".png"),
    plot = pa,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  {
    if (k == 17) {
      rm(k, x, x1, x2, Herb_community.mds, MDS_xy, pa)
    }
  }
}

# werner ID MDS family ####

{
  fam_carabid_modern <- data.table::fread(
    "./family_groups/carabids_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  fam_tenebrionid_modern <- data.table::fread(
    "./family_groups/tenebrionids_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  fam_cicindelid_modern <- data.table::fread(
    "./family_groups/cicindelids_modern.csv"
  ) %>%
    mutate(across(c(site, replicate, year), factor))

  {
    fam_carabid_modern_sums <- fam_carabid_modern %>%
      mutate("carabids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, carabids)
    fam_carabid_old_sums <- fam_carabid_old %>%
      mutate("carabids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, carabids)
    fam_carabid_sums <- rbind(fam_carabid_old_sums, fam_carabid_modern_sums)
    rm(fam_carabid_modern_sums, fam_carabid_old_sums)

    fam_tenebrionid_modern_sums <- fam_tenebrionid_modern %>%
      mutate("tenebrionids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, tenebrionids)
    fam_tenebrionid_old_sums <- fam_tenebrionid_old %>%
      mutate("tenebrionids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, tenebrionids)
    fam_tenebrionid_sums <- rbind(
      fam_tenebrionid_old_sums,
      fam_tenebrionid_modern_sums
    )
    rm(fam_tenebrionid_modern_sums, fam_tenebrionid_old_sums)

    fam_cicindelid_modern_sums <- fam_cicindelid_modern %>%
      mutate("cicindelids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, cicindelids)
    fam_cicindelid_old_sums <- fam_cicindelid_old %>%
      mutate("cicindelids" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, cicindelids)
    fam_cicindelid_sums <- rbind(
      fam_cicindelid_old_sums,
      fam_cicindelid_modern_sums
    )
    rm(fam_cicindelid_modern_sums, fam_cicindelid_old_sums)

    fam_sums <- inner_join(
      fam_carabid_sums,
      fam_cicindelid_sums,
      by = c("site", "replicate", "year")
    )
    fam_sums <- inner_join(
      fam_sums,
      fam_tenebrionid_sums,
      by = c("site", "replicate", "year")
    )
    rm(fam_cicindelid_sums, fam_carabid_sums, fam_tenebrionid_sums)
  }

  {
    x1 <- fam_carabid_modern %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("carabid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, carabid_diversity)
    x2 <- fam_carabid_old %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("carabid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, carabid_diversity)
    x <- rbind(x2, x1)

    y1 <- fam_cicindelid_modern %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("cicindelid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, cicindelid_diversity)
    y2 <- fam_cicindelid_old %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("cicindelid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, cicindelid_diversity)
    y <- rbind(y2, y1)

    z1 <- fam_tenebrionid_modern %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("tenebrionid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, tenebrionid_diversity)
    z2 <- fam_tenebrionid_old %>%
      mutate(decostand(across(where(is.numeric)), "pa")) %>%
      mutate("tenebrionid_diversity" = rowSums(across(where(is.numeric)))) %>%
      select(site, replicate, year, tenebrionid_diversity)
    z <- rbind(z2, z1)

    fam_presence <- inner_join(x, z, by = c("site", "replicate", "year"))
    fam_presence <- inner_join(
      fam_presence,
      y,
      by = c("site", "replicate", "year")
    )

    rm(x1, x2, x, y1, y2, y, z1, z2, z)
  }

  fam_sums <- fam_sums %>%
    mutate(
      period = factor(case_when(
        year %in% c(2002, 2003) ~ "old",
        year %in% c(2022, 2023) ~ "modern"
      ))
    ) %>%
    relocate(site, replicate, year, period)

  fam_presence <- fam_presence %>%
    mutate(
      period = factor(case_when(
        year %in% c(2002, 2003) ~ "old",
        year %in% c(2022, 2023) ~ "modern"
      ))
    ) %>%
    relocate(site, replicate, year, period)

  rm(
    fam_carabid_modern,
    fam_carabid_old,
    fam_cicindelid_modern,
    fam_cicindelid_old,
    fam_tenebrionid_modern,
    fam_tenebrionid_old
  )
}

# mds trials ###

# presence-absence ##
# all sites #

{
  fam_presence_values <- fam_presence[, 5:7] %>%
    mutate(decostand(across(where(is.numeric)), "pa")) %>%
    filter(rowSums(across(where(is.numeric))) != 0)
  fam_presence_cats <- fam_presence %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_presence_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress

  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family presence by year") +
    scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family presence by sampling period") +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family presence by year with sites, overlaps jittered") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family presence by year with sites") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family presence by sampling period") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023"),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family presence by site") +
    scale_color_discrete(
      name = "Site",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
}

{
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/year_all.png",
    plot = p1,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/period_all.png",
    plot = p2,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/year_all_sites_jitter.png",
    plot = p3,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/year_all_sites.png",
    plot = p4,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/period_all_sites.png",
    plot = p5,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/presence_absence/sites_jitter.png",
    plot = p6,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )

  rm(
    fam_presence_values,
    fam_presence_cats,
    Herb_community.mds,
    MDS_xy,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6
  )
}

dist <- vegdist(fam_presence_values, method = "bray")
anosim_type <- anosim(dist, fam_presence_cats$site)
summary(anosim_type)

# site-specific #

for (k in 1:17) {
  x <- fam_presence %>%
    filter(site == k) %>%
    mutate(decostand(across(where(is.numeric)), "pa"))
  fam_presence_values <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(5:7)
  fam_presence_cats <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_presence_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress

  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_jitter() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family presence by year")) +
    scale_color_discrete(name = "Year")
  (pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_jitter(width = 0.02, height = 0.02) +
    stat_ellipse() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family presence by sampling period")) +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    ))

  ggsave(
    paste0("./MDS_plots/family_grouped/presence_absence/year/site_", k, ".png"),
    plot = pa,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./MDS_plots/family_grouped/presence_absence/period/site_",
      k,
      ".png"
    ),
    plot = pb,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  {
    if (k == 17) {
      rm(
        k,
        x,
        fam_presence_values,
        fam_presence_cats,
        Herb_community.mds,
        MDS_xy,
        pa,
        pb
      )
    }
  }
}

# random

MDS_xy %>%
  filter(year == "2022") %>%
  ggplot(., aes(MDS1, MDS2, color = site)) +
  theme_bw() +
  geom_point(alpha = 0) +
  ggtitle("Family presence by site") +
  scale_color_discrete(
    name = "Year",
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  geom_text(
    data = MDS_xy,
    mapping = aes(MDS1, MDS2, label = site),
    show.legend = FALSE,
    max.overlaps = 1000
  ) +
  stat_ellipse()


# diversity ##
# all sites #

{
  fam_presence_values <- fam_presence[, 5:7] %>%
    filter(rowSums(across(where(is.numeric))) != 0)
  fam_presence_cats <- fam_presence %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_presence_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress

  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family diversity by year") +
    scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family diversity by sampling period") +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family diversity by year with sites, overlaps jittered") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family diversity by year with sites") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family diversity by sampling period") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023"),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family diversity by site") +
    scale_color_discrete(
      name = "Site",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
}

{
  ggsave(
    "./MDS_plots/family_grouped/diversity/year_all.png",
    plot = p1,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/diversity/period_all.png",
    plot = p2,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/diversity/year_all_sites_jitter.png",
    plot = p3,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/diversity/year_all_sites.png",
    plot = p4,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/diversity/period_all_sites.png",
    plot = p5,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/diversity/sites_jitter.png",
    plot = p6,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )

  rm(
    fam_presence_values,
    fam_presence_cats,
    Herb_community.mds,
    MDS_xy,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6
  )
}

# dist <- vegdist(fam_presence_values, method = "bray")
# anosim_type <- anosim(dist,fam_presence_cats$site)
# summary(anosim_type)

# site-specific #

for (k in 1:17) {
  x <- fam_presence %>% filter(site == k)
  fam_presence_values <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(5:7)
  fam_presence_cats <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_presence_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_presence_cats$site
  MDS_xy$replicate <- fam_presence_cats$replicate
  MDS_xy$year <- fam_presence_cats$year
  MDS_xy$period <- fam_presence_cats$period
  Herb_community.mds$stress

  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_jitter() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family diversity by year")) +
    scale_color_discrete(name = "Year")
  pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_jitter(width = 0.02, height = 0.02) +
    stat_ellipse() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family diversity by sampling period")) +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )

  ggsave(
    paste0("./MDS_plots/family_grouped/diversity/year/site_", k, ".png"),
    plot = pa,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0("./MDS_plots/family_grouped/diversity/period/site_", k, ".png"),
    plot = pb,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  {
    if (k == 17) {
      rm(
        k,
        x,
        fam_presence_values,
        fam_presence_cats,
        Herb_community.mds,
        MDS_xy,
        pa,
        pb
      )
    }
  }
}

# abundance ##
# all sites #

{
  fam_sums_values <- fam_sums[, 5:7] %>%
    filter(rowSums(across(where(is.numeric))) != 0)
  fam_sums_cats <- fam_sums %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_sums_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_sums_cats$site
  MDS_xy$replicate <- fam_sums_cats$replicate
  MDS_xy$year <- fam_sums_cats$year
  MDS_xy$period <- fam_sums_cats$period
  Herb_community.mds$stress

  p1 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family abundance by year") +
    scale_color_discrete(name = "Year")
  p2 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_point() +
    stat_ellipse() +
    theme_bw() +
    ggtitle("Family abundance by sampling period") +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
  p3 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family abundance by year with sites, overlaps jittered") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
  p4 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family abundance by year with sites") +
    scale_color_discrete(
      name = "Year",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p5 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family abundance by sampling period") +
    stat_ellipse() +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023"),
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    geom_text(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE
    )
  p6 <- ggplot(MDS_xy, aes(MDS1, MDS2, color = site)) +
    theme_bw() +
    geom_point(alpha = 0) +
    ggtitle("Family abundance by site") +
    scale_color_discrete(
      name = "Site",
      guide = guide_legend(override.aes = list(alpha = 1))
    ) +
    ggrepel::geom_text_repel(
      data = MDS_xy,
      mapping = aes(MDS1, MDS2, label = site),
      show.legend = FALSE,
      max.overlaps = 1000
    )
}

{
  ggsave(
    "./MDS_plots/family_grouped/abundance/year_all.png",
    plot = p1,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/abundance/period_all.png",
    plot = p2,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/abundance/year_all_sites_jitter.png",
    plot = p3,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/abundance/year_all_sites.png",
    plot = p4,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/abundance/period_all_sites.png",
    plot = p5,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )
  ggsave(
    "./MDS_plots/family_grouped/abundance/sites_jitter.png",
    plot = p6,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white"
  )

  rm(
    fam_sums_values,
    fam_sums_cats,
    Herb_community.mds,
    MDS_xy,
    p1,
    p2,
    p3,
    p4,
    p5,
    p6
  )
}

dist <- vegdist(fam_sums_values, method = "bray")
anosim_type <- anosim(dist, fam_sums_cats$site)
summary(anosim_type)

# site-specific #

for (k in 1:17) {
  x <- fam_sums %>% filter(site == k)
  fam_sums_values <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(5:7)
  fam_sums_cats <- x %>%
    filter(rowSums(across(where(is.numeric))) != 0) %>%
    select(1:4)
  Herb_community.mds <- fam_sums_values %>%
    metaMDS(comm = ., distance = "bray", trace = FALSE, autotransform = FALSE)
  MDS_xy <- data.frame(Herb_community.mds$points)
  MDS_xy$site <- fam_sums_cats$site
  MDS_xy$replicate <- fam_sums_cats$replicate
  MDS_xy$year <- fam_sums_cats$year
  MDS_xy$period <- fam_sums_cats$period
  Herb_community.mds$stress

  pa <- ggplot(MDS_xy, aes(MDS1, MDS2, color = year)) +
    geom_jitter() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family abundance by year")) +
    scale_color_discrete(name = "Year")
  (pb <- ggplot(MDS_xy, aes(MDS1, MDS2, color = period)) +
    geom_jitter(width = 0.02, height = 0.02) +
    stat_ellipse() +
    theme_bw() +
    ggtitle(paste0("Site ", k, ", family abundance by sampling period")) +
    scale_color_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    ))

  ggsave(
    paste0("./MDS_plots/family_grouped/abundance/year/site_", k, ".png"),
    plot = pa,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0("./MDS_plots/family_grouped/abundance/period/site_", k, ".png"),
    plot = pb,
    width = 4000,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  {
    if (k == 17) {
      rm(
        k,
        x,
        fam_sumse_values,
        fam_sums_cats,
        Herb_community.mds,
        MDS_xy,
        pa,
        pb
      )
    }
  }
}

# boxplots of all years ####

# anthia ###

{
  abund_box_data <- abund_all %>%
    mutate(
      anthia_decemguttata = anthia_decemguttata / 10,
      stenocara_dentata = stenocara_dentata / 10,
      zophosis_gracilicornis = zophosis_gracilicornis / 10
    ) %>%
    filter(paste(year, season) != "2023 March")

  p1 <- abund_box_data %>%
    ggplot(aes(x = year, y = anthia_decemguttata)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Mean Abundance") +
    ggtitle("Anthia decemguttata mean abundance per year") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$anthia_decemguttata)))

  p2 <- abund_box_data %>%
    ggplot(aes(x = site, y = anthia_decemguttata, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Site", y = "Mean Abundance") +
    ggtitle("Anthia decemguttata mean abundance per year") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$anthia_decemguttata))) +
    scale_fill_discrete(name = "Year")

  p3 <- abund_box_data %>%
    mutate(period = factor(period, levels = c("old", "modern"))) %>%
    ggplot(aes(x = period, y = anthia_decemguttata)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Mean Abundance") +
    ggtitle("Anthia decemguttata mean abundance per sampling period") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$anthia_decemguttata)))

  p4 <- abund_box_data %>%
    mutate(period = factor(period, levels = c("old", "modern"))) %>%
    ggplot(aes(x = site, y = anthia_decemguttata, fill = period)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Site", y = "Mean Abundance") +
    ggtitle("Anthia decemguttata mean abundance per samplling period") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$anthia_decemguttata))) +
    scale_fill_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
}

{
  ggsave(
    paste0("./beetle_boxplots/historical_comparison/anthia/by_year_all.png"),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0("./beetle_boxplots/historical_comparison/anthia/by_year_sites.png"),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0("./beetle_boxplots/historical_comparison/anthia/by_period_all.png"),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/anthia/by_period_sites.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
}

# stenocara ###

{
  abund_box_data <- abund_all %>%
    mutate(
      anthia_decemguttata = anthia_decemguttata / 10,
      stenocara_dentata = stenocara_dentata / 10,
      zophosis_gracilicornis = zophosis_gracilicornis / 10
    ) %>%
    filter(paste(year, season) != "2023 March")

  p1 <- abund_box_data %>%
    ggplot(aes(x = year, y = stenocara_dentata)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Mean Abundance") +
    ggtitle("Stenocara dentata mean abundance per year") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$stenocara_dentata)))

  p2 <- abund_box_data %>%
    ggplot(aes(x = site, y = stenocara_dentata, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Site", y = "Mean Abundance") +
    ggtitle("Stenocara dentata mean abundance per year") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$stenocara_dentata))) +
    scale_fill_discrete(name = "Year")

  p3 <- abund_box_data %>%
    mutate(period = factor(period, levels = c("old", "modern"))) %>%
    ggplot(aes(x = period, y = stenocara_dentata)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Mean Abundance") +
    ggtitle("Stenocara dentata mean abundance per sampling period") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$stenocara_dentata)))

  p4 <- abund_box_data %>%
    mutate(period = factor(period, levels = c("old", "modern"))) %>%
    ggplot(aes(x = site, y = stenocara_dentata, fill = period)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Site", y = "Mean Abundance") +
    ggtitle("Stenocara dentata mean abundance per samplling period") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$stenocara_dentata))) +
    scale_fill_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
}

{
  ggsave(
    paste0("./beetle_boxplots/historical_comparison/stenocara/by_year_all.png"),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/stenocara/by_year_sites.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/stenocara/by_period_all.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/stenocara/by_period_sites.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
}

# zophosis ###

{
  abund_box_data <- abund_all %>%
    mutate(
      anthia_decemguttata = anthia_decemguttata / 10,
      stenocara_dentata = stenocara_dentata / 10,
      zophosis_gracilicornis = zophosis_gracilicornis / 10
    ) %>%
    filter(paste(year, season) != "2023 March")

  p1 <- abund_box_data %>%
    ggplot(aes(x = year, y = zophosis_gracilicornis)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Mean Abundance") +
    ggtitle("Zophosis gracilicornis mean abundance per year") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(
      limits = c(0, max(abund_box_data$zophosis_gracilicornis))
    )

  p2 <- abund_box_data %>%
    ggplot(aes(x = site, y = zophosis_gracilicornis, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Site", y = "Mean Abundance") +
    ggtitle("Zophosis gracilicornis mean abundance per year") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(
      limits = c(0, max(abund_box_data$zophosis_gracilicornis))
    ) +
    scale_fill_discrete(name = "Year")

  p3 <- abund_box_data %>%
    mutate(period = factor(period, levels = c("old", "modern"))) %>%
    ggplot(aes(x = period, y = zophosis_gracilicornis)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Mean Abundance") +
    ggtitle("Zophosis gracilicornis mean abundance per sampling period") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(
      limits = c(0, max(abund_box_data$zophosis_gracilicornis))
    )

  p4 <- abund_box_data %>%
    mutate(period = factor(period, levels = c("old", "modern"))) %>%
    ggplot(aes(x = site, y = zophosis_gracilicornis, fill = period)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Site", y = "Mean Abundance") +
    ggtitle("Zophosis gracilicornis mean abundance per samplling period") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(
      limits = c(0, max(abund_box_data$zophosis_gracilicornis))
    ) +
    scale_fill_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
}

{
  ggsave(
    paste0("./beetle_boxplots/historical_comparison/zophosis/by_year_all.png"),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/zophosis/by_year_sites.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/zophosis/by_period_all.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/zophosis/by_period_sites.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
}

# all species ###

{
  abund_box_data <- abund_all %>%
    mutate(
      anthia_decemguttata = anthia_decemguttata / 10,
      stenocara_dentata = stenocara_dentata / 10,
      zophosis_gracilicornis = zophosis_gracilicornis / 10
    ) %>%
    filter(paste(year, season) != "2023 March") %>%
    mutate(
      abundance = anthia_decemguttata +
        stenocara_dentata +
        zophosis_gracilicornis
    )

  p1 <- abund_box_data %>%
    ggplot(aes(x = year, y = abundance)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Mean Abundance") +
    ggtitle("All species mean abundance per year") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$abundance)))

  p2 <- abund_box_data %>%
    ggplot(aes(x = site, y = abundance, fill = year)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Site", y = "Mean Abundance") +
    ggtitle("All species mean abundance per year") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$abundance))) +
    scale_fill_discrete(name = "Year")

  p3 <- abund_box_data %>%
    mutate(period = factor(period, levels = c("old", "modern"))) %>%
    ggplot(aes(x = period, y = abundance)) +
    geom_boxplot() +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Mean Abundance") +
    ggtitle("All species mean abundance per sampling period") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$abundance)))

  p4 <- abund_box_data %>%
    mutate(period = factor(period, levels = c("old", "modern"))) %>%
    ggplot(aes(x = site, y = abundance, fill = period)) +
    geom_boxplot(position = "dodge") +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Site", y = "Mean Abundance") +
    ggtitle("All species mean abundance per samplling period") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_y_continuous(limits = c(0, max(abund_box_data$abundance))) +
    scale_fill_discrete(
      name = "Sampling period",
      labels = c("old" = "2002/2003", "modern" = "2022/2023")
    )
}

{
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/all_species/by_year_all.png"
    ),
    plot = p1,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/all_species/by_year_sites.png"
    ),
    plot = p2,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/all_species/by_period_all.png"
    ),
    plot = p3,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  ggsave(
    paste0(
      "./beetle_boxplots/historical_comparison/all_species/by_period_sites.png"
    ),
    plot = p4,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )

  rm(p1, p2, p3, p4, abund_box_data)
}

# mixed models for difference significance in historical comparisons ####

# whole transect ##

abund_all <- abund_all %>%
  mutate(
    veg_type = factor(case_when(
      site %in% c(1) ~ "strandveld",
      site %in% c(2, 6, 16) ~ "restioid",
      site %in% c(3, 4, 5) ~ "proteoid",
      site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
      site %in% c(11) ~ "alpine",
      site %in% c(17) ~ "succulent_karoo"
    ))
  ) %>%
  mutate(
    site_rep = as.factor(paste(site, replicate)),
    abundance = as.integer(rowSums(across(
      .cols = c(anthia_decemguttata, stenocara_dentata, zophosis_gracilicornis)
    )))
  )
abund_all$period <- relevel(abund_all$period, ref = 2)

model_all <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = abundance ~ year,
    random = ~site,
    family = "poisson"
  )
summary(model_all)

model_anthia <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = anthia_decemguttata ~ year,
    random = ~site,
    family = "poisson"
  )
summary(model_anthia)

model_stenocara <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = stenocara_dentata ~ year,
    random = ~site,
    family = "poisson"
  )
summary(model_stenocara)

model_zophosis <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = zophosis_gracilicornis ~ year,
    random = ~site,
    family = "poisson"
  )
summary(model_zophosis)

model_all_period <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = abundance ~ period,
    random = ~site,
    family = "poisson"
  )
summary(model_all_period)

model_anthia_period <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = anthia_decemguttata ~ period,
    random = ~site,
    family = "poisson"
  )
summary(model_anthia_period)

model_stenocara_period <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = stenocara_dentata ~ period,
    random = ~site,
    family = "poisson"
  )
summary(model_stenocara_period)

model_zophosis_period <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = zophosis_gracilicornis ~ period,
    random = ~site,
    family = "poisson"
  )
summary(model_zophosis_period)

# veg type ##

model_veg <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = anthia_decemguttata ~ (year * veg_type),
    random = ~site,
    family = "poisson"
  )
summary(model_veg)

# site-by-site ##

model_anthia_sites <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = anthia_decemguttata ~ (period * veg_type),
    random = ~site,
    family = "poisson"
  )
summary(model_anthia_sites)

model_stenocara_sites <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = stenocara_dentata ~ (period * veg_type),
    random = ~site,
    family = "poisson"
  )
summary(model_stenocara_sites)

model_zophosis_sites <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = zophosis_gracilicornis ~ (period * veg_type),
    random = ~site,
    family = "poisson"
  )
summary(model_zophosis_sites)

abund_all %>% group_by(year) %>% summarise(sum(anthia_decemguttata))
abund_all %>% group_by(year) %>% summarise(sum(stenocara_dentata))
abund_all %>% group_by(year) %>% summarise(sum(zophosis_gracilicornis))

# results writing ####

beet_wide_results <- beet_wide %>%
  group_by(year) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE)

x <- beet_wide %>%
  filter(site == 11) %>%
  group_by(year) %>%
  summarise_if(is.numeric, sum, na.rm = TRUE)

beet_moder_summer_exclusive <- beet_wide_results %>% pivot_longer(cols = 2:43)
beet_wide_results %>%
  group_by(year) %>%
  mutate(abund = rowSums(across(where(is.numeric)))) %>%
  select(year, abund)
