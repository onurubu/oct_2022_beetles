#### Set-up ####
# packages

{
  library(tidyverse)
  library(weights)
  library(emmeans)
  library(EnvStats)
  library(sf)
  library(hms)
  library(vegan)
}

{
  comm <- abund_all %>%
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
    group_by(year, season, period, veg_type, site) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop") %>%
    relocate(year, season, period, veg_type, site) %>%
    arrange(site, year)
  comm <- comm %>% mutate(site = alt_labels[comm$site])
  comm[1:6, 1:6]
  comm$aspect <- sapply(strsplit(comm$site, "m", fixed = TRUE), "[[", 2)
  comm$aspect <- gsub("(Peak)", "W", comm$aspect)
  comm$aspect <- gsub("\\(", "", comm$aspect)
  comm$aspect <- gsub("\\)", "", comm$aspect)
  comm$site <- as.numeric(sapply(
    strsplit(comm$site, "m", fixed = TRUE),
    "[[",
    1
  ))
  siteinfo <- comm[, c(1:5, ncol(comm))]
  comm <- comm[, !colnames(comm) %in% colnames(siteinfo)]
}

#### Calculate elevation range in each year using a loop ----------------

spp.elev <- NULL
for (i in colnames(comm)) {
  dat <- siteinfo
  x <- pull(comm[, i])
  dat$abund <- x
  dat$abund <- ifelse(dat$abund > 1, 1, dat$abund)
  dat <- dat[dat$abund > 0, ]

  sampling <- length(unique(dat$period))

  if (sampling == 1) {
    tmp <- data.frame(
      species = i,
      max_E_02 = NA,
      max_W_02 = NA,
      min_E_02 = NA,
      min_W_02 = NA,
      max_E_22 = NA,
      max_W_22 = NA,
      min_E_22 = NA,
      min_W_22 = NA
    )
  } else {
    ## old
    dat02 <- dat[dat$period == "old", ]
    dat02 <- dat02[dat02$abund > 0, ]
    maxes <- tapply(dat02$site, dat02$aspect, max)
    mins <- tapply(dat02$site, dat02$aspect, min)

    asp.names <- names(maxes)
    if (length(asp.names) == 1) {
      if (names(maxes) == "W") {
        maxes <- c(NA, maxes)
        mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA)
        mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_02", "max_W_02")
    names(mins) <- c("min_E_02", "min_W_02")
    tmp02 <- t(data.frame(c(maxes, mins)))

    ## modern
    dat22 <- dat[dat$period == "modern", ]
    dat22 <- dat22[dat22$abund > 0, ]
    maxes <- tapply(dat22$site, dat22$aspect, max)
    mins <- tapply(dat22$site, dat22$aspect, min)

    asp.names <- names(maxes)
    if (length(asp.names) == 1) {
      if (names(maxes) == "W") {
        maxes <- c(NA, maxes)
        mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA)
        mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_22", "max_W_22")
    names(mins) <- c("min_E_22", "min_W_22")
    tmp22 <- t(data.frame(c(maxes, mins)))

    tmp <- data.frame(species = i, cbind(tmp02, tmp22))
  }

  spp.elev <- rbind(spp.elev, tmp)
  print(i)
}

rownames(spp.elev) <- 1:nrow(spp.elev)
spp.elev

#### Shifts + summaries + plot start here --------------------------------

### 1. Compute shifts ----
spp.elev$shift_max_E <- spp.elev$max_E_22 - spp.elev$max_E_02
spp.elev$shift_max_W <- spp.elev$max_W_22 - spp.elev$max_W_02
spp.elev$shift_min_E <- spp.elev$min_E_22 - spp.elev$min_E_02
spp.elev$shift_min_W <- spp.elev$min_W_22 - spp.elev$min_W_02

### 2. CLEAN SUMMARY TABLE (per species) ----
summary_table <- spp.elev %>%
  select(species, shift_max_E, shift_max_W, shift_min_E, shift_min_W) %>%
  mutate(
    max_E_dir = case_when(
      shift_max_E > 0 ~ "Upslope",
      shift_max_E < 0 ~ "Downslope",
      shift_max_E == 0 ~ "No change",
      TRUE ~ NA
    ),
    max_W_dir = case_when(
      shift_max_W > 0 ~ "Upslope",
      shift_max_W < 0 ~ "Downslope",
      shift_max_W == 0 ~ "No change",
      TRUE ~ NA
    ),
    min_E_dir = case_when(
      shift_min_E > 0 ~ "Upslope",
      shift_min_E < 0 ~ "Downslope",
      shift_min_E == 0 ~ "No change",
      TRUE ~ NA
    ),
    min_W_dir = case_when(
      shift_min_W > 0 ~ "Upslope",
      shift_min_W < 0 ~ "Downslope",
      shift_min_W == 0 ~ "No change",
      TRUE ~ NA
    )
  )

write.csv(
  summary_table,
  "./range_shifts/elevational_shift_summary.csv",
  row.names = FALSE
)

### 3. PROPORTIONS upslope / downslope / no change ----
count_dir <- function(x) {
  c(
    Upslope = sum(x > 0, na.rm = TRUE),
    Downslope = sum(x < 0, na.rm = TRUE),
    No_change = sum(x == 0, na.rm = TRUE)
  )
}

shift_summary <- data.frame(
  Upper_East = count_dir(spp.elev$shift_max_E),
  Upper_West = count_dir(spp.elev$shift_max_W),
  Lower_East = count_dir(spp.elev$shift_min_E),
  Lower_West = count_dir(spp.elev$shift_min_W)
)

shift_summary
round(prop.table(as.matrix(shift_summary), 2) * 100, 1) # % per column

### 4. Reshape to long format for ONE plot ----
spp_long <- spp.elev %>%
  select(species, shift_max_E, shift_max_W, shift_min_E, shift_min_W) %>%
  pivot_longer(
    cols = starts_with("shift_"),
    names_to = c("type", "aspect"),
    names_pattern = "shift_(max|min)_(E|W)",
    values_to = "shift"
  ) %>%
  mutate(
    type = ifelse(type == "max", "Upper limit", "Lower limit"),
    aspect = ifelse(aspect == "E", "East", "West")
  )

### 5. ONE combined plot ----
p_all <- ggplot(
  spp_long,
  aes(x = reorder(species, shift), y = shift, fill = aspect)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_wrap(~type, ncol = 1, scales = "free_x") +
  labs(
    title = "Beetle elevational range shifts (2002–2023)",
    x = "Species",
    y = "Shift in elevation (m)",
    fill = "Aspect"
  ) +
  theme_classic(base_size = 12) +
  scale_x_discrete(
    labels = c(
      "anthia_decemguttata" = "Anthia decemguttata",
      "stenocara_dentata" = "Stenocara dentata",
      "zophosis_gracilicornis" = "Zophosis gracilicornis"
    )
  )

ggsave(
  paste0("./range_shifts/all_slopes_shifts.png"),
  plot = p_all,
  width = 3840,
  height = 2160,
  units = "px",
  bg = "white",
  create.dir = TRUE
)

rm(
  dat,
  dat02,
  dat22,
  tmp,
  tmp02,
  tmp22,
  asp.names,
  i,
  maxes,
  mins,
  x,
  sampling,
  spp_long,
  spp.elev,
  count_dir
)
rm(shift_summary, summary_table)

# rm(dat, dat02, dat22, tmp, tmp02, tmp22, asp.names, i, maxes, mins, x, sampling, no_shift_both, not_plotted_both, not_plotted_east, not_plotted_west, shift, spp_long, spp_long_filt, spp.elev, count_dir)

#### Use this code for slopes in separate plots and no shift species removed
#### Calculate elevation range in each year using a loop ----------------

spp.elev <- NULL
for (i in colnames(comm)) {
  dat <- siteinfo
  x <- pull(comm[, i])
  dat$abund <- x
  dat$abund <- ifelse(dat$abund > 1, 1, dat$abund)
  dat <- dat[dat$abund > 0, ]

  sampling <- length(unique(dat$period))

  if (sampling == 1) {
    tmp <- data.frame(
      species = i,
      max_E_02 = NA,
      max_W_02 = NA,
      min_E_02 = NA,
      min_W_02 = NA,
      max_E_22 = NA,
      max_W_22 = NA,
      min_E_22 = NA,
      min_W_22 = NA
    )
  } else {
    ## old
    dat02 <- dat[dat$period == "old", ]
    dat02 <- dat02[dat02$abund > 0, ]
    maxes <- tapply(dat02$site, dat02$aspect, max)
    mins <- tapply(dat02$site, dat02$aspect, min)

    asp.names <- names(maxes)
    if (length(asp.names) == 1) {
      if (names(maxes) == "W") {
        maxes <- c(NA, maxes)
        mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA)
        mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_02", "max_W_02")
    names(mins) <- c("min_E_02", "min_W_02")
    tmp02 <- t(data.frame(c(maxes, mins)))

    ## modern
    dat22 <- dat[dat$period == "modern", ]
    dat22 <- dat22[dat22$abund > 0, ]
    maxes <- tapply(dat22$site, dat22$aspect, max)
    mins <- tapply(dat22$site, dat22$aspect, min)

    asp.names <- names(maxes)
    if (length(asp.names) == 1) {
      if (names(maxes) == "W") {
        maxes <- c(NA, maxes)
        mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA)
        mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_22", "max_W_22")
    names(mins) <- c("min_E_22", "min_W_22")
    tmp22 <- t(data.frame(c(maxes, mins)))

    tmp <- data.frame(species = i, cbind(tmp02, tmp22))
  }

  spp.elev <- rbind(spp.elev, tmp)
  print(i)
}

rownames(spp.elev) <- 1:nrow(spp.elev)
spp.elev

#### Shifts + summaries + plot start here --------------------------------
### 1. Compute shifts ----
spp.elev$shift_max_E <- spp.elev$max_E_22 - spp.elev$max_E_02
spp.elev$shift_max_W <- spp.elev$max_W_22 - spp.elev$max_W_02
spp.elev$shift_min_E <- spp.elev$min_E_22 - spp.elev$min_E_02
spp.elev$shift_min_W <- spp.elev$min_W_22 - spp.elev$min_W_02

# Species not plotted
not_plotted_east <- spp.elev %>%
  filter(
    (is.na(shift_min_E) | shift_min_E == 0) &
      (is.na(shift_max_E) | shift_max_E == 0)
  ) %>%
  select(species)

not_plotted_east

not_plotted_west <- spp.elev %>%
  filter(
    (is.na(shift_min_W) | shift_min_W == 0) &
      (is.na(shift_max_W) | shift_max_W == 0)
  ) %>%
  select(species)

not_plotted_west

not_plotted_both <- spp.elev %>%
  filter(
    (is.na(shift_min_E) | shift_min_E == 0) &
      (is.na(shift_max_E) | shift_max_E == 0) &
      (is.na(shift_min_W) | shift_min_W == 0) &
      (is.na(shift_max_W) | shift_max_W == 0)
  ) %>%
  select(species)

not_plotted_both

write.csv(
  not_plotted_east,
  "./range_shifts/not_plotted_east.csv",
  row.names = FALSE
)
write.csv(
  not_plotted_west,
  "./range_shifts/not_plotted_west.csv",
  row.names = FALSE
)
write.csv(
  not_plotted_both,
  "./range_shifts/not_plotted_both_slopes.csv",
  row.names = FALSE
)

### 2. CLEAN SUMMARY TABLE (per species) ----
summary_table <- spp.elev %>%
  select(species, shift_max_E, shift_max_W, shift_min_E, shift_min_W) %>%
  mutate(
    max_E_dir = case_when(
      shift_max_E > 0 ~ "Upslope",
      shift_max_E < 0 ~ "Downslope",
      shift_max_E == 0 ~ "No change",
      TRUE ~ NA
    ),
    max_W_dir = case_when(
      shift_max_W > 0 ~ "Upslope",
      shift_max_W < 0 ~ "Downslope",
      shift_max_W == 0 ~ "No change",
      TRUE ~ NA
    ),
    min_E_dir = case_when(
      shift_min_E > 0 ~ "Upslope",
      shift_min_E < 0 ~ "Downslope",
      shift_min_E == 0 ~ "No change",
      TRUE ~ NA
    ),
    min_W_dir = case_when(
      shift_min_W > 0 ~ "Upslope",
      shift_min_W < 0 ~ "Downslope",
      shift_min_W == 0 ~ "No change",
      TRUE ~ NA
    )
  )

View(summary_table)
write.csv(
  summary_table,
  "./range_shifts/elevational_shift_summary.csv",
  row.names = FALSE
)

### 3. PROPORTIONS upslope / downslope / no change ----
count_dir <- function(x) {
  c(
    Upslope = sum(x > 0, na.rm = TRUE),
    Downslope = sum(x < 0, na.rm = TRUE),
    No_change = sum(x == 0, na.rm = TRUE)
  )
}

shift_summary <- data.frame(
  Upper_East = count_dir(spp.elev$shift_max_E),
  Upper_West = count_dir(spp.elev$shift_max_W),
  Lower_East = count_dir(spp.elev$shift_min_E),
  Lower_West = count_dir(spp.elev$shift_min_W)
)

shift_summary
round(prop.table(as.matrix(shift_summary), 2) * 100, 1) # % per column

### 4. Reshape to long format ----
spp_long <- spp.elev %>%
  select(species, shift_max_E, shift_max_W, shift_min_E, shift_min_W) %>%
  pivot_longer(
    cols = starts_with("shift_"),
    names_to = c("type", "aspect"),
    names_pattern = "shift_(max|min)_(E|W)",
    values_to = "shift"
  ) %>%
  mutate(
    type = ifelse(type == "max", "Upper limit", "Lower limit"),
    aspect = ifelse(aspect == "E", "East", "West")
  )

## 4b. Exclude species with no change (shift = 0 or NA)
spp_long_filt <- spp_long %>%
  filter(!is.na(shift), shift != 0)

#no shift
no_shift_east <- spp.elev %>%
  filter(!is.na(shift_max_E), !is.na(shift_min_E)) %>%
  filter(shift_max_E == 0 & shift_min_E == 0) %>%
  pull(species)

no_shift_east

no_shift_west <- spp.elev %>%
  filter(!is.na(shift_max_W), !is.na(shift_min_W)) %>%
  filter(shift_max_W == 0 & shift_min_W == 0) %>%
  pull(species)

no_shift_west

no_shift_both <- spp.elev %>%
  filter(
    shift_max_E == 0,
    shift_min_E == 0,
    shift_max_W == 0,
    shift_min_W == 0
  ) %>%
  select(species)

no_shift_both

write.csv(
  data.frame(species = no_shift_east),
  "./range_shifts/no_shift_east_species.csv",
  row.names = FALSE
)

write.csv(
  data.frame(species = no_shift_west),
  "./range_shifts/no_shift_west_species.csv",
  row.names = FALSE
)

write.csv(
  no_shift_both,
  "./range_shifts/species_no_shift_both_slopes.csv",
  row.names = FALSE
)

### 5. Separate plots for East and West ----

## East slope only
p_east <- ggplot(
  spp_long_filt %>% filter(aspect == "East"),
  aes(x = reorder(species, shift), y = shift, fill = type) # colour by upper vs lower, or keep "aspect" if you prefer
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_wrap(~type, ncol = 1, scales = "free_x") +
  labs(
    title = "Eastern slope (old and modern)",
    x = "Species",
    y = "Shift in elevation (m)",
    fill = "Limit"
  ) +
  theme_classic(base_size = 12)

p_east # print East plot


## West slope only
p_west <- ggplot(
  spp_long_filt %>% filter(aspect == "West"),
  aes(x = reorder(species, shift), y = shift, fill = type)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_wrap(~type, ncol = 1, scales = "free_x") +
  labs(
    title = "West slope (old and modern)",
    x = "Species",
    y = "Shift in elevation (m)",
    fill = "Limit"
  ) +
  theme_classic(base_size = 12)

p_west # print West plot
ggsave(
  "./range_shifts/West_slope_shift.png",
  p_west,
  width = 3840,
  height = 2160,
  units = "px",
  bg = "white",
  create.dir = TRUE
)
ggsave(
  "./range_shifts/East_slope_shift.png",
  p_east,
  width = 3840,
  height = 2160,
  units = "px",
  bg = "white",
  create.dir = TRUE
)

rm(
  dat,
  dat02,
  dat22,
  tmp,
  tmp02,
  tmp22,
  asp.names,
  i,
  maxes,
  mins,
  x,
  sampling,
  spp_long,
  spp.elev,
  count_dir,
  p_all,
  p_east,
  p_west,
  no_shift_both,
  no_shift_east,
  not_plotted_west,
  no_shift_west,
  not_plotted_both,
  not_plotted_east,
  spp_long_filt
)
rm(shift_summary, summary_table)

# modern seasonal comparison ####

{
  comm <- beet_wide %>%
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
    group_by(year, season, veg_type, site) %>%
    summarise(across(where(is.numeric), sum), .groups = "drop") %>%
    relocate(year, season, veg_type, site) %>%
    arrange(site, year)
  comm <- comm %>% mutate(site = alt_labels[comm$site])
  comm[1:6, 1:6]
  comm$aspect <- sapply(strsplit(comm$site, "m", fixed = TRUE), "[[", 2)
  comm$aspect <- gsub("(Peak)", "W", comm$aspect)
  comm$aspect <- gsub("\\(", "", comm$aspect)
  comm$aspect <- gsub("\\)", "", comm$aspect)
  comm$site <- as.numeric(sapply(
    strsplit(comm$site, "m", fixed = TRUE),
    "[[",
    1
  ))
  siteinfo <- comm[, c(1:4, ncol(comm))]
  comm <- comm[, !colnames(comm) %in% colnames(siteinfo)]
}

#### Calculate elevation range in each year using a loop ----------------

spp.elev <- NULL
for (i in colnames(comm)) {
  dat <- siteinfo
  x <- pull(comm[, i])
  dat$abund <- x
  dat$abund <- ifelse(dat$abund > 1, 1, dat$abund)
  dat <- dat[dat$abund > 0, ]

  sampling <- length(unique(dat$season))

  if (sampling == 1) {
    tmp <- data.frame(
      species = i,
      max_E_02 = NA,
      max_W_02 = NA,
      min_E_02 = NA,
      min_W_02 = NA,
      max_E_22 = NA,
      max_W_22 = NA,
      min_E_22 = NA,
      min_W_22 = NA
    )
  } else {
    ## old
    dat02 <- dat[dat$season == "October", ]
    dat02 <- dat02[dat02$abund > 0, ]
    maxes <- tapply(dat02$site, dat02$aspect, max)
    mins <- tapply(dat02$site, dat02$aspect, min)

    asp.names <- names(maxes)
    if (length(asp.names) == 1) {
      if (names(maxes) == "W") {
        maxes <- c(NA, maxes)
        mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA)
        mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_02", "max_W_02")
    names(mins) <- c("min_E_02", "min_W_02")
    tmp02 <- t(data.frame(c(maxes, mins)))

    ## modern
    dat22 <- dat[dat$season == "March", ]
    dat22 <- dat22[dat22$abund > 0, ]
    maxes <- tapply(dat22$site, dat22$aspect, max)
    mins <- tapply(dat22$site, dat22$aspect, min)

    asp.names <- names(maxes)
    if (length(asp.names) == 1) {
      if (names(maxes) == "W") {
        maxes <- c(NA, maxes)
        mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA)
        mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_22", "max_W_22")
    names(mins) <- c("min_E_22", "min_W_22")
    tmp22 <- t(data.frame(c(maxes, mins)))

    tmp <- data.frame(species = i, cbind(tmp02, tmp22))
  }

  spp.elev <- rbind(spp.elev, tmp)
  print(i)
}

rownames(spp.elev) <- 1:nrow(spp.elev)
spp.elev

#### Shifts + summaries + plot start here --------------------------------

### 1. Compute shifts ----
spp.elev$shift_max_E <- spp.elev$max_E_22 - spp.elev$max_E_02
spp.elev$shift_max_W <- spp.elev$max_W_22 - spp.elev$max_W_02
spp.elev$shift_min_E <- spp.elev$min_E_22 - spp.elev$min_E_02
spp.elev$shift_min_W <- spp.elev$min_W_22 - spp.elev$min_W_02

### 2. CLEAN SUMMARY TABLE (per species) ----
summary_table <- spp.elev %>%
  select(species, shift_max_E, shift_max_W, shift_min_E, shift_min_W) %>%
  mutate(
    max_E_dir = case_when(
      shift_max_E > 0 ~ "Upslope",
      shift_max_E < 0 ~ "Downslope",
      shift_max_E == 0 ~ "No change",
      TRUE ~ NA
    ),
    max_W_dir = case_when(
      shift_max_W > 0 ~ "Upslope",
      shift_max_W < 0 ~ "Downslope",
      shift_max_W == 0 ~ "No change",
      TRUE ~ NA
    ),
    min_E_dir = case_when(
      shift_min_E > 0 ~ "Upslope",
      shift_min_E < 0 ~ "Downslope",
      shift_min_E == 0 ~ "No change",
      TRUE ~ NA
    ),
    min_W_dir = case_when(
      shift_min_W > 0 ~ "Upslope",
      shift_min_W < 0 ~ "Downslope",
      shift_min_W == 0 ~ "No change",
      TRUE ~ NA
    )
  )

write.csv(
  summary_table,
  "./range_shifts/modern_seasonal/elevational_shift_summary.csv",
  row.names = FALSE
)

### 3. PROPORTIONS upslope / downslope / no change ----
count_dir <- function(x) {
  c(
    Upslope = sum(x > 0, na.rm = TRUE),
    Downslope = sum(x < 0, na.rm = TRUE),
    No_change = sum(x == 0, na.rm = TRUE)
  )
}

shift_summary <- data.frame(
  Upper_East = count_dir(spp.elev$shift_max_E),
  Upper_West = count_dir(spp.elev$shift_max_W),
  Lower_East = count_dir(spp.elev$shift_min_E),
  Lower_West = count_dir(spp.elev$shift_min_W)
)

shift_summary
round(prop.table(as.matrix(shift_summary), 2) * 100, 1) # % per column

### 4. Reshape to long format for ONE plot ----
spp_long <- spp.elev %>%
  select(species, shift_max_E, shift_max_W, shift_min_E, shift_min_W) %>%
  pivot_longer(
    cols = starts_with("shift_"),
    names_to = c("type", "aspect"),
    names_pattern = "shift_(max|min)_(E|W)",
    values_to = "shift"
  ) %>%
  mutate(
    type = ifelse(type == "max", "Upper limit", "Lower limit"),
    aspect = ifelse(aspect == "E", "East", "West")
  )

### 5. ONE combined plot ----
p_all <- ggplot(
  spp_long,
  aes(x = reorder(species, shift), y = shift, fill = aspect)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_wrap(~type, ncol = 1, scales = "free_x") +
  labs(
    title = "Beetle elevational range shifts (2002–2023)",
    x = "Species",
    y = "Shift in elevation (m)",
    fill = "Aspect"
  ) +
  theme_classic(base_size = 12)

ggsave(
  paste0("./range_shifts/modern_seasonal/all_slopes_shifts.png"),
  plot = p_all,
  width = 3840,
  height = 2160,
  units = "px",
  bg = "white",
  create.dir = TRUE
)

rm(
  dat,
  dat02,
  dat22,
  tmp,
  tmp02,
  tmp22,
  asp.names,
  i,
  maxes,
  mins,
  x,
  sampling,
  spp_long,
  spp.elev,
  count_dir
)
rm(shift_summary, summary_table)

# rm(dat, dat02, dat22, tmp, tmp02, tmp22, asp.names, i, maxes, mins, x, sampling, no_shift_both, not_plotted_both, not_plotted_east, not_plotted_west, shift, spp_long, spp_long_filt, spp.elev, count_dir)

#### Use this code for slopes in separate plots and no shift species removed
#### Calculate elevation range in each year using a loop ----------------

spp.elev <- NULL
for (i in colnames(comm)) {
  dat <- siteinfo
  x <- pull(comm[, i])
  dat$abund <- x
  dat$abund <- ifelse(dat$abund > 1, 1, dat$abund)
  dat <- dat[dat$abund > 0, ]

  sampling <- length(unique(dat$season))

  if (sampling == 1) {
    tmp <- data.frame(
      species = i,
      max_E_02 = NA,
      max_W_02 = NA,
      min_E_02 = NA,
      min_W_02 = NA,
      max_E_22 = NA,
      max_W_22 = NA,
      min_E_22 = NA,
      min_W_22 = NA
    )
  } else {
    ## old
    dat02 <- dat[dat$season == "October", ]
    dat02 <- dat02[dat02$abund > 0, ]
    maxes <- tapply(dat02$site, dat02$aspect, max)
    mins <- tapply(dat02$site, dat02$aspect, min)

    asp.names <- names(maxes)
    if (length(asp.names) == 1) {
      if (names(maxes) == "W") {
        maxes <- c(NA, maxes)
        mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA)
        mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_02", "max_W_02")
    names(mins) <- c("min_E_02", "min_W_02")
    tmp02 <- t(data.frame(c(maxes, mins)))

    ## modern
    dat22 <- dat[dat$season == "March", ]
    dat22 <- dat22[dat22$abund > 0, ]
    maxes <- tapply(dat22$site, dat22$aspect, max)
    mins <- tapply(dat22$site, dat22$aspect, min)

    asp.names <- names(maxes)
    if (length(asp.names) == 1) {
      if (names(maxes) == "W") {
        maxes <- c(NA, maxes)
        mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA)
        mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_22", "max_W_22")
    names(mins) <- c("min_E_22", "min_W_22")
    tmp22 <- t(data.frame(c(maxes, mins)))

    tmp <- data.frame(species = i, cbind(tmp02, tmp22))
  }

  spp.elev <- rbind(spp.elev, tmp)
  print(i)
}

rownames(spp.elev) <- 1:nrow(spp.elev)
spp.elev

#### Shifts + summaries + plot start here --------------------------------
### 1. Compute shifts ----
spp.elev$shift_max_E <- spp.elev$max_E_22 - spp.elev$max_E_02
spp.elev$shift_max_W <- spp.elev$max_W_22 - spp.elev$max_W_02
spp.elev$shift_min_E <- spp.elev$min_E_22 - spp.elev$min_E_02
spp.elev$shift_min_W <- spp.elev$min_W_22 - spp.elev$min_W_02

# Species not plotted
not_plotted_east <- spp.elev %>%
  filter(
    (is.na(shift_min_E) | shift_min_E == 0) &
      (is.na(shift_max_E) | shift_max_E == 0)
  ) %>%
  select(species)

not_plotted_east

not_plotted_west <- spp.elev %>%
  filter(
    (is.na(shift_min_W) | shift_min_W == 0) &
      (is.na(shift_max_W) | shift_max_W == 0)
  ) %>%
  select(species)

not_plotted_west

not_plotted_both <- spp.elev %>%
  filter(
    (is.na(shift_min_E) | shift_min_E == 0) &
      (is.na(shift_max_E) | shift_max_E == 0) &
      (is.na(shift_min_W) | shift_min_W == 0) &
      (is.na(shift_max_W) | shift_max_W == 0)
  ) %>%
  select(species)

not_plotted_both

write.csv(
  not_plotted_east,
  "./range_shifts/modern_seasonal/not_plotted_east.csv",
  row.names = FALSE
)
write.csv(
  not_plotted_west,
  "./range_shifts/modern_seasonal/not_plotted_west.csv",
  row.names = FALSE
)
write.csv(
  not_plotted_both,
  "./range_shifts/modern_seasonal/not_plotted_both_slopes.csv",
  row.names = FALSE
)

### 2. CLEAN SUMMARY TABLE (per species) ----
summary_table <- spp.elev %>%
  select(species, shift_max_E, shift_max_W, shift_min_E, shift_min_W) %>%
  mutate(
    max_E_dir = case_when(
      shift_max_E > 0 ~ "Upslope",
      shift_max_E < 0 ~ "Downslope",
      shift_max_E == 0 ~ "No change",
      TRUE ~ NA
    ),
    max_W_dir = case_when(
      shift_max_W > 0 ~ "Upslope",
      shift_max_W < 0 ~ "Downslope",
      shift_max_W == 0 ~ "No change",
      TRUE ~ NA
    ),
    min_E_dir = case_when(
      shift_min_E > 0 ~ "Upslope",
      shift_min_E < 0 ~ "Downslope",
      shift_min_E == 0 ~ "No change",
      TRUE ~ NA
    ),
    min_W_dir = case_when(
      shift_min_W > 0 ~ "Upslope",
      shift_min_W < 0 ~ "Downslope",
      shift_min_W == 0 ~ "No change",
      TRUE ~ NA
    )
  )

View(summary_table)
write.csv(
  summary_table,
  "./range_shifts/modern_seasonal/elevational_shift_summary.csv",
  row.names = FALSE
)

### 3. PROPORTIONS upslope / downslope / no change ----
count_dir <- function(x) {
  c(
    Upslope = sum(x > 0, na.rm = TRUE),
    Downslope = sum(x < 0, na.rm = TRUE),
    No_change = sum(x == 0, na.rm = TRUE)
  )
}

shift_summary <- data.frame(
  Upper_East = count_dir(spp.elev$shift_max_E),
  Upper_West = count_dir(spp.elev$shift_max_W),
  Lower_East = count_dir(spp.elev$shift_min_E),
  Lower_West = count_dir(spp.elev$shift_min_W)
)

shift_summary
round(prop.table(as.matrix(shift_summary), 2) * 100, 1) # % per column

### 4. Reshape to long format ----
spp_long <- spp.elev %>%
  select(species, shift_max_E, shift_max_W, shift_min_E, shift_min_W) %>%
  pivot_longer(
    cols = starts_with("shift_"),
    names_to = c("type", "aspect"),
    names_pattern = "shift_(max|min)_(E|W)",
    values_to = "shift"
  ) %>%
  mutate(
    type = ifelse(type == "max", "Upper limit", "Lower limit"),
    aspect = ifelse(aspect == "E", "East", "West")
  )

## 4b. Exclude species with no change (shift = 0 or NA)
spp_long_filt <- spp_long %>%
  filter(!is.na(shift), shift != 0)

#no shift
no_shift_east <- spp.elev %>%
  filter(!is.na(shift_max_E), !is.na(shift_min_E)) %>%
  filter(shift_max_E == 0 & shift_min_E == 0) %>%
  pull(species)

no_shift_east

no_shift_west <- spp.elev %>%
  filter(!is.na(shift_max_W), !is.na(shift_min_W)) %>%
  filter(shift_max_W == 0 & shift_min_W == 0) %>%
  pull(species)

no_shift_west

no_shift_both <- spp.elev %>%
  filter(
    shift_max_E == 0,
    shift_min_E == 0,
    shift_max_W == 0,
    shift_min_W == 0
  ) %>%
  select(species)

no_shift_both

write.csv(
  data.frame(species = no_shift_east),
  "./range_shifts/modern_seasonal/no_shift_east_species.csv",
  row.names = FALSE
)

write.csv(
  data.frame(species = no_shift_west),
  "./range_shifts/modern_seasonal/no_shift_west_species.csv",
  row.names = FALSE
)

write.csv(
  no_shift_both,
  "./range_shifts/modern_seasonal/species_no_shift_both_slopes.csv",
  row.names = FALSE
)

### 5. Separate plots for East and West ----

## East slope only
p_east <- ggplot(
  spp_long_filt %>% filter(aspect == "East"),
  aes(x = reorder(species, shift), y = shift, fill = type) # colour by upper vs lower, or keep "aspect" if you prefer
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_wrap(~type, ncol = 1, scales = "free_x") +
  labs(
    title = "Eastern slope (October 2022 vs March 2023)",
    x = "Species",
    y = "Shift in elevation (m)",
    fill = "Limit"
  ) +
  scale_x_discrete(
    labels = c(
      "anthia_decemguttata" = "Anthia decemguttata",
      "stenocara_dentata" = "Stenocara dentata",
      "zophosis_gracilicornis" = "Zophosis gracilicornis"
    )
  ) +
  theme_classic(base_size = 15) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 25)
  ) +
  labs(tag = "(b)", size = 20) +
  theme(plot.tag.position = c(0.05, 0.98), plot.tag = element_text(size = 25))

p_east # print East plot


## West slope only
p_west <- ggplot(
  spp_long_filt %>% filter(aspect == "West"),
  aes(x = reorder(species, shift), y = shift, fill = type)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_wrap(~type, ncol = 1, scales = "free_x") +
  labs(
    title = "West slope (October 2022 vs March 2023)",
    x = "Species",
    y = "Shift in elevation (m)",
    fill = "Limit"
  ) +
  scale_x_discrete(
    labels = c(
      "anthia_decemguttata" = expression(paste(italic("Anthia decemguttata"))),
      "zophosis_gracilicornis" = expression(paste(italic(
        "Zophosis gracilicornis"
      ))),
      "Cryptochile.sp.1" = expression(paste(italic("Cryptochile"), " sp. 1")),
      "Gonogenia.sp..1" = expression(paste(italic("Gonogenia"), " sp. 1"))
    )
  ) +
  theme_classic(base_size = 15) +
  theme(axis.text.y = element_text(face = "italic")) +
  theme(
    legend.text = element_text(size = 25),
    legend.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.title = element_text(size = 25)
  ) +
  labs(tag = "(a)") +
  theme(plot.tag.position = c(0.05, 0.98), plot.tag = element_text(size = 25))

p_west # print West plot

gridExtra::grid.arrange(p_west, p_east)

### test
df <- data.frame(x = 1:10, y = 1:10)
(base <- ggplot(df, aes(x, y)) +
  geom_blank() +
  theme_bw())

# Full panel annotation
base +
  annotation_custom(
    grob = geom_text(aes("lol")),
    xmin = -Inf,
    xmax = Inf,
    ymin = -Inf,
    ymax = Inf
  )


ggsave(
  "./range_shifts/modern_seasonal/West_slope_shift.png",
  p_west,
  width = 3840,
  height = 2160,
  units = "px",
  bg = "white",
  create.dir = TRUE
)
ggsave(
  "./range_shifts/modern_seasonal/East_slope_shift.png",
  p_east,
  width = 3840,
  height = 2160,
  units = "px",
  bg = "white",
  create.dir = TRUE
)


rm(
  dat,
  dat02,
  dat22,
  tmp,
  tmp02,
  tmp22,
  asp.names,
  i,
  maxes,
  mins,
  x,
  sampling,
  spp_long,
  spp.elev,
  count_dir,
  p_all,
  p_east,
  p_west,
  no_shift_both,
  no_shift_east,
  not_plotted_west,
  no_shift_west,
  not_plotted_both,
  not_plotted_east,
  spp_long_filt
)
rm(shift_summary, summary_table)
rm(comm, siteinfo)
