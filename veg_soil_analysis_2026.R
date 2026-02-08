##### Coskun Kucukkaragoz
#### January 2026
### Veg and soil data analysis

# veg cover ####

veg_cover <- data.table::fread("./veg_cover_data/veg_cover_all.csv") %>%
  mutate(across(
    c(site, replicate, trap, year, season),
    factor
  ))


veg_cover_comp <- veg_cover %>%
  filter(
    paste(year, season) != "2023 March",
    paste(year, season) != "2003 March"
  ) %>%
  group_by(year, season, site, replicate) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  mutate(period = abund_all$period) %>%
  relocate(site, replicate, year, season, period)

veg_cover_season <- veg_cover %>%
  filter(
    paste(year, season) != "2002 October",
    paste(year, season) != "2023 October",
    paste(year, season) != "2003 October",
    paste(year, season) != "2003 March"
  )

veg_cover_comp_long <- veg_cover_comp %>%
  pivot_longer(cols = c(bare_ground, litter, rock, veg)) %>%
  mutate(
    period = factor(case_when(
      year %in% c(2002, 2003) ~ "old",
      year %in% c(2022, 2023) ~ "modern"
    ))
  ) %>%
  mutate(period = factor(period, levels = c("old", "modern"))) %>%
  relocate(year, season, site, replicate, period)

# year plots ###

{
  for (k in 1:17) {
    p <- veg_cover_comp_long %>%
      filter(site == k) %>%
      ggplot(aes(x = year, y = (value / 4), fill = name)) +
      geom_bar(position = "stack", stat = "identity") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(x = "Year", y = "Ground Composition (%)") +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(alt_labels[k]) +
      scale_fill_manual(
        values = c("#674533", "#B47A3B", "#9B9F9E", "#409B44"),
        name = "Ground cover",
        labels = c("Bare Ground", "Litter", "Exposed Rock", "Live Vegetation")
      )
    ggsave(
      paste0("./veg_plots/cover/years/", k, "_veg_cover_years.png"),
      plot = p,
      width = 3840,
      height = 2160,
      units = "px",
      bg = "white",
      create.dir = TRUE
    )
  }

  p <- veg_cover_comp_long %>%
    ggplot(aes(x = year, y = (value / 68), fill = name)) +
    geom_bar(position = "stack", stat = "identity") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Ground Composition (%)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("All Sites") +
    scale_fill_manual(
      values = c("#674533", "#B47A3B", "#9B9F9E", "#409B44"),
      name = "Ground cover",
      labels = c("Bare Ground", "Litter", "Exposed Rock", "Live Vegetation")
    )
  ggsave(
    paste0("./veg_plots/cover/years/all_veg_cover_years.png"),
    plot = p,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(k, p)
}

# period plots ###

{
  for (k in 1:17) {
    p <- veg_cover_comp_long %>%
      filter(site == k) %>%
      ggplot(aes(x = period, y = (value / 4), fill = name)) +
      geom_bar(position = "stack", stat = "identity") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(x = "Sampling Period", y = "Ground Composition (%)") +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(alt_labels[k]) +
      scale_fill_manual(
        values = c("#674533", "#B47A3B", "#9B9F9E", "#409B44"),
        name = "Ground cover",
        labels = c("Bare Ground", "Litter", "Exposed Rock", "Live Vegetation")
      ) +
      scale_x_discrete(labels = c("2002/2003", "2022/2023"))
    ggsave(
      paste0("./veg_plots/cover/period/", k, "_veg_cover_period.png"),
      plot = p,
      width = 3840,
      height = 2160,
      units = "px",
      bg = "white",
      create.dir = TRUE
    )
  }

  p <- veg_cover_comp_long %>%
    ggplot(aes(x = period, y = (value / 68), fill = name)) +
    geom_bar(position = "stack", stat = "identity") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Sampling Period", y = "Ground Composition (%)") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("All Sites") +
    scale_fill_manual(
      values = c("#674533", "#B47A3B", "#9B9F9E", "#409B44"),
      name = "Ground cover",
      labels = c("Bare Ground", "Litter", "Exposed Rock", "Live Vegetation")
    ) +
    scale_x_discrete(labels = c("2002/2003", "2022/2023"))
  ggsave(
    paste0("./veg_plots/cover/period/all_veg_cover_years.png"),
    plot = p,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(k, p)
}

model_veg_years <- veg_cover_comp_long %>%
  MCMCglmm(
    data = .,
    fixed = value ~ (year * name),
    random = ~site,
  )
summary(model_veg_years)

model_veg_period <- veg_cover_comp_long %>%
  MCMCglmm(
    data = .,
    fixed = value ~ (period * name),
    random = ~site,
  )
summary(model_veg_period)


Sys.sleep(0.1)


# soil ####

{
  soil_dat_02 <- data.table::fread("./soil/soil_data_2002_working.csv") %>%
    mutate(across(c(site, replicate, year), factor)) %>%
    rename("rock_soil" = "rock")
  soil_dat_22 <- data.table::fread("./soil/soil_data_2022_working.csv") %>%
    mutate(across(c(site, replicate, year), factor)) %>%
    rename("rock_soil" = "rock")

  soil_dat_02 <- rbind(soil_dat_02, soil_dat_02)
  soil_dat_22 <- rbind(soil_dat_22, soil_dat_22)

  soil_dat_02 <- soil_dat_02 %>%
    mutate(
      year = as.factor(rep(c(2002, 2003), each = (nrow(soil_dat_02) / 2))),
      season = as.factor(rep("October", nrow(.))),
      period = as.factor(rep("old", nrow(.)))
    ) %>%
    relocate(site, replicate, year, season, period)

  soil_dat_22 <- soil_dat_22 %>%
    mutate(
      year = as.factor(rep(c(2022, 2023), each = (nrow(soil_dat_22) / 2))),
      season = as.factor(rep("October", nrow(.))),
      period = as.factor(rep("modern", nrow(.)))
    ) %>%
    relocate(site, replicate, year, season, period)

  soil_dat <- rbind(soil_dat_02, soil_dat_22)
  rm(soil_dat_02, soil_dat_22)
}

# veg height ####

veg_height <- data.table::fread("./veg_height_data/veg_height_tothits.csv") %>%
  mutate(across(
    c(site, replicate, year, season),
    factor
  ))


veg_height_comp <- veg_height %>%
  filter(
    paste(year, season) != "2023 March",
    paste(year, season) != "2003 March"
  ) %>%
  group_by(year, season, site, replicate) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop") %>%
  mutate(period = abund_all$period) %>%
  relocate(site, replicate, year, season, period) %>%
  mutate(period = factor(period, levels = c("old", "modern")))

veg_height_season <- veg_height %>%
  filter(
    paste(year, season) != "2002 October",
    paste(year, season) != "2023 October",
    paste(year, season) != "2003 October",
    paste(year, season) != "2003 March"
  )

# tothits ###
# year plots ##

{
  for (k in 1:17) {
    p <- veg_height_comp %>%
      filter(site == k) %>%
      ggplot(aes(x = year, y = tothits)) +
      geom_boxplot() +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(x = "Year", y = "Vertical vegetation complexity") +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(alt_labels[k])
    ggsave(
      paste0("./veg_plots/height/tothits/years/", k, "_veg_height_years.png"),
      plot = p,
      width = 3840,
      height = 2160,
      units = "px",
      bg = "white",
      create.dir = TRUE
    )
  }

  p <- veg_height_comp %>%
    ggplot(aes(x = year, y = tothits)) +
    geom_boxplot() +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Vertical vegetation complexity") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("All Sites")
  ggsave(
    paste0("./veg_plots/height/tothits/years/all_veg_height_years.png"),
    plot = p,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(k, p)
}

# period plots ##

{
  for (k in 1:17) {
    p <- veg_height_comp %>%
      filter(site == k) %>%
      ggplot(aes(x = period, y = tothits)) +
      geom_boxplot() +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(x = "Sampling Period", y = "Vertical vegetation complexity") +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(alt_labels[k]) +
      scale_x_discrete(labels = c("2002/2003", "2022/2023"))
    ggsave(
      paste0("./veg_plots/height/tothits/period/", k, "_veg_height_period.png"),
      plot = p,
      width = 3840,
      height = 2160,
      units = "px",
      bg = "white",
      create.dir = TRUE
    )
  }

  p <- veg_height_comp %>%
    ggplot(aes(x = period, y = tothits)) +
    geom_boxplot() +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Sampling Period", y = "Vertical vegetation complexity") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("All Sites") +
    scale_x_discrete(labels = c("2002/2003", "2022/2023"))
  ggsave(
    paste0("./veg_plots/height/tothits/period/all_veg_height_period.png"),
    plot = p,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(k, p)
}

# maxhgt ###
# year plots ##

{
  for (k in 1:17) {
    p <- veg_height_comp %>%
      filter(site == k) %>%
      ggplot(aes(x = year, y = maxhgt)) +
      geom_boxplot() +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(x = "Year", y = "Maximum vegetation height") +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(alt_labels[k])
    ggsave(
      paste0("./veg_plots/height/maxhgt/years/", k, "_veg_height_years.png"),
      plot = p,
      width = 3840,
      height = 2160,
      units = "px",
      bg = "white",
      create.dir = TRUE
    )
  }

  p <- veg_height_comp %>%
    ggplot(aes(x = year, y = maxhgt)) +
    geom_boxplot() +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Maximum vegetation height") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("All Sites")
  ggsave(
    paste0("./veg_plots/height/maxhgt/years/all_veg_height_years.png"),
    plot = p,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(k, p)
}

# period plots ##

{
  for (k in 1:17) {
    p <- veg_height_comp %>%
      filter(site == k) %>%
      ggplot(aes(x = period, y = maxhgt)) +
      geom_boxplot() +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(x = "Sampling Period", y = "Maximum vegetation height") +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(alt_labels[k]) +
      scale_x_discrete(labels = c("2002/2003", "2022/2023"))
    ggsave(
      paste0("./veg_plots/height/maxhgt/period/", k, "_veg_height_period.png"),
      plot = p,
      width = 3840,
      height = 2160,
      units = "px",
      bg = "white",
      create.dir = TRUE
    )
  }

  p <- veg_height_comp %>%
    ggplot(aes(x = period, y = maxhgt)) +
    geom_boxplot() +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Sampling Period", y = "Maximum vegetation height") +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("All Sites") +
    scale_x_discrete(labels = c("2002/2003", "2022/2023"))
  ggsave(
    paste0("./veg_plots/height/maxhgt/period/all_veg_height_period.png"),
    plot = p,
    width = 3840,
    height = 2160,
    units = "px",
    bg = "white",
    create.dir = TRUE
  )
  rm(k, p)
}

model_tothits_years <- veg_height_comp %>%
  MCMCglmm(
    data = .,
    fixed = tothits ~ (year),
    random = ~site,
  )
summary(model_tothits_years)

model_tothits_years_sites <- veg_height_comp %>%
  MCMCglmm(
    data = .,
    fixed = tothits ~ (year * site),
    random = ~replicate,
  )
summary(model_tothits_years_sites)

model_tothits_period <- veg_height_comp %>%
  MCMCglmm(
    data = .,
    fixed = tothits ~ (period),
    random = ~site,
  )
summary(model_tothits_period)

model_tothits_period_sites <- veg_height_comp %>%
  MCMCglmm(
    data = .,
    fixed = tothits ~ (period * site),
    random = ~replicate,
  )
summary(model_tothits_period_sites)

model_maxhgt_years <- veg_height_comp %>%
  MCMCglmm(
    data = .,
    fixed = maxhgt ~ (year),
    random = ~site,
  )
summary(model_maxhgt_years)

model_maxhgt_years_sites <- veg_height_comp %>%
  MCMCglmm(
    data = .,
    fixed = maxhgt ~ (year * site),
    random = ~replicate,
  )
summary(model_maxhgt_years_sites)

model_maxhgt_period <- veg_height_comp %>%
  MCMCglmm(
    data = .,
    fixed = maxhgt ~ (period),
    random = ~site,
  )
summary(model_maxhgt_period)

model_maxhgt_period_sites <- veg_height_comp %>%
  MCMCglmm(
    data = .,
    fixed = maxhgt ~ (period * site),
    random = ~replicate,
  )
summary(model_maxhgt_period_sites)
