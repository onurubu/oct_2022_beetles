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
  relocate(year, season, site, replicate, period) %>%
  mutate(
    veg_type = factor(case_when(
      site %in% c(1) ~ "strandveld",
      site %in% c(2, 6, 16) ~ "restioid",
      site %in% c(3, 4, 5) ~ "proteoid",
      site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
      site %in% c(11) ~ "alpine",
      site %in% c(17) ~ "succulent_karoo"
    ))
  )

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
      ggplot(aes(x = period, y = (value / 8), fill = name)) +
      geom_bar(position = "stack", stat = "identity") +
      theme_minimal(base_size = 15) +
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
    ggplot(aes(x = period, y = (value / 136), fill = name)) +
    geom_bar(position = "stack", stat = "identity") +
    theme_minimal(base_size = 15) +
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

model_veg_period <- data.frame(veg_cover_comp_long) %>%
  MCMCglmm(
    data = .,
    fixed = value ~ (period * name),
    random = ~site,
  )
summary(model_veg_period)
library(easystats)
model_dashboard(model_veg_period)

model_veg_period_type <- veg_cover_comp_long %>%
  MCMCglmm(
    data = .,
    fixed = value ~ (period * name * veg_type),
    random = ~site,
  )
summary(model_veg_period_type)

x <- veg_cover_comp %>%
  group_by(period) %>%
  summarise(across(where(is.numeric), mean))
y <- x %>%
  group_by(period) %>%
  summarise(across(where(is.numeric), mean)) %>%
  summarise(
    bare_ground = (bare_ground[1] - bare_ground[2]) / bare_ground[2] * 100,
    litter = (litter[1] - litter[2]) / litter[2] * 100,
    rock = (rock[1] - rock[2]) / rock[2] * 100,
    veg = (veg[1] - veg[2]) / veg[2] * 100
  )

Sys.sleep(0.1)
x
y
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
      theme_minimal(base_size = 15) +
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
    theme_minimal(base_size = 15) +
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
      theme_minimal(base_size = 15) +
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
    theme_minimal(base_size = 15) +
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
      ggplot(aes(x = year, y = ((maxhgt - 1) * 25))) +
      geom_boxplot() +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(x = "Year", y = "Maximum vegetation height (cm)") +
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
    ggplot(aes(x = year, y = ((maxhgt - 1) * 25))) +
    geom_boxplot() +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Year", y = "Maximum vegetation height (cm)") +
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
      ggplot(aes(x = period, y = ((maxhgt - 1) * 25))) +
      geom_boxplot() +
      theme_minimal(base_size = 15) +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(x = "Sampling Period", y = "Maximum vegetation height (cm)") +
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
    ggplot(aes(x = period, y = ((maxhgt - 1) * 25))) +
    geom_boxplot() +
    theme_minimal(base_size = 15) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "Sampling Period", y = "Maximum vegetation height (cm)") +
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

(x <- veg_height_comp %>%
  group_by(period) %>%
  summarise(across(where(is.numeric), mean)))

(y <- x %>%
  group_by(period) %>%
  summarise(across(where(is.numeric), mean)) %>%
  summarise(
    tothits = (tothits[2] - tothits[1]) / tothits[1] * 100,
    maxhgt = (maxhgt[2] - maxhgt[1]) / maxhgt[1] * 100
  ))


# synthesis ####

{
  {
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
        relocate(site, replicate, year, season, period) %>%
        ungroup()

      soil_dat <- rbind(soil_dat_02, soil_dat_22)
      rm(soil_dat_02)
    }

    veg_height <- data.table::fread(
      "./veg_height_data/veg_height_tothits.csv"
    ) %>%
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

    site_data <- data.table::fread("./site_data/site_attributes.csv") %>%
      select(site, replicate, elevation, slope) %>%
      # select(site, replicate, elevation, slope, longitude, latitude) %>%
      # rename(x = longitude, y = latitude) %>%
      mutate(
        site = as.factor(site),
        replicate = as.factor(replicate)
      ) #,
    #   x2 = (x)^2,
    #   xy = (x * y),
    #   y2 = (y)^2,
    #   x3 = (x)^3,
    #   x2y = ((x)^2 * y),
    #   xy2 = (x * (y)^2),
    #   y3 = (y)^3
    # )

    fire_data <- readRDS("./fire_data/fire_history.rds") %>%
      select(-fire_years)
    fire_data <- fire_data %>%
      mutate(
        d_fire = case_when(d_fire >= 43 ~ d_fire + 11, d_fire < 43 ~ d_fire)
      )

    clim_data <- readRDS("./climate_data/clim_data.rds")
    clim_data_comp <- clim_data %>%
      filter(paste(year, month) != "2023 3") %>%
      select(-month)

    clim_data_ex_comp <- readRDS("./climate_data/extra_clim_data.rds") %>%
      group_by("year" = year(obs_time), "month" = month(obs_time)) %>%
      filter(month == 10, year %in% c(2002, 2003, 2022, 2023)) %>%
      ungroup() %>%
      group_by(year, site, replicate) %>%
      summarise(across(where(is.numeric), mean), .groups = "drop") %>%
      mutate(year = as.factor(year), month = as.factor(month)) %>%
      relocate(site, replicate, year, month)

    altitudes <- c(
      0,
      200,
      300,
      500,
      700,
      900,
      1100,
      1300,
      1500,
      1700,
      1900,
      1700,
      1500,
      1300,
      1100,
      900,
      500
    )

    site_areas <- c(
      3082.58,
      1531.61,
      1531.61,
      765.43,
      447.53,
      370.14,
      282.66,
      206.59,
      106.73,
      21.47,
      3.40,
      6.76,
      39.32,
      76.94,
      232.54,
      607.88,
      1353.37
    )

    {
      dat <- full_join(veg_cover_comp, soil_dat)
      dat <- full_join(dat, veg_height_comp)
      dat <- full_join(dat, fire_data)
      dat <- full_join(dat, clim_data_comp)
      dat <- full_join(dat, clim_data_ex_comp)

      env_all <- left_join(dat, site_data) %>%
        mutate(
          veg_type = factor(case_when(
            site %in% c(1) ~ "strandveld",
            site %in% c(2, 6, 16) ~ "restioid",
            site %in% c(3, 4, 5) ~ "proteoid",
            site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
            site %in% c(11) ~ "alpine",
            site %in% c(17) ~ "succulent_karoo"
          )),
          slope = (case_when(
            site %in% c(1:11) ~ 1,
            site %in% c(12:17) ~ 2,
          )),
          area = site_areas[site]
        ) %>%
        relocate(veg_type, .after = period)
      rm(dat)
    }
  }

  # seasonal ###

  abund_seasonal <- beet_wide %>%
    select(-ID) %>%
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
    arrange(year) %>%
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

  {
    xcover <- veg_cover_season %>%
      group_by(year, site, replicate) %>%
      summarise(across(where(is.numeric), mean), .groups = "drop")

    xfire <- fire_data %>% filter(year %in% c(2022, 2023))

    clim_data_season <- clim_data %>%
      filter(year != 2002, year != 2003, paste(year, month) != "2023 10") %>%
      select(-month)

    clim_data_ex_season <- readRDS("./climate_data/extra_clim_data.rds") %>%
      group_by("year" = year(obs_time), "month" = month(obs_time)) %>%
      filter(paste(year, month) %in% c("2022 10", "2023 3")) %>%
      ungroup() %>%
      group_by(year, site, replicate) %>%
      summarise(across(where(is.numeric), mean), .groups = "drop") %>%
      mutate(year = as.factor(year), month = as.factor(month)) %>%
      relocate(site, replicate, year, month)

    soil_dat_22 <- soil_dat_22 %>% select(-c(season, period))

    dat <- full_join(xcover, veg_height_season)
    dat <- full_join(dat, xfire)
    dat <- full_join(dat, clim_data_season)
    dat <- full_join(dat, clim_data_ex_season)
    dat <- dat %>% ungroup()
    dat <- full_join(dat, soil_dat_22)

    env_season <- left_join(dat, site_data) %>%
      mutate(
        veg_type = factor(case_when(
          site %in% c(1) ~ "strandveld",
          site %in% c(2, 6, 16) ~ "restioid",
          site %in% c(3, 4, 5) ~ "proteoid",
          site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
          site %in% c(11) ~ "alpine",
          site %in% c(17) ~ "succulent_karoo"
        )),
        slope = (case_when(
          site %in% c(1:11) ~ 1,
          site %in% c(12:17) ~ 2,
        )),
        area = site_areas[site]
      ) %>%
      relocate(veg_type, .after = replicate)

    rm(
      dat,
      xcover,
      xfire,
      clim_data_season,
      clim_data_ex_season,
      clim_data_comp,
      clim_data_ex_comp,
      clim_data,
      fire_data,
      soil_dat,
      veg_cover,
      veg_cover_comp_long,
      veg_cover_season,
      veg_height,
      veg_height_season,
      site_data,
      veg_cover_comp,
      veg_height_comp,
      soil_dat_22,
      site_areas
    )
  }
}
