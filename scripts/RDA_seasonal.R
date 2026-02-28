###################RDA analysis###################
# ==========================================
# Load libraries
# ==========================================
library(vegan)
library(grDevices)

# setwd("~/Documents/Analysis/PhD RDA data")

# ==========================================
# Import data
# ==========================================
env_data_all <- env_season
spp_data_all <- seasonal_data

# ==========================================
# Clean environmental data
# ==========================================
env_data <- env_data_all

# # Keep categorical columns as factors
# env_data$Season <- as.factor(env_data$Season)
# env_data$Sampling_event <- as.factor(env_data$Sampling_event)
# env_data$Slope <- as.factor(env_data$Slope)
# env_data$Site <- as.factor(env_data$Site)
# env_data$Replicate <- as.factor(env_data$Replicate)

# Convert all other columns to numeric (including Fire)
# numeric_cols <- setdiff(
#   names(env_data),
#   c("Season", "Vegetation_type", "Sampling_event", "Site", "Replicate")
# )
# env_data[numeric_cols] <- lapply(env_data[numeric_cols], function(x) {
#   as.numeric(gsub(",", ".", x))
# })

# Check structure to confirm
# str(env_data)

# ==========================================
# Prepare species data
# ==========================================
spp_data <- spp_data_all[, -c(1:5)] # remove non-species columns if needed

# ==========================================
# Standardize (z-score transform)
# ==========================================
Zspp_data <- as.data.frame(scale(spp_data, center = TRUE, scale = TRUE))
Zenv_data <- as.data.frame(scale(
  env_data[, sapply(env_data, is.numeric)],
  center = TRUE,
  scale = TRUE
))

# Add site labels for plotting
Zenv_data$Name <- spp_data_all$site
Zenv_data <- Zenv_data %>%
  mutate(
    slope = as.factor(case_when(
      slope < 0 ~ 1,
      slope > 0 ~ 2,
    )),
    veg_type = factor(case_when(
      Name %in% c(1) ~ "strandveld",
      Name %in% c(2, 6, 16) ~ "restioid",
      Name %in% c(3, 4, 5) ~ "proteoid",
      Name %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous",
      Name %in% c(11) ~ "alpine",
      Name %in% c(17) ~ "succulent_karoo"
    )),
    # period = factor(rep(c("old", "modern"), each = 136)),
    # year = factor(rep(c("2002", "2003", "2022", "2023"), each = 68))
  ) %>%
  select(-AMin_precip)

# Coskun trials ####

# cef <- envfit(tbRDA3, Zenv_data)
# cef #Keep significant or very important variables
model_full <- rda(Zspp_data ~ ., data = Zenv_data)

# Forward selection of variables:
fwd.sel <- ordistep(
  rda(Zspp_data ~ 1, data = Zenv_data), # lower model limit (simple!)
  scope = formula(model_full), # upper model limit (the "full" model)
  direction = "forward",
  R2scope = TRUE, # can't surpass the "full" model's R2
  pstep = 1000,
  trace = FALSE
)

fwd.sel$call

summary(fwd.sel)
summary(model_full)

# spe.rda.signif <- rda(
#   formula = as.formula(fwd.sel$call),
#   data = Zenv_data
# )

spe.rda.signif <- rda(formula(fwd.sel$call), data = Zenv_data)

summary(spe.rda.signif)
RsquareAdj(spe.rda.signif)

RsquareAdj(model_full)

anova.cca(spe.rda.signif, step = 1000)
anova.cca(spe.rda.signif, step = 1000, by = "term")
anova.cca(spe.rda.signif, step = 1000, by = "axis")
ordiplot(spe.rda.signif, scaling = 1, type = "text")

plot(
  spe.rda.signif,
  display = "bp",
) #Show sites and species + env. variables
# text(spe.rda.signif, "sites", col = "black", cex = 0.8) #Show sites names
# ordihull(spe.rda.signif, group = spp_data_all$period, label = TRUE, ) #Link sites with polygones representing groups
text(spe.rda.signif, "species", col = "red", cex = 0.8) #Show species names
text(spe.rda.signif, "species", col = "red", cex = 0.8) #Show species names
# points(spe.rda.signif, display = "sites") #Add points for sites
points(spe.rda.signif, display = "species", pch = "+")
x <- scores(spe.rda.signif)
anova(spe.rda.signif, by = "axis", perm.max = 500) #Test significance of axis
anova(spe.rda.signif) #Test significance of RDA


smry <- scores(spe.rda.signif)
df1 <- cbind(data.frame(smry$sites[, 1:2]), Zenv_data)
df2 <- data.frame(smry$species[, 1:2]) # loadings for PC1 and PC2

{
  df3 <- data.frame(smry$regression)
  df3 <- df3[!rownames(df3) %like% "Name", ]
}

rda.plot <- ggplot(df1, aes(x = RDA1, y = RDA2)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted")

rda.plot

rda.biplot <- rda.plot +
  geom_segment(
    data = df3,
    aes(x = 0, xend = RDA1, y = 0, yend = RDA2),
    color = "red",
    arrow = arrow(length = unit(0.02, "npc"))
  ) +
  geom_text(
    data = df3,
    aes(
      x = RDA1,
      y = RDA2,
      label = rownames(df3),
      hjust = 0.5 * (1 - sign(RDA1)),
      vjust = 0.5 * (1 - sign(RDA2))
    ),
    color = "red",
    size = 4
  ) +
  scale_x_continuous(limits = c(-1.1, 2.1)) +
  geom_point(
    data = df1,
    aes(x = RDA1, y = RDA2, shape = veg_type, col = veg_type),
    cex = 3,
    position = position_jitter()
  )

rda.biplot

# glm trials ####

####PLOTS####
#Plot results "tbRDA1" = selected variables + all significant and almost significant variables
#X11()
plot(
  tbRDA1,
  pch = as.character(Zenv_data$Name),
  col = 1:6,
) #Show sites and species + env. variables
# text(tbRDA1, "sites", col="black", cex=0.8) #Show sites names
ordihull(tbRDA1, group = spp_data_all$period, label = TRUE, ) #Link sites with polygones representing groups
text(tbRDA1, "species", col = "red", cex = 0.8) #Show species names
# points(tbRDA1, display = "sites") #Add points for sites
# points(tbRDA1, display = "species", pch = "+")
anova(tbRDA1, by = "axis", perm.max = 500) #Test significance of axis
anova(tbRDA1) #Test significance of RDA

#Plot results "tbRDA2" = all significant and almost significant variables
plot(
  tbRDA2,
  pch = as.character(Zenv_data$Name),
  col = "black",
) # affichage de tes sites et de tes especes (pas les noms, juste leurs points) + variables enviro
ordihull(tbRDA2, group = spp_data_all$period, label = TRUE) #Link sites with polygones
text(tbRDA2, "species", col = "red", cex = 0.8) #Show species names
#points(tbRDA2, display = "sites") #Add points for sites
# points(tbRDA2, display = "species", pch = "+")
anova(tbRDA2, by = "axis", perm.max = 500) #Test significance of axis
anova(tbRDA2) #Test significance of RDA

#Plot results "tbRDA3" = only significant variables
plot(
  tbRDA3,
  pch = as.character(Zenv_data$Name),
  col = "black",
) # affichage de tes sites et de tes especes (pas les noms, juste leurs points) + variables enviro
ordihull(tbRDA3, group = spp_data_all$period, label = TRUE) #Link sites with polygones
text(tbRDA3, "species", col = "red", cex = 0.8) #Show species names
#points(tbRDA3, display = "sites") #Add points for sites
#points(tbRDA3, display = "species", pch = "+")
anova(tbRDA3, by = "axis", perm.max = 500) #Test significance of axis
anova(tbRDA3) #Test significance of RDA


##Last test - all significant variables + sites as a condition - tbRDA4
tbRDA4 = rda(
  Zspp_data ~ bare_ground +
    litter +
    rock +
    veg +
    rock_soil +
    clay +
    sand +
    silt +
    `NO3(mg/kg)` +
    `C(%)` +
    maxhgt +
    m_tground +
    mMax_tground +
    mMin_tground +
    AMax_tground +
    AMin_tground +
    range_tground +
    m_tair +
    mMax_tair +
    mMin_tair +
    AMin_tair +
    m_precip +
    mMax_precip +
    mMin_precip +
    AMax_precip +
    range_precip +
    elevation +
    slope +
    Name,
  data = Zenv_data
)
summary(tbRDA4)

plot(
  tbRDA4,
  pch = as.character(Zenv_data$Name),
  col = "black",
)
ordihull(tbRDA4, group = spp_data_all$period, label = TRUE)
anova(tbRDA4, by = "axis", perm.max = 500)
anova(tbRDA4)


##############################################
# Indicator Species Analysis by Vegetation Type
# (2002–2003 vs 2022–2023)
##############################################

# =====================================================
# Load required packages
# =====================================================
library(indicspecies)
library(vegan)
library(dplyr)
library(ggplot2)

# =====================================================
# Confirm matching rows between environment and species data
# =====================================================
stopifnot(nrow(spp_data_all) == nrow(env_data_all))

# =====================================================
# Define helper function to run IndVal by vegetation type per sampling event
# =====================================================
run_indval_by_veg <- function(Sampling_event_label) {
  message("Running IndVal for: ", Sampling_event_label)

  # Subset environmental and species data
  env_sub <- env_data_all %>% filter(year == Sampling_event_label)
  spp_sub <- spp_data_all[env_data_all$year == Sampling_event_label, ]

  # Extract species abundance matrix
  spp_matrix <- spp_sub[, 6:ncol(spp_sub)]

  # Define vegetation grouping
  group <- as.factor(env_sub$veg_type)

  message("Unique vegetation types: ", length(unique(group)))

  if (length(unique(group)) < 2) {
    message(
      "⚠️ Skipping ",
      Sampling_event_label,
      ": only one vegetation type found."
    )
    return(NULL)
  }

  # =====================================================
  # Run IndVal
  # =====================================================
  indval_res <- multipatt(
    spp_matrix,
    cluster = group,
    func = "IndVal.g",
    duleg = TRUE,
    control = how(nperm = 999)
  )

  # Extract and filter significant species
  sig_tab <- data.frame(indval_res$sign)
  sig_tab <- sig_tab %>%
    mutate(Species = rownames(indval_res$sign)) %>%
    filter(p.value <= 0.05) %>%
    arrange(desc(stat))

  sig_tab$Year <- Sampling_event_label

  # Save individual file
  out_name <- paste0(
    "./indval/IndVal_",
    gsub("-", "_", Sampling_event_label),
    "_VegetationType.csv"
  )
  write.csv(sig_tab, out_name, row.names = FALSE)
  message("✅ Saved: ", out_name)

  return(sig_tab)
}

# =====================================================
# Run IndVal for both sampling events
# =====================================================
indval_2002_2003 <- run_indval_by_veg("2022")
indval_2022_2023 <- run_indval_by_veg("2023")

# =====================================================
# Combine all results
# =====================================================
indval_all <- bind_rows(indval_2002_2003, indval_2022_2023)
write.csv(
  indval_all,
  "./indval/IndVal_VegetationType_Combined.csv",
  row.names = FALSE
)
message("✅ Combined file saved: IndVal_VegetationType_Combined.csv")

# =====================================================
# Extract dominant vegetation type per species for each year
# =====================================================
extract_best_group <- function(df) {
  s_cols <- grep("^s\\.", colnames(df), value = TRUE)
  df$BestGroup <- apply(df[, s_cols], 1, function(x) s_cols[which.max(x)])
  df$BestGroup <- gsub("^s\\.", "", df$BestGroup)

  df_out <- df %>%
    select(Species, BestGroup, stat, Year) %>%
    rename(VegType = BestGroup, IndVal = stat)

  return(df_out)
}

# Apply to both years
ind_2002_2003_clean <- extract_best_group(indval_2002_2003)
ind_2022_2023_clean <- extract_best_group(indval_2022_2023)

# =====================================================
# Merge and classify indicator species
# =====================================================
ind_2002_2003 <- ind_2002_2003_clean %>%
  select(Species, VegType, IndVal) %>%
  rename(Veg_2002_2003 = VegType, Stat_2002_2003 = IndVal)

ind_2022_2023 <- ind_2022_2023_clean %>%
  select(Species, VegType, IndVal) %>%
  rename(Veg_2022_2023 = VegType, Stat_2022_2023 = IndVal)

compare_tbl <- full_join(ind_2002_2003, ind_2022_2023, by = "Species")

compare_tbl <- compare_tbl %>%
  mutate(
    Status = case_when(
      !is.na(Veg_2002_2003) &
        !is.na(Veg_2022_2023) &
        Veg_2002_2003 == Veg_2022_2023 ~ "Shared indicator",
      !is.na(Veg_2002_2003) &
        is.na(Veg_2022_2023) ~ "Lost indicator (2022 only)",
      is.na(Veg_2002_2003) &
        !is.na(Veg_2022_2023) ~ "New indicator (2023 only)",
      !is.na(Veg_2002_2003) &
        !is.na(Veg_2022_2023) &
        Veg_2002_2003 !=
          Veg_2022_2023 ~ "Shifted indicator (changed vegetation type)",
      TRUE ~ "Non-indicator"
    )
  )

write.csv(
  compare_tbl,
  "./indval/IndVal_VegetationType_Comparison_2002_2003_vs_2022_2023.csv",
  row.names = FALSE
)
message(
  "✅ Comparison table saved: IndVal_VegetationType_Comparison_2002_2003_vs_2022_2023.csv"
)

# =====================================================
# Summary + Visualization
# =====================================================
# =====================================================
# Improved plot showing both old and new indicator positions
# =====================================================

# Create a unified 'Vegetation Type' column for plotting
compare_tbl_plot <- compare_tbl %>%
  mutate(
    VegType_plot = case_when(
      Status == "New indicator (2023 only)" ~ Veg_2022_2023,
      TRUE ~ Veg_2002_2003
    )
  )

# Remove NAs (optional)
compare_tbl_plot <- compare_tbl_plot %>% filter(!is.na(VegType_plot))

# Plot
ggplot(
  compare_tbl_plot %>% filter(Status != "Non-indicator"),
  aes(x = VegType_plot, fill = Status)
) +
  geom_bar(position = "dodge") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Vegetation Type",
    y = "Number of indicator species"
  ) +
  scale_fill_manual(
    values = c(
      "Shared indicator" = "forestgreen",
      "New indicator (2022–2023 only)" = "steelblue3",
      "Lost indicator (2002–2003 only)" = "darkorange2",
      "Shifted indicator (changed vegetation type)" = "purple3"
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


############################################## IndVal Elevation
##############################################
# Indicator Species Analysis by Elevation × Slope
# (Treats 1700_1 ≠ 1700_2 as distinct ecological sites)
##############################################
# =====================================================
# Load required packages
# =====================================================
library(indicspecies)
library(vegan)
library(dplyr)
library(ggplot2)

# =====================================================
# Confirm matching rows between environment and species data
# =====================================================
stopifnot(nrow(spp_data_all) == nrow(env_data_all))

# =====================================================
# ADDED: Convert numeric slope codes to W/E labels
# =====================================================
env_data_all <- env_data_all %>%
  mutate(Slope = ifelse(slope == 1, "W", ifelse(slope == 2, "E", slope))) %>%
  select(-slope) %>%
  mutate(elevation = altitudes[site])

# If Slope is a factor, ensure it stays that way
env_data_all$Slope <- factor(env_data_all$Slope, levels = c("W", "E"))

# =====================================================
# Define helper function to run IndVal per Sampling Event
# =====================================================
run_indval_by_elevslope <- function(Sampling_event_label) {
  message("Running IndVal for: ", Sampling_event_label)

  # Subset environmental and species data
  env_sub <- env_data_all %>% filter(year == Sampling_event_label)
  spp_sub <- spp_data_all[env_data_all$year == Sampling_event_label, ]

  # Extract species abundance matrix
  spp_matrix <- spp_sub[, 6:ncol(spp_sub)]

  # =====================================================
  # Combine Elevation + Slope as grouping factor
  # =====================================================
  env_sub$ElevSlope <- paste0(env_sub$elevation, "_", env_sub$Slope)
  group <- as.factor(env_sub$ElevSlope)

  message("Unique Elevation×Slope groups: ", length(unique(group)))

  # Skip if only one group (edge case)
  if (length(unique(group)) < 2) {
    message(
      "⚠️ Skipping ",
      Sampling_event_label,
      ": only one elevation-slope level found."
    )
    return(NULL)
  }

  # =====================================================
  # Run IndVal
  # =====================================================
  indval_res <- multipatt(
    spp_matrix,
    cluster = group,
    func = "IndVal.g",
    duleg = TRUE,
    control = how(nperm = 999)
  )

  # Extract and filter significant species
  sig_tab <- data.frame(indval_res$sign)
  sig_tab <- sig_tab %>%
    mutate(Species = rownames(indval_res$sign)) %>%
    filter(p.value <= 0.05) %>%
    arrange(desc(stat))

  sig_tab$Year <- Sampling_event_label

  # Save results per sampling event
  out_name <- paste0(
    "IndVal_",
    gsub("-", "_", Sampling_event_label),
    "_ElevSlope.csv"
  )
  write.csv(sig_tab, out_name, row.names = FALSE)
  message("✅ Saved: ", out_name)

  return(sig_tab)
}

# =====================================================
# Run for both sampling events
# =====================================================
indval_2002_2003 <- run_indval_by_elevslope("2022")
indval_2022_2023 <- run_indval_by_elevslope("2023")

# Combine both
indval_all <- bind_rows(indval_2002_2003, indval_2022_2023)
write.csv(indval_all, "IndVal_ElevSlope_Combined.csv", row.names = FALSE)
message("✅ Combined file saved: IndVal_ElevSlope_Combined.csv")

# =====================================================
# Identify dominant ElevSlope per species per year
# =====================================================
extract_best_elevslope <- function(df) {
  s_cols <- grep("^s\\.", colnames(df), value = TRUE)

  df$BestGroup <- apply(df[, s_cols], 1, function(x) {
    s_cols[which.max(x)]
  })
  df$BestGroup <- gsub("^s\\.", "", df$BestGroup)

  df_out <- df %>%
    select(Species, BestGroup, stat, Year) %>%
    rename(ElevSlope = BestGroup, IndVal = stat)

  return(df_out)
}

# Apply to both
ind_2002_2003_clean <- extract_best_elevslope(indval_2002_2003)
ind_2022_2023_clean <- extract_best_elevslope(indval_2022_2023)

# =====================================================
# Merge and classify indicator species
# =====================================================
ind_2002_2003 <- ind_2002_2003_clean %>%
  select(Species, ElevSlope, IndVal) %>%
  rename(ElevSlope_2002_2003 = ElevSlope, Stat_2002_2003 = IndVal)

ind_2022_2023 <- ind_2022_2023_clean %>%
  select(Species, ElevSlope, IndVal) %>%
  rename(ElevSlope_2022_2023 = ElevSlope, Stat_2022_2023 = IndVal)

compare_tbl <- full_join(ind_2002_2003, ind_2022_2023, by = "Species")

compare_tbl <- compare_tbl %>%
  mutate(
    Status = case_when(
      !is.na(ElevSlope_2002_2003) &
        !is.na(ElevSlope_2022_2023) &
        ElevSlope_2002_2003 == ElevSlope_2022_2023 ~ "Shared indicator",
      !is.na(ElevSlope_2002_2003) &
        is.na(ElevSlope_2022_2023) ~ "Lost indicator (2022 only)",
      is.na(ElevSlope_2002_2003) &
        !is.na(ElevSlope_2022_2023) ~ "New indicator (2023 only)",
      !is.na(ElevSlope_2002_2003) &
        !is.na(ElevSlope_2022_2023) &
        ElevSlope_2002_2003 !=
          ElevSlope_2022_2023 ~ "Shifted indicator (changed elevation)",
      TRUE ~ "Non-indicator"
    )
  )

write.csv(
  compare_tbl,
  "IndVal_ElevSlope_Comparison_2002_2003_vs_2022_2023.csv",
  row.names = FALSE
)
message(
  "✅ Comparison table saved: IndVal_ElevSlope_Comparison_2002_2003_vs_2022_2023.csv"
)

# =====================================================
# Summary + Visualization
# =====================================================
# =====================================================
# Improved plot showing both old and new indicator positions
# =====================================================

compare_tbl_plot <- compare_tbl %>%
  mutate(
    ElevSlope_plot = case_when(
      Status == "New indicator (2023 only)" ~ ElevSlope_2022_2023,
      TRUE ~ ElevSlope_2002_2003
    )
  )

# Remove NAs (optional)
compare_tbl_plot <- compare_tbl_plot %>% filter(!is.na(ElevSlope_plot))

# Plot
ggplot(
  compare_tbl_plot %>% filter(Status != "Non-indicator"),
  aes(x = ElevSlope_plot, fill = Status)
) +
  geom_bar(position = "dodge") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Elevation", # ← Changed here
    y = "Number of indicator species"
  ) +
  scale_fill_manual(
    values = c(
      "Shared indicator" = "forestgreen",
      "New indicator (2022–2023 only)" = "steelblue3",
      "Lost indicator (2002–2003 only)" = "darkorange2",
      "Shifted indicator (changed elevation)" = "purple3"
    )
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Plot
ggplot(
  compare_tbl_plot %>% filter(Status != "Non-indicator"),
  aes(x = ElevSlope_plot, fill = Status)
) +
  geom_bar(position = "dodge") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major = element_blank(), # remove major gridlines
    panel.grid.minor = element_blank(), # remove minor gridlines
    panel.border = element_blank(), # keep clean panel (optional)
    axis.line = element_line(colour = "black"), # keep axis lines
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    x = "Elevation",
    y = "Number of indicator species"
  ) +
  scale_fill_manual(
    values = c(
      "Shared indicator" = "forestgreen",
      "New indicator (2022–2023 only)" = "steelblue3",
      "Lost indicator (2002–2003 only)" = "darkorange2",
      "Shifted indicator (changed elevation)" = "purple3"
    )
  )


################################ Shifts in all species_ VEG TYPES
library(dplyr)
library(tidyr)

# 1️⃣ Make sure metadata matches
stopifnot(nrow(spp_data_all) == nrow(env_data_all))

# 2️⃣ Combine environmental and species data
dat_all <- bind_cols(
  env_data_all[, c("year", "veg_type")],
  spp_data_all[, 6:ncol(spp_data_all)]
) # assuming species start at col 7

# 3️⃣ Reshape to long format
dat_long <- dat_all %>%
  pivot_longer(
    cols = -c(year, veg_type),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  mutate(Presence = ifelse(Abundance > 0, 1, 0))

# 4️⃣ Summarise presence by vegetation type and year
species_presence <- dat_long %>%
  group_by(year, veg_type, Species) %>%
  summarise(Present = as.numeric(sum(Presence) > 0), .groups = "drop")

# 5️⃣ Spread into two columns for comparison
species_wide <- species_presence %>%
  pivot_wider(
    names_from = year,
    values_from = Present,
    values_fill = 0
  ) %>%
  rename(Present_2002_2003 = "2022", Present_2022_2023 = "2023")

# 6️⃣ Identify species that changed vegetation type
species_shift <- species_wide %>%
  group_by(Species) %>%
  summarise(
    VegTypes_2002 = paste(
      veg_type[Present_2002_2003 == 1],
      collapse = ", "
    ),
    VegTypes_2022 = paste(
      veg_type[Present_2022_2023 == 1],
      collapse = ", "
    )
  ) %>%
  mutate(
    Status = case_when(
      VegTypes_2002 == "" & VegTypes_2022 != "" ~ "New species (2023 only)",
      VegTypes_2002 != "" & VegTypes_2022 == "" ~ "Lost species (2022 only)",
      VegTypes_2002 != VegTypes_2022 ~ "Shifted vegetation type",
      VegTypes_2002 == VegTypes_2022 ~ "Stable across vegetation types"
    )
  )

# ✅ View results
table(species_shift$Status)
head(species_shift %>% filter(Status == "Shifted vegetation type"))

# 💾 7️⃣ Save table to CSV
write.csv(
  species_shift,
  "Species_Shifts_All_VegetationTypes_2002_2003_vs_2022_2023.csv",
  row.names = FALSE
)

message(
  "✅ Table saved as: Species_Shifts_All_VegetationTypes_2002_2003_vs_2022_2023.csv"
)


################################ Shifts in all species (Elevation × Slope)
library(dplyr)
library(tidyr)

# 1️⃣ Make sure metadata matches
stopifnot(nrow(spp_data_all) == nrow(env_data_all))
x <- env_data_all #backup
# 2️⃣ Convert numeric slope codes (1/2) to W/E
env_data_all <- env_data_all %>%
  mutate(Slope = ifelse(Slope == 1, "W", ifelse(Slope == 2, "E", Slope)))

# Ensure slope is a factor for consistency
env_data_all$Slope <- factor(env_data_all$Slope, levels = c("W", "E"))

# 3️⃣ Combine Elevation × Slope first — clean and clear
env_data_all <- env_data_all %>%
  mutate(ElevSlope = paste0(elevation, "_", Slope)) %>%
  select(year, ElevSlope) # keep only what's needed for this analysis

# 4️⃣ Combine environmental and species data
dat_all <- bind_cols(env_data_all, spp_data_all[, 6:ncol(spp_data_all)]) # assuming species start at col 7

# 5️⃣ Reshape to long format
dat_long <- dat_all %>%
  pivot_longer(
    cols = -c(year, ElevSlope),
    names_to = "Species",
    values_to = "Abundance"
  ) %>%
  mutate(Presence = ifelse(Abundance > 0, 1, 0))

# 6️⃣ Summarise presence by Elevation × Slope and year
species_presence <- dat_long %>%
  group_by(year, ElevSlope, Species) %>%
  summarise(Present = as.numeric(sum(Presence) > 0), .groups = "drop")

# 7️⃣ Spread into two columns for comparison
species_wide <- species_presence %>%
  pivot_wider(
    names_from = year,
    values_from = Present,
    values_fill = 0
  ) %>%
  rename(Present_2002_2003 = "2022", Present_2022_2023 = "2023")

# 8️⃣ Identify species that changed Elevation × Slope
species_shift <- species_wide %>%
  group_by(Species) %>%
  summarise(
    ElevSlope_2002 = paste(ElevSlope[Present_2002_2003 == 1], collapse = ", "),
    ElevSlope_2022 = paste(ElevSlope[Present_2022_2023 == 1], collapse = ", ")
  ) %>%
  mutate(
    Status = case_when(
      ElevSlope_2002 == "" & ElevSlope_2022 != "" ~ "New species (2023 only)",
      ElevSlope_2002 != "" & ElevSlope_2022 == "" ~ "Lost species (2022 only)",
      ElevSlope_2002 != ElevSlope_2022 ~ "Shifted elevation",
      ElevSlope_2002 == ElevSlope_2022 ~ "Stable across elevations"
    )
  )

# ✅ View results
table(species_shift$Status)
head(species_shift %>% filter(Status == "Shifted elevation"))

# 💾 9️⃣ Save table to CSV
write.csv(
  species_shift,
  "Species_Shifts_All_ElevSlope_2002_2003_vs_2022_2023.csv",
  row.names = FALSE
)

message(
  "✅ Table saved as: Species_Shifts_All_ElevSlope_2002_2003_vs_2022_2023.csv"
)
