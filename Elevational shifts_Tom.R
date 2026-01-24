#### Import data
# comm <- read.csv("~/Documents/Analysis/Similarities/sampling event and seasons_elevation.csv", sep = ";") # Busi, remember to change this for your file path

comm <- abund_all %>% mutate(veg_type = factor(case_when(site %in% c(1) ~ "strandveld", site %in% c(2, 6, 16) ~ "restioid", site %in% c(3, 4, 5) ~ "proteoid", site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous", site %in% c(11) ~ "alpine", site %in% c(17) ~ "succulent_karoo"))) %>% relocate(site, replicate, year, veg_type) %>% filter(paste(year, season) != "2023 March") %>% group_by(year, season, period, veg_type, site) %>% summarise(across(where(is.numeric), sum), .groups = "drop") %>% relocate(year, season, period, veg_type, site) %>% arrange(site, year)

##

comm[1:6, 1:6]
comm$aspect <- sapply(strsplit(comm$Elevation, " "), "[[", 4)
comm$aspect <- gsub("(Peak)", "W", comm$aspect)
comm$aspect <- gsub("\\(", "", comm$aspect)
comm$aspect <- gsub("\\)", "", comm$aspect)
comm$Elevation <- as.numeric(sapply(strsplit(comm$Elevation, " "), "[[", 1))

siteinfo <- comm[, c(1:5, ncol(comm))]
comm <- comm[, !colnames(comm) %in% colnames(siteinfo)]

#### Calculate elevation range in each year using a loop

# Loop through each species and calculate range shifts. 
# Do this separately for each aspect - they have different thermal regimes. 
# For species that are not present in both time periods - ignore or provide NAs. 


i <- "Acropyga.arnoldi"
i <- "Anoplolepis.steingroeveri"

spp.elev <- NULL
for(i in colnames(comm)){
  
  # Get the right data
  dat <- siteinfo
  dat$abund <- comm[, i]
  dat$abund <- ifelse(dat$abund > 1, 1, dat$abund)
  
  # Remove records with no abundance
  dat <- dat[dat$abund > 0, ]
  
  # Is the species sampled in both time periods?
  sampling <- length(unique(dat$Sampling_event))
  
  # Calculate range shifts
  if(sampling == 1){
    tmp <- data.frame(species = i, 
                      max_E_02 = NA,
                      max_W_02 = NA, 
                      min_E_02 = NA, 
                      min_W_02 = NA,
                      max_E_22 = NA,
                      max_W_22 = NA, 
                      min_E_22 = NA, 
                      min_W_22 = NA)
  } else{
      dat02 <- dat[dat$Sampling_event == "2002-2003", ]
      dat02 <- dat02[dat02$abund > 0, ]
      maxes <- tapply(dat02$Elevation, dat02$aspect, max)
      mins <- tapply(dat02$Elevation, dat02$aspect, min)
      
      # Fix missing aspects
      asp.names <- names(maxes)
      if(length(asp.names) == 1){
        
        if(names(maxes) == "W"){
          
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

      # For 2022...
      dat22 <- dat[dat$Sampling_event == "2022-2023", ]
      dat22 <- dat22[dat22$abund > 0, ]
      maxes <- tapply(dat22$Elevation, dat22$aspect, max)
      mins <- tapply(dat22$Elevation, dat22$aspect, min)
      
      # Fix missing aspects
      asp.names <- names(maxes)
      if(length(asp.names) == 1){
        
        if(names(maxes) == "W"){
          
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

      # Tidy up
      tmp <- data.frame(species = i, cbind(tmp02, tmp22))
      
    
  }
 
    spp.elev <- rbind(spp.elev, tmp)
    print(i)

}

rownames(spp.elev) <- 1:nrow(spp.elev)
spp.elev

library(dplyr)
library(tidyr)
library(ggplot2)

### 1. Compute shifts ----
spp.elev$shift_max_E <- spp.elev$max_E_22 - spp.elev$max_E_02
spp.elev$shift_max_W <- spp.elev$max_W_22 - spp.elev$max_W_02
spp.elev$shift_min_E <- spp.elev$min_E_22 - spp.elev$min_E_02
spp.elev$shift_min_W <- spp.elev$min_W_22 - spp.elev$min_W_02

### 2. Reshape to long format for ONE plot ----
spp_long <- spp.elev %>%
  select(species,
         shift_max_E, shift_max_W,
         shift_min_E, shift_min_W) %>%
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

### 3. ONE combined plot ----
ggplot(spp_long,
       aes(x = reorder(species, shift),
           y = shift,
           fill = aspect)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_wrap(~ type, ncol = 1, scales = "free_x") +
  labs(
    title = "Ant elevational range shifts (2002–2023)",
    x = "Species",
    y = "Shift in elevation (m)",
    fill = "Aspect"
  ) +
  theme_classic(base_size = 12)







#### Use this code 
#### Import data ---------------------------------------------------------
comm <- read.csv("~/Documents/Analysis/Similarities/sampling event and seasons_elevation.csv",
                 sep = ";") # Busi, remember to change this for your file path
comm[1:6, 1:6]

comm$aspect <- sapply(strsplit(comm$Elevation, " "), "[[", 4)
comm$aspect <- gsub("(Peak)", "W", comm$aspect)
comm$aspect <- gsub("\\(", "", comm$aspect)
comm$aspect <- gsub("\\)", "", comm$aspect)
comm$Elevation <- as.numeric(sapply(strsplit(comm$Elevation, " "), "[[", 1))

siteinfo <- comm[, c(1:5, ncol(comm))]
comm <- comm[, !colnames(comm) %in% colnames(siteinfo)]

#### Calculate elevation range in each year using a loop ----------------

spp.elev <- NULL
for(i in colnames(comm)){
  
  dat <- siteinfo
  dat$abund <- comm[, i]
  dat$abund <- ifelse(dat$abund > 1, 1, dat$abund)
  dat <- dat[dat$abund > 0, ]
  
  sampling <- length(unique(dat$Sampling_event))
  
  if(sampling == 1){
    tmp <- data.frame(species = i, 
                      max_E_02 = NA,
                      max_W_02 = NA, 
                      min_E_02 = NA, 
                      min_W_02 = NA,
                      max_E_22 = NA,
                      max_W_22 = NA, 
                      min_E_22 = NA, 
                      min_W_22 = NA)
  } else {
    ## 2002–2003
    dat02 <- dat[dat$Sampling_event == "2002-2003", ]
    dat02 <- dat02[dat02$abund > 0, ]
    maxes <- tapply(dat02$Elevation, dat02$aspect, max)
    mins  <- tapply(dat02$Elevation, dat02$aspect, min)
    
    asp.names <- names(maxes)
    if(length(asp.names) == 1){
      if(names(maxes) == "W"){
        maxes <- c(NA, maxes); mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA); mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_02", "max_W_02")
    names(mins)  <- c("min_E_02", "min_W_02")
    tmp02 <- t(data.frame(c(maxes, mins)))
    
    ## 2022–2023
    dat22 <- dat[dat$Sampling_event == "2022-2023", ]
    dat22 <- dat22[dat22$abund > 0, ]
    maxes <- tapply(dat22$Elevation, dat22$aspect, max)
    mins  <- tapply(dat22$Elevation, dat22$aspect, min)
    
    asp.names <- names(maxes)
    if(length(asp.names) == 1){
      if(names(maxes) == "W"){
        maxes <- c(NA, maxes); mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA); mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_22", "max_W_22")
    names(mins)  <- c("min_E_22", "min_W_22")
    tmp22 <- t(data.frame(c(maxes, mins)))
    
    tmp <- data.frame(species = i, cbind(tmp02, tmp22))
  }
  
  spp.elev <- rbind(spp.elev, tmp)
  print(i)
}

rownames(spp.elev) <- 1:nrow(spp.elev)
spp.elev

#### Shifts + summaries + plot start here --------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

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

 View(summary_table)
 write.csv(summary_table, "elevational_shift_summary.csv", row.names = FALSE)

### 3. PROPORTIONS upslope / downslope / no change ----
count_dir <- function(x){
  c(
    Upslope   = sum(x > 0, na.rm = TRUE),
    Downslope = sum(x < 0, na.rm = TRUE),
    No_change = sum(x == 0, na.rm = TRUE)
  )
}

shift_summary <- data.frame(
  Upper_East  = count_dir(spp.elev$shift_max_E),
  Upper_West  = count_dir(spp.elev$shift_max_W),
  Lower_East  = count_dir(spp.elev$shift_min_E),
  Lower_West  = count_dir(spp.elev$shift_min_W)
)

shift_summary
round(prop.table(as.matrix(shift_summary), 2) * 100, 1)  # % per column

### 4. Reshape to long format for ONE plot ----
spp_long <- spp.elev %>%
  select(species,
         shift_max_E, shift_max_W,
         shift_min_E, shift_min_W) %>%
  pivot_longer(
    cols = starts_with("shift_"),
    names_to = c("type", "aspect"),
    names_pattern = "shift_(max|min)_(E|W)",
    values_to = "shift"
  ) %>%
  mutate(
    type   = ifelse(type == "max", "Upper limit", "Lower limit"),
    aspect = ifelse(aspect == "E", "East", "West")
  )

### 5. ONE combined plot ----
ggplot(spp_long,
       aes(x = reorder(species, shift),
           y = shift,
           fill = aspect)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_wrap(~ type, ncol = 1, scales = "free_x") +
  labs(
    title = "Ant elevational range shifts (2002–2023)",
    x = "Species",
    y = "Shift in elevation (m)",
    fill = "Aspect"
  ) +
  theme_classic(base_size = 12)






#### Use this code for slopes in separate plots and no shift species removed
#### Import data ---------------------------------------------------------
comm <- read.csv("~/Documents/Analysis/Similarities/sampling event and seasons_elevation.csv",
                 sep = ";") # Busi, remember to change this for your file path
comm[1:6, 1:6]

comm$aspect <- sapply(strsplit(comm$Elevation, " "), "[[", 4)
comm$aspect <- gsub("(Peak)", "W", comm$aspect)
comm$aspect <- gsub("\\(", "", comm$aspect)
comm$aspect <- gsub("\\)", "", comm$aspect)
comm$Elevation <- as.numeric(sapply(strsplit(comm$Elevation, " "), "[[", 1))

siteinfo <- comm[, c(1:5, ncol(comm))]
comm <- comm[, !colnames(comm) %in% colnames(siteinfo)]

#### Calculate elevation range in each year using a loop ----------------

spp.elev <- NULL
for(i in colnames(comm)){
  
  dat <- siteinfo
  dat$abund <- comm[, i]
  dat$abund <- ifelse(dat$abund > 1, 1, dat$abund)
  dat <- dat[dat$abund > 0, ]
  
  sampling <- length(unique(dat$Sampling_event))
  
  if(sampling == 1){
    tmp <- data.frame(species = i, 
                      max_E_02 = NA,
                      max_W_02 = NA, 
                      min_E_02 = NA, 
                      min_W_02 = NA,
                      max_E_22 = NA,
                      max_W_22 = NA, 
                      min_E_22 = NA, 
                      min_W_22 = NA)
  } else {
    ## 2002–2003
    dat02 <- dat[dat$Sampling_event == "2002-2003", ]
    dat02 <- dat02[dat02$abund > 0, ]
    maxes <- tapply(dat02$Elevation, dat02$aspect, max)
    mins  <- tapply(dat02$Elevation, dat02$aspect, min)
    
    asp.names <- names(maxes)
    if(length(asp.names) == 1){
      if(names(maxes) == "W"){
        maxes <- c(NA, maxes); mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA); mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_02", "max_W_02")
    names(mins)  <- c("min_E_02", "min_W_02")
    tmp02 <- t(data.frame(c(maxes, mins)))
    
    ## 2022–2023
    dat22 <- dat[dat$Sampling_event == "2022-2023", ]
    dat22 <- dat22[dat22$abund > 0, ]
    maxes <- tapply(dat22$Elevation, dat22$aspect, max)
    mins  <- tapply(dat22$Elevation, dat22$aspect, min)
    
    asp.names <- names(maxes)
    if(length(asp.names) == 1){
      if(names(maxes) == "W"){
        maxes <- c(NA, maxes); mins <- c(NA, mins)
      } else {
        maxes <- c(maxes, NA); mins <- c(mins, NA)
      }
    }
    names(maxes) <- c("max_E_22", "max_W_22")
    names(mins)  <- c("min_E_22", "min_W_22")
    tmp22 <- t(data.frame(c(maxes, mins)))
    
    tmp <- data.frame(species = i, cbind(tmp02, tmp22))
  }
  
  spp.elev <- rbind(spp.elev, tmp)
  print(i)
}

rownames(spp.elev) <- 1:nrow(spp.elev)
spp.elev

#### Shifts + summaries + plot start here --------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

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

write.csv(not_plotted_east, "~/Documents/Analysis/not_plotted_east.csv", row.names = FALSE)
write.csv(not_plotted_west, "~/Documents/Analysis/not_plotted_west.csv", row.names = FALSE)
write.csv(not_plotted_both, "~/Documents/Analysis/not_plotted_both_slopes.csv", row.names = FALSE)

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
write.csv(summary_table, "elevational_shift_summary.csv", row.names = FALSE)

### 3. PROPORTIONS upslope / downslope / no change ----
count_dir <- function(x){
  c(
    Upslope   = sum(x > 0, na.rm = TRUE),
    Downslope = sum(x < 0, na.rm = TRUE),
    No_change = sum(x == 0, na.rm = TRUE)
  )
}

shift_summary <- data.frame(
  Upper_East  = count_dir(spp.elev$shift_max_E),
  Upper_West  = count_dir(spp.elev$shift_max_W),
  Lower_East  = count_dir(spp.elev$shift_min_E),
  Lower_West  = count_dir(spp.elev$shift_min_W)
)

shift_summary
round(prop.table(as.matrix(shift_summary), 2) * 100, 1)  # % per column

### 4. Reshape to long format ----
spp_long <- spp.elev %>%
  select(species,
         shift_max_E, shift_max_W,
         shift_min_E, shift_min_W) %>%
  pivot_longer(
    cols = starts_with("shift_"),
    names_to = c("type", "aspect"),
    names_pattern = "shift_(max|min)_(E|W)",
    values_to = "shift"
  ) %>%
  mutate(
    type   = ifelse(type == "max", "Upper limit", "Lower limit"),
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

write.csv(data.frame(species = no_shift_east),
          "~/Documents/Analysis/no_shift_east_species.csv",
          row.names = FALSE)

write.csv(data.frame(species = no_shift_west),
          "~/Documents/Analysis/no_shift_west_species.csv",
          row.names = FALSE)

write.csv(no_shift_both,
          "~/Documents/Analysis/species_no_shift_both_slopes.csv",
          row.names = FALSE)

### 5. Separate plots for East and West ----

## East slope only
p_east <- ggplot(
  spp_long_filt %>% filter(aspect == "East"),
  aes(x = reorder(species, shift),
      y = shift,
      fill = type)   # colour by upper vs lower, or keep "aspect" if you prefer
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_wrap(~ type, ncol = 1, scales = "free_x") +
  labs(
    title = "Eastern slope (2002-2003 and 2022–2023)",
    x = "Species",
    y = "Shift in elevation (m)",
    fill = "Limit"
  ) +
  theme_classic(base_size = 12)

p_east   # print East plot


## West slope only
p_west <- ggplot(
  spp_long_filt %>% filter(aspect == "West"),
  aes(x = reorder(species, shift),
      y = shift,
      fill = type)
) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "red") +
  coord_flip() +
  facet_wrap(~ type, ncol = 1, scales = "free_x") +
  labs(
    title = "West slope (2002-2003 and 2022–2023)",
    x = "Species",
    y = "Shift in elevation (m)",
    fill = "Limit"
  ) +
  theme_classic(base_size = 12)

p_west   # print West plot
ggsave("~/Documents/Analysis/West_slope_shift.png", p_west, width = 7, height = 12, dpi = 300)
