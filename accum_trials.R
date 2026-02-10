############################################
# Sample-based rarefaction curves
############################################

##################for all elevations, years (2002-2005 and 2022-2023) separate#####################################
############################################
# Sample-based accumulation curves
# Using Elevation instead of Site
# Years: 2002–2003 vs 2022–2023
############################################

# Load package
library(vegan)
library(tidyverse)

# 1. Read CSV
beetles <- beet_wide %>%
  select(-ID, -season) %>%
  rename(
    "Year" = "year",
    "Elevation" = "site",
    "Replicate" = "replicate",
    "Pitfall.trap" = "trap"
  ) %>%
  mutate(Elevation = as.factor(alt_labels[Elevation])) %>%
  relocate(Year, Elevation, Replicate, Pitfall.trap)

# Drop empty column if present
beetles <- beetles[, !grepl("^Unnamed", names(beetles))]

# 2. Standardise column names
names(beetles) <- make.names(names(beetles))

# 3. Identify metadata and species columns
meta_cols <- c("Year", "Elevation", "Replicate", "Pitfall.trap")
species_cols <- setdiff(names(beetles), meta_cols)

# 4. Convert species columns to numeric
beetles[, species_cols] <- lapply(beetles[, species_cols], function(x) {
  suppressWarnings(as.numeric(x))
})

# 5. Check the years available
cat("Year present:\n")
print(unique(beetles$Year))

cat("Elevations present:\n")
print(unique(beetles$Elevation))

# 6. Define elevation groups
group1 <- c(alt_labels[1:6])

group2 <- c(alt_labels[7:10])

group3 <- c(alt_labels[12:17])

# 7. Helper: subset by year
subset_by_year <- function(year_tag) {
  subset(beetles, Year == year_tag)
}

# 8. Compute accumulation curve for one elevation × year
get_accum <- function(data_year, elev_label) {
  sub <- data_year[data_year$Elevation == elev_label, ]

  # Need at least 2 samples
  if (nrow(sub) < 2) {
    cat("Skipping:", elev_label, "| too few samples\n")
    return(NULL)
  }

  comm <- sub[, species_cols, drop = FALSE]
  comm <- comm[, colSums(comm, na.rm = TRUE) > 0, drop = FALSE]

  # Need at least 2 species
  if (ncol(comm) < 2) {
    cat("Skipping:", elev_label, "| too few species\n")
    return(NULL)
  }

  cat(
    "Year =",
    unique(sub$Year),
    "| Elevation =",
    elev_label,
    "| Samples =",
    nrow(sub),
    "| Species cols =",
    ncol(comm),
    "\n"
  )

  specaccum(comm, method = "random", permutations = 1000)
}


# 9. Plot one group for a given year
plot_group <- function(data_year, group_elevs, year_label, title_text) {
  acc_list <- lapply(group_elevs, function(e) get_accum(data_year, e))
  keep <- !sapply(acc_list, is.null)
  acc_list <- acc_list[keep]
  good_elevs <- group_elevs[keep]

  if (length(acc_list) == 0) {
    plot(
      1,
      1,
      type = "n",
      main = paste(title_text, year_label),
      xlab = "",
      ylab = "",
      sub = "No data available"
    )
    return()
  }

  xmax <- max(sapply(acc_list, function(a) max(a$sites, na.rm = TRUE)))
  ymax <- max(sapply(acc_list, function(a) max(a$richness, na.rm = TRUE)))

  plot(
    0,
    0,
    type = "n",
    xlim = c(1, xmax),
    ylim = c(0, ymax),
    xlab = "Samples",
    ylab = "Species richness",
    main = paste(title_text, "–", year_label)
  )

  cols <- seq_along(good_elevs)

  for (i in seq_along(good_elevs)) {
    acc <- acc_list[[i]]
    lines(acc$sites, acc$richness, col = cols[i], lty = i, lwd = 2) # thicker lines
  }

  legend(
    "bottomright",
    legend = good_elevs,
    col = cols,
    lty = seq_along(good_elevs),
    lwd = 2, # match curve thickness
    bty = "n"
  )
}

# 10. Make six plots
data_2002 <- subset_by_year(2002)
data_2003 <- subset_by_year(2003)
data_2022 <- subset_by_year(2022)
data_2023 <- subset_by_year(2023)

# 2002
plot_group(data_2002, group1, "2002", "0-900 m a.s.l. W")
plot_group(data_2002, group2, "2002", "1100-1900 m a.s.l. W")
plot_group(data_2002, group3, "2002", "1700-500 m a.s.l. E")

# 2003
plot_group(data_2002, group1, "2003", "0-900 m a.s.l. W")
plot_group(data_2002, group2, "2003", "1100-1900 m a.s.l. W")
plot_group(data_2002, group3, "2003", "1700-500 m a.s.l. E")

# # 2004
# plot_group(data_2002, group1, "2004", "0-900 m a.s.l. W")
# plot_group(data_2002, group2, "2004", "1100-1900 m a.s.l. W")
# plot_group(data_2002, group3, "2004", "1700-500 m a.s.l. E")

# # 2005
# plot_group(data_2002, group1, "2005", "0-900 m a.s.l. W")
# plot_group(data_2002, group2, "2005", "1100-1900 m a.s.l. W")
# plot_group(data_2002, group3, "2005", "1700-500 m a.s.l. E")

# 2022
plot_group(data_2022, group1, "2022", "0-900 m a.s.l. W")
plot_group(data_2022, group2, "2022", "1100-1900 m a.s.l. W")
plot_group(data_2022, group3, "2022", "1700-500 m a.s.l. E")

# 2023
plot_group(data_2023, group1, "2023", "0-900 m a.s.l. W")
plot_group(data_2023, group2, "2023", "1100-1900 m a.s.l. W")
plot_group(data_2023, group3, "2023", "1700-500 m a.s.l. E")


########Years pulled together##########

# 1. Read CSV
beetles <- beet_wide %>%
  select(-ID, -season) %>%
  rename(
    "Year" = "year",
    "Elevation" = "site",
    "Replicate" = "replicate",
    "Pitfall.trap" = "trap"
  ) %>%
  mutate(Elevation = as.factor(alt_labels[Elevation])) %>%
  relocate(Year, Elevation, Replicate, Pitfall.trap)
# Drop empty column if present
beetles <- beetles[, !grepl("^Unnamed", names(beetles))]

# 2. Standardise column names
names(beetles) <- make.names(names(beetles))

# 3. Identify metadata and species columns
meta_cols <- c("Year", "Elevation", "Replicate", "Pitfall.trap")
species_cols <- setdiff(names(beetles), meta_cols)

# 4. Convert species columns to numeric
beetles[, species_cols] <- lapply(beetles[, species_cols], function(x) {
  suppressWarnings(as.numeric(x))
})

# 5. Check contents
cat("Years present:\n")
print(sort(unique(beetles$Year)))

cat("Elevations present:\n")
print(unique(beetles$Elevation))

# 6. Define elevation groups (MATCH DATA EXACTLY)
group1 <- c(alt_labels[1:6])

group2 <- c(alt_labels[7:10])

group3 <- c(alt_labels[12:17])

# 7. Compute accumulation curve for one elevation (ALL YEARS POOLED)
get_accum <- function(data_all, elev_label) {
  sub <- data_all[data_all$Elevation == elev_label, ]

  # Need at least 2 samples
  if (nrow(sub) < 2) {
    cat("Skipping:", elev_label, "| too few samples\n")
    return(NULL)
  }

  comm <- sub[, species_cols, drop = FALSE]
  comm <- comm[, colSums(comm, na.rm = TRUE) > 0, drop = FALSE]

  # Need at least 2 species
  if (ncol(comm) < 2) {
    cat("Skipping:", elev_label, "| too few species\n")
    return(NULL)
  }

  cat(
    "Elevation =",
    elev_label,
    "| Samples =",
    nrow(sub),
    "| Species cols =",
    ncol(comm),
    "\n"
  )

  specaccum(comm, method = "random", permutations = 1000)
}

# 8. Plot one elevation group (ALL YEARS)
plot_group <- function(data_all, group_elevs, title_text) {
  acc_list <- lapply(group_elevs, function(e) get_accum(data_all, e))
  keep <- !sapply(acc_list, is.null)
  acc_list <- acc_list[keep]
  good_elevs <- group_elevs[keep]

  if (length(acc_list) == 0) {
    plot(
      1,
      1,
      type = "n",
      main = title_text,
      xlab = "",
      ylab = "",
      sub = "No data available"
    )
    return()
  }

  xmax <- max(sapply(acc_list, function(a) max(a$sites)))
  ymax <- max(sapply(acc_list, function(a) max(a$richness)))

  plot(
    0,
    0,
    type = "n",
    xlim = c(1, xmax),
    ylim = c(0, ymax),
    xlab = "Samples",
    ylab = "Species richness",
    main = title_text
  )

  cols <- seq_along(good_elevs)

  for (i in seq_along(good_elevs)) {
    acc <- acc_list[[i]]
    lines(acc$sites, acc$richness, col = cols[i], lty = i, lwd = 2)
  }

  legend(
    "bottomright",
    legend = good_elevs,
    col = cols,
    lty = seq_along(good_elevs),
    lwd = 2,
    bty = "n"
  )
}

# 9. FINAL: THREE PLOTS (YEARS POOLED)

# West – low elevations
plot_group(beetles, group1, "0–900 m a.s.l. W (2002-2005 and 2022-2023)")

# West – high elevations
plot_group(beetles, group2, "1100–1900 m a.s.l. W (2002-2005 and 2022-2023)")

# East – elevations
plot_group(beetles, group3, "1700–500 m a.s.l. E (2002-2005 and 2022-2023)")
