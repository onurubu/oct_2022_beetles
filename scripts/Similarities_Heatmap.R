
###################################### Heatmap with singletons ###################

# ==========================================
# 1. Libraries
# ==========================================
{library(vegan)
library(tidyverse)     # includes ggplot2, dplyr, tidyr
library(ggvenn)        # ggplot-friendly Venn diagrams
library(openxlsx)}

# ==========================================
# 2. Set working directory & read data
# ==========================================
# setwd("~/Documents/Analysis/Similarities")

# dat <- read.csv("Sampling events_elevation.csv", 
                # sep = ";", header = TRUE, 
                # check.names = FALSE, stringsAsFactors = FALSE)

dat <- abund_all %>% mutate(veg_type = factor(case_when(site %in% c(1) ~ "strandveld", site %in% c(2, 6, 16) ~ "restioid", site %in% c(3, 4, 5) ~ "proteoid", site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous", site %in% c(11) ~ "alpine", site %in% c(17) ~ "succulent_karoo"))) %>% relocate(site, replicate, year, veg_type) %>% group_by(period, year, season, site, replicate, veg_type) %>% summarise(across(where(is.numeric), sum), .groups = "drop") %>% filter(paste(year, season) != "2023 March") %>% group_by(period, site) %>% summarise(across(where(is.numeric), sum), .groups = "drop") %>% mutate(veg_type = factor(case_when(site %in% c(1) ~ "strandveld", site %in% c(2, 6, 16) ~ "restioid", site %in% c(3, 4, 5) ~ "proteoid", site %in% c(7, 8, 9, 10, 12, 13, 14, 15) ~ "ericaceous", site %in% c(11) ~ "alpine", site %in% c(17) ~ "succulent_karoo"))) %>% relocate(period, site, veg_type) %>% select(-veg_type)

cat("Data loaded with", nrow(dat), "rows and", ncol(dat), "columns\n")

# ==========================================
# 3. Convert abundance to presence/absence
# ==========================================
species_cols <- setdiff(names(dat), c("period", "veg_type", "site"))
dat[, species_cols] <- lapply(dat[, species_cols], function(x) ifelse(x > 0, 1, 0))

# ==========================================
# 4. Aggregate by period × Elevation
# ==========================================
dat_agg <- aggregate(. ~ period + `site`,
                     data = dat, FUN = max)

# Split into the two sampling events
dat_2002 <- subset(dat_agg, period == "old")
dat_2022 <- subset(dat_agg, period == "modern")

# Drop the period column
dat_2002 <- dat_2002[, -1]
dat_2022 <- dat_2022[, -1]

# ==========================================
# 5. Define elevation order
# ==========================================
elev_levels <- c(
  "0 m a.s.l. W (Strandveld)",
  "200 m a.s.l. W (Restioid)",
  "300 m a.s.l. W (Proteoid)",
  "500 m a.s.l. W (Proteoid)",
  "700 m a.s.l. W (Proteoid)",
  "900 m a.s.l. W (Restioid)",
  "1100 m a.s.l. W (Ericaceous)",
  "1300 m a.s.l. W (Ericaceous)",
  "1500 m a.s.l. W (Ericaceous)",
  "1700 m a.s.l. W (Ericaceous)",
  "1900 m a.s.l. (Peak) (Alpine)",
  "1700 m a.s.l. E (Ericaceous)",
  "1500 m a.s.l. E (Ericaceous)",
  "1300 m a.s.l. E (Ericaceous)",
  "1100 m a.s.l. E (Ericaceous)",
  "900 m a.s.l. E (Restioid)",
  "500 m a.s.l. E (Succulent_Karoo)"
)

dat_2002$`site` <- factor(elev_levels)
dat_2022$`site` <- factor(elev_levels)

# Sort by elevation
dat_2002 <- dat_2002[order(dat_2002$`site`), ]
dat_2022 <- dat_2022[order(dat_2022$`site`), ]

# Set row names and drop factor column
rownames(dat_2002) <- dat_2002$`site`
rownames(dat_2022) <- dat_2022$`site`
dat_2002 <- dat_2002[, -1]
dat_2022 <- dat_2022[, -1]

# ==========================================
# 6. Compute Jaccard similarity
# ==========================================
mat_2002 <- as.matrix(dat_2002)
mat_2022 <- as.matrix(dat_2022)

sim_mat <- matrix(NA, nrow = nrow(mat_2002), ncol = nrow(mat_2022),
                  dimnames = list(rownames(mat_2002), rownames(mat_2022)))

for (i in seq_len(nrow(mat_2002))) {
  for (j in seq_len(nrow(mat_2022))) {
    sim_mat[i, j] <- 1 - vegdist(rbind(mat_2002[i, ], mat_2022[j, ]), method = "jaccard")
  }
}

# ==========================================
# 7. Reshape similarity matrix for ggplot
# ==========================================
heat_df <- as.data.frame(as.table(sim_mat))
colnames(heat_df) <- c("Site_2002", "Site_2022", "Similarity")

# Reverse rows so top = highest elevation
heat_df$Site_2002 <- factor(heat_df$Site_2002, levels = rev(elev_levels))
heat_df$Site_2022 <- factor(heat_df$Site_2022, levels = elev_levels)

# ==========================================
# 8. Create heatmap with ggplot
# ==========================================
heatmap_plot <- ggplot(heat_df, aes(x = Site_2022, y = Site_2002, fill = Similarity)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Jaccard\nSimilarity") +
  labs(x = "2022–2023", y = "2002–2003", 
       title = "Jaccard Similarity Heatmap") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# ==========================================
# 9. Prepare Venn diagram data
# ==========================================
venn_cols <- c("blue", "red")   # 2002-2003 = blue, 2022-2023 = red

set_2002_all <- colnames(dat_2002)[colSums(dat_2002) > 0]
set_2022_all <- colnames(dat_2022)[colSums(dat_2022) > 0]

overall_species <- list(
  "old" = set_2002_all,
  "modern" = set_2022_all
)

overall_venn <- ggvenn(
  overall_species, 
  fill_color = venn_cols,
  stroke_size = 0.8,
  set_name_size = 5
) + labs(title = "Overall Species Overlap")

# ==========================================
# 10. Output to PDF
# ==========================================
library(openxlsx)    # if not installed, run: install.packages("openxlsx")
write.xlsx(heat_df, "Heatmap_results.xlsx", rowNames = FALSE)

pdf("Heatmap_and_Venn.pdf", width = 10, height = 8)

# ---- Heatmap
print(heatmap_plot)

# ---- Dendrograms (added here, after the heatmap)
dist_2002 <- vegdist(dat_2002[-10,], method = "jaccard")
clust_2002 <- hclust(dist_2002, method = "average")

dist_2022 <- vegdist(dat_2022, method = "jaccard")
clust_2022 <- hclust(dist_2022, method = "average")

plot(clust_2002,
     main = "Dendrogram – 2002–2003",
     xlab = "",
     ylab = "Distance",      # <-- added label
     sub = "")

plot(clust_2022,
     main = "Dendrogram – 2022–2023",
     xlab = "",
     ylab = "Distance",      # <-- added label
     sub = "")

# ---- Venn diagrams
print(overall_venn)

for (site in elev_levels) {
  if (!(site %in% rownames(dat_2002)) || !(site %in% rownames(dat_2022))) next
  
  sp_2002 <- colnames(dat_2002)[dat_2002[site, ] > 0]
  sp_2022 <- colnames(dat_2022)[dat_2022[site, ] > 0]
  
  if (length(sp_2002) == 0 & length(sp_2022) == 0) next
  
  site_species <- list(
    "2002-2003" = sp_2002,
    "2022-2023" = sp_2022
  )
  
  site_venn <- ggvenn(
    site_species, 
    fill_color = venn_cols,
    stroke_size = 0.8,
    set_name_size = 5
  ) + labs(title = paste("Species Overlap –", site))
  
  print(site_venn)
}

dev.off()

cat("All plots (heatmap, dendrograms, and Venn diagrams) saved to Heatmap_and_Venn.pdf\n")






################# Without singletons ############################################################
# ==========================================
# 1. Libraries
# ==========================================
library(vegan)
library(tidyverse)     # includes ggplot2, dplyr, tidyr
library(ggvenn)        # ggplot-friendly Venn diagrams
library(openxlsx)

# ==========================================
# 2. Set working directory & read data
# ==========================================
setwd("~/Documents/Analysis/Similarities")

dat <- read.csv("Sampling events_elevation.csv", 
                sep = ";", header = TRUE, 
                check.names = FALSE, stringsAsFactors = FALSE)

cat("Data loaded with", nrow(dat), "rows and", ncol(dat), "columns\n")

# ==========================================
# 3. Remove true singletons (total abundance = 1) then convert to presence/absence
# ==========================================
species_cols <- setdiff(names(dat), c("period", "site"))

# remove species whose total abundance = 1
singleton_counts <- colSums(dat[, species_cols])
true_singletons <- names(singleton_counts[singleton_counts == 1])
cat("Removing", length(true_singletons), "true singletons\n")
dat <- dat[, !(names(dat) %in% true_singletons)]

# convert remaining species to presence/absence
species_cols <- setdiff(names(dat), c("period", "site"))
dat[, species_cols] <- lapply(dat[, species_cols], function(x) ifelse(x > 0, 1, 0))

# ==========================================
# 4. Aggregate by period × Elevation
# ==========================================
dat_agg <- aggregate(. ~ period + `site`,
                     data = dat, FUN = max)

# Split into the two sampling events
dat_2002 <- subset(dat_agg, period == "2002-2003")
dat_2022 <- subset(dat_agg, period == "2022-2023")

# Drop the period column
dat_2002 <- dat_2002[, -1]
dat_2022 <- dat_2022[, -1]

# ==========================================
# 5. Define elevation order
# ==========================================
elev_levels <- c(
  "0 m a.s.l. W (Strandveld)",
  "200 m a.s.l. W (Restioid)",
  "300 m a.s.l. W (Proteoid)",
  "500 m a.s.l. W (Proteoid)",
  "700 m a.s.l. W (Proteoid)",
  "900 m a.s.l. W (Restioid)",
  "1100 m a.s.l. W (Ericaceous)",
  "1300 m a.s.l. W (Ericaceous)",
  "1500 m a.s.l. W (Ericaceous)",
  "1700 m a.s.l. W (Ericaceous)",
  "1900 m a.s.l. (Peak) (Alpine)",
  "1700 m a.s.l. E (Ericaceous)",
  "1500 m a.s.l. E (Ericaceous)",
  "1300 m a.s.l. E (Ericaceous)",
  "1100 m a.s.l. E (Ericaceous)",
  "900 m a.s.l. E (Restioid)",
  "500 m a.s.l. E (Succulent_Karoo)"
)

dat_2002$`site` <- factor(dat_2002$`site`, levels = elev_levels)
dat_2022$`site` <- factor(dat_2022$`site`, levels = elev_levels)

# Sort by elevation
dat_2002 <- dat_2002[order(dat_2002$`site`), ]
dat_2022 <- dat_2022[order(dat_2022$`site`), ]

# Set row names and drop factor column
rownames(dat_2002) <- dat_2002$`site`
rownames(dat_2022) <- dat_2022$`site`
dat_2002 <- dat_2002[, -1]
dat_2022 <- dat_2022[, -1]

# ==========================================
# 6. Compute Jaccard similarity
# ==========================================
mat_2002 <- as.matrix(dat_2002)
mat_2022 <- as.matrix(dat_2022)

sim_mat <- matrix(NA, nrow = nrow(mat_2002), ncol = nrow(mat_2022),
                  dimnames = list(rownames(mat_2002), rownames(mat_2022)))

for (i in seq_len(nrow(mat_2002))) {
  for (j in seq_len(nrow(mat_2022))) {
    sim_mat[i, j] <- 1 - vegdist(rbind(mat_2002[i, ], mat_2022[j, ]), method = "jaccard")
  }
}

# ==========================================
# 7. Reshape similarity matrix for ggplot
# ==========================================
heat_df <- as.data.frame(as.table(sim_mat))
colnames(heat_df) <- c("Site_2002", "Site_2022", "Similarity")

# Reverse rows so top = highest elevation
heat_df$Site_2002 <- factor(heat_df$Site_2002, levels = rev(elev_levels))
heat_df$Site_2022 <- factor(heat_df$Site_2022, levels = elev_levels)

# ==========================================
# 8. Create heatmap with ggplot
# ==========================================
heatmap_plot <- ggplot(heat_df, aes(x = Site_2022, y = Site_2002, fill = Similarity)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Jaccard\nSimilarity") +
  labs(x = "2022–2023", y = "2002–2003", 
       title = "Jaccard Similarity Heatmap") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5)
  )

# ==========================================
# 9. Prepare Venn diagram data
# ==========================================
venn_cols <- c("blue", "red")   # 2002-2003 = blue, 2022-2023 = red

set_2002_all <- colnames(dat_2002)[colSums(dat_2002) > 0]
set_2022_all <- colnames(dat_2022)[colSums(dat_2022) > 0]

overall_species <- list(
  "2002-2003" = set_2002_all,
  "2022-2023" = set_2022_all
)

overall_venn <- ggvenn(
  overall_species, 
  fill_color = venn_cols,
  stroke_size = 0.8,
  set_name_size = 5
) + labs(title = "Overall Species Overlap")

# ==========================================
# 10. Output to PDF
# ==========================================
library(openxlsx)    # if not installed, run: install.packages("openxlsx")
write.xlsx(heat_df, "Heatmap_results.xlsx", rowNames = FALSE)

pdf("Heatmap_and_Venn.pdf", width = 10, height = 8)

# ---- Heatmap
print(heatmap_plot)

# ---- Dendrograms (placed right after heatmap)
dist_2002 <- vegdist(dat_2002, method = "jaccard")
clust_2002 <- hclust(dist_2002, method = "average")

dist_2022 <- vegdist(dat_2022, method = "jaccard")
clust_2022 <- hclust(dist_2022, method = "average")

plot(clust_2002,
     main = "Dendrogram – 2002–2003",
     xlab = "",
     ylab = "Distance",      # <-- added label
     sub = "")

plot(clust_2022,
     main = "Dendrogram – 2022–2023",
     xlab = "",
     ylab = "Distance",      # <-- added label
     sub = "")

# ---- Venn diagrams
print(overall_venn)

for (site in elev_levels) {
  if (!(site %in% rownames(dat_2002)) || !(site %in% rownames(dat_2022))) next
  
  sp_2002 <- colnames(dat_2002)[dat_2002[site, ] > 0]
  sp_2022 <- colnames(dat_2022)[dat_2022[site, ] > 0]
  
  if (length(sp_2002) == 0 & length(sp_2022) == 0) next
  
  site_species <- list(
    "2002-2003" = sp_2002,
    "2022-2023" = sp_2022
  )
  
  site_venn <- ggvenn(
    site_species, 
    fill_color = venn_cols,
    stroke_size = 0.8,
    set_name_size = 5
  ) + labs(title = paste("Species Overlap –", site))
  
  print(site_venn)
}

dev.off()

cat("All plots (heatmap, dendrograms, and Venn diagrams) saved to Heatmap_and_Venn.pdf\n")

