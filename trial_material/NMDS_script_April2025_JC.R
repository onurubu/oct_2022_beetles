#### Coskun Kucukkaragoz: I have edited this script provided by Abusisiwe Ndaba to suit my needs

# Load necessary libraries
library(tidyverse)
library(vegan)
library(ggplot2)
library(data.table)

#####Run this before to start####
# Read in your data

data_all <- abund_all %>% filter(rowSums(across(where(is.numeric)))!=0)

# data_all <- fread("complete NMDS data_JC.csv") %>% mutate(across(where(is.character), factor)) %>% mutate(across(where(is.integer), factor))
str(data_all)
summary(data_all)

#####To be reported for NMDS results#####
#Plot + stress value (last line of all NMDS runs) + ANOSIM R and significance

#####Graph to compare all sampling events (all data 2002, 2003, 2022 and 2023)#####
# Select variables
data_species <- data_all[, 5:8] # exclude Replicate, year, group and site columns

# Convert all columns to numeric (handle factors/characters properly)
data_species <- data.frame(lapply(data_species, function(x) as.numeric(as.character(x))))
# Check for NAs
if (anyNA(data_species)) {
  cat("NAs found - removing rows with NAs.\n")
  data_species <- na.omit(data_species)
}
# Apply log(x + 1) transformation
#data_species <- log1p(data_species)

# Convert to Presence-absence
data_species<-decostand(data_species,"pa") #switch to presence-absence data

# Convert species abundance data to a matrix
species_data <- data_species %>% as.matrix(.)

# Perform NMDS analysis
nmds <- metaMDS(species_data, distance = "bray", k = 2, trymax = 100)

# Create a data frame for plotting
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event = data_all$year)

# Plot the NMDS results - Community in 2002 & 2022
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sampling_event)) +
  geom_point(aes(shape= Sampling_event, fill=Sampling_event)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,21,24,24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Sampling_event")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_all$year)
summary(anosim_type)

##################Graph to compare years (all data 2002/2003 and 2022/2023)#####
# Select variables
data_species<-data_all[,5:8] # exclude Replicate, year, group and site columns

# Convert all columns to numeric (handle factors/characters properly)
data_species <- data.frame(lapply(data_species, function(x) as.numeric(as.character(x))))
# Check for NAs
if (anyNA(data_species)) {
  cat("NAs found - removing rows with NAs.\n")
  data_species <- na.omit(data_species)
}
# Apply log(x + 1) transformation
#data_species <- log1p(data_species)

# Convert to Presence-absence
data_species<-decostand(data_species,"pa") #switch to presence-absence data

# Convert species abundance data to a matrix
species_data <- as.matrix(data_species)  

# Perform NMDS analysis
nmds <- metaMDS(species_data, distance = "bray", k = 2, trymax = 100)

# Create a data frame for plotting
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event = data_all$period)

# Plot the NMDS results - Community in 2002 & 2022
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sampling_event)) +
  geom_point(aes(shape= Sampling_event, fill=Sampling_event)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,21,24,24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Sampling_event")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_all$year)
summary(anosim_type)

##########ALTITUDE##################--------------------------------------------
#####Graph for all altitude one year (all sites for the selected year)#####
####Year 1####
data_species<-data_all[data_all$period=="old",5:8] #Select year number here !!
data_factor<-data_all[data_all$period=="old",1:4] #Select year number here !!
title<-"Year 2002" #Select year number here !!

# Convert to Presence-absence
data_species<-decostand(data_species,"pa") #switch to presence-absence data

# Apply log(x + 1) transformation
#data_species <- log1p(data_species)

species_data <- as.matrix(data_species)  

nmds <- metaMDS(species_data, distance = "bray", k = 2, trymax = 100)
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event=data_factor$Date, Season = data_factor$Season, Year = data_factor$Year, Site = data_factor$Site)
# Plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Site)) +
  ggtitle(title)+
  geom_point(aes(shape = Site, fill= Site)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,22,23,24,25,21,22,23,24,25,21,22,23,24,25,21,22)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Site")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_factor$Site)
summary(anosim_type)

####Year 2####
data_species<-data_all[data_all$Year=="2022-2023", 6:120] #Select year number here !!
data_factor<-data_all[data_all$Year=="2022-2023",1:5] #Select year number here !!
title<-"Year 2022-2023" #Select year number here !!

# Convert to Presence-absence
data_species<-decostand(data_species,"pa") #switch to presence-absence data

# Apply log(x + 1) transformation
#data_species <- log1p(data_species)

species_data <- as.matrix(data_species)  

nmds <- metaMDS(species_data, distance = "bray", k = 2, trymax = 100)
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event=data_factor$Date, Season = data_factor$Season, Year = data_factor$Year, Site = data_factor$Site)
# Plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Site)) +
  ggtitle(title)+
  geom_point(aes(shape = Site, fill= Site)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,22,23,24,25,21,22,23,24,25,21,22,23,24,25,21,22)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Site")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_factor$Site)
summary(anosim_type)


#####Graph for each altitude separated years & grouped seasons (all replicates for the selected site for both years)#####
#Here site 1 as an exemple - change site number for other sites
data_species<-data_all[data_all$Site=="17",6:120] #Select site number here !!
data_factor<-data_all[data_all$Site=="17",1:5] #Select site number here !!
title<-"Site 17" #Select site number here !!

# Convert to Presence-absence
data_species<-decostand(data_species,"pa") #switch to presence-absence data

# Apply log(x + 1) transformation
#data_species <- log1p(data_species)

species_data <- as.matrix(data_species)

nmds <- metaMDS(species_data, distance = "bray", k = 2)
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event=data_factor$Date, Season = data_factor$Season, Year = data_factor$Year, Site = data_factor$Site)
# Plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Year)) +
  ggtitle(title)+
  geom_point(aes(shape= Year, fill=Year)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,21,24,24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Year")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_factor$Year)
summary(anosim_type)


#####Graph for each altitude separated years & separated seasons (all replicates for the selected site for both years)#####
#Here site 1 as an exemple - change site number for other sites
data_species<-data_all[data_all$Site=="17",6:120] #Select site number here !!
data_factor<-data_all[data_all$Site=="17",1:5] #Select site number here !!
title<-"Site 17" #Select site number here !!

# Convert to Presence-absence
data_species<-decostand(data_species,"pa") #switch to presence-absence data

#Apply log(x + 1) transformation
#data_species <- log1p(data_species)

species_data <- as.matrix(data_species)

nmds <- metaMDS(species_data, distance = "bray", k = 2)
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event=data_factor$Date, Season = data_factor$Season, Year = data_factor$Year, Site = data_factor$Site)
# Plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sampling_event)) +
  ggtitle(title)+
  geom_point(aes(shape= Sampling_event, fill=Sampling_event)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,21,24,24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Sampling_event")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_factor$Sampling_event)
summary(anosim_type)




######################################################--------------------------------------------

#################VEGETATION TYPES#####################--------------------------------------------

######################################################--------------------------------------------

# Load necessary libraries
library(vegan)
library(ggplot2)

# Set working directory (both seasons csv.csv)
setwd("~/Documents/Analysis/PhD NMDS data")

#####Run this before to start####
# Read in your data
# data_all <- read.csv("complete NMDS data_JC.csv", sep=";", dec=",")
# data_all[sapply(data_all, is.character)] <- lapply(data_all[sapply(data_all, is.character)], as.factor)
# data_all[sapply(data_all, is.integer)] <- lapply(data_all[sapply(data_all, is.integer)], as.factor)
# str(data_all)
# summary(data_all)

#####To be reported for NMDS results#####
#Plot + stress value (last line of all NMDS runs) + ANOSIM R and significance

#####Graph to compare all seasons (all data 2002, 2003, 2022 and 2023)#####
# Select variables
data_species<-data_all[, 6:120] # exclude Replicate, year, group and site columns

# Convert all columns to numeric (handle factors/characters properly)
data_species <- data.frame(lapply(data_species, function(x) as.numeric(as.character(x))))
# Check for NAs
if (anyNA(data_species)) {
  cat("NAs found - removing rows with NAs.\n")
  data_species <- na.omit(data_species)
}
# Apply log(x + 1) transformation
data_species <- log1p(data_species)

# Convert to Presence-absence
#data_species<-decostand(data_species,"pa") #switch to presence-absence data

# Convert species abundance data to a matrix
species_data <- as.matrix(data_species)  

# Perform NMDS analysis
nmds <- metaMDS(species_data, distance = "bray", k = 3, trymax = 100)

# Create a data frame for plotting
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event = data_all$Date)

# Plot the NMDS results - Community in 2002 & 2022
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sampling_event)) +
  geom_point(aes(shape= Sampling_event, fill=Sampling_event)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,21,24,24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Sampling_event")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_all$Date)
summary(anosim_type)

########### without singletons #######################################################
# ==============================
# 1. Load necessary libraries
# ==============================
library(vegan)
library(ggplot2)

# ==============================
# 2. Set working directory
# ==============================
setwd("~/Documents/Analysis/PhD NMDS data")

##### Run this before to start #####
# Read in your data
data_all <- read.csv("complete NMDS data_JC_no singletons.csv", sep=";", dec=",")

# Convert character/integer columns to factor
data_all[sapply(data_all, is.character)] <- lapply(data_all[sapply(data_all, is.character)], as.factor)
data_all[sapply(data_all, is.integer)] <- lapply(data_all[sapply(data_all, is.integer)], as.factor)

str(data_all)
summary(data_all)

##### Graph to compare years (all data 2002 and 2022) #####
# Select only species columns
data_species <- data_all[, 7:109]  # adjust indices if your species columns differ

# Convert all columns to numeric
data_species <- data.frame(lapply(data_species, function(x) as.numeric(as.character(x))))

# Apply log(x + 1) transformation
data_species <- log1p(data_species)

# ==============================
# 5. Convert to presence/absence if needed (optional)
# Uncomment below if you want PA data
# data_species <- decostand(data_species, "pa")

# ==============================
# 6. Convert to matrix
# ==============================
species_data <- as.matrix(data_species)

# ==============================
# 7. Perform NMDS
# ==============================
nmds <- metaMDS(species_data, distance = "bray", k = 3, trymax = 100)

# ==============================
# 8. Prepare data for plotting
# ==============================
nmds_df <- data.frame(NMDS1 = nmds$points[, 1],
                      NMDS2 = nmds$points[, 2],
                      Sampling_event = data_all$Sampling_event)

# ==============================
# 9. Plot NMDS results
# ==============================
# Create NMDS plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sampling_event)) +
  geom_point(aes(shape= Sampling_event, fill=Sampling_event)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,21,24,24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Sampling_event")

# ==============================
# 10. ANOSIM test for group differences
# ==============================
dist <- vegdist(species_data, method = "bray")
anosim_type <- anosim(dist, data_all$Sampling_event)
summary(anosim_type)








##################Graph to compare years (all data 2002/2003 and 2022/2023)#####
# Select variables
library(vegan)
library(ggplot2)

# Set working directory (both seasons csv.csv)
setwd("~/Documents/Analysis/PhD NMDS data")

#####Run this before to start####
# Read in your data
data_all <- read.csv("complete NMDS data_JC.csv", sep=";", dec=",")
data_all[sapply(data_all, is.character)] <- lapply(data_all[sapply(data_all, is.character)], as.factor)
data_all[sapply(data_all, is.integer)] <- lapply(data_all[sapply(data_all, is.integer)], as.factor)
str(data_all)
summary(data_all)

data_species<-data_all[,7:120] # exclude Replicate, year, group and site columns

# Convert all columns to numeric (handle factors/characters properly)
data_species <- data.frame(lapply(data_species, function(x) as.numeric(as.character(x))))
# Check for NAs
if (anyNA(data_species)) {
  cat("NAs found - removing rows with NAs.\n")
  data_species <- na.omit(data_species)
}

# Apply log(x + 1) transformation
data_species <- log1p(data_species)

# Convert to Presence-absence
#data_species<-decostand(data_species,"pa") #switch to presence-absence data

# Convert species abundance data to a matrix
species_data <- as.matrix(data_species)  

# Perform NMDS analysis
nmds <- metaMDS(species_data, distance = "bray", k = 3, trymax = 100)

# Create a data frame for plotting
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event = data_all$Sampling_event)

# Plot the NMDS results - Community in 2002 & 2022
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sampling_event)) +
  geom_point(aes(shape= Sampling_event, fill=Sampling_event)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,21,24,24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Sampling_event")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_all$Sampling_event)
summary(anosim_type)


################################# without singletons ###########################
library(vegan)
library(ggplot2)

# Set working directory (both seasons csv.csv)
setwd("~/Documents/Analysis/PhD NMDS data")

#####Run this before to start####
# Read in your data
data_all <- read.csv("complete NMDS data_JC.csv", sep=";", dec=",")
data_all[sapply(data_all, is.character)] <- lapply(data_all[sapply(data_all, is.character)], as.factor)
data_all[sapply(data_all, is.integer)] <- lapply(data_all[sapply(data_all, is.integer)], as.factor)
str(data_all)
summary(data_all)

# ==============================
# Select variables
# ==============================
data_species <- data_all[, 7:109]   # adjust if your species columns differ

# Convert all columns to numeric (handle factors/characters properly)
data_species <- data.frame(lapply(data_species, function(x) as.numeric(as.character(x))))


# Apply log(x + 1) transformation
data_species <- log1p(data_species)

# ==============================
# (Optional) Convert to Presence–Absence
# ==============================
# data_species <- decostand(data_species, "pa")  # uncomment if you want PA data instead of abundance

# ==============================
# Convert to matrix
# ==============================
species_data <- as.matrix(data_species)

# ==============================
# Perform NMDS analysis
# ==============================
nmds <- metaMDS(species_data, distance = "bray", k = 3, trymax = 100)

# ==============================
# Prepare data for plotting
# ==============================
nmds_df <- data.frame(
  NMDS1 = nmds$points[, 1],
  NMDS2 = nmds$points[, 2],
  Sampling_event = data_all$Sampling_event
)

# ==============================
# Plot the NMDS results
# ==============================
# Plot the NMDS results - Community in 2002 & 2022
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sampling_event)) +
  geom_point(aes(shape= Sampling_event, fill=Sampling_event)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,21,24,24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Sampling_event")

# ==============================
# Test group differences with ANOSIM
# ==============================
dist <- vegdist(species_data, method = "bray")
anosim_type <- anosim(dist, data_all$Sampling_event)
summary(anosim_type)



##########VEGETATION TYPE##################--------------------------------------------
#####Graph for all vegetation type one year (all sites for the selected year)#####
####Year 1####
library(vegan)
library(ggplot2)

# Set working directory (both seasons csv.csv)
setwd("~/Documents/Analysis/PhD NMDS data")

#####Run this before to start####
# Read in your data
data_all <- read.csv("complete NMDS data_JC_no singletons.csv", sep=";", dec=",")
data_all[sapply(data_all, is.character)] <- lapply(data_all[sapply(data_all, is.character)], as.factor)
data_all[sapply(data_all, is.integer)] <- lapply(data_all[sapply(data_all, is.integer)], as.factor)
str(data_all)
summary(data_all)

data_species<-data_all[data_all$Year=="2002-2003",7:109] #Select year number here !!
data_factor<-data_all[data_all$Year=="2002-2003",1:6] #Select year number here !!
title<-"Year 2002-2003" #Select year number here !!

# Convert to Presence-absence
#data_species<-decostand(data_species,"pa") #switch to presence-absence data

# Apply log(x + 1) transformation
data_species <- log1p(data_species)

species_data <- as.matrix(data_species)  

nmds <- metaMDS(species_data, distance = "bray", k = 3, trymax = 100)
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event=data_factor$Sampling_event, Season = data_factor$Season, Year = data_factor$Year, Vegetation_type = data_factor$Vegetation_type)
# Plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Vegetation_type)) +
  ggtitle(title)+
  geom_point(aes(shape = Vegetation_type, fill= Vegetation_type)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,22,23,24,25,21,22,23,24,25,21,22,23,24,25,21,22)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Vegetation_type")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_factor$Vegetation_type)
summary(anosim_type)



####Year 2####
library(vegan)
library(ggplot2)

# Set working directory (both seasons csv.csv)
setwd("~/Documents/Analysis/PhD NMDS data")

#####Run this before to start####
# Read in your data
data_all <- read.csv("complete NMDS data_JC_no singletons.csv", sep=";", dec=",")
data_all[sapply(data_all, is.character)] <- lapply(data_all[sapply(data_all, is.character)], as.factor)
data_all[sapply(data_all, is.integer)] <- lapply(data_all[sapply(data_all, is.integer)], as.factor)
str(data_all)
summary(data_all)

data_species<-data_all[data_all$Year=="2022-2023",7:109] #Select year number here !!
data_factor<-data_all[data_all$Year=="2022-2023",1:6] #Select year number here !!
title<-"Year 2022-2023" #Select year number here !!

# Convert to Presence-absence
#data_species<-decostand(data_species,"pa") #switch to presence-absence data

# Apply log(x + 1) transformation
data_species <- log1p(data_species)

species_data <- as.matrix(data_species)  

nmds <- metaMDS(species_data, distance = "bray", k = 3, trymax = 100)
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event=data_factor$Sampling_event, Season = data_factor$Season, Year = data_factor$Year,Vegetation_type = data_factor$Vegetation_type)
# Plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Vegetation_type)) +
  ggtitle(title)+
  geom_point(aes(shape = Vegetation_type, fill= Vegetation_type)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,22,23,24,25,21,22,23,24,25,21,22,23,24,25,21,22)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Vegetation_type")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_factor$Vegetation_type)
summary(anosim_type)


#####Graph for each vegetation type separated years & grouped seasons (all replicates for the selected site for both years)#####
library(vegan)
library(ggplot2)

# Set working directory (both seasons csv.csv)
setwd("~/Documents/Analysis/PhD NMDS data")

#####Run this before to start####
# Read in your data
data_all <- read.csv("complete NMDS data_JC.csv", sep=";", dec=",")
data_all[sapply(data_all, is.character)] <- lapply(data_all[sapply(data_all, is.character)], as.factor)
data_all[sapply(data_all, is.integer)] <- lapply(data_all[sapply(data_all, is.integer)], as.factor)
str(data_all)
summary(data_all)

#Here Alpine as an example - change veg types for other veg types
data_species<-data_all[data_all$Vegetation_type=="Alpine",7:120] #Select vegetation type here !!
data_factor<-data_all[data_all$Vegetation_type=="Alpine",1:6] #Select vegetation type here !!
title<-"Alpine" #Select vegetation type here !!

# Convert to Presence-absence
data_species<-decostand(data_species,"pa") #switch to presence-absence data

# Apply log(x + 1) transformation
#data_species <- log1p(data_species)

species_data <- as.matrix(data_species)

nmds <- metaMDS(species_data, distance = "bray", k = 3)
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Date=data_factor$Date, Season = data_factor$Season, Sampling_event = data_factor$Sampling_event, Vegetation_type = data_factor$Vegetation_type)
# Plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sampling_event)) +
  ggtitle(title)+
  geom_point(aes(shape= Sampling_event, fill=Sampling_event)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,21,24,24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Sampling_event")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_factor$Sampling_event)
summary(anosim_type)


#####Graph for each altitude separated years & separated seasons (all replicates for the selected site for both years)#####
library(vegan)
library(ggplot2)

# Set working directory (both seasons csv.csv)
setwd("~/Documents/Analysis/PhD NMDS data")

#####Run this before to start####
# Read in your data
data_all <- read.csv("complete NMDS data_JC.csv", sep=";", dec=",")
data_all[sapply(data_all, is.character)] <- lapply(data_all[sapply(data_all, is.character)], as.factor)
data_all[sapply(data_all, is.integer)] <- lapply(data_all[sapply(data_all, is.integer)], as.factor)
str(data_all)
summary(data_all)

#Here is Alpine as an example - change veg type for other veg types
data_species<-data_all[data_all$Vegetation_type=="Alpine",7:120] #Select vegetation type here !!
data_factor<-data_all[data_all$Vegetation_type=="Alpine",1:6] #Select vegetation type here !!
title<-"Alpine" #Select vegetation type here !!

# Convert to Presence-absence
#data_species<-decostand(data_species,"pa") #switch to presence-absence data

#Apply log(x + 1) transformation
data_species <- log1p(data_species)

species_data <- as.matrix(data_species)

nmds <- metaMDS(species_data, distance = "bray", k = 3)
nmds_df <- data.frame(NMDS1 = nmds$points[, 1], NMDS2 = nmds$points[, 2], 
                      Sampling_event=data_factor$Sampling_event, Season = data_factor$Season, Year = data_factor$Year, Vegetation_type = data_factor$Vegetation_type)
# Plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sampling_event)) +
  ggtitle(title)+
  geom_point(aes(shape= Sampling_event, fill=Sampling_event)) +
  stat_ellipse() +
  theme_classic(base_size=14) +
  scale_shape_manual(values=c(21,21,24,24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Sampling_event")
#Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type<-anosim(dist,data_factor$Sampling_event)
summary(anosim_type)


###### Without singletons ######
library(vegan)
library(ggplot2)

# Set working directory (both seasons csv.csv)
setwd("~/Documents/Analysis/PhD NMDS data")

#####Run this before to start####
# Read in your data
data_all <- read.csv("complete NMDS data_JC_no singletons.csv", sep=";", dec=",")
data_all[sapply(data_all, is.character)] <- lapply(data_all[sapply(data_all, is.character)], as.factor)
data_all[sapply(data_all, is.integer)] <- lapply(data_all[sapply(data_all, is.integer)], as.factor)
str(data_all)
summary(data_all)

#Here is Alpine as an example - change veg type for other veg types
data_species <- data_all[data_all$Vegetation_type == "Succulent Karoo", 7:109]  # Select vegetation type here !!
data_factor  <- data_all[data_all$Vegetation_type == "Succulent Karoo", 1:6]    # Select vegetation type here !!
title <- "Succulent Karoo"  # Select vegetation type here !!

# Convert to Presence-absence if desired
 data_species <- decostand(data_species, "pa")  # switch to presence-absence data

# Apply log(x + 1) transformation
#data_species <- log1p(data_species)

species_data <- as.matrix(data_species)

# NMDS
nmds <- metaMDS(species_data, distance = "bray", k = 3)

nmds_df <- data.frame(
  NMDS1 = nmds$points[, 1],
  NMDS2 = nmds$points[, 2],
  Sampling_event = data_factor$Sampling_event,
  Season = data_factor$Season,
  Year = data_factor$Year,
  Vegetation_type = data_factor$Vegetation_type
)

# Plot
ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = Sampling_event)) +
  ggtitle(title) +
  geom_point(aes(shape = Sampling_event, fill = Sampling_event)) +
  stat_ellipse() +
  theme_classic(base_size = 14) +
  scale_shape_manual(values = c(21, 21, 24, 24)) +
  labs(x = "NMDS Axis 1", y = "NMDS Axis 2", color = "Sampling_event")

# Check significance of the difference between groups
dist <- vegdist(species_data, method = "bray")
anosim_type <- anosim(dist, data_factor$Sampling_event)
summary(anosim_type)

