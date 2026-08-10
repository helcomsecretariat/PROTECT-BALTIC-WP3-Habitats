########## Setup ##########
library(terra)
library(mdendro)
library(dendextend)
library(cluster)
library(dplyr)
library(geodata)
library(randomcoloR)
library(colorspace)
library(dynamicTreeCut)
library(fpc)
library(ape)
library(vegan)

run_id <- readRDS("current_run_id.RDS")
#run_id <- "run1_4manSalinityZones_10perc_binStrictLVL"

fp <- paste0("outputs/", run_id, "/dendrograms/")
l <- list.files(fp, recursive = TRUE, full.names = TRUE)
l <- l[grepl("_species_df.csv", l)]

zones <- basename(dirname(l))

for(i in 1:length(zones)){
  temp <- read.csv(l[i])
  temp$zone <- zones[i]
  temp$cluster_name <- paste0(zones[i], "_", temp$cluster)
  if(i == 1){df <- temp}
  if(i > 1){df <- rbind(df, temp)}
}

# Cluter zones
df_zones <- df %>%
  group_by(cluster_name) %>%
  summarise(zone = first(zone))

k <- length(levels(as.factor(df$cluster_name)))
clust_names <- levels(as.factor(df$cluster_name))

## Specify SDMs
#spec_df <- read.csv("biotope_species_list.csv")
#
## Filter based on mobility
#spec_df <- filter(spec_df, mobility == "non-mobile")
#lf <- paste0("inputs/sdms/", spec_df$scientific_name, ".tif")
#spec_n <- spec_df$scientific_name
#
# Import rasters
#r <- rast(lf)
#names(r) <- spec_n

for(i in 1:k){
 
  ########## Import matrix ###########
  zone <- df_zones$zone[i]
  clust <- df_zones$cluster_name[i]
  #out_dir <- paste0("outputs/", run_id, "/dendrograms/", zone)
  #if(!file.exists(out_dir)){dir.create(out_dir)}
  fp <- paste0("outputs/", run_id, "/matrices/", zone, ".rds")
  m_diss <- readRDS(fp)
  m <- as.matrix(m_diss)
  
  m <- 1-m
  #m <- round(m, digits = 3)
  
  
  # Identify species in cluster
  com_spec <- df$species[which(df$cluster_name == clust)]
  m <- m[com_spec,com_spec]

  # Calculate the mean Sorensen-Dice (or Bray-Curtis) coefficients across
  # species in the group
  spec_cor <- colMeans(m)
  
  clust_cors <- data.frame(cluster_name = clust_names[i], 
                           scientific_name = names(spec_cor),
                           mean_coefficient = spec_cor)
  
  row.names(clust_cors) <- NULL
  
  if(i == 1){df_cors <- clust_cors}
  if(i > 1){df_cors <- rbind(df_cors, clust_cors)}
  
}

# Export
write.csv(df_cors, paste0("outputs/", run_id, "/species_cluster_coefficients.csv"))


