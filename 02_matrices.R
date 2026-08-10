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
source("_functions/annotate_clusters.R")
source("_functions/save_biotopes.R")
source("_functions/medoid_ttest.R")

# Specify the proportion of the species' total distribution that must be in the
# zone to qualify for inclusion. E.g. if 0.01, then species with less than 1 % of their
# total distribution in the zone relative to the whole region are excluded
sdm_prop_cutoff <- 0.1 # Use 0 to always include all species

run_id <- readRDS("current_run_id.RDS")
dir.create(paste0("outputs/", run_id, "/matrices/"), recursive = TRUE)
#file.copy("zones.png", paste0("outputs/", run_id, "/zones.png"))

grid <- rast("inputs/HELCOM_WP3_grid.tif")
grid_n <- sum(values(grid, na.rm = T))

# Specify SDMs
spec_df <- read.csv("biotope_species_list.csv")

# Filter based on mobility
spec_df <- filter(spec_df, mobility == "non-mobile")
lf <- paste0("inputs/sdms/", spec_df$scientific_name, ".tif")
spec_n <- spec_df$scientific_name

# Import rasters
r <- rast(lf)
names(r) <- spec_n

# Import zones
zones <- vect(paste0("outputs/", run_id, "/zones/zones.shp"))
zones_n <- zones$zone

# Import species total coverage sums
spec_areas <- spec_df %>% select(scientific_name, presence_sum)
#idx <- which(!is.na(values(grid)))
#rs <- r[idx]
#spec_sums <- colSums(rs, na.rm = TRUE)
#
#spec_areas <- data.frame("scientific_name" = names(spec_sums), "sum" = spec_sums)
#spec_areas <- read.csv(paste0("inputs/sdm_data/", sdm_type, "/sdm_area.csv"))

########## Caclulate Bray-Curtis dissimilarity matrix ##########

for(i in 1:length(zones_n)){
  # Specify zone
  temp <- zones[zones$zone == zones_n[i],]
  tgrid <- terra::mask(grid, temp, touches = TRUE)
  tgrid <- crop(tgrid, ext(temp))
  r_crop <- crop(r, ext(tgrid))
  r_crop <- r_crop*tgrid
  
  
  # Sample
  #n_cells <- ncell_zone
  # Get indices of valid cells
  idx <- which(!is.na(values(tgrid)) & !is.na(values(r_crop[[1]])))
  rs <- r_crop[idx]
  #rs <- terra::extract(r_crop, idx) # Other option
  ncell_zone <- length(idx)
  
  #rs <- terra::spatSample(r_crop, size = n_cells, method = "regular", na.rm = TRUE)
  # Remove species with 0 presences in zone
  spec_nonzero <- colSums(rs)
  spec_nonzero <- spec_nonzero != 0
  rs <- rs[,spec_nonzero]
  
  spec_sums <- colSums(rs)
  spec_props <- spec_sums / ncell_zone
  spec_zone <- data.frame("scientific_name" = names(spec_sums), "zone_sum" = spec_sums)
  spec_zone <- left_join(spec_zone, spec_areas, by = "scientific_name")
  spec_zone$zone_prop <- spec_zone$zone_sum /spec_zone$presence_sum
  
  # Remove species with only a small fraction of their distribution in the subbasin?
  spec_zone$zone_prop <- round(spec_zone$zone_prop, digits = 5)
  write.csv(spec_zone, paste0("outputs/", run_id, "/matrices/", zones_n[i], "_zone_prop.csv"), row.names = FALSE)
  # Remove species?
  spec_select <- spec_zone$scientific_name[spec_zone$zone_prop > sdm_prop_cutoff]
  rs <- rs[,spec_select]
  
  ### Calculate pairwise bray-curtis dissimilarity
  #v <- terra::values(rs, mat = TRUE)
  v  <- as.matrix(rs)
  #v[is.na(v)] <- 0
  m <- t(v)
  
  keep <- colSums(m) > 0
  m <- m[, keep, drop = FALSE]
  
  m_diss <- vegdist(m, method = "bray")
  saveRDS(m_diss, paste0("outputs/", run_id, "/matrices/", zones_n[i], ".rds"))
}



