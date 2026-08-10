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

# Specify SDMs
spec_df <- read.csv("biotope_species_list.csv")

# Filter based on mobility
spec_df <- filter(spec_df, mobility == "non-mobile")
lf <- paste0("inputs/sdms/", spec_df$scientific_name, ".tif")
spec_n <- spec_df$scientific_name

# Import rasters
r <- rast(lf)
names(r) <- spec_n

k <- length(levels(as.factor(df$cluster_name)))
clust_names <- levels(as.factor(df$cluster_name))

habs <- rast(lf)
names(habs) <- spec_n
habs_all <- app(habs, fun = "sum", na.rm = TRUE)

dir.create(paste0("outputs/", run_id, "/_tifs/individ_cluster_specprop_tif/"), recursive = TRUE)
dir.create(paste0("outputs/", run_id, "/_tifs/individ_cluster_binary_tif/"), recursive = TRUE)
dir.create(paste0("outputs/", run_id, "/_tifs/individ_cluster_sorensen_tif/"), recursive = TRUE)


for(i in 1:k){
  
  com_spec <- df$species[which(df$cluster_name == clust_names[i])]
  habs_com <- habs[[com_spec]]
  n <- nlyr(habs_com)
  
  intersection <- app(habs_com, fun = "sum", na.rm = TRUE)
  union <- habs_all + n - intersection
  
  # Jaccard similarity
  #comr <- intersection / union
  
  #Sorensen–Dice coefficient
  comr <- (2 * intersection) / (habs_all + n)
  
  #-------- Habitat maps by species proportion
  hab_single <- intersection
  hab_single <- hab_single / n
  writeRaster(hab_single, 
              paste0("outputs/", run_id, "/_tifs/individ_cluster_specprop_tif/", clust_names[i], ".tif"),
              overwrite = TRUE)
  
  #-------- Habitat maps by Sorensen-Dice coefficient
  hab_single <- comr
  writeRaster(hab_single, 
              paste0("outputs/", run_id, "/_tifs/individ_cluster_sorensen_tif/", clust_names[i], ".tif"),
              overwrite = TRUE)
  
  #-------- Binary habitat maps classified as Sorensen-Dice coefficient greater than X
  hab_single <- comr
  hab_single <- ifel(hab_single > 0.2, 1, 0)
  writeRaster(hab_single, 
              paste0("outputs/", run_id, "/_tifs/individ_cluster_binary_tif/", clust_names[i], ".tif"),
              overwrite = TRUE)
  
  names(comr) <- clust_names[i]
  if(i == 1){coms <- comr}
  if(i > 1){coms <- c(coms, comr)}
  cat("Iteration", i, "complete - ")
  
}

#plot(coms, background = "white", maxnl = 32)

# Colours
clrs <- qualitative_hcl(k, palette = "Dark 3")
clrs <- clrs[sample(1:length(clrs))] # Randomize colours?

# Winning cluster
com_maj <- terra::app(coms, fun = "which.max", na.rm = FALSE)
#com_maj[is.na(grid)] <- NA
#com_maj <- project(com_maj, "EPSG:4326", method = "near")
#com_maj[com_maj == 0] <- NA
cls <- data.frame(id=1:k, com = clust_names)
levels(com_maj) <- cls
#plot(com_maj, col = clrs, reverse = TRUE)

# Check for cells with poor match
winning_score  <- max(coms)
score_vals <- values(winning_score, na.rm = TRUE)

hist(score_vals, breaks = 50, 
     main = "Distribution of best-cluster Sørensen scores",
     xlab = "Sørensen similarity to winning cluster")
abline(v = c(0.05, 0.10, 0.15, 0.20), col = "red", lty = 2)

com_maj[winning_score < 0.03] <- 0

no_hab <- data.frame(id=0, com = "Habitat undefined")
cls <- rbind(no_hab, cls)
levels(com_maj) <- cls
clrs <- c("grey", clrs)

#temp <- habs_all
#temp[temp == 0] <- 9999


### Export
dir.create(paste0("outputs/", run_id, "/app_data/categorical_all"), recursive = TRUE)
out_dir <- paste0("outputs/", run_id, "/app_data/categorical_all/")

# Export map
png(paste0(out_dir, "/biotopes_categorical_map.png"), 
    width = 22000, height = 18000, res = 2400)
plot(com_maj, col = clrs, maxcell = ncell(com_maj), axes = FALSE, frame.plot = FALSE, reverse = TRUE)
dev.off()

# Species df
species_df <- select(df, cluster_name, species)
colnames(species_df) <- c("habitat_name", "species")
write.csv(species_df, paste0(out_dir, "species_df.csv"), row.names = FALSE)

# Habitat df
habitat_df <- cls
colnames(habitat_df) <- c("habitat_id", "habitat_name")
write.csv(habitat_df, paste0(out_dir, "habitat_df.csv"), row.names = FALSE)

# Habitat raster
names(com_maj) <- "habitat_id"
writeRaster(com_maj, paste0(out_dir, "cluster_majority.tif"), overwrite = TRUE)

