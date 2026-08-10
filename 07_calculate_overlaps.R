########## Setup ##########
library(terra)
library(mdendro)
library(dendextend)
library(cluster)
library(dplyr)
library(tidyr)
library(geodata)
library(randomcoloR)
library(colorspace)
library(dynamicTreeCut)
library(fpc)
library(ape)
library(vegan)

#run_id <- "run3_4manSalinityZones_10perc_binConf" # Specify run id
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

k <- length(levels(as.factor(df$cluster_name)))
clust_names <- levels(as.factor(df$cluster_name))

#---------------------------------------------------#
# Binary habitat stack
l <- paste0("outputs/", run_id, "/_tifs/individ_cluster_binary_tif/", clust_names, ".tif")
habs_bin <- rast(l)
names(habs_bin) <- clust_names

#-------------------------------------------------------------------------
### Calculate % of substrate per habitat
l <- list.files("inputs/substrate", full.names = TRUE)
substrate <- rast(l)
names(substrate) <- c("Coarse", "Hard", "Sand", "Soft")

for(i in 1:nlyr(habs_bin)){
  
  r <- habs_bin[[i]]
  rsub <- substrate*r
  
  sub_vals <- values(rsub, mat = TRUE, na.rm = TRUE)
  sub_vals <- colSums(sub_vals)
  sub_vals <- (sub_vals/sum(sub_vals))*100
  
  tdf <- t(as.data.frame(sub_vals))
  
  if(i == 1){sub_df <- tdf}
  if(i > 1){sub_df <- rbind(sub_df, tdf)}
  
  cat("Iteration", i, "complete - ")
}

sub_df <- as.data.frame(sub_df)
tdf <- data.frame("cluster_names" = names(habs_bin))
sub_df <- cbind(tdf, sub_df)
row.names(sub_df) <- NULL

write.csv(sub_df, paste0("outputs/", run_id, "/substrate_overlap.csv"), row.names = FALSE)

#---------------------------------------------------#
# Calculate % of habitat per mobile species
# Specify SDMs
spec_df <- read.csv("biotope_species_list.csv")

# Filter based on mobility
spec_df <- filter(spec_df, mobility == "mobile")
lf <- paste0("inputs/sdms/", spec_df$scientific_name, ".tif")
spec_n <- spec_df$scientific_name

# Import rasters
r <- rast(lf)
names(r) <- spec_n
mob_species <- r

# Cells × layers matrices
sp_mat  <- values(mob_species)
hab_mat <- values(habs_bin)

# Treat NA as 0 (absence/no habitat)
sp_mat[is.na(sp_mat)]   <- 0
hab_mat[is.na(hab_mat)] <- 0

# Overlap counts: species (rows) × habitats (cols)
overlap <- crossprod(sp_mat, hab_mat)   # = t(sp_mat) %*% hab_mat

# Per-habitat denominator 
# (This gives percent of each habitat covered by species (e.g. if 100% then 100% of the habitat is covered by species X))
hab_totals <- colSums(hab_mat)
pct <- 100 * sweep(overlap, 2, hab_totals, "/")
pct[is.nan(pct)] <- 0   # species with no presence cells

mob_species_df <- pct
mob_species_df <- round(mob_species_df, digits = 2)
mob_species_df <- as.data.frame(mob_species_df)

s <- data.frame("scientific_name" = row.names(mob_species_df))
sg <- select(spec_df, scientific_name, species_group)
s <- left_join(s, sg, by = "scientific_name")
mob_species_df <- cbind(s, mob_species_df)
write.csv(mob_species_df, paste0("outputs/", run_id, "/mobile_overlap_habitat.csv"), row.names = FALSE)

#---------------------------------------------------#
# Calculate % of habitat per non-mobile species
# Specify SDMs
spec_df <- read.csv("biotope_species_list.csv")

# Filter based on mobility
spec_df <- filter(spec_df, mobility == "non-mobile")
lf <- paste0("inputs/sdms/", spec_df$scientific_name, ".tif")
spec_n <- spec_df$scientific_name

# Import rasters
r <- rast(lf)
names(r) <- spec_n
mob_species <- r

# Cells × layers matrices
sp_mat  <- values(mob_species)
hab_mat <- values(habs_bin)

# Treat NA as 0 (absence/no habitat)
sp_mat[is.na(sp_mat)]   <- 0
hab_mat[is.na(hab_mat)] <- 0

# Overlap counts: species (rows) × habitats (cols)
overlap <- crossprod(sp_mat, hab_mat)   # = t(sp_mat) %*% hab_mat

# Per-habitat denominator 
# (This gives percent of each habitat covered by species (e.g. if 100% then 100% of the habitat is covered by species X))
hab_totals <- colSums(hab_mat)
pct <- 100 * sweep(overlap, 2, hab_totals, "/")
pct[is.nan(pct)] <- 0   # species with no presence cells

mob_species_df <- pct
mob_species_df <- round(mob_species_df, digits = 2)
mob_species_df <- as.data.frame(mob_species_df)

s <- data.frame("scientific_name" = row.names(mob_species_df))
sg <- select(spec_df, scientific_name, species_group)
s <- left_join(s, sg, by = "scientific_name")
mob_species_df <- cbind(s, mob_species_df)
write.csv(mob_species_df, paste0("outputs/", run_id, "/non_mobile_overlap_habitat.csv"), row.names = FALSE)

#---------------------------------------------------#
# Calculate % of habitat in different depth ranges
depth <- rast("inputs/predictors/depth.tif")
#breaks <- seq(0, minmax(depth)[2], by = 30)
breaks <- c(0, 10, 30, 60, 100, 200, 300)
bin_labels <- paste0(head(breaks, -1), "-", tail(breaks, -1), " m")

breaks <- c(breaks, minmax(depth)[2])
bin_labels <- c(bin_labels, "> 300 m")

# Extract aligned values
depth_vals <- values(depth, mat = FALSE)
hab_mat    <- values(habs_bin)
hab_mat[is.na(hab_mat)] <- 0

# Classify depth into bins
depth_bin <- cut(depth_vals, breaks = breaks,
                 include.lowest = TRUE, right = FALSE,
                 labels = bin_labels)

# Drop cells with NA or out-of-range depth
keep      <- !is.na(depth_bin)
depth_bin <- depth_bin[keep]
hab_mat   <- hab_mat[keep, , drop = FALSE]

# Count bin membership per habitat
counts <- vapply(seq_len(ncol(hab_mat)), function(j) {
  tabulate(as.integer(depth_bin[hab_mat[, j] == 1]),
           nbins = length(bin_labels))
}, integer(length(bin_labels)))

rownames(counts) <- bin_labels
colnames(counts) <- colnames(hab_mat)

# Percentages per habitat (columns sum to 100)
hab_totals <- colSums(counts)
pct <- 100 * sweep(counts, 2, hab_totals, "/")
pct[is.nan(pct)] <- 0

depth_df <- pct
depth_df  <- round(depth_df , digits = 2)
depth_df  <- as.data.frame(depth_df)

s <- data.frame("depth_range" = row.names(depth_df))
depth_df  <- cbind(s, depth_df )
write.csv(depth_df , paste0("outputs/", run_id, "/depth_overlap_habitat.csv"), row.names = FALSE)

