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
library(ggplot2)
library(Polychrome)

run_id <- readRDS("current_run_id.RDS")
land <- rast("inputs/WP3_land_mask.tif")

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

dir.create(paste0("outputs/", run_id, "/_pngs/individ_cluster_specprop_png/"), recursive = TRUE)
dir.create(paste0("outputs/", run_id, "/_pngs/individ_cluster_binary_png/"), recursive = TRUE)
dir.create(paste0("outputs/", run_id, "/_pngs/individ_cluster_sorensen_png/"), recursive = TRUE)

#---------------------------------------------------#
# Binary clusters
l <- paste0("outputs/", run_id, "/_tifs/individ_cluster_binary_tif/", clust_names, ".tif")
habs <- rast(l)

for(i in 1:k){
  r <- habs[[i]] + 1 # Add 1 to specify land as zero
  dest_fp <- paste0("outputs/", run_id, "/_pngs/individ_cluster_binary_png/", clust_names[i], ".png")
  r[land == 1] <- 0
  
  
  # Export plot
  dfr <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  labs <- c("Land", "Low similarity", "Moderate to high similarity")
  cols <- c('grey80', 'white', '#d7301f')
  names(cols) <- labs
  names(dfr)[3] <- "class"
  
  dfr$category <- factor(dfr$class, levels = 0:2, labels = labs)
  
  temp_plot <- ggplot(dfr, aes(x, y, fill = category)) +
    geom_raster() +
    coord_equal() +
    scale_fill_manual(
      values = cols,
      name = "",
      breaks = labs[2:5]  # exclude "Land" from the legend
    ) +
    theme(panel.background = element_rect(fill = "grey90"), panel.grid = element_blank(),
          axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())
  
  png(filename = dest_fp,
      width = 9.5,
      height = 8.5,
      units = "in",
      res = 600)
  print(temp_plot)
  dev.off()
  
  cat("Iteration", i, "complete - ")
}

