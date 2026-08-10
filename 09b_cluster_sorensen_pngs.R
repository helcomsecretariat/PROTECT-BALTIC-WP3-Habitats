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
library(scico)
library(RColorBrewer)

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
# Sorensen clusters
l <- paste0("outputs/", run_id, "/_tifs/individ_cluster_sorensen_tif/", clust_names, ".tif")
habs <- rast(l)

# Build a land-only data frame once per species group (land doesn't change)
land_df <- as.data.frame(land, xy = TRUE, na.rm = TRUE)
names(land_df)[3] <- "land"
land_df <- land_df[land_df$land == 1, ]

for(i in 1:k){
  r <- habs[[i]]
  dest_fp <- paste0("outputs/", run_id, "/_pngs/individ_cluster_sorensen_png/", clust_names[i], ".png")
  r[land == 1] <- NA
  r[r==0] <- NA
  
  dfr <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(dfr)[3] <- "value"
  
  temp_plot <- ggplot() +
    # Land underneath
    geom_raster(data = land_df, aes(x = x, y = y), fill = "grey80") +
    # Continuous data on top
    geom_raster(data = dfr, aes(x = x, y = y, fill = value)) +
    coord_equal() +
    #scale_fill_scico(name = "Relative probability\nabove presence\nthreshold",
    #                 na.value = "transparent",
    #                 palette = "berlin",
    #                 direction = 1) +
    scale_fill_distiller(name = "Sorensen-Dice\ncoefficient",
                         na.value = "transparent",
                         palette = "YlOrRd",
                         direction = 1) +
    theme(panel.background = element_rect(fill = "white"), 
          plot.background = element_rect(fill = "white", color = NA),
          legend.background = element_rect(fill = "white", color = NA),
          legend.key = element_rect(fill = "white", color = NA),
          panel.grid = element_blank(),
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  png(filename = dest_fp,
      width = 9.5,
      height = 8.5,
      units = "in",
      res = 600)
  print(temp_plot)
  dev.off()
  
  cat("Iteration", i, "complete - ")
}