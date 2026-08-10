### Setup
library(terra)

run_id <- "WP3-Habitats"
saveRDS(run_id, "current_run_id.RDS")
dir.create(paste0("outputs/", run_id, "/zones"), recursive = TRUE)
file.copy("biotope_species_list.csv", paste0("outputs/", run_id, "/biotope_species_list.csv"))

zone_opt <- 3

######## Option 1 - Defining by salinity zones
if(zone_opt == 1){
  n_zones <- 4
  
  salinity <- rast("inputs/predictors/salinity_surface_mean.tif")
  valid_idx <- which(!is.na(values(salinity)))
  sal_vals <- values(salinity)[valid_idx]
  
  breaks <- quantile(sal_vals, probs = seq(0, 1, length.out = n_zones + 1))
  breaks[1] <- breaks[1] - 1e-6
  breaks[n_zones + 1] <- breaks[n_zones + 1] + 1e-6
  clusters <- as.integer(cut(sal_vals, breaks = breaks, labels = FALSE))
  
  zone_raster <- rast(salinity)
  values(zone_raster) <- NA_integer_
  values(zone_raster)[valid_idx] <- clusters
  names(zone_raster) <- paste0("ecoregion_k", n_zones)
  
  subs <- as.polygons(zone_raster)
  names(subs)[1] <- "zone"
  subs$zone <- paste0("Zone", 1:n_zones)
  
  writeVector(subs, "inputs/zones/zones.shp", filetype = "ESRI Shapefile", overwrite = TRUE)
}

######## Option 2 - Defining salinity zones with equal interval
if(zone_opt == 2){
  n_zones <- 6
  
  salinity <- rast("inputs/predictors/salinity_surface_mean.tif")
  valid_idx <- which(!is.na(values(salinity)))
  sal_vals <- values(salinity)[valid_idx]
  
  breaks <- seq(0, 33, length.out = n_zones + 1)
  #breaks <- quantile(sal_vals, probs = seq(0, 1, length.out = n_zones + 1))
  #breaks[1] <- breaks[1] - 1e-6
  #breaks[n_zones + 1] <- breaks[n_zones + 1] + 1e-6
  clusters <- as.integer(cut(sal_vals, breaks = breaks, labels = FALSE))
  
  zone_raster <- rast(salinity)
  values(zone_raster) <- NA_integer_
  values(zone_raster)[valid_idx] <- clusters
  names(zone_raster) <- paste0("ecoregion_k", n_zones)
  
  subs <- as.polygons(zone_raster)
  names(subs)[1] <- "zone"
  subs$zone <- paste0("Zone", 1:n_zones)
  
  writeVector(subs, "inputs/zones/zones.shp", filetype = "ESRI Shapefile", overwrite = TRUE)
}

######## Option 3 - Defining salinity zones manually
# 0-3, 3-6, 6-12, 12-31

if(zone_opt == 3){
  n_zones <- 4
  
  salinity <- rast("inputs/predictors/salinity_surface_mean.tif")
  valid_idx <- which(!is.na(values(salinity)))
  sal_vals <- values(salinity)[valid_idx]
  
  #breaks <- seq(0, 33, length.out = n_zones + 1)
  breaks <- c(0, 3, 6, 12, 31)
  #breaks <- quantile(sal_vals, probs = seq(0, 1, length.out = n_zones + 1))
  #breaks[1] <- breaks[1] - 1e-6
  #breaks[n_zones + 1] <- breaks[n_zones + 1] + 1e-6
  clusters <- as.integer(cut(sal_vals, breaks = breaks, labels = FALSE))
  
  zone_raster <- rast(salinity)
  values(zone_raster) <- NA_integer_
  values(zone_raster)[valid_idx] <- clusters
  names(zone_raster) <- paste0("ecoregion_k", n_zones)
  
  subs <- as.polygons(zone_raster)
  names(subs)[1] <- "zone"
  subs$zone <- paste0("Zone", 1:n_zones)
  
  writeVector(subs, paste0("outputs/", run_id, "/zones/zones.shp"), filetype = "ESRI Shapefile", overwrite = TRUE)
}

########## Plot ##########
n <- nrow(subs)
pal <- hcl.colors(n, palette = "Set 2")

png(paste0("outputs/", run_id, "/zones.png"), width = 1600, height = 1200, res = 150)
plot(subs, col = pal, border = "grey30", main = "Zones")
cents <- centroids(subs)
coords <- crds(cents)
text(coords[, 1], coords[, 2], labels = subs$zone, cex = 0.7, font = 2, col = "black")
dev.off()
