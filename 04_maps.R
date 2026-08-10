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
source("_functions/save_biotopes_poly.R")

run_id <- readRDS("current_run_id.RDS")
#sdm_type <- "binary"

if(exists("run_id_master")){run_id <- run_id_master}

grid <- rast("inputs/HELCOM_WP3_grid.tif")
grid_n <- sum(values(grid, na.rm = T))

# Import zones
zones <- vect(paste0("outputs/", run_id, "/zones/zones.shp"))
zones_n <- zones$zone

########## SDMs ###########
# Specify SDMs
spec_df <- read.csv("biotope_species_list.csv")

# Filter based on mobility
spec_df <- filter(spec_df, mobility == "non-mobile")
lf <- paste0("inputs/sdms/", spec_df$scientific_name, ".tif")
spec_n <- spec_df$scientific_name

# Import rasters
r <- rast(lf)
names(r) <- spec_n

########## k details ###########
df_k <- read.csv(paste0("outputs/", run_id, "/df_optimal_k.csv"))

########## Loop by zone ###########
for(z in 1:length(zones_n)){

  zone <- zones_n[z]
  zone_poly <- zones[zones$zone == zones_n[z],]
  out_dir <- paste0("outputs/", run_id, "/dendrograms/", zone)
  if(!file.exists(out_dir)){dir.create(out_dir)}
  k <- df_k$k[df_k$zone == zone]
  df <- read.csv(paste0("outputs/", run_id, "/dendrograms/", zone, "/", zone, "_species_df.csv"))
  clust_names <- rev(c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[1:k])
  
  # Map clusters
  for(i in 1:k){
    habs <- rast(lf)
    names(habs) <- spec_n
    
    com_spec <- df$species[which(df$cluster == clust_names[i])]
    habs_com <- habs[[com_spec]]
    
    n <- nlyr(habs_com)
    comr <- app(habs_com, fun = "sum", na.rm = TRUE)
    #comr <- comr/n # proportions instead of species number
    names(comr) <- clust_names[i]
    
    if(i == 1){coms <- comr}
    if(i > 1){coms <- c(coms, comr)}
  }
  
  while (dev.cur() > 1) dev.off()
  
  save_biotopes_poly(coms, 
                     paste0(out_dir, "/", zone, "_clusters_separate.png"), 
                     spat_vector = zone_poly,
                     res = 300,
                     rast_col = hcl.colors(10, "Turku"),
                     vec_col = "cyan",
                     vec_fill = NULL,  
                     panel_factor = 1200,
                     vec_lwd = 0.1
                     )
  
}
