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
  temp$cluster_name <- paste0(zones[i], "_", temp$cluster)
  if(i == 1){df <- temp}
  if(i > 1){df <- rbind(df, temp)}
}

write.csv(df, paste0("outputs/", run_id, "/full_cluster_df.csv"), row.names = FALSE)

# Get unique cluster names
clusters <- unique(df$cluster_name)

# Build species sets per cluster
species_by_cluster <- split(df$species, df$cluster_name)

# Generate all pairwise combinations
pairs <- combn(clusters, 2, simplify = FALSE)

# Calculate overlap for each pair
overlap_df <- do.call(rbind, lapply(pairs, function(p) {
  sp_a <- species_by_cluster[[p[1]]]
  sp_b <- species_by_cluster[[p[2]]]
  
  n_shared    <- length(intersect(sp_a, sp_b))
  n_a         <- length(sp_a)
  n_b         <- length(sp_b)
  n_union     <- length(union(sp_a, sp_b))
  
  data.frame(
    cluster_1        = p[1],
    cluster_2        = p[2],
    n_species_1      = n_a,
    n_species_2      = n_b,
    n_shared         = n_shared,
    pct_of_cluster_1 = round(100 * n_shared / n_a, 1),
    pct_of_cluster_2 = round(100 * n_shared / n_b, 1),
    total_pct      = round(100 * n_shared / n_union, 1),
    stringsAsFactors  = FALSE
  )
}))

# Order by Jaccard similarity (highest overlap first)
overlap_df <- overlap_df %>% arrange(desc(total_pct))

write.csv(overlap_df, paste0("outputs/", run_id, "/cluster_overlap_df.csv"), row.names = FALSE)