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

run_id <- readRDS("current_run_id.RDS")

ttest_alpha <- 0.05
dir.create(paste0("outputs/", run_id, "/dendrograms/"), recursive = TRUE)

grid <- rast("inputs/HELCOM_WP3_grid.tif")
grid_n <- sum(values(grid, na.rm = T))

# Import zones
zones <- vect(paste0("outputs/", run_id, "/zones/zones.shp"))
zones_n <- zones$zone

# Setup k data frame
df_k <- data.frame("zone" = zones_n, "k" = NA)

for(i in 1:length(zones_n)){

  ########## Import matrix ###########
  zone <- zones_n[i]
  out_dir <- paste0("outputs/", run_id, "/dendrograms/", zone)
  if(!file.exists(out_dir)){dir.create(out_dir)}
  fp <- paste0("outputs/", run_id, "/matrices/", zone, ".rds")
  m_diss <- readRDS(fp)
  
  ########## Determine optimal k ###########
  # Correct to Euclidian distances (if using Ward method)
  pcoa_result <- cmdscale(m_diss, k = nrow(as.matrix(m_diss)) - 1, eig = TRUE)
  m_diss <- dist(pcoa_result$points)
  
  # If proportion of variance is high in negative eigenvalues, e.g. > 0.05, perform correction
  eig <- pcoa_result$eig
  neg_var <- sum(eig[eig < 0]) / sum(abs(eig))
  
  if(abs(neg_var) > 0.05){
    pcoa_result <- ape::pcoa(m_diss, correction = "cailliez")
    m_diss <- dist(pcoa_result$vectors)
  }
  
  # Determine k
  k <- optimal_k_medoid_ttest(m_diss,
                              k_max_rule = "n/2",
                              k_min = 2,
                              alpha = ttest_alpha) ###!!! NOTE: here the function is called for determining optimal k
  
  k <- k$optimal_k
  
  ########## Create dendrogram ###########
  hc <- hclust(m_diss, method = "ward.D2")
  cl <- cutree(hc, k = k)
  lnk.dend <- as.dendrogram(hc)
  
  clrs <- qualitative_hcl(k, palette = "Dark 3")
  clrs <- clrs[sample(1:length(clrs))] # Randomize colours?
  
  dend <- lnk.dend %>% 
    dendextend::set("branches_k_col", clrs, k = k) %>% 
    dendextend::set("labels_colors", clrs, k = k) %>%
    dendextend::set("labels_cex", 0.7) %>%
    dendextend::set("branches_lwd", 2)
  
  df <- data.frame("species" = labels(dend), "colour" = get_leaves_branches_attr(dend, "col"))
  df$colour <- factor(df$colour, levels = clrs)
  temp <- df %>%
    group_by(colour) %>%
    summarise(colour = first(colour))
  clust_names <- rev(c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))[1:k])
  temp$cluster <- clust_names
  temp$cluster_int <- rev(1:k)
  df <- left_join(df, temp, by = "colour")
  
  # Labels for dendrogram
  cluster_map <- clust_names
  names(cluster_map) <- rev(1:k)
  h <- attr(dend, "height"); if (is.null(h)) h <- max(get_branches_heights(dend))
  
  ### Export dendrogram
  jpeg(filename = paste0(out_dir, "/", zone, "_dendrogram.jpg"), 
       width=5000, height=nrow(df)*90, res=600)
  par(mar=c(2,2,2,10))
  plot(dend, nodePar = list(cex = 0.5, lab.cex = 0.5), horiz = TRUE)
  annotate_clusters(dend, df = df, k = k, clrs = rev(clrs), labels_by_cluster = cluster_map, offset = -0.25)
  dev.off()
  
  ### Export bifurcation t-test plot
  jpeg(filename = paste0(out_dir, "/", zone, "_ttest.jpg"), 
       width=3000, height=2500, res=400)
  par(mar=c(2,2,2,10))
  optimal_k_medoid_ttest(m_diss,
                         k_max_rule = "n/2",
                         k_min = 2)
  dev.off()
  
  # Export / save
  df_k$k[df_k$zone == zone] <- k
  write.csv(df, paste0(out_dir, "/", zone, "_species_df.csv"), row.names = FALSE)

}

# Export k data frame
write.csv(df_k, paste0("outputs/", run_id, "/df_optimal_k.csv"), row.names = FALSE)
