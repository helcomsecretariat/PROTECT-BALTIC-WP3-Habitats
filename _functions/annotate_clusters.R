annotate_clusters <- function(dend, df, k, clrs, labels_by_cluster = NULL, offset = 0.04, cex = 0.9, ...) {
  # cut into k clusters
  cl <- df$cluster_int # named by leaf labels
  names(cl) <- df$species
  labs_in_order <- labels(dend)           # plotting order
  leaf_pos <- setNames(seq_along(labs_in_order), labs_in_order)  # 1..n
  
  # midpoint for each cluster = middle of the range of its leaf positions
  midpos <- tapply(names(cl), cl, function(members) {
    p <- leaf_pos[members]
    (min(p) + max(p)) / 2
  })
  
  #ord <- order(as.numeric(midpos), decreasing = FALSE)
  
  # nice labels: use provided names (e.g. A/B/C) or fallback to numeric cluster IDs
  labels_vec <- labels_by_cluster[names(midpos)]
  #midpos <- as.numeric(midpos)[ord]
  #labels_vec <- labels_vec[ord]
  
  # figure out where to place the text, slightly outside the tree
  h <- attr(dend, "height")
  if (is.null(h)) h <- max(get_branches_heights(dend))
  pad <- h * offset
  
  usr <- par("usr")
  x_target <- usr[2] - diff(usr[1:2]) * offset
  text(x = x_target, 
       y = as.numeric(midpos), 
       labels = labels_vec,
       col = clrs, font = 2,
       xpd = NA, pos = 4, cex = cex, ...)
}