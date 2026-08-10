# medoid_ttest.R
#
# Optimal k via Bifurcation Paired T-test (medoid method)
#
# Method (Pang et al. 2023, Appendix S1, Fig. S3):
#   Build a ward.D2 dendrogram. At each step k -> k+1, the dendrogram
#   prescribes exactly one cluster to split. Work directly in dissimilarity
#   space: the cluster medoid is the object with the lowest total dissimilarity
#   to all other members. Measure each object's dissimilarity to its cluster
#   medoid before and after the split, then test for a significant decrease
#   using a one-sided paired t-test.
#   Optimal k = the first k at which the prescribed split is non-significant
#   (i.e., no further bifurcation is needed).
#
# Usage:
#   mat <- readRDS("data/matrices/Arkona_Basin.rds")
#   if (!inherits(mat, "dist")) mat <- as.dist(mat)
#   result <- optimal_k_medoid_ttest(mat)
#   result$optimal_k

library(ggplot2)


optimal_k_medoid_ttest <- function(
    dist_matrix,          # Bray-Curtis dissimilarity matrix, class 'dist'
    k_range     = NULL,   # explicit integer vector, overrides k_min/k_max_rule if supplied
    k_min       = 2,      # minimum k to test (>= 2)
    k_max_rule  = "n/2",  # upper bound rule: "sqrt(n)","n/10","n/6","n/4","n/3","n/2","n-1"
    alpha       = 0.05,   # significance threshold for the paired t-test
    plot        = TRUE,   # produce a criterion plot?
    matrix_name = NULL    # optional label for the plot title
) {
  
  # ---- Input checks --------------------------------------------------------
  if (!inherits(dist_matrix, "dist"))
    stop("dist_matrix must be a 'dist' object.")
  
  n <- attr(dist_matrix, "Size")
  
  # Resolve k_range: explicit vector takes priority; otherwise apply k_min/k_max_rule.
  if (is.null(k_range)) {
    if (!is.numeric(k_min) || k_min < 2)
      stop("k_min must be an integer >= 2.")
    k_max <- switch(k_max_rule,
                    "sqrt(n)" = floor(sqrt(n)),
                    "n/10"    = floor(n / 10),
                    "n/6"     = floor(n /  6),
                    "n/4"     = floor(n /  4),
                    "n/3"     = floor(n /  3),
                    "n/2"     = floor(n /  2),
                    "n-1"     = n - 1L,
                    stop('k_max_rule must be one of: "sqrt(n)", "n/10", "n/6", "n/4", "n/3", "n/2", "n-1"')
    )
    k_range <- as.integer(k_min):max(as.integer(k_min), k_max)
  }
  
  if (length(k_range) < 2)
    stop("k_range must span at least 2 values.")
  
  
  # ---- Build dendrogram and full dissimilarity matrix ----------------------
  # The dendrogram prescribes which cluster splits at each step (Fig. S3a).
  hc <- hclust(dist_matrix, method = "ward.D2")
  
  # Full n x n matrix for direct distance look-up (avoids repeated as.matrix calls).
  dm <- as.matrix(dist_matrix)
  
  
  # ---- Bifurcation paired t-test across k_range ----------------------------
  t_stats     <- numeric(length(k_range))   # t-statistic at each k (for plot)
  significant <- logical(length(k_range))   # TRUE if split at k is significant
  
  for (i in seq_along(k_range)) {
    
    k         <- k_range[i]
    labels_k  <- cutree(hc, k = k)
    labels_k1 <- cutree(hc, k = k + 1)   # one cluster in labels_k will split
    
    # Identify the cluster that splits: objects in it map to 2 groups at k+1.
    split_cluster <- NA
    for (cl in unique(labels_k)) {
      if (length(unique(labels_k1[labels_k == cl])) >= 2) {
        split_cluster <- cl
        break
      }
    }
    
    if (is.na(split_cluster)) next
    
    idx <- which(labels_k == split_cluster)
    
    if (length(idx) < 4) next   # need at least 4 objects to form two sub-groups
    
    sub_dm <- dm[idx, idx, drop = FALSE]
    
    # Medoid = the object with the smallest total dissimilarity to all others
    # in the cluster (i.e., the most central object).
    medoid_before <- idx[which.min(colSums(sub_dm))]   # global index
    
    # Distance of each object to the pre-split medoid (before).
    d_before <- dm[idx, medoid_before]
    
    # Distance of each object to its new sub-cluster medoid (after).
    sub_labels <- labels_k1[idx]
    d_after    <- numeric(length(idx))
    
    for (sub_cl in unique(sub_labels)) {
      members      <- which(sub_labels == sub_cl)
      sub2         <- sub_dm[members, members, drop = FALSE]
      medoid_after <- idx[members[which.min(colSums(sub2))]]   # global index
      d_after[members] <- dm[idx[members], medoid_after]
    }
    
    # One-sided paired t-test: H1 = d_before > d_after (significant WCV decrease).
    test <- tryCatch(
      t.test(d_before, d_after, paired = TRUE, alternative = "greater"),
      error = function(e) list(statistic = 0, p.value = 1)
    )
    
    t_stats[i]     <- as.numeric(test$statistic)
    significant[i] <- test$p.value < alpha
  }
  
  
  # ---- Select optimal k ----------------------------------------------------
  # Optimal k = first k where the prescribed split is non-significant,
  # meaning objects are already optimally clustered (Appendix S1, Fig. S3c).
  non_sig   <- which(!significant)
  optimal_k <- if (length(non_sig) == 0) max(k_range) else k_range[min(non_sig)]
  
  
  # ---- Plot ----------------------------------------------------------------
  if (plot) {
    title_str <- paste0(
      "Medoid: Bifurcation Paired T-test",
      if (!is.null(matrix_name)) paste0(" \u2014 ", matrix_name) else ""
    )
    
    df <- data.frame(k = k_range, t_stat = t_stats,
                     sig = factor(significant, levels = c(TRUE, FALSE),
                                  labels = c("Significant", "Non-significant")))
    
    p <- ggplot(df, aes(x = k, y = t_stat, colour = sig)) +
      geom_point(size = 2.5) +
      geom_line(aes(group = 1), colour = "grey60", linewidth = 0.5) +
      geom_hline(yintercept = qt(1 - alpha, df = length(k_range) - 1),
                 linetype = "dotted", colour = "grey40") +
      geom_vline(xintercept = optimal_k, linetype = "dashed",
                 colour = "firebrick", linewidth = 0.8) +
      annotate("text", x = optimal_k, y = max(t_stats),
               label = paste0("k = ", optimal_k),
               hjust = -0.15, vjust = 1, colour = "firebrick", size = 3.5) +
      scale_colour_manual(values = c("Significant"     = "steelblue",
                                     "Non-significant" = "coral"),
                          name = paste0("Split (alpha=", alpha, ")")) +
      scale_x_continuous(breaks = pretty(k_range, n = 8)) +
      labs(title = title_str,
           x     = "Number of clusters (k)",
           y     = "Paired t-statistic") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom")
    
    print(p)
  }
  
  
  # ---- Return --------------------------------------------------------------
  list(
    optimal_k   = optimal_k,
    t_stats     = t_stats,
    significant = significant,
    k_range     = k_range
  )
}
