library(terra)

save_biotopes_poly <- function(rast_stack, 
                               output_file  = "output.png",
                               spat_vector  = NULL,
                               title        = NULL, 
                               res          = 150, 
                               panel_factor = 600,
                               # Raster colour scheme
                               rast_col     = hcl.colors(100, "viridis"),
                               # SpatVector styling
                               vec_col      = "#FF4500",   # border colour
                               vec_fill     = "#FF450030", # fill (transparent)
                               vec_lwd      = 1.2) {
  
  n <- nlyr(rast_stack)
  
  ncols <- ceiling(sqrt(n))
  nrows <- ceiling(n / ncols)
  
  panel_size <- round(panel_factor * (res / 150))
  
  png(output_file, 
      width  = ncols * panel_size, 
      height = nrows * panel_size, 
      res    = res)
  
  if (!is.null(title)) {
    layout(matrix(c(rep(0, ncols),
                    seq_len(nrows * ncols)),
                  nrow  = nrows + 1, 
                  byrow = TRUE),
           heights = c(0.08, rep(1, nrows)))
    par(mar = c(0.5, 0.5, 1.5, 2.5))
    
    plot.new()
    text(0.5, 0.5, title, cex = 1, font = 2, xpd = NA)
  } else {
    par(mfrow = c(nrows, ncols), mar = c(0.5, 0.5, 1.5, 2.5))
  }
  
  # Reproject vector once if provided, to match raster CRS
  if (!is.null(spat_vector)) {
    if (!same.crs(spat_vector, rast_stack)) {
      message("Reprojecting SpatVector to match raster CRS.")
      spat_vector <- project(spat_vector, crs(rast_stack))
    }
  }
  
  for (i in seq_len(n)) {
    plot(rast_stack[[i]], 
         main     = names(rast_stack)[i],
         cex.main = 0.9,
         axes     = FALSE,
         legend   = TRUE,
         col      = rast_col)
    
    if (!is.null(spat_vector)) {
      plot(spat_vector, 
           add    = TRUE, 
           border = vec_col, 
           col    = vec_fill, 
           lwd    = vec_lwd)
    }
  }
  
  if (n < nrows * ncols) {
    for (i in seq_len(nrows * ncols - n)) plot.new()
  }
  
  dev.off()
  message("Saved: ", output_file)
}