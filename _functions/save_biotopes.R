library(terra)

save_biotopes <- function(rast_stack, output_file = "output.png", 
                                  title = NULL, res = 150, panel_factor = 600) {
  n <- nlyr(rast_stack)
  
  ncols <- ceiling(sqrt(n))
  nrows <- ceiling(n / ncols)
  
  # Scale panel size proportionally to res so margins never overwhelm the plot
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
  
  for (i in seq_len(n)) {
    plot(rast_stack[[i]], 
         main     = names(rast_stack)[i],
         cex.main = 0.9,
         axes     = FALSE,
         legend   = TRUE)
  }
  
  if (n < nrows * ncols) {
    for (i in seq_len(nrows * ncols - n)) plot.new()
  }
  
  dev.off()
  message("Saved: ", output_file)
}