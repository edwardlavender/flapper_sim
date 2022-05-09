######################################
######################################
#### helper.R

#### This code:
# 1) Defines helper functions used by multiple scripts.

#### Steps preceding this code:
# 1) NA


######################################
######################################
#### Function definitions

#### Scale raster (to a maximum value of one)
scale_raster <- function(x) {
  sm <- raster::cellStats(x, "sum")
  print(paste0("sum = ", sm))
  x <- x/sm
  x <- x/raster::cellStats(x, "max")
  # print(raster::cellStats(x, "max"))
  return(x)
}

#### Plot raster using interpolation
plot_raster_img <-
  function(x,...) terra::plot(terra::rast(x), smooth = TRUE,
                              axes = FALSE, legend = FALSE,...)
plot_raster <-
  function(x,...) raster::plot(x, interpolate = TRUE,...)

#### Add contour to a plot
add_contour <-
  function(x, p = 0.5, ext = NULL){
    if(!is.null(ext)) x <- raster::crop(x, ext)
    quant <- raster::quantile(x, p)
    raster::contour(x, levels = quant, lwd = 0.5, labels = "", add = TRUE)
    return(invisible())
  }


#### End of code.
######################################
######################################
