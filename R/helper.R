######################################
######################################
#### helper.R

#### This code:
# 1) Defines helper functions used by multiple scripts.

#### Steps preceding this code:
# 1) NA


######################################
######################################
#### Register suggested dependencies with renv

# Flag commonmark/(r)markdown packages (for README documentation)
if(!requireNamespace("commonmark", quietly = TRUE))
  renv::install("commonmark")
if(!requireNamespace("markdown", quietly = TRUE))
  renv::install("markdown")
if(!requireNamespace("rmarkdown", quietly = TRUE))
  renv::install("rmarkdown")


######################################
######################################
#### Function definitions

#### Read raster (if exists)
read_raster_if_exists <- function(file, read_all = TRUE){
  if(file.exists(file)){
    raster::readAll(raster::raster(file))
  } else NULL
}

#### Maximum raster value
cellStatsMx <- function(x){
  if(!is.null(x)){
    raster::cellStats(x, "max")
  } else NULL
}

#### Scale raster (to a maximum value of one)
scale_raster <- function(x, scale){
  x/scale
}

#### White out rasters
white_out <- function(x){
  if(is.null(x)) return(NULL)
  x[x == 0] <- NA
  return(x)
}

#### Plot raster using interpolation
plot_raster_img <-
  function(x,...) terra::plot(terra::rast(x), smooth = TRUE,
                              axes = FALSE, legend = FALSE,...)
plot_raster <-
  function(x,...) raster::plot(x, interpolate = TRUE,...)

#### Add volume (home range) contour to a plot
add_contour <-
  function(x, p = 0.5, ext = NULL, lwd = 0.5,...){
    if(is.null(x)) return(NULL)
    if(!is.null(ext)) x <- raster::crop(x, ext)
    x <- spatialEco::raster.vol(terra::rast(x), p = p, sample = FALSE)
    x <- raster::raster(x)
    raster::contour(x, nlevels = 1, drawlabels = FALSE, add = TRUE, lwd = lwd,...)
    return(invisible())
  }

#### Count the number of decimal places (ndp)
ndp <- function(x){
  xc <- format(x, scientific = FALSE, trim = TRUE)
  dp <- stringr::str_split_fixed(xc, "[.]", 2)[, 2]
  nchar(dp)
}


#### End of code.
######################################
######################################
