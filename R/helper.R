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

#### pf_dens()
pf_dens <- function(xpf, coord = NULL, plot = TRUE,
                    use_tryCatch = TRUE,
                    ...) {

  # Process raster
  # spatstat assumes planar coordinates
  raster::crs(xpf) <- NA

  # Set up window
  rim  <- maptools::as.im.RasterLayer(xpf)
  rwin <- spatstat.geom::as.owin(rim)

  # Get ppp
  if (is.null(coord)) {
    coord <- raster::rasterToPoints(xpf) |> as.data.frame()
    colnames(coord) <- c("x", "y", "mark")
    marks <- coord[, 3]
  } else {
    # Check coord
    if (inherits(coord, "matrix")) {
      coord <- as.data.frame(coord[, 1:2])
      colnames(coord) <- c("x", "y")
    }
    if (!inherits(coord, "data.frame")) stop("`coord` should a data.frame.")
    if (!all(colnames(coord) %in% c("x", "y"))) stop("`coord` should contain x and y columns.")
    # Drop duplicate coordinates & adjust marks (weights) accordingly
    # ... This assumes that coordinates that are uniquely define in `coord`
    # ... are uniquely defined on xpf (see definition of weights below).
    marks <- rep(1/nrow(coord), nrow(coord))
    nrw_0 <- nrow(coord)
    coord <-
      coord |>
      dplyr::mutate(mark = marks,
                    key = paste(.data$x, .data$y)) |>
      dplyr::group_by(.data$key) |>
      dplyr::mutate(mark = sum(.data$mark)) |>
      dplyr::slice(1L)
    nrw_1 <- nrow(coord)
    stopifnot(all.equal(1, sum(coord$mark)))
    if (nrw_0 != nrw_1) warning("`coords` contains duplicate coordinates that have been processed.")
    # Define marks/weights
    marks <- coord$mark
    rim <- raster::rasterize(x = coord[, c("x", "y")], y = xpf, field = coord$mark)
    rim <- maptools::as.im.RasterLayer(rim)
  }
  rppp <- spatstat.geom::ppp(x = coord$x, y = coord$y, window = rwin, marks = marks)

  # Get intensity (expected number of points PER UNIT AREA)
  dens <- tryCatch(spatstat.explore::density.ppp(rppp, weights = rim, ...),
                   error = function(e) e)
  if (inherits(dens, "error")) {
    if (!use_tryCatch) {
      stop(dens)
    } else {
      warning(paste("\n", paste(dens, collapse = "\n ")), immediate. = TRUE, call. = FALSE)
      return(NULL)
    }
  }
  # Translate intensity into expected number of points PER PIXEL
  dens <- raster::raster(dens)
  dens <- dens * raster::area(dens)
  # Translate expect counts into proportion of points per pixel
  dens <- dens/raster::cellStats(dens, "sum")
  stopifnot(all.equal(1, raster::cellStats(dens, "sum")))

  # Plot
  if (plot) {
    ext <- raster::extent(xpf)
    prettyGraphics::pretty_map(add_rasters = list(x = dens),
                               xlim = ext[1:2], ylim = ext[3:4])
  }

  # Return map
  dens

}


#### End of code.
######################################
######################################
