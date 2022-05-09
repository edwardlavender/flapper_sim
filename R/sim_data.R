######################################
######################################
#### sim_data.R

#### This code
# Simulates arrays, a movement paths and movement
# ... for evaluating the performance of algorithms
# ... for inferring patterns of space use.

#### Steps preceding this code
# NA.


######################################
######################################
#### Set up

#### Wipe workspace
rm(list = ls())

#### Essential packages
library(prettyGraphics)
library(flapper)
source("./R/helper.R")

#### Global parameters
seed <- 20210828
proj_wgs84 <- sp::CRS("+init=epsg:4326")


######################################
######################################
#### Simulate arrays

#### Define general area
# We will consider a small area, without coastline, for speed.
gebco <- flapper::dat_gebco
raster::plot(gebco)
# centre <- locator()
centre <- sp::SpatialPoints(matrix(c(709411.6, 6258554), nrow = 1),
                            proj4string = raster::crs(gebco))
area <- rgeos::gBuffer(centre, width = 2000, capStyle = "SQUARE") # width is radius, diameter is width *2, etc.
raster::crs(area) <- raster::crs(gebco)
raster::lines(area)
area_size <- raster::area(area) # m2
area_size_km <- area_size * 1/(1000^2)
area_size_km

# crop_from_click(gebco)
ext <- raster::extent(area)
gebco <- raster::crop(gebco, ext)
raster::plot(gebco)

#### Define study area grid
# It is essential that this grid is used for simulation and modelling
blank <- raster::raster(raster::extent(gebco), res = c(75, 75), crs = raster::crs(gebco))
grid <- raster::resample(gebco, blank)
rm(gebco, blank)
raster::plot(grid)

#### Define array properties
n_receivers <- c(25, 15, 5)
arrangement <- c("regular", "random", "clustered")
clustering  <- c(0.4, 0.6)
# seed <- seed
dat_sim_array_info <-
  expand.grid(arrangement = arrangement,
              n_receivers = n_receivers,
              clustering = clustering,
              seed = seed,
              stringsAsFactors = FALSE) %>%
  dplyr::arrange(arrangement, n_receivers)
dat_sim_array_info$clustering[dat_sim_array_info$arrangement != "clustered"] <- NA
dat_sim_array_info <- dat_sim_array_info[!duplicated(dat_sim_array_info), ]
dat_sim_array_info$n_clusters <- ceiling(dat_sim_array_info$n_receivers * dat_sim_array_info$clustering)
dat_sim_array_info$name <- paste0("dat_sim_array_",
                                  dat_sim_array_info$arrangement, "_",
                                  dat_sim_array_info$n_receivers, "_",
                                  dat_sim_array_info$clustering, ".rds")
dat_sim_array_info <- dat_sim_array_info[, c("name", "arrangement", "n_receivers", "clustering", "n_clusters", "seed")]
rownames(dat_sim_array_info) <- 1:nrow(dat_sim_array_info)
dat_sim_array_info

#### Simulate arrays from properties
simulate_array <- FALSE
if(simulate_array){

  pp <- par(mfrow = prettyGraphics::par_mf(nrow(dat_sim_array_info)))
  dat_sim_arrays <-
    lapply(split(dat_sim_array_info, 1:nrow(dat_sim_array_info)), function(array_info){

      #### Simulate array
      # array_info <- dat_sim_array_info[12, , drop = FALSE]
      print(array_info)
      dat_sim_array <- sim_array(boundaries = ext,
                                 n_receivers = array_info$n_receivers,
                                 arrangement = array_info$arrangement,
                                 nclusters = array_info$n_receivers * array_info$clustering,
                                 seed = array_info$seed,
                                 verbose = FALSE
      )
      # dist_btw_clicks(lonlat = FALSE)

      #### Re-express coordinates on grid
      dat_sim_array$array$xy <- raster::xyFromCell(grid, raster::cellFromXY(grid, sp::coordinates(dat_sim_array$array$xy)))
      dat_sim_array$array$xy <- dat_sim_array$array$xy[complete.cases(dat_sim_array$array$xy), ]
      dat_sim_array$array$xy <- sp::SpatialPoints(dat_sim_array$array$xy, raster::crs(grid))

      #### Save
      prompt <- FALSE
      if(prompt) readline("Press [Enter] to continue...")
      return(dat_sim_array)
    })
  par(pp)
  saveRDS(dat_sim_arrays, "./data/arrays/dat_sim_arrays.rds")
} else {
  dat_sim_arrays <- readRDS("./data/arrays/dat_sim_arrays.rds")
}

#### Visualise arrays with detection centroids given some detection range
pp <- par(mfrow = prettyGraphics::par_mf(nrow(dat_sim_array_info)))
lapply(dat_sim_arrays, function(array){
  # array <- dat_sim_arrays[[1]]
  centroids <- get_detection_centroids(xy = array$array$xy, detection_range = 300,
                          boundaries = area,
                          byid = TRUE, plot = FALSE)
  raster::plot(centroids, xlim = ext[1:2], ylim = ext[3:4])
  raster::lines(ext, col = "royalblue")
}) %>% invisible()
par(pp)

#### Define study duration
# We will keep this as small as possible for simulation speed.
study_period <- seq(as.POSIXct("2021-01-01 12:00:00", tz = "UTC"),
                    length.out = 500, by = "2 mins")
range(study_period)
length(study_period)

#### Define moorings dataframes
dat_sim_moorings <- lapply(dat_sim_arrays, function(array){
  moorings <- data.frame(receiver_id = 1:length(array$array$xy),
                         receiver_start_date = as.Date(min(study_period)),
                         receiver_end_date = as.Date(max(study_period))
  )
  return(moorings)
})

#### Save tidy table of array properties
dat_sim_array_info_tidy <- dat_sim_array_info
# Add area coverage assuming a detection radius
dat_sim_array_info_tidy$cov_m2 <- NA
dat_sim_array_info_tidy$cov_pc  <- NA
for(i in 1:length(dat_sim_arrays)){
  array <- dat_sim_arrays[[i]]
  dat_sim_array_info_tidy$n_receivers[i] <- length(array$array$xy)
  dat_sim_array_info_tidy$cov_m2[i] <-
    get_detection_area_sum(array$array$xy, detection_range = 300, boundaries = ext)
  dat_sim_array_info_tidy$cov_pc[i] <- (dat_sim_array_info_tidy$cov_m2[i]/area_size_km) * 100
}
# Tidy existing columns
dat_sim_array_info_tidy$arrangement <- factor(dat_sim_array_info_tidy$arrangement,
                                              levels = c("regular", "random", "clustered"),
                                              labels = c("Regular", "Random", "Clustered"))
dat_sim_array_info_tidy <-
  dat_sim_array_info_tidy %>%
  dplyr::arrange(arrangement, n_receivers, clustering, n_clusters) %>%
  dplyr::mutate(ID = 1:dplyr::n()) %>%
  dplyr::select(ID = ID,
                Arrangement = arrangement,
                `Receivers (N)` = n_receivers,
                `Clusters (Pr)` = clustering,
                `Clusters (N)` = n_clusters,
                `Coverage(m2)` = cov_m2,
                `Coverage(%)` = cov_pc)
# Round area coverage estimates
dat_sim_array_info_tidy$`Coverage(m2)` <-
  prettyGraphics::add_lagging_point_zero(
    round(dat_sim_array_info_tidy$`Coverage(m2)`, 2), 2)
dat_sim_array_info_tidy$`Coverage(%)` <-
  prettyGraphics::add_lagging_point_zero(
  round(dat_sim_array_info_tidy$`Coverage(%)`, 2), 2)
# Save tidy table
ncol(dat_sim_array_info_tidy)
write.table(dat_sim_array_info_tidy, "./fig/dat_sim_array_info_tidy.txt", na = "",
            quote = FALSE, row.names = FALSE)


######################################
######################################
#### Simulate movement

#### Define a list of model parameters
dat_sim_mvt_info <- list("1" = list())

#### Path sim (1)

## Define movement model
mob_sim <- 500
calc_mpr <- function(distance,...) {
  pr <- stats::plogis(7.5 + distance * -0.025)
  pr[distance > mob_sim] <- 0
  return(pr)
}
plot(1:mob_sim, calc_mpr(1:mob_sim), ylim = c(0, 1))

## Simulate step lengths from movement model
steps <- data.frame(distance = seq(0, mob_sim, length.out = 1e4))
steps$pr <- calc_mpr(steps$distance)
sim_step_every_2_mins <- function(...,data = steps, size = 1) {
  sample(x = data$distance, size = size, prob = data$pr)
}
prettyGraphics::pretty_hist(sim_step_every_2_mins(size = 1e3))

## Add details to list
dat_sim_mvt_info[[1]]$sim_step_every_2_mins <- sim_step_every_2_mins
dat_sim_mvt_info[[1]]$seed <- seed
dat_sim_mvt_info[[1]]$name <- "dat_sim_path_1.rds"

simulate_movement <- FALSE
if(simulate_movement){

  #### Simulate paths
  dat_sim_paths <- lapply(dat_sim_mvt_info, function(mvt_info){

    ## Get parameters
    # mvt_info <- dat_sim_mvt_info[[1]]
    seed <- mvt_info$seed
    sim_step_every_2_mins <- mvt_info$sim_step_every_2_mins

    ## Simulate starting location
    set.seed(seed)
    # p_1 <- sp::coordinates(sp::spsample(area, n = 1, type = "random"))
    p_1 <- sp::coordinates(rgeos::gCentroid(area))
    raster::plot(area)
    points(p_1, col = "red")

    ## Simulate movement in area
    dat_sim_path <- sim_path_sa(n = length(study_period),
                                p_1 = p_1,
                                area = area,
                                sim_step = sim_step_every_2_mins,
                                seed = seed)

    #### Re-express sampled locations on the grid
    dat_sim_path$xy_mat_on_grid <- raster::xyFromCell(grid,
                                                      raster::cellFromXY(grid, dat_sim_path$xy_mat))

    ## Return path
    prompt <- FALSE
    if(prompt) readline("Press [Enter] to continue...")
    return(dat_sim_path)
  })
  saveRDS(dat_sim_paths, "./data/paths/dat_sim_paths.rds")

} else {
  dat_sim_paths <- readRDS("./data/paths/dat_sim_paths.rds")
}

#### Check distances along the simulated path
lapply(dat_sim_paths, function(dat_sim_path){
  ## Original path
  print(utils.add::basic_stats(
    raster::pointDistance(
      dat_sim_path$xy_mat[1:(nrow(dat_sim_path$xy_mat)- 1), ],
      dat_sim_path$xy_mat[2:nrow(dat_sim_path$xy_mat), ],
      lonlat = FALSE)
  ))
  ## Path re-expressed on grid
  dat_sim_path$dist_on_grid <-   rep(NA, nrow(dat_sim_path$xy_mat_on_grid))
  dat_sim_path$dist_on_grid[2:length(dat_sim_path$dist_on_grid)] <-
    raster::pointDistance(
      dat_sim_path$xy_mat_on_grid[1:(nrow(dat_sim_path$xy_mat_on_grid)- 1), ],
      dat_sim_path$xy_mat_on_grid[2:nrow(dat_sim_path$xy_mat_on_grid), ],
      lonlat = FALSE)
  print(utils.add::basic_stats(dat_sim_path$dist_on_grid, na.rm = TRUE))
  invisible()
}) %>% invisible()



######################################
######################################
#### Simulate observed movement time series

######################################
#### Depth time series

#### Generate archival time series
simulate_archival <- FALSE
if(simulate_archival){
  dat_sim_archival_by_path <- lapply(1:length(dat_sim_paths), function(path_id){
    # Get path
    dat_sim_path <- dat_sim_paths[[path_id]]
    # Get depth corresponding to locations
    dat_sim_archival <- data.frame(indivdual_id = 1,
                                   timestamp = study_period,
                                   depth_sim = raster::extract(grid, dat_sim_path$xy_mat_on_grid) # on grid
                                   )
    stopifnot(all(!is.na(dat_sim_archival$depth_sim)))
    # Simulate depths using a depth error model
    dat_sim_archival$depth <- runif(1:nrow(dat_sim_archival),
                                    dat_sim_archival$depth_sim - 5,
                                    dat_sim_archival$depth_sim + 5)
    return(dat_sim_archival)
  })
  saveRDS(dat_sim_archival_by_path, "./data/movement/dat_sim_archival_by_path.rds")

} else {
  dat_sim_archival_by_path <-  readRDS("./data/movement/dat_sim_archival_by_path.rds")
}


######################################
#### Acoustic time series

#### Define detection probability function based on distance
detection_range <- 300
calc_dpr <-
  function(x){
    ifelse(x <= detection_range, stats::plogis(3 + -0.03 * x), 0)
  }
plot(1:detection_range, calc_dpr(1:detection_range), type = "l")

#### Simulate acoustic time series for each path and array design
simulate_detections  <- FALSE
if(simulate_detections){
  dat_sim_detections_by_path <- lapply(1:length(dat_sim_paths), function(path_id){

    #### Extract path
    # path_id <- 1L
    dat_sim_path <- dat_sim_paths[[path_id]]

    #### For the selected path, simulate detections for each array design
    dat_sim_detections_by_array <-
      lapply(dat_sim_arrays, function(array){

        ## Simulate detections at receivers (across grid)
        # array <- dat_sim_arrays[[1]]
        dat_sim_detections <- sim_detections(path = dat_sim_path$xy_mat_on_grid,   # expressed on grid
                                             xy = sp::coordinates(array$array$xy), # expressed on grid
                                             calc_detection_pr = calc_dpr,
                                             seed = seed,
                                             plot = FALSE,
                                             verbose = FALSE)
        rownames(dat_sim_detections$det_mat) <- as.character(study_period)
        colnames(dat_sim_detections$det_mat) <- as.character(1:ncol(dat_sim_detections$det_mat))
        # Define 'acoustics' dataframe
        dat_sim_acoustics <- make_df_detections(dat_sim_detections$det_mat,
                                                only_keep_detections = TRUE,
                                                set_names = TRUE,
                                                as_POSIXct = function(x) as.POSIXct(x, tz = "UTC"))
        dat_sim_acoustics$receiver_id <- as.integer(as.character(dat_sim_acoustics$receiver_id))
        return(dat_sim_acoustics)
      })
    return(dat_sim_detections_by_array)
  })
  print(summary(dat_sim_detections_by_path))
  saveRDS(dat_sim_detections_by_path, "./data/movement/dat_sim_detections_by_path.rds")

} else{
  dat_sim_detections_by_path <- readRDS("./data/movement/dat_sim_detections_by_path.rds")
}


######################################
######################################
#### Examine simulated patterns of space use

#### Estimate UD for simulated data
pp <- par(mfrow = prettyGraphics::par_mf(length(dat_sim_paths)))
dat_sim_paths_ud <- lapply(dat_sim_paths, function(dat_sim_path){
  # Get UD
  # dat_sim_path <- dat_sim_paths[[1]]
  dat_sim_path_spdf <- sp::SpatialPointsDataFrame(
    dat_sim_path$xy_mat_on_grid,
    data = data.frame(ID = factor(rep(1, nrow(dat_sim_path$xy_mat)))),
    proj4string = raster::crs(grid))
  dat_sim_path_ud <- adehabitatHR::kernelUD(xy = dat_sim_path_spdf, grid = 500)
  dat_sim_path_ud <- raster::raster(dat_sim_path_ud[[1]])
  # Plot UD for simulated data
  raster::plot(dat_sim_path_ud)
  prettyGraphics::add_sp_path(dat_sim_path$xy_mat, length = 0.01, lwd = 0.1)
  return(dat_sim_path_ud)
})


######################################
######################################
#### Set up algorithms (shared param)

#### Define study area
# Define above.

#### Set numeric algorithm parameters
step             <- 120
det_rng          <- 300      # as defined previously because on the scale of the grid
clock_drift      <- 5
mob_on_grid      <- mob_sim + raster::res(grid)[1] # mobility is higher than simulated when expressed at the resolution of the grid

#### Define movement model over grid
# [with a relaxation of the maximum mobility to allow for the effects of grid resolution]
calc_mpr_on_grid <- function(distance,...) {
  pr <- stats::plogis(7.5 + distance * -0.025)
  pr[distance > mob_on_grid] <- 0 # relaxation of maximum mobility
  return(pr)
}
plot(1:mob_on_grid, calc_mpr_on_grid(1:mob_on_grid), type = "l")
lines(1:mob_on_grid, calc_mpr(1:mob_on_grid), col = "red")

#### Prepare movement time series
## Examine frequency of detections
pp <- par(mfrow = c(3, 5), mar = c(1, 1, 1, 1))
lapply(dat_sim_detections_by_path, function(dat_sim_detections_by_array){
  lapply(dat_sim_detections_by_array, function(dat_sim_detections){
    # dat_sim_detections <- dat_sim_detections_by_path[[1]][[1]]
    prettyGraphics::pretty_line(dat_sim_detections$timestamp,
                                pretty_axis_args =
                                  list(side = 1,
                                       axis = list(at = utils.add::seq_range(range(study_period), length.out = 5)))
                                )
  })
}) %>% invisible()
par(pp)

#### Check mobility estimates for each path/array

#### Make folders to store outputs
lapply(1:length(dat_sim_detections_by_path), function(path_id){
  root_to_path <- paste0("./data/estimates/path", "_", path_id)
  dir.create(root_to_path)
  lapply(1:length(dat_sim_arrays), function(array_id){
    dir.create(paste0(root_to_path, "/array_", array_id))
    lapply(c("ac", "dc", "acdc", "acpf", "dcpf", "acdcpf"), function(alg){
      dir.create(paste0(root_to_path, "/array_", array_id, "/", alg, "/"))
    })
    dir.create(paste0(root_to_path, "/array_", array_id, "/ac/record/"))
    dir.create(paste0(root_to_path, "/array_", array_id, "/dc/record/"))
    dir.create(paste0(root_to_path, "/array_", array_id, "/acdc/record/"))
  })
}) %>% invisible()

#### Now proceed to implement 'flapper' algorithms.

#### End of code.
######################################
######################################
