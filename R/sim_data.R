######################################
######################################
#### sim_data.R

#### This code
# Simulates array(s), movement path(s) and movement tracks
# ... for evaluating the performance of the flapper algorithms
# ... for inferring patterns of space use
# ... and exploring their sensitivity to different parameters.

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
seed       <- 20210828
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

#### Define the detection range
true_detection_range <- 300

#### Visualise arrays with detection containers given some detection range
if (FALSE) {
  pp <- par(mfrow = prettyGraphics::par_mf(nrow(dat_sim_array_info)))
  lapply(dat_sim_arrays, function(array){
    # array <- dat_sim_arrays[[1]]
    containers <- get_detection_containers(xy = array$array$xy, detection_range = true_detection_range,
                                           boundaries = area,
                                           byid = TRUE, plot = FALSE)
    raster::plot(containers, xlim = ext[1:2], ylim = ext[3:4])
    raster::lines(ext, col = "royalblue")
  }) %>% invisible()
  par(pp)
}

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
    get_detection_area_sum(array$array$xy, detection_range = true_detection_range, boundaries = ext)
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
true_mobility <- 500
calc_mpr <- function(distance,...) {
  pr <- stats::plogis(7.5 + distance * -0.025)
  pr[distance > true_mobility] <- 0
  return(pr)
}
plot(1:true_mobility, calc_mpr(1:true_mobility), ylim = c(0, 1))

## Simulate step lengths from movement model
steps <- data.frame(distance = seq(0, true_mobility, length.out = 1e4))
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
calc_dpr <-
  function(x){
    ifelse(x <= true_detection_range, stats::plogis(3 + -0.03 * x), 0)
  }
plot(1:true_detection_range, calc_dpr(1:true_detection_range), type = "l")

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
#### Examine simulated data and emergent properties

if (FALSE) {

  #### Examine frequency of detections
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

  #### Estimate UD for simulated data
  pp <- par(mfrow = prettyGraphics::par_mf(length(dat_sim_paths)))
  dat_sim_paths_ud <- lapply(dat_sim_paths, function(dat_sim_path){
    # Get UD
    # dat_sim_path <- dat_sim_paths[[1]]
    dat_sim_path_spdf <- sp::SpatialPointsDataFrame(
      dat_sim_path$xy_mat_on_grid,
      data = data.frame(ID = factor(rep(1, nrow(dat_sim_path$xy_mat)))),
      proj4string = raster::crs(grid))
    dat_sim_path_ud <- adehabitatHR::kernelUD(xy = dat_sim_path_spdf, grid = true_mobility)
    dat_sim_path_ud <- raster::raster(dat_sim_path_ud[[1]])
    # Plot UD for simulated data
    raster::plot(dat_sim_path_ud)
    prettyGraphics::add_sp_path(dat_sim_path$xy_mat, length = 0.01, lwd = 0.1)
    return(dat_sim_path_ud)
  })

}


######################################
######################################
#### Define algorithm parameters

######################################
#### Method

# We will implement the algorithms with the 'correct' parameters as well as
# 'incorrect' ones use to examine algorithm sensitivity.
# Here, we define the list of algorithm parameters. For each simulated dataset
# we'll define a set of folders in which to store the results of the algorithms
# as applied with a particular set of parameters.


######################################
#### Define 'simple' parameters

#### Study area
# Define above.

#### Step length
step             <- 120

#### Clock drift
clock_drift      <- 5

#### Helper function
adjust_beta <- function(true_beta, true_param, param){
  true_beta * true_param/param
}


######################################
#### Define detection probability parameters

#### Define detection probability calculator wrapper
# This function takes the true detection range as 300 m
# And, if this is supplied, returns the detection probability function
# used for the simulations above. If another detection range is inputted,
# the detection range and beta parameters in the detection probability
# function are adjusted accordingly and an updated function is returned
true_beta_detection <- -0.03
calc_dprw <- function(distance, detection_range, verbose = FALSE){
  beta <- adjust_beta(true_beta_detection, true_detection_range, detection_range)
  if (verbose)
    cat(paste0("detection range: ", detection_range, "; beta: ", beta))
  f <- function(distance)
    ifelse(distance <= detection_range, stats::plogis(3 + beta * distance), 0)
  if (!missing(distance)) {
    f(distance)
  } else f
}

# Visualise examples
# Selected ('true') detection probability parameters
# Underestimation: half detection range (half), increase gradient (double beta)
# Overestimation: double detection range (double), reduce gradient (half beta)
dist <- 1:1000
plot(dist, calc_dprw(dist, 300, verbose = TRUE), type = "l")
lines(dist, calc_dprw(dist, 150, verbose = TRUE), type = "l")
lines(dist, calc_dprw(dist, 600, verbose = TRUE), type = "l")


######################################
#### Define mobility parameters

#### Define movement probability calculator wrapper
# Note that a suitable choice for mobility is higher than simulated
# when expressed at the resolution of the grid
# And the movement model needs to be relaxed accordingly
true_beta_mobility <- -0.025
calc_mprw <- function(distance, mobility, verbose = FALSE) {
  mobility_on_grid      <- mobility + 75 # raster::res(grid)[1]
  beta <- adjust_beta(true_beta_mobility, true_mobility, mobility)
  if (verbose)
    cat(paste0("mobility: ", mobility, "; beta: ", beta))
  f <- function(distance, ...) {
    pr <- stats::plogis(7.5 + beta * distance)
    # relax movement model by mobility_on_grid:
    pr[distance > mobility_on_grid] <- 0
    return(pr)
  }
  if (!missing(distance)) {
    f(distance)
  } else f
}

# Visualise examples
# Selected ('true') detection probability parameters
# Underestimation: half detection range (half), increase gradient (double beta)
# Overestimation: double detection range (double), reduce gradient (half beta)
dist <- 1:1000
plot(dist, calc_mprw(dist, 500, verbose = TRUE), type = "l")
lines(dist, calc_mprw(dist, 250, verbose = TRUE), type = "l")
lines(dist, calc_mprw(dist, 1000, verbose = TRUE), type = "l")


######################################
#### Collect simulation parameters

#### Define parameter combinations of interest
# * Correct/under/overestimate detection range with 'correct' mobility
# * Correct/under/overestimate mobility with 'correct' detection range
# (* Mis-specification of both parameters)

#### Define list of algorithm parameters
## Define detection ranges and mobilities
# We will implement the algorithms for different:
# a) detection ranges (and associated parameters)
# b) mobility values  (and associated models)
# We only need to specify the ranges of these parameters
# ... and other parameters in the model are changed proportionally
# ... via the calc_dprw() and calc_mprw() 'wrapper' functions
do2d <- TRUE
if(do2d){
  s <- c(seq(0.05, 0.25, by = 0.05), seq(0.25, 2.5, by = 0.25))
  try_detection_ranges <- c(true_detection_range, true_detection_range * s)
  try_detection_ranges <- try_detection_ranges[!duplicated(try_detection_ranges)]
  try_mobilities       <- c(true_mobility, true_mobility * s)
  try_mobilities       <- try_mobilities[!duplicated(try_mobilities)]
  combs <-
    rbind(
      data.frame(detection_range = try_detection_ranges,
                 mobility = true_mobility,
                 is_mpr = FALSE,
                 is_dpr = TRUE),
      data.frame(detection_range = true_detection_range,
                 mobility = try_mobilities,
                 is_mpr = TRUE,
                 is_dpr = FALSE)
    )
} else {
  s <- seq(0.25, 2, by = 0.25)
  try_detection_ranges <- c(true_detection_range * s)
  try_detection_ranges <- try_detection_ranges[!duplicated(try_detection_ranges)]
  try_mobilities       <- c(true_mobility, true_mobility * s)
  try_mobilities       <- try_mobilities[!duplicated(try_mobilities)]
  combs <- expand.grid(detection_range = try_detection_ranges,
                       mobility = try_mobilities,
                       is_mpr = NA,
                       is_dpr = NA)
  nrow(combs)
}
# Check selected combinations of parameters
alg_param <-
  combs |>
  dplyr::mutate(id = paste0(detection_range, "_", mobility),
                mobility_on_grid = mobility + raster::res(grid)[1],
                correct = detection_range == true_detection_range & mobility == true_mobility,
                is_dpr = ifelse(correct, TRUE, is_dpr),
                is_mpr = ifelse(correct, TRUE, is_mpr))  |>
  dplyr::filter(!duplicated(id)) |>
  dplyr::mutate(index = dplyr::row_number()) |>
  dplyr::select(index, id, is_mpr, is_dpr, correct, detection_range, mobility, mobility_on_grid)

# Write tidy table
if(do2d){
  alg_param |>
    dplyr::slice(-1L) |>
    dplyr::mutate(ID = paste0("S3 (", index - 1, ")"),
                  dpr_beta = adjust_beta(true_beta_detection, true_detection_range, detection_range) |>
                    round(3) |> add_lagging_point_zero(n = 3),
                  dpr_gamma = detection_range,
                  mpr_beta = adjust_beta(true_beta_mobility, true_mobility, mobility)
                  |> round(3) |> add_lagging_point_zero(n = 3),
                  mpr_delta = mobility
    ) |>
    dplyr::select(ID,
                  dpr_beta, dpr_gamma,
                  mpr_beta, mpr_delta
    ) |>
    write.table("./fig/dat_sim_param.txt",
                quote = FALSE, row.names = FALSE, na = "", sep = ",")
  # Visualise algorithm parameters
  png("./fig/algorithm_parameters.png",
      height = 6, width = 11, units = "in", res = 600 )
  pp <- par(mfrow = c(1, 2), oma = c(2, 2, 2, 2))
  pretty_plot(range(alg_param$detection_range), c(0, 1),
              type = "n", xlab = "", ylab = "")
  lapply(try_detection_ranges, function(rng) {
    x <- seq(0, rng, length.out = 100)
    lines(x, calc_dprw(x, rng, verbose = TRUE))
    points(rng, 0, pch = 21, bg = "black", cex = 0.75)
  })
  mtext(side = 1, "Distance (m) from receiver location", line = 2)
  mtext(side = 2, "Detection probability", line = 2.5)
  mtext(side = 3, "A", adj = 0, font = 2, cex = 1.5)
  x <- seq(0, true_detection_range, length.out = 100)
  lines(x, calc_dprw(x, true_detection_range), lwd = 3)
  legend("topright",
         legend = c("'correct' model", expression(gamma)),
         lwd = c(3, NA),
         pch = c(NA, 21),
         pt.bg = c(NA, "black"),
         pt.cex = 0.75,
         y.intersp = 1.25,
         bty = "n")
  pretty_plot(range(alg_param$mobility), c(0, 1),
              type = "n", xlab = "", ylab = "")
  lapply(try_mobilities, function(rng) {
    x <- seq(0, rng, length.out = 100)
    lines(x, calc_mprw(x, rng, verbose = TRUE))
    points(rng, 0, pch = 21, bg = "black", cex = 0.75)
  })
  x <- seq(0, true_mobility, length.out = 100)
  lines(x, calc_mprw(x, true_mobility), lwd = 3)
  legend("topright",
         legend = c("'correct' model", expression(Delta(T[1], T[2]))),
         lwd = c(3, NA),
         pch = c(NA, 21),
         pt.bg = c(NA, "black"),
         pt.cex = 0.75,
         y.intersp = 1.25,
         bty = "n")
  mtext(side = 1, "Distance (m) from previous location", line = 2)
  mtext(side = 2, "Movement probability", line = 2.5)
  mtext(side = 3, "B", adj = 0, font = 2, cex = 1.5)
  par(pp)
  dev.off()
}

# Define list
alg_param <- split(alg_param, alg_param$index)
names(alg_param) <- sapply(alg_param, \(elm) elm$id)
names(alg_param)


######################################
#### Set up directory system

#### Folder structure
# {dataset}/          {algorithm implementation}/ {algorithm} ->
# {path_*}/{array_*}/ {alg_imp_*}/ {algorithm}

#### Make folders to store outputs
clean <- FALSE
if (clean) unlink(file.path("data", "estimates", "path_1"), recursive = TRUE)
lapply(1:length(dat_sim_detections_by_path), function(path_id){
  root_to_path <- file.path("data", "estimates", paste0("path_", path_id))
  dir.create(root_to_path)
  lapply(1:length(dat_sim_arrays), function(array_id){
    root_to_array <- file.path(root_to_path, paste0("array_", array_id))
    dir.create(root_to_array)
    lapply(names(alg_param), function(alg_id){
      root_to_alg <- file.path(root_to_array, paste0("alg_imp_", alg_id))
      dir.create(root_to_alg)
      lapply(c("ac", "dc", "acdc", "acpf", "dcpf", "acdcpf"), function(alg){
        dir.create(file.path(root_to_alg, alg))
        })
      dir.create(file.path(root_to_alg, "ac", "record"))
      dir.create(file.path(root_to_alg, "dc", "record"))
      dir.create(file.path(root_to_alg, "acdc", "record"))
    })
  })
}) %>% invisible()

#### Define folder name with the 'correct' outputs (for any given array)
alg_true     <- paste0(true_detection_range, "_", true_mobility)
alg_imp_true <- paste0("alg_imp_", true_detection_range, "_", true_mobility)

#### Now proceed to implement 'flapper' algorithms.

#### End of code.
######################################
######################################
