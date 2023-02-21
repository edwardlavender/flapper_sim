######################################
######################################
#### sim_workhorse.R

#### This code:
# 1) Implements the COA-KUD approach and the flapper family of algorithms
# ... to reconstruct movement patterns and patterns of space use from simulated data,
# ... given a specific movement path and array ID. The main purpose is to generate
# ... outputs from the following algorithms:
# ... * ac, acdc
# ... * acpf, acdcpf

#### Steps preceding this code:
# 1) Implement sim_data.R


######################################
######################################
#### Define data and algorithm parameters

#### Define input controls
# If the script is implemented manually ('manual = TRUE'), we will run all code
# Otherwise, plots and 'extra' code will be suppressed (to improve speed) when
# ... applied iteratively (over multiple paths/arrays/algorithm implementations)
manual <- FALSE
if(manual){
  path_id  <- 1; array_id <- 1  # dataset controls
  alg <- alg_param[["300_50"]] # algorithm parameters
}

#### Get array/movement path and associated data
con_root   <- paste0("./data/estimates/path_", path_id, "/array_", array_id, "/",
                     "alg_imp_", alg$id, "/")
array      <- dat_sim_arrays[[array_id]]
moorings   <- dat_sim_moorings[[array_id]]
archival   <- dat_sim_archival_by_path[[path_id]]
acoustics  <- dat_sim_detections_by_path[[path_id]][[array_id]]
path       <- dat_sim_paths[[path_id]]

#### Define-algorithm-specific parameters
det_rng          <- alg$detection_range
calc_dpr         <- calc_dprw(detection_range = alg$detection_range)
mob_on_grid      <- alg$mobility_on_grid
calc_mpr_on_grid <- calc_mprw(mobility = alg$mobility)

#### Define 'global' param
kud_grid_resolution <- 120


######################################
######################################
#### Implement COA approach

if (manual){

  #### Guess a suitable delta_t interval
  delta_t <- "2 hours"

  pp <- par(mfrow = c(1, 2))
  coa_setup_delta_t(acoustics, delta_t = delta_t, method = 1:2L)
  par(pp)

  #### Calculate COAs
  acoustics$individual_id <- 1L
  acoustics$receiver_id   <- factor(acoustics$receiver_id, moorings$receiver_id)
  acoustics_mat           <- make_matrix_detections(acoustics, delta_t = delta_t)
  out_coa                 <- coa(acoustics_mat, sp::coordinates(array$array$xy))
  head(out_coa)

  #### Examine COAs
  raster::plot(area)
  points(out_coa, col ="red")

  #### Estimate KUD from COAs
  if(nrow(out_coa) >= 5){
    out_coa_spdf <- sp::SpatialPointsDataFrame(
      out_coa[, c("x", "y")],
      data = data.frame(ID = factor(rep(1, nrow(out_coa)))),
      proj4string = raster::crs(grid))
    out_coa_ud <- adehabitatHR::kernelUD(xy = out_coa_spdf, grid = kud_grid_resolution)
    out_coa_ud <- raster::raster(out_coa_ud[[1]])
    out_coa_ud <- raster::resample(out_coa_ud, grid)
    out_coa_ud <- out_coa_ud/raster::cellStats(out_coa_ud, "sum")
    raster::cellStats(out_coa_ud, sum)
    raster::plot(out_coa_ud)
    prettyGraphics::add_sp_path(path$xy_mat_on_grid, length = 0.01, lwd = 0.25)
  } else out_coa_ud <- NULL

}


######################################
######################################
#### Flapper algorithm: array/algorithm-specific set up

#### Define acoustic containers
detection_containers_file <- paste0(con_root, "detection_containers.rds")
if(!file.exists(detection_containers_file)){
  rxy <- array$array$xy
  rxy <- sp::SpatialPointsDataFrame(rxy, data.frame(receiver_id = 1:length(rxy)))
  detection_containers <- acs_setup_containers(xy = rxy,
                                               detection_range = det_rng,
                                               boundaries = area,
                                               resolution = 1000,
                                               plot = manual,
                                               verbose = TRUE
  )
  if(manual){
    coords <- sp::coordinates(array$array$xy)
    text(coords[, 1], coords[, 2], rxy$receiver_id, font = 2, cex = 2)
  }
  saveRDS(detection_containers, detection_containers_file)
} else {
  detection_containers <- readRDS(detection_containers_file)
}

#### Detection container overlaps
detection_centroids_overlaps_file <-  paste0(con_root, "detection_containers_overlaps.rds")
if(!file.exists(detection_centroids_overlaps_file)){
  containers_dets <- get_detection_containers(xy = array$array$xy,
                                              detection_range = det_rng,
                                              boundaries = area,
                                              byid = TRUE, resolution = 1000)
  containers_df <- moorings
  row.names(containers_df) <- names(containers_dets)
  containers_dets <- sp::SpatialPolygonsDataFrame(containers_dets, containers_df)
  overlaps <- get_detection_containers_overlap(containers =  containers_dets)
  saveRDS(overlaps, detection_centroids_overlaps_file)
} else {
  overlaps <- readRDS(detection_centroids_overlaps_file)
}

#### Define detection kernels
detection_kernels_file <- paste0(con_root, "detection_kernels.rds")
if(!file.exists(detection_kernels_file)){
  rxy <- array$array$xy
  rxy <- sp::SpatialPointsDataFrame(rxy, moorings)
  kernels <- acs_setup_detection_kernels(xy = rxy,
                                         containers = detection_containers,
                                         overlaps = overlaps,
                                         calc_detection_pr = calc_dpr,
                                         bathy = grid)
  saveRDS(kernels, detection_kernels_file)
} else {
  kernels <- readRDS(detection_kernels_file)
}


######################################
######################################
#### Prepare data

#### Prepare movement data
# For improved efficiency, this preparation could be implemented only once
# for each dataset (and not also for each algorithm implementation) but
# the current cost of the processing implemented here is small.
acoustics$receiver_id <- as.integer(as.character(acoustics$receiver_id))
acoustics$index      <- 1:nrow(acoustics)
archival$index <- 1:nrow(archival)
archival       <- archival[archival$timestamp >= min(acoustics$timestamp) &
                             archival$timestamp <= max(acoustics$timestamp), ]

#### Focus on parts of paths that detections enable us to resolve
path$xy_mat_on_grid_within_acoustics <- path$xy_mat_on_grid[archival$index, ]

#### Estimate KUD from path
if(manual){
  if(nrow(path$xy_mat_on_grid_within_acoustics) >= 5){
    path_spdf <- sp::SpatialPointsDataFrame(
      path$xy_mat_on_grid_within_acoustics,
      data = data.frame(ID = factor(rep(1, nrow(path$xy_mat_on_grid_within_acoustics)))),
      proj4string = raster::crs(grid))
    path_ud <- adehabitatHR::kernelUD(xy = path_spdf, grid = kud_grid_resolution)
    path_ud <- raster::raster(path_ud[[1]])
    path_ud <- raster::resample(path_ud, grid)
    path_ud <- path_ud/raster::cellStats(path_ud, sum)
    raster::cellStats(path_ud, sum)
    if(manual){
      raster::plot(path_ud)
      prettyGraphics::add_sp_path(path$xy_mat_on_grid_within_acoustics,
                                  length = 0.05)
    }
  } else path_ud <- NULL
}


######################################
######################################
#### Implement AC and ACDC algorithms

#### Implement AC algorithm
out_ac_file         <- paste0(con_root, "ac/out_ac.rds")
out_ac_success_file <- paste0(con_root, "ac/out_ac_success.rds")
if(!file.exists(out_ac_success_file)){
  out_ac <- tryCatch(
    ac(acoustics = acoustics,
       step = step,
       bathy = grid,
       detection_containers = detection_containers,
       detection_kernels = kernels, detection_kernels_overlap = overlaps, detection_time_window = clock_drift,
       normalise = TRUE,
       mobility = mob_on_grid,
       save_record_spatial = 0,
       con = paste0(con_root, "ac/"),
       write_record_spatial_for_pf = list(filename = paste0(con_root, "ac/record/"), format = "GTiff",
                                          overwrite = TRUE)
    ), error = function(e) e)
  if (inherits(out_ac, "error")) {
    out_ac_success <- FALSE
    print("ERROR")
    print(out_ac)
  } else {
    out_ac_success <- TRUE
    saveRDS(out_ac, out_ac_file)
  }
  saveRDS(out_ac_success, out_ac_success_file)
} else {
  out_ac_success <- readRDS(out_ac_success_file)
  if (out_ac_success) out_ac <- readRDS(out_ac_file)
}

#### Check normalisation
if(FALSE){
  out_ac_record_sums <-
    pbapply::pblapply(pf_setup_record(paste0(con_root, "ac/record/"))[1:10], function(x){
      raster::cellStats(raster::raster(x), "sum")
    }) %>% unlist()
  unique(out_ac_record_sums)
}

#### Implement ACDC algorithm
out_acdc_success <- out_ac_success
if(out_ac_success){
  out_acdc_file         <- paste0(con_root, "acdc/out_acdc.rds")
  out_acdc_success_file <- paste0(con_root, "acdc/out_acdc_success.rds")
  if(!file.exists(out_acdc_success_file)){
    out_acdc <-
      tryCatch(
        acdc(acoustics = acoustics,
             archival = archival,
             bathy = grid,
             detection_containers = detection_containers,
             detection_kernels = kernels, detection_kernels_overlap = overlaps, detection_time_window = clock_drift,
             normalise = TRUE,
             mobility = mob_on_grid,
             calc_depth_error = function(...) matrix(c(-5, 5), nrow = 2),
             save_record_spatial = 0,
             con = paste0(con_root, "acdc/"),
             write_record_spatial_for_pf = list(filename = paste0(con_root, "acdc/record/"), format = "GTiff",
                                                overwrite = TRUE)
        ), error = function(e) e)
    if (inherits(out_acdc, "error")) {
      out_acdc_success <- FALSE
      print("ERROR")
      print(out_acdc)
    } else {
      out_acdc_success <- TRUE
      saveRDS(out_acdc, out_acdc_file)
    }
    saveRDS(out_acdc_success, out_acdc_success_file)
  } else {
    out_acdc_success <- readRDS(out_acdc_success_file)
    if (out_acdc_success) out_acdc <- readRDS(out_acdc_file)
  }
}

#### Check normalisation
if(FALSE){
  out_acdc_record_sums <-
    pbapply::pblapply(pf_setup_record(paste0(con_root, "acdc/record/"))[1:10], function(x){
      raster::cellStats(raster::raster(x), "sum")
    }) %>% unlist()
  unique(out_acdc_record_sums)
}


######################################
######################################
#### Examine AC/ACDC results

#### Collate AC/ACDC outputs
if(manual){

  if(out_ac_success) out_ac_s     <- acdc_simplify(out_ac)
  if(out_acdc_success) out_acdc_s <- acdc_simplify(out_acdc)

  #### Plotting window
  pp <- par(mfrow = c(1, 2))

  #### AC
  if(out_ac_success){
    raster::plot(out_ac_s$map)
    prettyGraphics::add_sp_path(path$xy_mat_on_grid, length = 0.01)
  }

  #### ACDC
  if(out_acdc_success){
    raster::plot(out_acdc_s$map)
    prettyGraphics::add_sp_path(path$xy_mat_on_grid, length = 0.01)
  }
  par(pp)

}


######################################
######################################
#### Implement particle filtering

#### ACPF
out_acpf_success <- out_ac_success
if(out_ac_success) {
  ## Implement algorithm
  cat(paste0(con_root, "ac/record/"))
  print(list.files(paste0(con_root, "ac/record/")))
  ## Define record
  out_ac_record <- pf_setup_record(paste0(con_root, "ac/record/"), pattern = ".tif")
  out_acpf_file <- paste0(con_root, "acpf/out_acpf.rds")
  if(!file.exists(out_acpf_file)){
    out_acpf <- pf(record = out_ac_record,
                   calc_movement_pr = calc_mpr_on_grid, # note relaxed movement model for grid
                   mobility = mob_on_grid,              # note relaxed mobility
                   n = 100L,
                   con = paste0(con_root, "acpf/acpf_log.txt"),
                   seed = seed
    )
    saveRDS(out_acpf, out_acpf_file)
  } else {
    out_acpf <- readRDS(out_acpf_file)
  }
  # Check convergence
  # If the algorithm fails to converge, pf() will return the outputs up to the
  # moment in time at which the convergence failure was realised.
  out_acpf_success <- ifelse(length(out_acpf$history) != length(out_ac_record), FALSE, TRUE)
}

#### ACDCPF
out_acdcpf_success <- out_acdc_success
if(out_acdc_success){
  ## Define record
  out_acdc_record <- pf_setup_record(paste0(con_root, "acdc/record/"))
  ## Implement algorithm
  out_acdcpf_file <- paste0(con_root, "acdcpf/out_acdcpf.rds")
  if(!file.exists(out_acdcpf_file)){
    out_acdcpf <- pf(record = out_acdc_record,
                     calc_movement_pr = calc_mpr_on_grid, # note relaxed movement model for grid
                     mobility = mob_on_grid,              # note relaxed mobility
                     n = 100L,
                     con = paste0(con_root, "acdcpf/acdcpf_log.txt"),
                     seed = seed
    )
    saveRDS(out_acdcpf, out_acdcpf_file)
  } else {
    out_acdcpf <- readRDS(out_acdcpf_file)
  }
  # Check convergence
  out_acdcpf_success <- ifelse(length(out_acdcpf$history) != length(out_acdc_record), FALSE, TRUE)
}


######################################
######################################
#### Process particle samples

#### Examine overall histories for example time steps
if(manual){
  pp <- par(mfrow = c(1, 2))
  if(out_ac_success) pf_plot_history(out_acpf, time_steps = 1)
  if(out_acdc_success) pf_plot_history(out_acdcpf, time_steps = 1)
  par(pp)
}

#### Assemble particle histories for connected cell pairs
## ACPF [~2.5 minutes]
if(out_acpf_success){
  out_acpf_pairs_file <- paste0(con_root, "acpf/out_acpf_pairs.rds")
  if(!file.exists(out_acpf_pairs_file)){
    out_acpf$args$bathy <- grid
    out_acpf_pairs <- pf_simplify(out_acpf,
                                  cl = NULL,
                                  return = "archive"
    )
    saveRDS(out_acpf_pairs, out_acpf_pairs_file)
  } else {
    out_acpf_pairs <- readRDS(out_acpf_pairs_file)
  }
}
## ACDCPF
if(out_acdcpf_success){
  out_acdcpf_pairs_file <- paste0(con_root, "acdcpf/out_acdcpf_pairs.rds")
  if(!file.exists(out_acdcpf_pairs_file)){
    out_acdcpf$args$bathy <- grid
    out_acdcpf_pairs <- pf_simplify(out_acdcpf,
                                    cl = NULL,
                                    return = "archive"
    )
    saveRDS(out_acdcpf_pairs, out_acdcpf_file)
  } else {
    out_acdcpf_pairs <- readRDS(out_acdcpf_pairs_file)
  }
}

#### Simplify particle histories to retain unique particles
if(out_acpf_success){
  out_acpf_pairs_unq_file <- paste0(con_root, "acpf/out_acpf_pairs_unq.rds")
  if (!file.exists(out_acpf_pairs_unq_file)) {
    out_acpf_pairs_unq   <- pf_simplify(out_acpf_pairs, summarise_pr = TRUE, return = "archive")
    saveRDS(out_acpf_pairs_unq, out_acpf_pairs_unq_file)
  } else {
    out_acpf_pairs_unq <- readRDS(out_acpf_pairs_unq_file)
  }
  # Check normalisation for an example time step
  sum(out_acpf_pairs_unq[["history"]][[2]][["pr_current"]])
}
if(out_acdcpf_success){
  out_acdcpf_pairs_unq_file <- paste0(con_root, "acdcpf/out_acdcpf_pairs_unq.rds")
  if (!file.exists(out_acdcpf_pairs_unq_file)) {
    out_acdcpf_pairs_unq   <- pf_simplify(out_acdcpf_pairs, summarise_pr = TRUE, return = "archive")
    saveRDS(out_acdcpf_pairs_unq, out_acdcpf_pairs_unq_file)
  } else {
    out_acdcpf_pairs_unq <- readRDS(out_acdcpf_pairs_unq_file)
  }
  # Check normalisation for an example time step
  sum(out_acdcpf_pairs_unq[["history"]][[2]][["pr_current"]])
}


######################################
######################################
#### Generate maps

#### Define scale parameters
# Scale maps between 0 and 1 for comparability
scale_pou <- "max"
scale_kud <- "max"
kud_grid  <- kud_habitat(grid)
kud_size  <- 100L

#### POU maps
if(out_acpf_success){
  out_acpf_pou_file <- paste0(con_root, "acpf/out_acpf_pou.tif")
  if(!file.exists(out_acpf_pou_file)){
    out_acpf_pou_raw  <- pf_plot_map(out_acpf_pairs_unq, grid, scale = "original")
    out_acpf_pou      <- out_acpf_pou_raw/raster::cellStats(out_acpf_pou_raw, "max")
    raster::writeRaster(out_acpf_pou_raw, paste0(con_root, "acpf/out_acpf_pou_raw.tif"), overwrite = TRUE)
    raster::writeRaster(out_acpf_pou, out_acpf_pou_file, overwrite = TRUE)
  } else {
    out_acpf_pou <- raster::raster(out_acpf_pou_file)
  }
}
if(out_acdcpf_success){
  out_acdcpf_pou_file <- paste0(con_root, "acdcpf/out_acdcpf_pou.tif")
  if(!file.exists(out_acdcpf_pou_file)){
    out_acdcpf_pou_raw  <- pf_plot_map(out_acdcpf_pairs_unq, grid, scale = "original")
    out_acdcpf_pou      <- out_acdcpf_pou_raw/raster::cellStats(out_acdcpf_pou_raw, "max")
    raster::writeRaster(out_acdcpf_pou_raw, paste0(con_root, "acdcpf/out_acdcpf_pou_raw.tif"), overwrite = TRUE)
    raster::writeRaster(out_acdcpf_pou, out_acdcpf_pou_file, overwrite = TRUE)
  } else {
    out_acdcpf_pou <- raster::raster(out_acdcpf_pou_file)
  }
}

#### Get KUDs
if(out_acpf_success){
  out_acpf_kud_file <- paste0(con_root, "acpf/out_acpf_kud.tif")
  if(!file.exists(out_acpf_kud_file)){
    out_acpf_kud_raw   <- pf_kud(out_acpf_pou_raw,
                                 sample_size = kud_size,
                                 grid = kud_grid)
    out_acpf_kud   <- out_acpf_kud_raw/raster::cellStats(out_acpf_kud_raw, scale_kud)
    raster::writeRaster(out_acpf_kud_raw, paste0(con_root, "acpf/out_acpf_kud_raw.tif"), overwrite = TRUE)
    raster::writeRaster(out_acpf_kud, out_acpf_kud_file, overwrite = TRUE)
  } else {
    out_acpf_kud <- raster::raster(out_acpf_kud_file)
  }
}
if(out_acdcpf_success){
  out_acdcpf_kud_file <- paste0(con_root, "acdcpf/out_acdcpf_kud.tif")
  if(!file.exists(out_acdcpf_kud_file)){
    out_acdcpf_kud_raw   <- pf_kud(out_acdcpf_pou_raw,
                                 sample_size = kud_size,
                                 grid = kud_grid)
    out_acdcpf_kud   <- out_acdcpf_kud_raw/raster::cellStats(out_acdcpf_kud_raw, scale_kud)
    raster::writeRaster(out_acdcpf_kud_raw, paste0(con_root, "acdcpf/out_acdcpf_kud_raw.tif"), overwrite = TRUE)
    raster::writeRaster(out_acdcpf_kud, out_acdcpf_kud_file, overwrite = TRUE)
  } else {
    out_acdcpf_kud <- raster::raster(out_acdcpf_kud_file)
  }
}


#### End of code.
######################################
######################################
