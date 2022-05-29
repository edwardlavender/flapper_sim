######################################
######################################
#### sim_workhorse.R

#### This code:
# 1) Implements the COA-KUD approach and the flapper family of algorithms
# ... to reconstruct movement patterns and patterns of space use from simulated data,
# ... given a specific movement path and array ID.

#### Steps preceding this code:
# 1) Implement sim_data.R


######################################
######################################
#### Path and array details

#### Get array/path details
# path_id  <- 1; array_id <- 12
con_root <- paste0("./data/estimates/path_", path_id, "/array_", array_id, "/")

#### Get associated data
array      <- dat_sim_arrays[[array_id]]
moorings   <- dat_sim_moorings[[array_id]]
archival   <- dat_sim_archival_by_path[[path_id]]
acoustics  <- dat_sim_detections_by_path[[path_id]][[array_id]]
path       <- dat_sim_paths[[path_id]]

#### Global param
kud_grid_resolution <- 120


######################################
######################################
#### Implement COA approach

#### Guess a suitable delta_t interval
delta_t <- "2 hours"
check_delta_t <- FALSE
if(check_delta_t) {
  pp <- par(mfrow = c(1, 2))
  coa_setup_delta_t(acoustics, delta_t = delta_t, method = 1:2L)
  par(pp)
}

#### Calculate COAs
acoustics$individual_id <- 1L
acoustics$receiver_id <- factor(acoustics$receiver_id, moorings$receiver_id)
acoustics_mat <- make_matrix_detections(acoustics, delta_t = delta_t)
out_coa <- coa(acoustics_mat, sp::coordinates(array$array$xy))
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


######################################
######################################
#### Flapper algorithm: array-specific set up

#### Define acoustic containers
make_detection_containers <- TRUE
if(make_detection_containers){
  rxy <- array$array$xy
  rxy <- sp::SpatialPointsDataFrame(rxy, data.frame(receiver_id = 1:length(rxy)))
  detection_containers <- acs_setup_containers(xy = rxy,
                                            detection_range = det_rng,
                                            boundaries = area,
                                            resolution = 1000,
                                            plot = FALSE,
                                            verbose = TRUE
                                            )
  saveRDS(detection_containers, paste0(con_root, "detection_containers.rds"))
} else {
  detection_containers <- readRDS(paste0(con_root, "detection_containers.rds"))
}

#### Detection container overlaps
make_detection_centriods_overlaps <- TRUE
if(make_detection_centriods_overlaps){
  containers_dets <- get_detection_containers(xy = array$array$xy,
                                            detection_range = det_rng,
                                            boundaries = area,
                                            byid = TRUE, resolution = 1000)
  containers_df <- moorings
  row.names(containers_df) <- names(containers_dets)
  containers_dets <- sp::SpatialPolygonsDataFrame(containers_dets, containers_df)
  overlaps <- get_detection_containers_overlap(containers =  containers_dets)
  saveRDS(overlaps, paste0(con_root, "detection_containers_overlaps.rds"))
} else {
  overlaps <- readRDS(paste0(con_root, "detection_containers_overlaps.rds"))
}

#### Define detection kernels
make_detection_kernels <- TRUE
if(make_detection_kernels){
  rxy <- array$array$xy
  rxy <- sp::SpatialPointsDataFrame(rxy, moorings)
  kernels <- acs_setup_detection_kernels(xy = rxy,
                                         containers = detection_containers,
                                         overlaps = overlaps,
                                         calc_detection_pr = calc_dpr,
                                         bathy = grid)
  saveRDS(kernels, paste0(con_root, "detection_kernels.rds"))
} else {
  kernels <- readRDS(paste0(con_root, "detection_kernels.rds"))
}


######################################
######################################
#### Prepare data

#### Prepare movement data
acoustics$receiver_id <- as.integer(as.character(acoustics$receiver_id))
acoustics$index      <- 1:nrow(acoustics)
archival$index <- 1:nrow(archival)
archival       <- archival[archival$timestamp >= min(acoustics$timestamp) &
                             archival$timestamp <= max(acoustics$timestamp), ]

#### Focus on parts of paths that detections enable us to resolve
path$xy_mat_on_grid_within_acoustics <- path$xy_mat_on_grid[archival$index, ]

#### Estimate KUD from path
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
  raster::plot(path_ud)
  prettyGraphics::add_sp_path(path$xy_mat_on_grid_within_acoustics,
                              length = 0.05)
} else path_ud <- NULL


######################################
######################################
#### Implement AC and ACDC algorithms

#### Implement AC algorithm
run_ac <- TRUE
if(run_ac){
  out_ac <- ac(acoustics = acoustics,
               step = step,
               bathy = grid,
               detection_containers = detection_containers,
               detection_kernels = kernels, detection_kernels_overlap = overlaps, detection_time_window = clock_drift,
               normalise = TRUE,
               mobility = mob_on_grid,
               save_record_spatial = 0,
               con = paste0(con_root, "ac/"),
               write_record_spatial_for_pf = list(filename = paste0(con_root, "ac/record/"), format = "GTiff",
                                                  overwrite = TRUE),
               normalise = TRUE
               )
  saveRDS(out_ac, paste0(con_root, "ac/out_ac.rds"))
} else{
  out_ac <- readRDS(paste0(con_root, "ac/out_ac.rds"))
}

#### Check normalisation
check_normalisation <- FALSE
if(check_normalisation){
  out_ac_record_sums <-
    pbapply::pblapply(pf_setup_record(paste0(con_root, "ac/record/"))[1:10], function(x){
      raster::cellStats(raster::raster(x), "sum")
    }) %>% unlist()
  unique(out_ac_record_sums)
}

#### Implement ACDC algorithm
run_acdc <- TRUE
if(run_acdc){
  out_acdc <- acdc(acoustics = acoustics,
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
                   )
  saveRDS(out_acdc, paste0(con_root, "acdc/out_acdc.rds"))
} else{
  out_acdc <- readRDS(paste0(con_root, "acdc/out_acdc.rds"))
}

#### Check normalisation
check_normalisation <- FALSE
if(check_normalisation){
  out_acdc_record_sums <-
    pbapply::pblapply(pf_setup_record(paste0(con_root, "acdc/record/"))[1:10], function(x){
      raster::cellStats(raster::raster(x), "sum")
    }) %>% unlist()
  unique(out_acdc_record_sums)
}


######################################
######################################
#### Examine AC/ACDC results

#### Plotting window
pp <- par(mfrow = c(1, 2))

#### AC
## Simplify outputs
out_ac_s <- acdc_simplify(out_ac)
## Overall map
raster::plot(out_ac_s$map)
prettyGraphics::add_sp_path(path$xy_mat_on_grid, length = 0.01)

#### ACDC
## Simplify outputs
out_acdc_s <- acdc_simplify(out_acdc)
## Overall map
raster::plot(out_acdc_s$map)
prettyGraphics::add_sp_path(path$xy_mat_on_grid, length = 0.01)
par(pp)


######################################
######################################
#### Implement particle filtering

#### ACPF
## Define record
out_ac_record <- pf_setup_record(paste0(con_root, "ac/record/"))

## Implement algorithm [5 minutes]
run_acpf <- TRUE
if(run_acpf){
  out_acpf <- pf(record = out_ac_record,
                 calc_movement_pr = calc_mpr_on_grid, # note relaxed movement model for grid
                 mobility = mob_on_grid,              # note relaxed mobility
                 n = 100L,
                 con = paste0(con_root, "acpf/acpf_log.txt"),
                 seed = seed
                 )
  saveRDS(out_acpf, paste0(con_root, "acpf/out_acpf.rds"))

} else {
  out_acpf <- readRDS(paste0(con_root, "acpf/out_acpf.rds"))
}

#### ACDCPF
## Define record
out_acdc_record <- pf_setup_record(paste0(con_root, "acdc/record/"))
## Implement algorithm
run_acdcpf <- TRUE
if(run_acdcpf){
  out_acdcpf <- pf(record = out_acdc_record,
                   calc_movement_pr = calc_mpr_on_grid, # note relaxed movement model for grid
                   mobility = mob_on_grid,              # note relaxed mobility
                   n = 100L,
                   con = paste0(con_root, "acdcpf/acdcpf_log.txt"),
                   seed = seed
                   )
  saveRDS(out_acdcpf, paste0(con_root, "acdcpf/out_acdcpf.rds"))
} else {
  out_acdcpf <- readRDS(paste0(con_root, "acdcpf/out_acdcpf.rds"))
}


######################################
######################################
#### Process particle samples

#### Examine overall histories for example time steps
pp <- par(mfrow = c(1, 2))
pf_plot_history(out_acpf, time_steps = 1)
pf_plot_history(out_acdcpf, time_steps = 1)
par(pp)

#### Assemble particle histories for connected cell pairs
## ACPF [~2.5 minutes]
run_pf_simplify <- TRUE
if(run_pf_simplify){
  out_acpf$args$bathy <- grid
  out_acpf_pairs <- pf_simplify(out_acpf,
                                cl = NULL,
                                return = "archive"
                                )
  saveRDS(out_acpf_pairs, paste0(con_root, "acpf/out_acpf_pairs.rds"))
} else {
  out_acpf_pairs <- readRDS(paste0(con_root, "acpf/out_acpf_pairs.rds"))
}
## ACDCPF
run_pf_simplify <- TRUE
if(run_pf_simplify){
  out_acdcpf$args$bathy <- grid
  out_acdcpf_pairs <- pf_simplify(out_acdcpf,
                                  cl = NULL,
                                  return = "archive"
                                  )
  saveRDS(out_acdcpf_pairs, paste0(con_root, "acdcpf/out_acdcpf_pairs.rds"))
} else {
  out_acdcpf_pairs <- readRDS(paste0(con_root, "acdcpf/out_acdcpf_pairs.rds"))
}

#### Simplify particle histories to retain unique particles
out_acpf_pairs_unq   <- pf_simplify(out_acpf_pairs, summarise_pr = TRUE, return = "archive")
out_acdcpf_pairs_unq <- pf_simplify(out_acdcpf_pairs, summarise_pr = TRUE, return = "archive")
# Check normalisation for an example time step
sum(out_acpf_pairs_unq[["history"]][[2]][["pr_current"]])
sum(out_acdcpf_pairs_unq[["history"]][[2]][["pr_current"]])

#### Build a sample of paths
# To assemble a paths with max_n_copies = 5L and max_n_paths = 10000L this takes ~ 1 minute
build_paths <- TRUE
if(build_paths){
  ## ACPF paths
  set.seed(seed)
  out_acpf_paths <- pf_simplify(out_acpf_pairs,
                                bathy = grid,
                                max_n_copies = NULL,
                                max_n_paths = 1000L,
                                return = "path"
                                )
  max(out_acpf_paths$path_id); nrow(out_acpf_paths)
  saveRDS(out_acpf_paths, paste0(con_root, "acpf/out_acpf_paths.rds"))
  ## ACDCPF paths
  set.seed(seed)
  out_acdcpf_paths <- pf_simplify(out_acdcpf_pairs,
                                  bathy = grid,
                                  max_n_copies = NULL,
                                  max_n_paths = 1000L,
                                  return = "path"
                                  )
  max(out_acdcpf_paths$path_id); nrow(out_acdcpf_paths)
  saveRDS(out_acdcpf_paths, paste0(con_root, "acdcpf/out_acdcpf_paths.rds"))
} else {
  out_acpf_paths <- readRDS(paste0(con_root, "acpf/out_acpf_paths.rds"))
  out_acdcpf_paths <- readRDS(paste0(con_root, "acdcpf/out_acdcpf_paths.rds"))
}

#### Sub-sample a selection of paths for estimation speed
subset_paths <- FALSE
if(subset_paths){
  n_paths <- 50
  out_acpf_paths_ll   <- pf_loglik(out_acpf_paths)
  out_acpf_paths      <- out_acpf_paths[out_acpf_paths$path_id %in% out_acpf_paths_ll$path_id[1:n_paths], ]
  out_acdcpf_paths_ll <- pf_loglik(out_acdcpf_paths)
  out_acdcpf_paths    <- out_acdcpf_paths[out_acdcpf_paths$path_id %in% out_acdcpf_paths_ll$path_id[1:n_paths], ]
}

#### Examine overall maps
pp <- par(mfrow = c(2, 2))
pf_plot_map(out_acpf_pairs_unq, map = grid, scale = "sum")
pf_plot_map(out_acdcpf_pairs_unq, map = grid, scale = "sum")
pf_plot_map(out_acpf_paths, map = grid, scale = "sum")
pf_plot_map(out_acdcpf_paths, map = grid, scale = "sum")
par(pp)

#### Check paths
pf_plot_1d(out_acdcpf_paths, archival)
pf_plot_2d(out_acdcpf_paths[out_acdcpf_paths$path_id == 1, ],
           bathy = grid,
           add_paths = list(length = 0.01))

#### Examine alternative KUD approaches
trial_kuds <- FALSE
if(trial_kuds){

  ## Trial approaches based on particles
  # The truth
  pp <- par(mfrow = c(2, 3))
  prettyGraphics::pretty_map(add_rasters = list(x = path_ud), main = "Sim")
  # pf_kud_1() approach without resampling produces a highly skewed picture
  out_acpf_pairs_ud_a <- pf_kud_1(out_acpf_pairs_unq,
                                  bathy = grid, sample_size = NULL,
                                  estimate_ud = adehabitatHR::kernelUD, grid = kud_grid_resolution,
                                  chunks = 10L, cl = parallel::makeCluster(10L))
  # pf_kud_1() approach with re-sampling similarly produces a highly skewed picture
  out_acpf_pairs_ud_b <- pf_kud_1(out_acpf_pairs_unq,
                                  bathy = grid, sample_size = 5000L,
                                  estimate_ud = adehabitatHR::kernelUD, grid = kud_grid_resolution,
                                  chunks = 10L, cl = parallel::makeCluster(10L))
  add_contour(out_acpf_pairs_ud_b)
  par(pp)


  ## Trial approaches based on paths
  # The truth
  pp <- par(mfrow = c(2, 2))
  prettyGraphics::pretty_map(add_rasters = list(x = path_ud), main = "Sim")
  # pf_kud_1() approach produces a highly skewed picture
  # ... This requires more unique cell samples at each time step than assembled by current routines
  cl <- parallel::makeCluster(10L)
  parallel::clusterEvalQ(cl, library(raster))
  out_acpf_pairs_ud_a <- pf_kud_1(out_acpf_paths,
                                  bathy = grid, sample_size = 5000L,
                                  estimate_ud = adehabitatHR::kernelUD, grid = kud_grid_resolution,
                                  chunks = 10L, cl = cl)
  # pf_kud_2() approach over all paths looks better
  out_acpf_pairs_ud_c <- pf_kud_2(out_acpf_paths, sample_size = 100L,
                                  bathy = grid,
                                  estimate_ud = adehabitatHR::kernelUD, grid = kud_grid_resolution)
  par(pp)
}


#### Apply KUD to particle samples and paths
pp <- par(mfrow = c(2, 2))
out_acpf_pairs_ud   <- pf_kud_2(out_acpf_pairs_unq,
                                bathy = grid, sample_size = NULL,
                                estimate_ud = adehabitatHR::kernelUD, grid = kud_grid_resolution)
out_acdcpf_pairs_ud <- pf_kud_2(out_acdcpf_pairs_unq,
                                bathy = grid, sample_size = NULL,
                                estimate_ud = adehabitatHR::kernelUD, grid = kud_grid_resolution)
out_acpf_paths_ud   <- pf_kud_2(out_acpf_paths,
                                bathy = grid, sample_size = NULL,
                                estimate_ud = adehabitatHR::kernelUD, grid = kud_grid_resolution)
out_acdcpf_paths_ud <- pf_kud_2(out_acdcpf_paths,
                                bathy = grid, sample_size = NULL,
                                estimate_ud = adehabitatHR::kernelUD, grid = kud_grid_resolution)
par(pp)


######################################
######################################
#### Examine algorithm performance

#### Compare the true path, COA estimates and AC, ACPF and ACDCPF outputs
save_png <- TRUE
if(save_png) png(paste0("./fig/path_", path_id, "_array_", array_id, ".png"),
                 height = 10, width = 7, units = "in", res = 600)
## Graphical param
pp <- par(mfrow = c(5, 3), oma = c(1, 1, 2, 1), mar = c(1, 0, 1, 0))
xlim <- ext[1:2]; ylim <- ext[3:4]
paa <- list(side = 1:4, axis = list(labels = FALSE))
cex_main <- 1.25
adj_1 <- 0.125
adj_2 <- 0.125
spaces <- "        "
## Make plots
# Simulated array and path
array$array$xy_buf <- rgeos::gBuffer(sp::SpatialPoints(array$array$xy), width = 300, quadsegs = 1000)
prettyGraphics::pretty_map(add_rasters = list(x = grid, plot_method = raster::plot, legend = FALSE),
                           # add_paths = list(x = path$xy_mat, lwd = 1, length = 0.025),
                           add_polys = list(x = array$array$xy_buf, pch = 21, border = "royalblue", lwd = 1),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE
                           )
mtext(side = 3, "A", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(array & path)"), adj = adj_2)
# Simulate path (over the range of acoustic data) and KUD
if(!is.null(path_ud)) add_path_ud <- list(x = path_ud, plot_method = raster::plot, legend = FALSE) else add_path_ud <- NULL
prettyGraphics::pretty_map(add_rasters = add_path_ud,
                           add_paths = list(x = path$xy_mat_on_grid_within_acoustics, lwd = 0.75, length = 0.025),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "B", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(path [det] & KUD)"), adj = adj_2)
# COAs and COA KUD
if(!is.null(out_coa_ud)) add_coa_ud <- list(x = out_coa_ud, plot_method = raster::plot, legend = FALSE) else add_coa_ud <- NULL
prettyGraphics::pretty_map(add_rasters = add_coa_ud,
                           add_points = list(x = out_coa),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "C", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(COAs & KUD)"), adj = adj_2)
# AC
prettyGraphics::pretty_map(add_rasters = list(x = out_ac_s$map, plot_method = raster::plot, legend = FALSE),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "D", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(AC)"), adj = adj_2)
# ACPF map (from particles)
pf_plot_map(out_acpf_pairs_unq,
            add_rasters = list(plot_method = raster::plot, legend = FALSE),
            map = grid, scale = "sum",
            xlim = xlim, ylim = ylim,
            pretty_axis_args = paa,
            crop_spatial = TRUE)
mtext(side = 3, "E", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACPF)"), adj = adj_2)
# ACPF + KUD (from particles)
prettyGraphics::pretty_map(add_rasters = list(x = out_acpf_pairs_ud, plot_method = raster::plot, legend = FALSE),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "F", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACPF KUD)"), adj = adj_2)
# ACPF path
path_eg <-  out_acpf_paths[out_acpf_paths$path_id == out_acpf_paths$path_id[1], ]
prettyGraphics::pretty_map(add_paths = list(x = path_eg$cell_x, path_eg$cell_y,
                                            col = viridis::viridis(nrow(path_eg)),
                                            length = 0.025, lwd = 0.75),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "G", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACPF path)"), adj = adj_2)
# ACPF map (from paths)
pf_plot_map(out_acpf_paths,
            add_rasters = list(plot_method = raster::plot, legend = FALSE),
            map = grid, scale = "sum",
            xlim = xlim, ylim = ylim,
            pretty_axis_args = paa,
            crop_spatial = TRUE)
mtext(side = 3, "H", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACPF [2])"), adj = adj_2)
# ACPF + KUD (from paths)
prettyGraphics::pretty_map(add_rasters = list(x = out_acpf_paths_ud, plot_method = raster::plot, legend = FALSE),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "I", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACPF KUD [2])"), adj = adj_2)
## ACDC
prettyGraphics::pretty_map(add_rasters = list(x = out_acdc_s$map, plot_method = raster::plot, legend = FALSE),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "J", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDC)"), adj = adj_2)
## ACDCPF
pf_plot_map(out_acdcpf_pairs_unq,
            add_rasters = list(plot_method = raster::plot, legend = FALSE),
            map = grid, scale = "sum",
            xlim = xlim, ylim = ylim,
            pretty_axis_args = paa,
            crop_spatial = TRUE)
mtext(side = 3, "K", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDCPF)"), adj = adj_2)
## ACDCPF + KUD
prettyGraphics::pretty_map(add_rasters = list(x = out_acdcpf_pairs_ud, plot_method = raster::plot, legend = FALSE),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "L", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDCPF KUD)"), adj = adj_2)
# ACDCPF path
path_eg <-  out_acdcpf_paths[out_acdcpf_paths$path_id == out_acdcpf_paths$path_id[1], ]
prettyGraphics::pretty_map(add_paths = list(x = path_eg$cell_x, path_eg$cell_y,
                                            col = viridis::viridis(nrow(path_eg)),
                                            length = 0.025, lwd = 0.75),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "M", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDCPF path)"), adj = adj_2)
# ACDCPF map (from paths)
pf_plot_map(out_acdcpf_paths,
            add_rasters = list(plot_method = raster::plot, legend = FALSE),
            map = grid, scale = "sum",
            xlim = xlim, ylim = ylim,
            pretty_axis_args = paa,
            crop_spatial = TRUE)
mtext(side = 3, "N", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDCPF [2])"), adj = adj_2)
# ACDCPF + KUD (from paths)
prettyGraphics::pretty_map(add_rasters = list(x = out_acdcpf_paths_ud, plot_method = raster::plot, legend = FALSE),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
mtext(side = 3, "O", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDCPF KUD [2])"), adj = adj_2)
par(pp)
if(save_png) dev.off()


#### End of code.
######################################
######################################
