######################################
######################################
#### sim_evaluate.R

#### This code:
# 1) Evaluates outputs of the flapper family of algorithms for different array designs
# ... when implemented 'correctly' by comparing POU/KUD maps for the simulated path
# ... to those reconstructed by the COA & flapper family of algorithms

#### Steps preceding this code:
# 1) Implement sim_data.R (see below)
# 2) Implement sim_workhorse.R for each array design to generate outputs
# ... (e.g., via sim_implement.R)


######################################
######################################
#### Setup

#### Simulate data
source("./R/sim_data.R")

#### Define global param
path_id   <- 1
path      <- dat_sim_paths[[1]]
delta_t   <- "2 hours"
kud_grid  <- kud_habitat(grid)
# Define map scaling parameters: 'original' or 'max'
# (Deprecated) Use scale_* = "original" here to load 'raw' maps
# Scaling is implemented below in the plotting section instead
scale_pou <- "original"
scale_kud <- "original"
stopifnot(scale_pou == "original")
stopifnot(scale_kud == "original")

#### Check the number of acoustic detections for each array
sapply(dat_sim_detections_by_path[[path_id]], function(acc) nrow(acc))

#### POU and KUDs for each array
estimates_by_array_best_file <-
  paste0("./data/estimates/path_", path_id, "/estimates_by_array_best.rds")
if(!file.exists(estimates_by_array_best_file)){

  tictoc::tic()
  estimates_by_array <- lapply(1:length(dat_sim_arrays), function(array_id){

    #### Get array/path details
    # array_id <- 10
    con_root <- paste0("./data/estimates/path_", path_id, "/array_", array_id, "/", alg_imp_true, "/")
    print(array_id)

    #### Get associated data
    array      <- dat_sim_arrays[[array_id]]
    archival   <- dat_sim_archival_by_path[[path_id]]
    acoustics  <- dat_sim_detections_by_path[[path_id]][[array_id]]
    moorings   <- dat_sim_moorings[[array_id]]

    #### Process simulated time series
    acoustics$receiver_id <- as.integer(as.character(acoustics$receiver_id))
    archival$index <- 1:nrow(archival)
    archival <- archival[archival$timestamp >= min(acoustics$timestamp) &
                           archival$timestamp <= max(acoustics$timestamp), ]
    # Focus on parts of paths that detections enable us to resolve
    path$xy_mat_on_grid_within_acoustics <- path$xy_mat_on_grid[archival$index, ]

    #### POU for simulated path
    dat_sim_path_pf <-
      data.frame(path_id = 1,
                 timestep = 1:nrow(path$xy_mat_on_grid_within_acoustics),
                 cell_x = path$xy_mat_on_grid_within_acoustics[, 1],
                 cell_y = path$xy_mat_on_grid_within_acoustics[, 2],
                 cell_pr = 1)
    dat_sim_path_pf$cell_id <- raster::cellFromXY(grid, dat_sim_path_pf[, c("cell_x", "cell_y")])
    dat_sim_path_pf$cell_z  <- raster::extract(grid, dat_sim_path_pf$cell_id)
    dat_sim_path_pf <- dat_sim_path_pf[, c("path_id",
                                           "timestep",
                                           "cell_id",
                                           "cell_x", "cell_y", "cell_z",
                                           "cell_pr")]
    class(dat_sim_path_pf) <- c(class(dat_sim_path_pf), "pf_path")
    path_pou <- pf_plot_map(dat_sim_path_pf,
                            map = grid,
                            scale = scale_pou)

    #### KUD for simulated path
    path_kud <- NULL
    if(nrow(path$xy_mat_on_grid_within_acoustics) >= 5){
      path_spdf <- sp::SpatialPointsDataFrame(
        path$xy_mat_on_grid_within_acoustics,
        data = data.frame(ID = factor(rep(1, nrow(path$xy_mat_on_grid_within_acoustics)))),
        proj4string = raster::crs(grid))
      path_kud <- adehabitatHR::kernelUD(xy = path_spdf, grid = kud_grid)
      path_kud <- raster::raster(path_kud[[1]])
      path_kud <- raster::resample(path_kud, grid)
      path_kud <- path_kud/raster::cellStats(path_kud, "sum")
      if(scale_kud != "original") path_kud <- path_kud/raster::cellStats(path_kud, scale_kud)
      pretty_map(add_rasters = list(x = path_kud, interpolation = TRUE))
    }

    #### COA for simulated time series
    acoustics$individual_id <- 1L
    acoustics$receiver_id <- factor(acoustics$receiver_id, moorings$receiver_id)
    acoustics_mat <- make_matrix_detections(acoustics, delta_t = delta_t)
    out_coa <- coa(acoustics_mat, sp::coordinates(array$array$xy))
    out_coa <- data.frame(x = out_coa[, 1], y = out_coa[, 2])

    #### POU for COAs
    coa_path_pf <-
      data.frame(path_id = 1,
                 timestep = 1:nrow(out_coa),
                 cell_x = out_coa[, "x"],
                 cell_y = out_coa[, "y"],
                 cell_pr = 1)
    coa_path_pf$cell_id <- raster::cellFromXY(grid, coa_path_pf[, c("cell_x", "cell_y")])
    coa_path_pf$cell_z  <- raster::extract(grid, coa_path_pf$cell_id)
    coa_path_pf <- coa_path_pf[, c("path_id",
                                   "timestep",
                                   "cell_id",
                                   "cell_x", "cell_y", "cell_z",
                                   "cell_pr")]
    class(coa_path_pf) <- c(class(coa_path_pf), "pf_path")
    coa_pou <- pf_plot_map(coa_path_pf,
                           map = grid,
                           scale = scale_pou)

    #### KUD for COAs
    coa_kud <- NULL
    if(nrow(out_coa) >= 5){
      out_coa_spdf <- sp::SpatialPointsDataFrame(
        out_coa[, c("x", "y")],
        data = data.frame(ID = factor(rep(1, nrow(out_coa)))),
        proj4string = raster::crs(grid))
      coa_kud <- adehabitatHR::kernelUD(xy = out_coa_spdf, grid = kud_grid)
      coa_kud <- raster::raster(coa_kud[[1]])
      coa_kud <- raster::resample(coa_kud, grid)
      coa_kud <- coa_kud/raster::cellStats(coa_kud, "sum")
      if (scale_kud != "original") coa_kud <- coa_kud/raster::cellStats(coa_kud, scale_kud)
      if(all(is.na(raster::getValues(coa_kud)))){
        coa_kud <- NULL
      } else {
        pretty_map(add_rasters = list(x = coa_kud, interpolation = TRUE))
      }
    }

    #### ACPF and ACDCPF maps
    if(scale_pou == "original"){
      acpf_pou   <- raster::raster(paste0(con_root, "acpf/out_acpf_pou_raw.tif"))
      acdcpf_pou <- raster::raster(paste0(con_root, "acdcpf/out_acdcpf_pou_raw.tif"))
    } else if (scale_pou == "max"){
      acpf_pou   <- raster::raster(paste0(con_root, "acpf/out_acpf_pou.tif"))
      acdcpf_pou   <- raster::raster(paste0(con_root, "acdcpf/out_acdcpf_pou.tif"))
    }
    if(scale_kud == "original"){
      acpf_kud   <- raster::raster(paste0(con_root, "acpf/out_acpf_kud_raw.tif"))
      acdcpf_kud <- raster::raster(paste0(con_root, "acdcpf/out_acdcpf_kud_raw.tif"))
    } else if(scale_kud == "max"){
      acpf_kud   <- raster::raster(paste0(con_root, "acpf/out_acpf_kud.tif"))
      acdcpf_kud <- raster::raster(paste0(con_root, "acdcpf/out_acdcpf_kud.tif"))
    }

    #### Return list
    list(path = path, coa_xy = out_coa,
         pou = list(sim = path_pou, coa = coa_pou, acpf = acpf_pou, acdcpf = acdcpf_pou),
         kud = list(sim = path_kud, coa = coa_kud, acpf = acpf_kud, acdcpf = acdcpf_kud)
    )
  })
  tictoc::toc()
  saveRDS(estimates_by_array, estimates_by_array_best_file)
} else estimates_by_array <- readRDS(estimates_by_array_best_file)


######################################
######################################
#### Make plots

#### Choose whether or not to plot POUs or KUDs
type <- c("kud", "pou")
type <- type[1]

#### Choose whether or not to focus on a small selection of arrays
dat_sim_array_info_2 <- dat_sim_array_info
dat_sim_array_info_2$index <- 1:nrow(dat_sim_array_info_2)
dat_sim_array_info_2$arrangement <- factor(dat_sim_array_info_2$arrangement,
                                           levels = c("regular", "random", "clustered"),
                                           labels = c("Reg.", "Ran.", "Clu."))
dat_sim_array_info_2 <-
  dat_sim_array_info_2 %>%
  dplyr::arrange(arrangement, n_receivers, clustering, n_clusters)
dat_sim_array_info_2$label <- 1:nrow(dat_sim_array_info_2)
sbt <- c("sbt", "full")
sbt <- sbt[2]
if(sbt == "sbt") {
  array_ids <- dat_sim_array_info_2$index[which(dat_sim_array_info_2$n_receivers %in% c(5, 25) &
                                          (is.na(dat_sim_array_info_2$clustering) |
                                            dat_sim_array_info_2$clustering == 0.6))]
  dat_sim_array_info_2[dat_sim_array_info_2$index %in% array_ids, ]
  width <- 5; height <- length(array_ids)
} else {
  array_ids <- dat_sim_array_info_2$index
  array_ids <- array_ids[dat_sim_array_info_2$label]
  width <- 5; height <- 12
}

#### Set up image
save_png <- TRUE
if(save_png) png(paste0("./fig/evaluation/path_", path_id, "_array_evaluation_", type, "_", sbt, ".png"),
                 height = height, width = width, units = "in", res = 600)
pp <- par(mfrow = c(length(array_ids), width), oma = c(1, 1, 3, 1), mar = c(0, 0, 0, 0))
xlim <- ext[1:2]; ylim <- ext[3:4]
paa <- list(side = 1:4, axis = list(labels = FALSE, lwd.ticks = 0))
cex_main <- 1.25
cex_axis_title <- 0.8
add_title_1 <- list(side = 3, text = "", cex = cex_axis_title, line = 0)
add_title_2 <- list(side = 4, text = "", cex = cex_axis_title, line = -0.25)
atx <- c(xlim[1] + 70, xlim[1] + 1250)
atm <- c(0, 0, -200, 250)
add_textbox <- list(x = atx,
                    y = ylim[1] + 750,
                    textlist = "",
                    margin = atm, justify = "c", font = 2, cex = cex_main,
                    fill = scales::alpha("white", 0.80), border = scales::alpha("white", 0.80))

#### Define colour schemes
bathy_zlim <- c(0, 225)
bathy_col_param <- pretty_cols_brewer(bathy_zlim, scheme = "Blues", n_breaks = max(bathy_zlim))
bathy_cols <- bathy_col_param$col
add_bathy <- list(x = grid, plot_method = plot_raster_img, zlim = bathy_zlim, col = bathy_cols)
add_paths <- list(length = 0.01, lwd = 0.5)
if(type == "pou") add_paths$lwd <- 0.1

#### Make plots for each array
lapply(array_ids, function(array_id){

  #### Get array/path details
  # array_id <- 7
  print(array_id)
  array                  <- dat_sim_arrays[[array_id]]
  array_info             <- dat_sim_array_info_2[dat_sim_array_info_2$index == array_id, ]
  array_info$n_receivers <- length(array$array$xy)
  array_label            <- array_info$label
  estimates_for_array    <- estimates_by_array[[array_id]]
  path_for_array         <- estimates_for_array$path

  #### Get movement time series
  archival  <- dat_sim_archival_by_path[[path_id]]
  acoustics <- dat_sim_detections_by_path[[path_id]][[array_id]]
  archival  <- archival[archival$timestamp >= min(acoustics$timestamp) &
                          archival$timestamp <= max(acoustics$timestamp), ]

  #### Plot (1): array
  array$array$xy_buf <- rgeos::gBuffer(sp::SpatialPoints(array$array$xy), width = true_detection_range, quadsegs = 1000)
  prettyGraphics::pretty_map(add_rasters = add_bathy,
                             add_points = list(x = sp::coordinates(array$array$xy),
                                               pch = 21, col = "black", bg = "black", cex = 1),
                             add_polys = list(x = array$array$xy_buf, border = "black", lwd = 1.5, lty = 3),
                             xlim = xlim, ylim = ylim,
                             pretty_axis_args = paa,
                             crop_spatial = TRUE
                             )
  add_textbox$textlist <- paste0(array_label, "a")
  do.call(plotrix::textbox, add_textbox)
  if(is.na(array_info$clustering)){
    mtext(side = 2, text = paste0(array_info$arrangement, " (", array_info$n_receivers, ")"),
          cex = cex_axis_title, line = -0.2)
  } else {
    mtext(side = 2, text = paste0(array_info$arrangement, " (", array_info$n_receivers, " [", array_info$clustering, "])"),
          cex = cex_axis_title, line = -0.2)
  }
  if(array_id == array_ids[1]){
    add_title_1$text <- "Array"
    do.call(mtext, add_title_1)
  }

  #### Define algorithm plotting param
  # Get the maximum value across all rasters (for the current array) for scaling
  scale <- TRUE
  if(scale){
    maps <-
      list(white_out(estimates_for_array[[type]]$sim),
           white_out(estimates_for_array[[type]]$coa),
           white_out(estimates_for_array[[type]]$acpf),
           white_out(estimates_for_array[[type]]$acdcpf)) |>
      plyr::compact()
    mx <-
      sapply(maps, function(r){
        if(!is.null(r)) raster::cellStats(r, "max")
      }) |>
      max(na.rm = TRUE)
    zlim <- c(0, 1)
  } else {
    zlim <- NULL
  }
  map_param <- map_param_raw <- list(x = grid,
                                     add_rasters = list(x = NULL,
                                                        plot_method = raster::image,
                                                        zlim = zlim),
                                     xlim = xlim, ylim = ylim,
                                     pretty_axis_args = paa,
                                     crop_spatial = TRUE)

  #### Plot (2): path (between detections)
  rx <- white_out(estimates_for_array[[type]]$sim)
  if (!is.null(rx)) stopifnot(all.equal(1, raster::cellStats(rx, "sum")))
  if(scale) rx <- scale_raster(rx, mx)
  map_param$add_rasters$x <- rx
  add_paths$x             <- path_for_array$xy_mat_on_grid_within_acoustics
  map_param$add_paths     <- add_paths
  do.call(prettyGraphics::pretty_map, map_param)
  map_param$add_rasters$x <- NULL
  map_param$add_paths     <- NULL
  add_textbox$textlist    <- paste0(array_label, "b")
  do.call(plotrix::textbox, add_textbox)
  if(array_id == array_ids[1]){
    add_title_1$text <- "Path"
    do.call(mtext, add_title_1)
  }

  #### Plot (3): COA
  rr <- estimates_for_array[[type]]$coa
  rx <- white_out(rr)
  if (!is.null(rx)) stopifnot(all.equal(1, raster::cellStats(rx, "sum")))
  if(scale) rx <- scale_raster(rx, mx)
  map_param$add_rasters$x <- rx
  if(is.null(map_param$add_rasters$x)) map_param$add_rasters <- NULL
  map_param$add_points    <- list(x = estimates_for_array$coa_xy$x, estimates_for_array$coa_xy$y, pch = 17, bg = "black")
  do.call(prettyGraphics::pretty_map, map_param)
  if(type == "kud" && !is.null(rr)) add_contour(rr, ext = ext)
  map_param$add_rasters$x <- NULL
  map_param$add_points    <- NULL
  if(is.null(map_param$add_rasters)) map_param$add_rasters <- map_param_raw$add_rasters
  add_textbox$textlist <- paste0(array_label, "c")
  do.call(plotrix::textbox, add_textbox)
  if(array_id == array_ids[1]){
    add_title_1$text <- "COA"
    do.call(mtext, add_title_1)
  }

  #### Plot (4): ACPF
  map_param$add_rasters$plot_method <- plot_raster_img
  rr <- estimates_for_array[[type]]$acpf
  rx <- white_out(rr)
  if (!is.null(rx)) stopifnot(all.equal(1, raster::cellStats(rx, "sum")))
  if(scale) rx <- scale_raster(rx, mx)
  map_param$add_rasters$x <- rx
  do.call(prettyGraphics::pretty_map, map_param)
  add_contour(rr, ext = ext)
  map_param$add_rasters$x <- NULL
  add_textbox$textlist    <- paste0(array_label, "d")
  do.call(plotrix::textbox, add_textbox)
  if(array_id == array_ids[1]){
    add_title_1$text <- "ACPF"
    do.call(mtext, add_title_1)
  }

  #### Plot (5): ACDCPF
  rr <- estimates_for_array[[type]]$acdcpf
  rx <- white_out(rr)
  if (!is.null(rx))  stopifnot(all.equal(1, raster::cellStats(rx, "sum")))
  if(scale) rx <- scale_raster(rx, mx)
  map_param$add_rasters$x <- rx
  do.call(prettyGraphics::pretty_map, map_param)
  add_contour(rr, ext = ext)
  map_param$add_rasters$x <- NULL
  add_textbox$textlist    <- paste0(array_label, "e")
  do.call(plotrix::textbox, add_textbox)
  if(array_id == array_ids[1]){
    add_title_1$text <- "ACDCPF"
    do.call(mtext, add_title_1)
  }
  add_title_2$text <- paste0("(", nrow(acoustics), ", ", nrow(archival), ")")
  do.call(mtext, add_title_2)

  return(invisible())
}) %>% invisible()


#### Save figure
dev.off()


#### End of code.
######################################
######################################
