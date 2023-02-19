######################################
######################################
#### sim_illustrate.R

#### This code:
# 1) Illustrates example outputs of the flapper family of algorithms

#### Steps preceding this code:
# 1) Implement sim_data.R


######################################
######################################
#### Setup

#### Run sim_data.R
source("./R/sim_data.R")

run <- TRUE # suppress setup once it has been run once for quick editing of plots (below)
if(run){


  ######################################
  #### Get simulated data

  #### Get array/path details
  path_id  <- 1; array_id <- 12
  con_root <- paste0("./data/estimates/path_", path_id, "/array_", array_id, "/", alg_imp_true, "/")

  #### Get associated data
  array      <- dat_sim_arrays[[array_id]]
  moorings   <- dat_sim_moorings[[array_id]]
  archival   <- dat_sim_archival_by_path[[path_id]]
  acoustics  <- dat_sim_detections_by_path[[path_id]][[array_id]]
  path       <- dat_sim_paths[[path_id]]

  #### Global param
  kud_grid <- kud_habitat(grid)


  ######################################
  #### Get algorithm results

  #### Scale outputs to a maximum value of one
  # This means that all maps can be generated using one scale bar
  # This has no discernible influence on the maps.
  scale <- TRUE

  #### DC algorithm(s)
  ## DC algorithm (0.07 mins)
  out_dc <- dc(archival = archival,
               bathy = grid,
               normalise = TRUE,
               calc_depth_error = function(...) matrix(c(-5, 5), nrow = 2),
               save_record_spatial = 1:30L)
  out_dc_s <- acdc_simplify(out_dc, type = "dc")
  if(scale) out_dc_s$map <- scale_raster(out_dc_s$map)
  ## DCPF algorithm (0.02 mins)
  out_dc_record <- acdc_access_maps(out_dc_s)
  out_dcpf <- pf(out_dc_record,
                 bathy = grid,
                 origin = path$xy_mat_on_grid[1, , drop = FALSE],
                 mobility_from_origin = 1,
                 mobility = alg_param[[alg_true]]$mobility_on_grid,
                 calc_movement_pr = calc_mprw(mobility = alg_param[[alg_true]]$mobility))
  # Define paths (0.03 mins)
  out_dcpf_paths <- pf_simplify(out_dcpf, max_n_paths = 1000)
  # Define a subset of paths
  out_dcpf_paths_ll  <- pf_loglik(out_dcpf_paths)
  out_dcpf_paths_sbt <- out_dcpf_paths[out_dcpf_paths$path_id %in% out_dcpf_paths_ll$path_id[1], ]
  # Interpret LCPs (0.01 mins)
  out_dcpf_lcps <- lcp_interp(paths = out_dcpf_paths_sbt,
                              surface = grid,
                              calc_distance = FALSE)

  #### AC* algorithms
  # AC algorithm
  out_ac   <- readRDS(paste0(con_root, "ac/out_ac.rds"))
  out_ac_s <- acdc_simplify(out_ac)
  if(scale) out_ac_s$map <- scale_raster(out_ac_s$map)
  out_acpf_pairs <- readRDS(paste0(con_root, "acpf/out_acpf_pairs.rds"))
  out_acpf_pairs_unq <- pf_simplify(out_acpf_pairs, summarise_pr = TRUE, return = "archive")
  out_acpf_pairs_pou <- pf_plot_map(out_acpf_pairs_unq, map = grid)
  if(scale) out_acpf_pairs_pou <- scale_raster(out_acpf_pairs_pou)
  out_acpf_pairs_ud  <- pf_kud(pf_plot_map(out_acpf_pairs_unq, map = grid),
                               sample_size = 100,
                               estimate_ud = adehabitatHR::kernelUD,
                               grid = kud_grid)
  if(scale) out_acpf_pairs_ud <- scale_raster(out_acpf_pairs_ud)
  out_acpf_paths_file <- paste0(con_root, "acpf/out_acpf_paths.rds")
  if(!file.exists(out_acpf_paths_file)) {
    out_acpf_paths <- pf_simplify(out_acpf_pairs,
                                  bathy = grid,
                                  max_n_copies = NULL,
                                  max_n_paths = 1000L,
                                  return = "path"
    )
    saveRDS(out_acpf_paths, out_acpf_paths_file)
  } else out_acpf_paths <- readRDS(out_acpf_paths_file)
  out_acpf_paths_ll  <- pf_loglik(out_acpf_paths)
  out_acpf_paths_sbt <-  out_acpf_paths[out_acpf_paths$path_id %in% out_acpf_paths_ll$path_id[1], ]
  # ACDC algorithm(s)
  out_acdc   <- readRDS(paste0(con_root, "acdc/out_acdc.rds"))
  out_acdc_s <- acdc_simplify(out_acdc)
  if(scale) out_acdc_s$map <- scale_raster(out_acdc_s$map)
  out_acdcpf_pairs     <- readRDS(paste0(con_root, "acdcpf/out_acdcpf_pairs.rds"))
  out_acdcpf_pairs_unq <- pf_simplify(out_acdcpf_pairs, summarise_pr = TRUE, return = "archive")
  out_acdcpf_pairs_pou <- pf_plot_map(out_acdcpf_pairs_unq, map = grid)
  if(scale) out_acdcpf_pairs_pou <- scale_raster(out_acdcpf_pairs_pou)
  out_acdcpf_pairs_ud  <- pf_kud(pf_plot_map(out_acdcpf_pairs_unq, map = grid),
                                 sample_size = 100,
                                 estimate_ud = adehabitatHR::kernelUD,
                                 grid = kud_grid)
  if(scale) out_acdcpf_pairs_ud <- scale_raster(out_acdcpf_pairs_ud)
  out_acdcpf_paths_file <- paste0(con_root, "acdcpf/out_acdcpf_paths.rds")
  if(!file.exists(out_acdcpf_paths_file)) {
    out_acdcpf_paths <- pf_simplify(out_acdcpf_pairs,
                                  bathy = grid,
                                  max_n_copies = NULL,
                                  max_n_paths = 1000L,
                                  return = "path"
    )
    saveRDS(out_acdcpf_paths, out_acdcpf_paths_file)
  } else out_acdcpf_paths <- readRDS(out_acdcpf_paths_file)

  out_acdcpf_paths_ll  <- pf_loglik(out_acdcpf_paths)
  out_acdcpf_paths_sbt <-
    out_acdcpf_paths[out_acdcpf_paths$path_id %in% out_acdcpf_paths_ll$path_id[1], ]
}


######################################
######################################
#### Visualise algorithms

#### Comments
# The code below uses layout() to define a plot layout
# Rasters are plotted using raster::image() rather than raster::plot() to accommodate this implementation.

### Set up figure to save
save_png <- TRUE
if(save_png) png(paste0("./fig/path_", path_id, "_array_", array_id,  "_illustration.png"),
                 height = 10, width = 12, units = "in", res = 600)


### Define graphical param
mat <- matrix(c(1, 1, 3, 5, 9,
                1, 1, 4, 6, 10,
                2, 2, 2, 7, 11,
                2, 2, 2, 8, 12), ncol = 5, nrow = 4, byrow = TRUE)
layout(mat)
# layout.show(12)
pp <- par(oma = c(2, 5, 2, 8), mar = c(1, 0, 1, 0))
ext <- raster::extent(grid)
xlim <- ext[1:2]; ylim <- ext[3:4]
paa <- list(side = 1:4, axis = list(labels = FALSE))
cex_main <- 1.25
cex_axis <- cex_main
adj_1 <- 0.025
adj_2 <- 0.025
spaces <- "    "

#### Define colour schemes
bathy_zlim <- c(0, 225)
bathy_col_param <- pretty_cols_brewer(bathy_zlim, scheme = "Blues", n_breaks = max(bathy_zlim))
bathy_cols <- bathy_col_param$col
add_bathy <- list(x = grid, plot_method = plot_raster_img, zlim = bathy_zlim, col = bathy_cols)
add_paths <- list(x = NULL, col = viridis::viridis(nrow(path$xy_mat_on_grid)), lwd = 2, length = 0.05)

#### Plot the simulated array and path
array$array$xy_buf <- rgeos::gBuffer(sp::SpatialPoints(array$array$xy), width = 300, quadsegs = 1000)
add_paths$x <- path$xy_mat_on_grid
pretty_map(add_rasters = add_bathy,
           add_paths = add_paths,
           add_points = list(x = sp::coordinates(array$array$xy),
                             pch = 21, col = "black", bg = "black", cex = 1.5),
           add_polys = list(x = array$array$xy_buf, border = "black", lwd = 1.5, lty = 3),
           xlim = xlim, ylim = ylim,
           pretty_axis_args = paa,
           crop_spatial = TRUE
)
mtext(side = 3, "A", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(array & path)"), adj = adj_2, cex = cex_main)

#### Plot the 'observed' movement time series
## Depth time series
depth_scale <- data.frame(depth = bathy_col_param$breaks[1:(length(bathy_col_param$breaks) - 1)],
                          col = bathy_cols)
archival$breaks <- cut(archival$depth, depth_scale$depth, labels = FALSE)
archival$col <- depth_scale$col[archival$breaks]
axis_ls <- pretty_plot(archival$timestamp, archival$depth * - 1,
                       pretty_axis_args = list(side = 3:2,
                                               axis = list(list(labels = FALSE), list())),
                       cex.axis = cex_axis + 0.1,
                       xlab = "", ylab = "",
                       type = "n")
n <- nrow(archival)
x <- archival$timestamp
y <- archival$depth * -1
arrows(x0 = x[1:(n - 1)],
       x1 = x[2:n],
       y0 = y[1:(n - 1)],
       y1 =  y[2:n],
       col = archival$col,
       length = 0,
       lwd = 3)
points(x, y, pch = 21, col = archival$col, bg = archival$col, cex = 1)
## Acoustic time series
pretty_line(acoustics$timestamp, pretty_axis_args = list(axis_ls = axis_ls), inherit = 1L,
            pch = 21, bg = "black", cex = 2,
            add = TRUE)
## Labels
mtext(side = 2, "Depth (m)", cex = cex_main, line = 2.5)
mtext(side = 3, "B", adj = adj_1 - 0.01, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(movement time series)"), adj = adj_2 - 0.01, cex = cex_main)


#### Plot DC algorithm(s)
## DC
pretty_map(add_rasters = list(x = white_out(out_dc_s$map), plot_method = plot_raster_img),
           pretty_axis_args = paa)
add_contour(out_dc_s$map)
mtext(side = 3, "C", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(DC)"), adj = adj_2, cex = cex_main)
## DCPF
ext_dcpf <-
  raster::extent(
    rgeos::gBuffer(sp::SpatialPoints(out_dcpf_paths_sbt[, c("cell_x", "cell_y")]), width = 500))
pf_plot_2d(out_dcpf_paths_sbt,
           bathy = grid,
           add_bathy = add_bathy,
           add_paths = add_paths,
           xlim = ext_dcpf[1:2], ylim = ext_dcpf[3:4],
           pretty_axis_args = paa,
           crop_spatial = TRUE)
mtext(side = 3, "F", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(DCPF)"), adj = adj_2, cex = cex_main)

#### Plot AC algorithm(s)
## AC
prettyGraphics::pretty_map(add_rasters = list(x = white_out(out_ac_s$map),
                                              plot_method = plot_raster_img),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
add_contour(out_ac_s$map)
mtext(side = 3, "D", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(AC)"), adj = adj_2, cex = cex_main)
## AC path
pf_plot_2d(out_acpf_paths_sbt,
           bathy = grid,
           add_bathy = add_bathy,
           add_paths = add_paths,
           pretty_axis_args = paa,
           crop_spatial = TRUE)
mtext(side = 3, "G", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACPF [path])"), adj = adj_2, cex = cex_main)
## AC particles
prettyGraphics::pretty_map(add_rasters = list(x = white_out(out_acpf_pairs_pou),
                                              plot_method = plot_raster_img),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
add_contour(out_acpf_pairs_pou)
mtext(side = 3, "I", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACPF [POU])"), adj = adj_2, cex = cex_main)
## AC UD
prettyGraphics::pretty_map(add_rasters = list(x = white_out(out_acpf_pairs_ud),
                                              plot_method = plot_raster_img),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
add_contour(out_acpf_pairs_ud)
mtext(side = 3, "K", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACPF [KUD])"), adj = adj_2, cex = cex_main)

#### Plot ACDC algorithm(s)
## ACDC
prettyGraphics::pretty_map(add_rasters = list(x = white_out(out_acdc_s$map),
                                              plot_method = plot_raster_img),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
add_contour(out_acdc_s$map)
mtext(side = 3, "E", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDC)"), adj = adj_2, cex = cex_main)
## ACDC path
pf_plot_2d(out_acdcpf_paths_sbt,
           bathy = grid,
           add_bathy = add_bathy,
           add_paths = add_paths,
           pretty_axis_args = paa,
           crop_spatial = TRUE)
mtext(side = 3, "H", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDCPF [path])"), adj = adj_2, cex = cex_main)
## ACDC particles
prettyGraphics::pretty_map(add_rasters = list(x = white_out(out_acdcpf_pairs_pou),
                                              plot_method = plot_raster_img),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
add_contour(out_acdcpf_pairs_pou)
mtext(side = 3, "J", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDCPF [POU])"), adj = adj_2, cex = cex_main)
## ACDC UD
prettyGraphics::pretty_map(add_rasters = list(x = white_out(out_acdcpf_pairs_ud),
                                              plot_method = plot_raster_img),
                           xlim = xlim, ylim = ylim,
                           pretty_axis_args = paa,
                           crop_spatial = TRUE)
add_contour(out_acdcpf_pairs_ud)
mtext(side = 3, "L", adj = adj_1, font = 2, cex = cex_main)
mtext(side = 3, paste0(spaces, "(ACDCPF [KUD])"), adj = adj_2, cex = cex_main)

#### Close graphical device
par(pp)
dev.off()


#### Plot legends [combined with figures above manually using Preview on Mac]
## Bathymetry
png(paste0("./fig/path_", path_id, "_array_", array_id,  "_illustration_legend_1.png"),
    height = 5, width = 3, units = "in", res = 600)
pp <- par(oma = c(1, 1, 1, 4))
fields::image.plot(grid, zlim = sort(bathy_zlim*-1), legend.only = TRUE, col = rev(bathy_cols))
mtext(side = 4, "Depth (m)", cex = cex_main, line = 2.5)
dev.off()
## Time
png(paste0("./fig/path_", path_id, "_array_", array_id,  "_illustration_legend_2.png"),
    height = 5, width = 3, units = "in", res = 600)
pp <- par(oma = c(1, 1, 1, 4))
fields::image.plot(grid, zlim = c(0, 500), legend.only = TRUE, col = viridis::viridis(501))
mtext(side = 4, "Time (steps)", cex = cex_main, line = 2.5)
dev.off()
## Surfaces
png(paste0("./fig/path_", path_id, "_array_", array_id,  "_illustration_legend_3.png"),
    height = 5, width = 3, units = "in", res = 600)
pp <- par(oma = c(1, 1, 1, 4))
fields::image.plot(grid, zlim = c(0, 1), legend.only = TRUE, col = rev(grDevices::terrain.colors(255)))
mtext(side = 4, "Score", cex = cex_main, line = 2.5)
dev.off()


#### End of code.
######################################
######################################
