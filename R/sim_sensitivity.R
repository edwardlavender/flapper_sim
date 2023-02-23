######################################
######################################
#### sim_sensitivity.R

#### This code:
# 1) Examines the sensitivity of the flapper algorithms to the values
# ... of the detection probability and movement parameters
# ... by comparing the 'best' KUDs to 'mis-specified' KUDs
# ... for each array design and calculating supremum norm scores

#### Steps preceding this code:
# 1) Implement sim_data.R
# 2) Implement sim_workhorse.R for each array design to generate outputs
# 3) This script is based on sim_evaluate.R


######################################
######################################
#### Compute KUDs

# Suppress setup once run once for quick modification of figures
run_setup <- FALSE
if(run_setup){

  #### Simulate data
  source("./R/sim_data.R")

  #### Define global param
  path_id   <- 1

  #### Load data
  estimates_by_array <-
    readRDS(paste0("./data/estimates/path_", path_id, "/estimates_by_array_best.rds"))

  #### Maps for each array & algorithm implementation
  # This code defines a list, with one element for each array design
  # Each element is a list with one element for each algorithm implementation
  # ... that contains a list of POU/KUD maps for the ACPF and ACDCPF algorithms
  estimates_by_array_and_imp_file <-
    paste0("./data/estimates/path_", path_id, "/estimates_by_array_and_imp.rds")
  if(!file.exists(estimates_by_array_and_imp_file)){

    estimates_by_array_and_imp <- lapply(1:length(dat_sim_arrays), function(array_id){

      #### Define root
      # array_id <- 12L
      print(array_id)
      con_root <- paste0("./data/estimates/path_", path_id, "/array_", array_id, "/")

      #### Pull out the KUDs for the ACPF and ACDCPF algorithms for each algorithm implementation
      kud_by_alg <-
        pbapply::pblapply(alg_param, function(alg){
          # Define path
          # alg <- alg_param[[1]]
          con_root_alg <- paste0(con_root, "alg_imp_", alg$id, "/")
          # Read maps
          acpf_pou   <- read_raster_if_exists(paste0(con_root_alg, "acpf/out_acpf_pou_raw.tif"))
          acpf_kud   <- read_raster_if_exists(paste0(con_root_alg, "acpf/out_acpf_kud_raw.tif"))
          acdcpf_pou <- read_raster_if_exists(paste0(con_root_alg, "acdcpf/out_acdcpf_pou_raw.tif"))
          acdcpf_kud <- read_raster_if_exists(paste0(con_root_alg, "acdcpf/out_acdcpf_kud_raw.tif"))
          list(
            pou = list(acpf = acpf_pou, acdcpf = acdcpf_pou),
            kud = list(acpf = acpf_kud, acdcpf = acdcpf_kud)
          )
        })
      names(kud_by_alg) <- names(alg_param)
      kud_by_alg
    })
    saveRDS(estimates_by_array_and_imp, estimates_by_array_and_imp_file)
  } else {
    estimates_by_array_and_imp <- readRDS(estimates_by_array_and_imp_file)
  }

  #### Define a list of information for plotting
  # For each array, we will plot the the following:
  # * the array design
  # * the 'true' KUD
  # * the 'best' KUD
  # * a KUD based on an underestimated parameter value
  # * a KUD based on an overestimated parameter value

  #### Calculate supremum norms
  # We need a list, with one element for each array
  # And a subsequent list for each algorithm
  # And that list should contain a dataframe with the scores for each algorithm implementation
  # elm <- estimates_by_array_and_imp[[1]]
  # type <- "pou"
  # algorithm <- "acpf"
  alg_param <- do.call(rbind, alg_param)
  sup_by_array <-
    pbapply::pblapply(1:length(estimates_by_array_and_imp), function(i) {
      elm <- estimates_by_array_and_imp[[i]] # elm <- kud_by_alg
      best_maps <- elm[[alg_true]]
      by_type <-
        lapply(c("pou", "kud"), function(type){
          by_algorithm <-
            lapply(c("acpf", "acdcpf"), function(algorithm){
              sup <- sapply(elm, function(elm){
                best <- best_maps[[type]][[algorithm]]
                misp <- elm[[type]][[algorithm]]
                if (is.null(misp)) return(NA)
                raster::cellStats(abs(best - misp), "max")
              }, USE.NAMES = TRUE)
              d <- data.frame(array = i,
                              algorithm = algorithm,
                              id = names(sup),
                              sup = as.numeric(sup))
              ind      <- match(d$id, alg_param$id)
              d$dpr    <- alg_param$detection_range[ind]
              d$mpr    <- alg_param$mobility[ind]
              d$is_mpr <- alg_param$is_mpr[ind]
              d$is_dpr <- alg_param$is_dpr[ind]
              d
            })
          names(by_algorithm) <- c("acpf", "acdcpf")
          by_algorithm
        })
      names(by_type) <- c("pou", "kud")
      by_type
    })

}


######################################
######################################
#### Make plots

#### Set options
# POU versus KUDs
type <- c("kud", "pou")
type <- type[1]
# Subset/full selection of arrays
# ... set below
# Algorithm
algorithm <- c("acpf", "acdcpf")
algorithm <- algorithm[2]
# Variable
variable  <- c("dpr", "mpr")
variable  <- variable[2]
alg_param[, "dpr"] <- alg_param$detection_range
alg_param[, "mpr"] <- alg_param$mobility
if(variable == "dpr"){
  variable_under <- c("150_500")
  variable_over  <- c("600_500")
} else if(variable == "mpr"){
  variable_under <- c("300_250")
  variable_over  <- c("300_1000")
}

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
  width <- 6; height <- length(array_ids)
} else {
  array_ids <- dat_sim_array_info_2$index
  array_ids <- array_ids[dat_sim_array_info_2$label]
  width <- 6; height <- 12
}


#### Set up image
save_png <- TRUE
if(save_png) png(paste0("./fig/sensitivity/path_", path_id, "_", algorithm, "_", variable, "_", type, "_", sbt, ".png"),
                 height = height, width = width, units = "in", res = 600)
pp <- par(mfrow = c(length(array_ids), width), oma = c(1, 1, 3, 1), mar = c(0, 0, 0, 0))
xlim <- ext[1:2]; ylim <- ext[3:4]
paa <- list(side = 1:4, axis = list(labels = FALSE, lwd.ticks = 0))
cex_main <- 1.25
cex_axis_title <- 0.8
add_title_1 <- list(side = 3, text = "", cex = cex_axis_title, line = 0)
add_title_2 <- list(side = 4, text = "", cex = cex_axis_title, line = -0.1)
atx <- c(xlim[1] + 70, xlim[1] + 1250)
atp <- diff(atx)/diff(xlim)
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
  # array_id <- 5
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

  #### Define helper functions
  add_receivers_to_background <- function(r, id){
      do.call(graphics::points,
              list(x = sp::coordinates(array$array$xy),
                   pch = 21, cex = 0.3,
                   col = scales::alpha("dimgrey", 0.95), bg = scales::alpha("dimgrey", 0.95)))
  }
  add_scale_bar <- function(id){
    dist <- alg_param[alg_param$id == id, variable]
    raster::scalebar(d = dist,
                     xy = c(xlim[2] - 125 - dist, ylim[1] + 125),
                     label = "", lonlat = FALSE, lwd = 2)
  }

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
      list(estimates_for_array[[type]]$sim,
           estimates_for_array[[type]][[algorithm]],
           estimates_by_array_and_imp[[array_id]][[variable_under]][[type]][[algorithm]],
           estimates_by_array_and_imp[[array_id]][[variable_over]][[type]][[algorithm]]
      ) |>
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
  rr <- estimates_for_array[[type]]$sim
  rx <- white_out(rr)
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

  #### Plot (3): Algorithm (best)
  map_param$add_rasters$plot_method <- plot_raster_img
  rr <- estimates_for_array[[type]][[algorithm]]
  rx <- white_out(rr)
  if(scale) rx <- scale_raster(rx, mx)
  map_param$add_rasters$x <- rx
  do.call(prettyGraphics::pretty_map, map_param)
  add_contour(rr, ext = ext)
  add_receivers_to_background(map_param$add_rasters$x, id = alg_true)
  add_scale_bar(alg_true)
  map_param$add_rasters$x <- NULL
  add_textbox$textlist    <- paste0(array_label, "c")
  do.call(plotrix::textbox, add_textbox)
  if(array_id == array_ids[1]){
    add_title_1$text <- "Best"
    do.call(mtext, add_title_1)
  }

  #### Plot (4): Algorithm (underestimate)
  map_param$add_rasters$plot_method <- plot_raster_img
  rr <- estimates_by_array_and_imp[[array_id]][[variable_under]][[type]][[algorithm]]
  rx <- white_out(rr)
  if(scale) rx <- scale_raster(rx, mx)
  map_param$add_rasters$x <- rx
  do.call(prettyGraphics::pretty_map, map_param)
  add_contour(rr, ext = ext)
  add_receivers_to_background(map_param$add_rasters$x, variable_under)
  add_scale_bar(variable_under)
  map_param$add_rasters$x <- NULL
  add_textbox$textlist    <- paste0(array_label, "d")
  do.call(plotrix::textbox, add_textbox)
  if(array_id == array_ids[1]){
    add_title_1$text <- "Under"
    do.call(mtext, add_title_1)
  }

  #### Plot (5): Algorithm (overestimate)
  map_param$add_rasters$plot_method <- plot_raster_img
  rr <- estimates_by_array_and_imp[[array_id]][[variable_over]][[type]][[algorithm]]
  rx <- white_out(rr)
  if(scale) rx <- scale_raster(rx, mx)
  map_param$add_rasters$x <- rx
  do.call(prettyGraphics::pretty_map, map_param)
  add_contour(rr, ext = ext)
  add_receivers_to_background(map_param$add_rasters$x, variable_over)
  add_scale_bar(variable_over)
  map_param$add_rasters$x <- NULL
  add_textbox$textlist    <- paste0(array_label, "e")
  do.call(plotrix::textbox, add_textbox)
  if(array_id == array_ids[1]){
    add_title_1$text <- "Over"
    do.call(mtext, add_title_1)
  }

  #### Plot (6): Supremum norm
  # Pull out values
  sup <- sup_by_array[[array_id]][[type]][[algorithm]]
  if (variable == "mpr") {
    sup <- sup[sup$is_mpr, ]
    true_param <- true_mobility
  } else if (variable == "dpr") {
    sup <- sup[sup$is_dpr, ]
    true_param <- true_detection_range
  }
  sup <- sup[order(sup[, variable]), ]
  # Generate tidier tables (*100)
  sup$sup  <- sup$sup * 1e2
  sup$labels <- ""
  pos_labels <- sup[, variable] >= 0.25 * true_param
  sup[pos_labels, "labels"] <- sup[pos_labels, variable]
  sup$cex <- 0.35
  sup[pos_labels, "cex"] <- 0.7
  eps <- 1e-5
  yaxis <- pretty_seq(c(0, max(sup$sup, na.rm = TRUE) + eps), pretty_args = list(n = 5))
  if(length(yaxis$at) <= 4L){
    yaxis <- pretty_seq(c(0, max(sup$sup, na.rm = TRUE) + eps), pretty_args = list(n = 6))
  }
  stopifnot(length(yaxis$at) > 2L)
  ylab <- prettyGraphics:::pretty_labels(yaxis$at, yaxis$at)
  pos_labels <- seq(1, length(ylab), 2)
  ylab[pos_labels] <- ""
  ylab[length(ylab)] <- ""
  if(length(which(ylab != "")) < 2L) warning(paste("For array", array_id, "there are Less than two y axis labels!"))
  pretty_plot(sup[, variable], sup$sup,
              pretty_axis_args = list(x = list(x = range(c(0, sup[, variable])),
                                               y = range(yaxis$at)),
                                      axis = list(list(at = sup[, variable],
                                                       labels = sup$labels,
                                                       mgp = c(3, -0.3, 0)),
                                                  list(at = yaxis$at,
                                                       labels = ylab,
                                                       mgp = c(3, -0.5 - 0.3 * max(ndp(yaxis$at)) , 0))),
                                      control_axis = list(tck = 0.02, las = TRUE, cex.axis = 0.5)),
              type = "n",
              xlab = "", ylab = "")
  lines(sup[, variable], sup$sup, col = scales::alpha("grey", 0.7))
  points(sup[, variable], sup$sup, pch = 21,
         bg = scales::alpha("black", 0.75),
         col = NA,
         cex = sup$cex)
  add_textbox_scatter <- add_textbox
  legend_x <- "bottomright"
  if(algorithm == "acpf" & variable == "mpr") legend_x <- "topright"
  legend(legend_x, legend = paste0(array_label, "f"),
         bty = "n", text.font = 2, cex = add_textbox$cex, adj = c(0, 0.625)) # larger adj[2], lower label
  if(array_id == array_ids[1]){
    add_title_1$text <- "Sup"
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
