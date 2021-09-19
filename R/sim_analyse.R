######################################
######################################
#### sim_analyse.R

#### This code
# Analyses simulated movement patterns using multiple methods
# ... designed for particle samples or reconstructed paths to
# ... examine their performance for real data.

#### Steps preceding this code
# Implementation of sim_data.R


######################################
######################################
#### Set up

#### Define simulated path [[1]] as a pf_path object
# ... so that we can use flapper functions to analyse the path
dat_sim_path_pf <-
  data.frame(path_id = 1,
             timestep = 1:nrow(dat_sim_paths[[1]]$xy_mat_on_grid),
             cell_x = dat_sim_paths[[1]]$xy_mat_on_grid[, 1],
             cell_y = dat_sim_paths[[1]]$xy_mat_on_grid[, 2],
             cell_pr = 1)
dat_sim_path_pf$cell_id <- raster::cellFromXY(grid, dat_sim_path_pf[, c("cell_x", "cell_y")])
dat_sim_path_pf$cell_z  <- raster::extract(grid, dat_sim_path_pf$cell_id)
dat_sim_path_pf <- dat_sim_path_pf[, c("path_id",
                                       "timestep",
                                       "cell_id",
                                       "cell_x", "cell_y", "cell_z",
                                       "cell_pr")]
class(dat_sim_path_pf) <- c(class(dat_sim_path_pf), "pf_path")


######################################
######################################
#### Probability of use maps

#### Plot `probability of use' for the simulated path
pf_plot_map(xpf = dat_sim_path_pf,
            map = grid)
add_sp_path(dat_sim_path_pf$cell_x, dat_sim_path_pf$cell_y, length = 0.01, lwd = 0.1)


######################################
######################################
#### Kud maps

#### Plotting window
pp <- par(mfrow = c(2, 3))

#### Visualise the 'true' KUD based on pf_kud_1()
# This approach cannot be implemented for simulated paths because there is only
# ... one location at each time step.

#### Visualise the 'true' KUD based on a few sampled locations based on pf_kud_2()
pf_kud_2(xpf = dat_sim_path_pf,
         bathy = grid,
         sample_size = 50,
         grid = 120)

#### Visualise 'true' KUD via pf_kud_2 for all particles
pf_kud_2(xpf = dat_sim_path_pf,
         bathy = grid,
         sample_size = NULL,
         grid = 120)

#### Visualise 'true' KUD via pf_kud_2 for more particles
pf_kud_2(xpf = dat_sim_path_pf,
         bathy = grid,
         sample_size = 5000,
         grid = 120)
pf_kud_2(xpf = dat_sim_path_pf,
         bathy = grid,
         sample_size = 50000,
         grid = 120)

#### Visualise 'true' KUD with less smoothing
ud <- adehabitatHR::kernelUD(sp::SpatialPoints(dat_sim_path_pf[, c("cell_x", "cell_y")]))
ud@h$h
ud <- adehabitatHR::kernelUD(sp::SpatialPoints(dat_sim_path_pf[, c("cell_x", "cell_y")]),
                             h = 100)
raster::plot(raster::raster(ud))
par(pp)


#### Results
# The 'true' KUD based on more sampled locations becomes more pixelated
# ... and actually better represents the simulated movement path
# ... by increasing the weight applied to locations with data and
# ... decreasing the weight applied to locations without data
# ... This looks equivalent to decreasing the smoothing parameter
# ... However, unlike for the smoothing parameter for which, in theory, there are principled
# ... ways of deciding how much smoothing to do, there isn't an obvious, principled way to decide
# ... how much resampling to do.


#### End of code.
######################################
######################################
