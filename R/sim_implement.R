######################################
######################################
#### sim_implement.R

#### This code
# Implements the sim_workhorse.R script for all movement paths/array designs.

#### Steps preceding this code
source("./R/sim_data.R")
try(dev.off(), silent = TRUE)

#### Extensions
# 1) To improve speed, implement this script in parallel


######################################
######################################
#### Implementation

cat("\014")
cout <- paste0("./data/estimates/log_sim_implement_", as.numeric(Sys.time()), ".txt")
if (file.exists(cout)) unlink(cout)
sink(cout)
cat("log_sim_implement.txt")
ns     <- "\n \n"
spaces <- function() cat(ns)
breaker  <- function() cat(paste0(paste0(rep("-", 150), collapse = ""), "\n"))
t1 <- Sys.time()
lapply(1:length(dat_sim_paths), function(path_id){

  # Set up
  spaces()
  invisible(replicate(3, breaker()))
  cat(paste0("path ", path_id, ns))

  # (optional) Define paths to skip
  # if (!(path_id %in% 71:80)) return(NULL)

  # Loop over arrays
  array_ids <- 1:length(dat_sim_arrays)
  # array_ids <- 12L
  lapply(array_ids, function(array_id){

    # (optional) check array
    # raster::plot(dat_sim_arrays[[array_id]]$array$xy)

    # Skip simulations that failed to produce observations
    if(nrow(fails[fails$path == path_id & fails$array == array_id, ]) > 0L) {
      cat(paste("SKIP: path", path_id, "array", array_id, "!"))
      return(NULL)
    }

    # Implement simulations
    spaces()
    invisible(replicate(2, breaker()))
    cat(paste0("path ", path_id, ": array ", array_id, ns))
    lapply(alg_param, function(alg) {
      spaces()
      invisible(replicate(1, breaker()))
      cat(paste0("path ", path_id, ": array ", array_id, ": implementation ", alg$id, ns))
      tictoc::tic()
      source("./R/sim_workhorse.R", local = TRUE)
      spaces()
      tictoc::toc()
    })

  })
})
sink()
t2 <- Sys.time()
difftime(t2, t1, units = "mins")
beepr::beep(10)


#### End of code
######################################
######################################
