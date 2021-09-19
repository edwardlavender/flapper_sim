######################################
######################################
#### sim_implement.R

#### This code
# Implements the sim_workhorse script for all movement paths/array designs.

#### Steps preceding this code
# Implement sim_data.R


######################################
######################################
#### Implementation

sink("./data/estimates/log_sim_implement.txt")
cat("log_sim_implement.txt")
breaker <- "-----------------------------------------------------------------------\n"
t1 <- Sys.time()
lapply(1:length(dat_sim_paths), function(path_id){

  cat("\n \n \n")
  cat(breaker); cat(breaker); cat(breaker); cat(breaker); cat(breaker);
  cat(paste("path", path_id, "\n \n \n"))

  lapply(1:length(dat_sim_arrays), function(array_id){

    cat("\n \n \n")
    cat(breaker); cat(breaker);
    cat(paste("array", array_id, "\n \n \n"))
    source("./R/sim_workhorse.R", local = TRUE)
  })
})
sink()
t2 <- Sys.time()
difftime(t2, t1, units = "mins")

#### End of code
######################################
######################################
