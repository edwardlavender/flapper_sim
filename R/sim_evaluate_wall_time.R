######################################
######################################
#### sim_evaluate_wall_time.R

#### This code
# Calculates the wall time for simulations used for algorithm evaluation.

#### Steps preceding this code
# 1) Run sim_data.R
# 2) Run flapper algorithms via sim_workhorse.R


######################################
######################################
#### Estimate wall times

#### Define 'dat_sim_array_info_2' via sim_evaluate.R
# ... okay.

#### Define the path ID.
path_id <- 1L

#### Get a dataframe of wall times
sim_wall_time <- lapply(split(dat_sim_array_info_2, 1:nrow(dat_sim_array_info_2)), function(array){

  #### Get param
  # array <- dat_sim_array_info_2[1, , drop = FALSE]
  array_id <- array$index

  #### AC time
  con <- paste0("./data/estimates/path_", path_id, "/array_", array_id, "/", alg_imp_true, "/ac/acdc_log.txt")
  txt <- readLines(con)
  array$ac_time <- readr::parse_number(stringr::str_split_fixed(txt[length(txt)], "~", n = 2)[, 2])

  #### ACDC time
  con <- paste0("./data/estimates/path_", path_id, "/array_", array_id, "/", alg_imp_true, "/acdc/acdc_log.txt")
  txt <- readLines(con)
  array$acdc_time <- readr::parse_number(stringr::str_split_fixed(txt[length(txt)], "~", n = 2)[, 2])

  #### ACPF time
  con <- paste0("./data/estimates/path_", path_id, "/array_", array_id, "/", alg_imp_true, "/acpf/acpf_log.txt")
  txt <- readLines(con)
  array$acpf_time <- readr::parse_number(stringr::str_split_fixed(txt[length(txt)], "~", n = 2)[, 2])

  #### ACDCPF time
  con <- paste0("./data/estimates/path_", path_id, "/array_", array_id, "/", alg_imp_true, "/acdcpf/acdcpf_log.txt")
  txt <- readLines(con)
  array$acdcpf_time <- readr::parse_number(stringr::str_split_fixed(txt[length(txt)], "~", n = 2)[, 2])

  #### Return data
  return(array)
}) %>% dplyr::bind_rows()

#### Calculate summary statistics
utils.add::basic_stats(unlist(sim_wall_time[, c("ac_time", "acdc_time")]))
utils.add::basic_stats(unlist(sim_wall_time[, c("acpf_time", "acdcpf_time")]))

#### Create a tidy table of wall times
sim_wall_time <-
  sim_wall_time %>%
  dplyr::mutate(ac_time     = add_lagging_point_zero(ac_time, n = 2),
                acdc_time   = add_lagging_point_zero(acdc_time, n = 2),
                acpf_time   = add_lagging_point_zero(acpf_time, n = 2),
                acdcpf_time = add_lagging_point_zero(acdcpf_time, n = 2)) %>%
  dplyr::select(ID = label,
                AC = ac_time,
                ACDC = acdc_time,
                ACPF = acpf_time,
                ACDCPF = acdcpf_time) %>% dplyr::arrange(ID)

#### Save wall times to file.
write.table(sim_wall_time, "./fig/sim_evaluate_wall_time.txt", quote = FALSE, sep = ",", row.names = FALSE)


#### End of code.
######################################
######################################
