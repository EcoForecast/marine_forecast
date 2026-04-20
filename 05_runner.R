suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(rjags)
  library(coda)
  library(neon4cast)
})

source("01_download_data.R")
source("02_process_data.R")
source("03_make_visuals.R")
source("05_EnKF.R")

# Download and process targets
targets_raw <- download_targets()
targets_processed <- process_targets(targets_raw)

# Example run
res_daily <- run_iterative_enkf_daily(
  targets_processed = targets_processed,
  site_select       = 1,
  response_var      = "chlora_cci",
  forecast_date     = Sys.Date(),
  horizon_days      = 30L,
  n_members         = 31L,
  out_dir           = "outputs/advanced_ssm",
  model_id          = "bu4cast_coastal_ssm",
  force_refit       = FALSE,
  min_obs           = 100,
  use_last_forecast_day_for_restart = TRUE,
  do_submit         = TRUE,
  task              = "aquatics"
)
