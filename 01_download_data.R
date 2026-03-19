##' Download Coastal Targets
##' @return data.frame of coastal targets (raw)
download_targets <- function() {
  base_url <- "https://minio-s3.apps.shift.nerc.mghpcc.org/bu4cast-ci-read/challenges/project_id=bu4cast/targets/"
  coastal_url <- paste0(base_url, "coastal-targets.csv")
  readr::read_csv(coastal_url, show_col_types = FALSE)
}