##' Process Coastal Targets
##' @param target data.frame from download_targets()
##' @return list with:
##'   - target_long_daily: daily-aggregated long data (datetime, site_id, variable, observation)
##'   - target_wide_daily: daily-aggregated wide data (datetime, site_id, <variables...>)
process_targets <- function(target) {
  
  # ---- Validate required columns ----
  req <- c("site_id", "datetime", "variable", "observation")
  if (!all(req %in% names(target))) {
    stop(
      "target is missing required columns: ",
      paste(setdiff(req, names(target)), collapse = ", ")
    )
  }
  
  # ---- Aggregate to daily (robust even if input isn’t perfectly daily) ----
  target_long_daily <- target %>%
    dplyr::mutate(datetime = as.Date(datetime)) %>%
    dplyr::group_by(datetime, site_id, variable) %>%
    dplyr::summarise(observation = mean(observation, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(observation = ifelse(is.nan(observation), NA, observation)) %>%
    dplyr::arrange(site_id, datetime, variable)
  
  # ---- Pivot wide: one row per site/day ----
  target_wide_daily <- target_long_daily %>%
    tidyr::pivot_wider(names_from = variable, values_from = observation) %>%
    dplyr::arrange(site_id, datetime)
  
  list(
    target_long_daily = target_long_daily,
    target_wide_daily = target_wide_daily
  )
}


#### test