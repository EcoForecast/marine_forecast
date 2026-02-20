library(dplyr)
targets_raw <- download_targets()
targets_processed <- process_targets(targets_raw)
make_example_visualizations(
  target_long_daily = targets_processed$target_long_daily
)