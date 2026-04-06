# state-space model 
dlm_ssm_string <- "
model {

  for (t in 1:n) {
    y[t] ~ dnorm(x[t], tau_obs)
  }

  x[1] ~ dnorm(mu_init, tau_init)

  for (t in 2:n) {
    mu[t] <- beta0 +
             (1 + beta_x) * x[t-1] +
             beta_s1 * sin1[t] +
             beta_c1 * cos1[t] +
             beta_s2 * sin2[t] +
             beta_c2 * cos2[t]

    x[t] ~ dnorm(mu[t], tau_proc)
  }

  mu_init <- beta0 +
             beta_s1 * sin1[1] +
             beta_c1 * cos1[1] +
             beta_s2 * sin2[1] +
             beta_c2 * cos2[1]

  beta0   ~ dnorm(0, 1)
  beta_x  ~ dnorm(0, 1)
  beta_s1 ~ dnorm(0, 1)
  beta_c1 ~ dnorm(0, 1)
  beta_s2 ~ dnorm(0, 1)
  beta_c2 ~ dnorm(0, 1)

  sigma_obs  ~ dunif(0, 5)
  sigma_proc ~ dunif(0, 5)
  sigma_init ~ dunif(0, 5)

  tau_obs  <- pow(sigma_obs,  -2)
  tau_proc <- pow(sigma_proc, -2)
  tau_init <- pow(sigma_init, -2)
}
"

# Prepare data
prepare_dlm_data <- function(targets_processed,
                             site_select  = 1,
                             response_var = "chlora_cci",
                             min_obs      = 100) {
  
  df <- targets_processed$target_wide_daily %>%
    dplyr::filter(site_id == site_select) %>%
    dplyr::mutate(datetime = as.Date(datetime)) %>%
    dplyr::arrange(datetime) %>%
    tidyr::complete(datetime = seq(min(datetime), max(datetime), by = "day")) %>%
    dplyr::mutate(site_id = site_select)
  
  if (!response_var %in% names(df)) {
    stop("response_var '", response_var, "' not found in wide data.")
  }
  
  y_raw <- df[[response_var]]
  
  if (sum(!is.na(y_raw)) < min_obs) {
    stop("Too few observed values (", sum(!is.na(y_raw)),
         ") in ", response_var, " — need at least ", min_obs)
  }
  
  y_center <- mean(y_raw, na.rm = TRUE)
  y_scale  <- sd(y_raw, na.rm = TRUE)
  if (!is.finite(y_scale) || y_scale == 0) y_scale <- 1
  
  df <- df %>%
    dplyr::mutate(
      doy  = lubridate::yday(datetime),
      y    = (.data[[response_var]] - y_center) / y_scale,
      sin1 = sin(2 * pi * doy / 365.25),
      cos1 = cos(2 * pi * doy / 365.25),
      sin2 = sin(4 * pi * doy / 365.25),
      cos2 = cos(4 * pi * doy / 365.25)
    )
  
  message("Series length  : ", nrow(df), " days")
  message("Observed y     : ", sum(!is.na(df$y)),
          " (", round(100 * mean(!is.na(df$y)), 1), "%)")
  message("Missing y      : ", sum(is.na(df$y)),
          " — will be sampled by JAGS")
  
  list(
    data         = df,
    response_var = response_var,
    y_center     = y_center,
    y_scale      = y_scale
  )
}

# Fit SSM with JAGS
fit_dlm_ssm <- function(prepped,
                        n_chains = 3,
                        burn_in  = 15000,
                        n_iter   = 10000,
                        thin     = 5) {
  
  dat   <- prepped$data
  y_obs <- dat$y[!is.na(dat$y)]
  n     <- nrow(dat)
  
  inits <- replicate(n_chains, list(
    x          = ifelse(is.na(dat$y),
                        rnorm(n, mean(y_obs), sd(y_obs)),
                        dat$y + rnorm(n, 0, 0.05)),
    beta0      = rnorm(1, 0, 0.2),
    beta_x     = rnorm(1, 0, 0.2),
    beta_s1    = rnorm(1, 0, 0.2),
    beta_c1    = rnorm(1, 0, 0.2),
    beta_s2    = rnorm(1, 0, 0.2),
    beta_c2    = rnorm(1, 0, 0.2),
    sigma_obs  = runif(1, 0.05, 1),
    sigma_proc = runif(1, 0.05, 1),
    sigma_init = runif(1, 0.05, 1)
  ), simplify = FALSE)
  
  jags_data <- list(
    n    = n,
    y    = dat$y,
    sin1 = dat$sin1,
    cos1 = dat$cos1,
    sin2 = dat$sin2,
    cos2 = dat$cos2
  )
  
  message("Fitting JAGS model (", n_chains, " chains, ",
          burn_in, " burn-in + ", n_iter, " sampling)...")
  
  jm <- rjags::jags.model(
    file     = textConnection(dlm_ssm_string),
    data     = jags_data,
    inits    = inits,
    n.chains = n_chains,
    quiet    = TRUE
  )
  
  update(jm, n.iter = burn_in, progress.bar = "none")
  
  samp <- rjags::coda.samples(
    model          = jm,
    variable.names = c(
      "x",
      "beta0", "beta_x",
      "beta_s1", "beta_c1",
      "beta_s2", "beta_c2",
      "sigma_obs", "sigma_proc", "sigma_init"
    ),
    n.iter       = n_iter,
    thin         = thin,
    progress.bar = "none"
  )
  
  list(model = jm, samples = samp)
}

# 5. Posterior summary for fit
summarize_dlm_fit <- function(prepped, samples) {
  draws  <- as.matrix(samples)
  x_cols <- grep("^x\\[", colnames(draws))
  
  x_draws <- draws[, x_cols, drop = FALSE] * prepped$y_scale + prepped$y_center
  
  ci <- apply(
    x_draws, 2, quantile,
    probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
    na.rm = TRUE
  )
  
  prepped$data %>%
    dplyr::transmute(
      datetime,
      site_id,
      observation = .data[[prepped$response_var]]
    ) %>%
    dplyr::mutate(
      lower_95 = ci[1, ],
      lower_50 = ci[2, ],
      median   = ci[3, ],
      upper_50 = ci[4, ],
      upper_95 = ci[5, ]
    )
}

# Plot fit
plot_dlm_fit <- function(posterior_df, response_var = "chlora_cci") {
  ggplot2::ggplot(posterior_df, ggplot2::aes(x = datetime)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower_95, ymax = upper_95),
      fill = "skyblue", alpha = 0.25
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower_50, ymax = upper_50),
      fill = "steelblue", alpha = 0.40
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = median),
      color = "blue4", linewidth = 0.8
    ) +
    ggplot2::geom_point(
      ggplot2::aes(y = observation),
      color = "black", size = 0.9, alpha = 0.7, na.rm = TRUE
    ) +
    ggplot2::labs(
      title    = paste("State-space model fit:", response_var),
      subtitle = "50% (dark) and 95% (light) credible intervals + posterior median",
      x        = "Date",
      y        = response_var
    ) +
    ggplot2::theme_bw()
}

# Extract restart state from fit (first run only)
extract_restart_state_from_fit <- function(prepped,
                                           samples,
                                           forecast_date,
                                           n_members = 31L,
                                           seed = 123) {
  draws   <- as.matrix(samples)
  x_cols  <- grep("^x\\[", colnames(draws), value = TRUE)
  n_t     <- length(x_cols)
  n_draws <- nrow(draws)
  
  if (n_members > n_draws) {
    stop("n_members is larger than posterior draw count.")
  }
  
  set.seed(seed)
  idx <- sample(seq_len(n_draws), size = n_members, replace = FALSE)
  
  tibble::tibble(
    site_id            = as.character(prepped$data$site_id[1]),
    variable           = prepped$response_var,
    parameter          = seq_len(n_members),
    reference_datetime = as.POSIXct(forecast_date),
    
    x_last_scaled = draws[idx, x_cols[n_t]],
    
    beta0   = draws[idx, "beta0"],
    beta_x  = draws[idx, "beta_x"],
    beta_s1 = draws[idx, "beta_s1"],
    beta_c1 = draws[idx, "beta_c1"],
    beta_s2 = draws[idx, "beta_s2"],
    beta_c2 = draws[idx, "beta_c2"],
    
    sigma_proc = draws[idx, "sigma_proc"],
    sigma_obs  = draws[idx, "sigma_obs"],
    
    y_center = prepped$y_center,
    y_scale  = prepped$y_scale
  )
}

# Save / load restart state locally
save_restart_state <- function(restart_state,
                               out_dir = "outputs/advanced_ssm/restart_states") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  site_id  <- restart_state$site_id[1]
  variable <- restart_state$variable[1]
  ref_date <- as.Date(restart_state$reference_datetime[1])
  
  file_path <- file.path(
    out_dir,
    paste0("restart_", site_id, "_", variable, "_", ref_date, ".rds")
  )
  
  saveRDS(restart_state, file_path)
  message("Saved restart state to: ", file_path)
  invisible(file_path)
}

load_restart_state <- function(site_id,
                               variable,
                               out_dir = "outputs/advanced_ssm/restart_states") {
  patt <- paste0("^restart_", site_id, "_", variable, "_.*\\.rds$")
  files <- list.files(out_dir, pattern = patt, full.names = TRUE)
  
  if (length(files) == 0) return(NULL)
  
  latest <- files[order(file.info(files)$mtime, decreasing = TRUE)][1]
  message("Loading restart state from: ", latest)
  readRDS(latest)
}

# Forecast from saved restart state
forecast_from_restart_state <- function(restart_state,
                                        forecast_date,
                                        horizon_days = 14L,
                                        model_id = "bu4cast_coastal_ssm") {
  
  forecast_dates <- seq(
    as.Date(forecast_date) + 1,
    as.Date(forecast_date) + horizon_days,
    by = "day"
  )
  
  n_members <- nrow(restart_state)
  
  forecast <- purrr::map_dfr(seq_len(n_members), function(m) {
    
    x_t <- restart_state$x_last_scaled[m]
    
    beta0   <- restart_state$beta0[m]
    beta_x  <- restart_state$beta_x[m]
    beta_s1 <- restart_state$beta_s1[m]
    beta_c1 <- restart_state$beta_c1[m]
    beta_s2 <- restart_state$beta_s2[m]
    beta_c2 <- restart_state$beta_c2[m]
    
    sig_proc <- restart_state$sigma_proc[m]
    sig_obs  <- restart_state$sigma_obs[m]
    
    y_center <- restart_state$y_center[m]
    y_scale  <- restart_state$y_scale[m]
    
    purrr::map_dfr(seq_along(forecast_dates), function(i) {
      
      fdate <- forecast_dates[i]
      doy   <- lubridate::yday(fdate)
      
      sin1 <- sin(2 * pi * doy / 365.25)
      cos1 <- cos(2 * pi * doy / 365.25)
      sin2 <- sin(4 * pi * doy / 365.25)
      cos2 <- cos(4 * pi * doy / 365.25)
      
      mu_t <- beta0 +
        (1 + beta_x) * x_t +
        beta_s1 * sin1 + beta_c1 * cos1 +
        beta_s2 * sin2 + beta_c2 * cos2
      
      x_t <<- rnorm(1, mean = mu_t, sd = sig_proc)
      
      pred_scaled <- rnorm(1, mean = x_t, sd = sig_obs)
      pred <- pred_scaled * y_scale + y_center
      pred <- max(pred, 0)
      
      tibble::tibble(
        model_id           = model_id,
        datetime           = as.POSIXct(fdate),
        reference_datetime = as.POSIXct(forecast_date),
        site_id            = restart_state$site_id[m],
        family             = "ensemble",
        parameter          = restart_state$parameter[m],
        variable           = restart_state$variable[m],
        prediction         = pred,
        
        # internal state for next restart
        x_state_scaled_end = x_t
      )
    })
  })
  
  forecast
}

# Update restart state after forecast
update_restart_state_from_forecast <- function(restart_state,
                                               forecast_df,
                                               use_last_forecast_day = TRUE) {
  state_rows <- if (use_last_forecast_day) {
    forecast_df %>%
      dplyr::group_by(parameter) %>%
      dplyr::slice_tail(n = 1) %>%
      dplyr::ungroup()
  } else {
    forecast_df %>%
      dplyr::group_by(parameter) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::ungroup()
  }
  
  last_states <- state_rows %>%
    dplyr::select(parameter, datetime, x_state_scaled_end)
  
  restart_new <- restart_state %>%
    dplyr::select(-reference_datetime, -x_last_scaled) %>%
    dplyr::left_join(last_states, by = "parameter") %>%
    dplyr::rename(
      reference_datetime = datetime,
      x_last_scaled      = x_state_scaled_end
    )
  
  restart_new
}

# Summarize forecast ensemble
summarize_forecast <- function(forecast_df) {
  forecast_df %>%
    dplyr::mutate(datetime = as.Date(datetime)) %>%
    dplyr::group_by(datetime) %>%
    dplyr::summarise(
      p05 = quantile(prediction, 0.05, na.rm = TRUE),
      p25 = quantile(prediction, 0.25, na.rm = TRUE),
      p50 = median(prediction, na.rm = TRUE),
      p75 = quantile(prediction, 0.75, na.rm = TRUE),
      p95 = quantile(prediction, 0.95, na.rm = TRUE),
      .groups = "drop"
    )
}

plot_forecast_summary <- function(fc_summary, response_var, forecast_date) {
  ggplot2::ggplot(fc_summary, ggplot2::aes(x = datetime)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = p05, ymax = p95),
      fill = "skyblue", alpha = 0.25
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = p25, ymax = p75),
      fill = "steelblue", alpha = 0.40
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = p50),
      color = "blue4", linewidth = 0.8
    ) +
    ggplot2::labs(
      title    = paste("Iterative SSM forecast:", response_var, "—", forecast_date),
      subtitle = "50% and 95% posterior predictive intervals",
      x = "Date", y = response_var
    ) +
    ggplot2::theme_bw()
}

# Main daily iterative runner
run_iterative_ssm_daily <- function(targets_processed,
                                    site_select   = 1,
                                    response_var  = "chlora_cci",
                                    forecast_date = Sys.Date(),
                                    horizon_days  = 14L,
                                    n_members     = 31L,
                                    out_dir       = "outputs/advanced_ssm",
                                    model_id      = "bu4cast_coastal_ssm",
                                    force_refit   = FALSE,
                                    min_obs       = 100,
                                    use_last_forecast_day_for_restart = TRUE,
                                    do_submit     = FALSE,
                                    task          = "aquatics") {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  restart_dir <- file.path(out_dir, "restart_states")
  dir.create(restart_dir, recursive = TRUE, showWarnings = FALSE)
  
  tag <- paste0("site", site_select, "_", response_var, "_", as.Date(forecast_date))
  
  site_id_chr <- as.character(site_select)
  
  restart_state <- NULL
  fit <- NULL
  prepped <- NULL
  posterior_df <- NULL
  
  if (!force_refit) {
    restart_state <- load_restart_state(
      site_id  = site_id_chr,
      variable = response_var,
      out_dir  = restart_dir
    )
  }
  
  if (is.null(restart_state)) {
    message("No restart state found. Running initial fit...")
    
    prepped <- prepare_dlm_data(
      targets_processed = targets_processed,
      site_select       = site_select,
      response_var      = response_var,
      min_obs           = min_obs
    )
    
    fit <- fit_dlm_ssm(prepped)
    posterior_df <- summarize_dlm_fit(prepped, fit$samples)
    
    restart_state <- extract_restart_state_from_fit(
      prepped       = prepped,
      samples       = fit$samples,
      forecast_date = forecast_date,
      n_members     = n_members
    )
    
    # Save fit diagnostics on initial/refit run
    p_fit <- plot_dlm_fit(posterior_df, response_var)
    ggplot2::ggsave(
      file.path(out_dir, paste0("ssm_fit_", tag, ".png")),
      p_fit, width = 12, height = 5.5, dpi = 160
    )
    
    saveRDS(
      list(
        prepped      = prepped,
        posterior    = posterior_df,
        samples      = fit$samples,
        restart_init = restart_state
      ),
      file.path(out_dir, paste0("ssm_fit_object_", tag, ".rds"))
    )
    
    draws <- as.matrix(fit$samples)
    param_cols <- c("beta0", "beta_x", "beta_s1", "beta_c1",
                    "beta_s2", "beta_c2", "sigma_obs", "sigma_proc")
    
    param_table <- t(apply(draws[, param_cols, drop = FALSE], 2, function(x) {
      c(
        Mean    = round(mean(x), 4),
        SD      = round(sd(x), 4),
        `2.5%`  = round(quantile(x, 0.025), 4),
        `50%`   = round(quantile(x, 0.500), 4),
        `97.5%` = round(quantile(x, 0.975), 4)
      )
    })) %>%
      as.data.frame() %>%
      tibble::rownames_to_column("parameter_name")
    
    readr::write_csv(
      param_table,
      file.path(out_dir, paste0("ssm_param_table_", tag, ".csv"))
    )
    
  } else {
    message("Restart state found. Continue forecasting from saved state.")
  }
  
  forecast <- forecast_from_restart_state(
    restart_state = restart_state,
    forecast_date = forecast_date,
    horizon_days  = horizon_days,
    model_id      = model_id
  )
  
  restart_new <- update_restart_state_from_forecast(
    restart_state = restart_state,
    forecast_df   = forecast,
    use_last_forecast_day = use_last_forecast_day_for_restart
  )
  
  save_restart_state(
    restart_state = restart_new,
    out_dir       = restart_dir
  )
  
  forecast_out <- forecast %>%
    dplyr::select(
      model_id, datetime, reference_datetime,
      site_id, family, parameter, variable, prediction
    )
  
  # local forecast
  forecast_file_local <- file.path(
    out_dir,
    paste0("forecast_", tag, ".csv")
  )
  readr::write_csv(forecast_out, forecast_file_local)
  
  # summary and plot
  fc_summary <- summarize_forecast(forecast_out)
  readr::write_csv(
    fc_summary,
    file.path(out_dir, paste0("forecast_summary_", tag, ".csv"))
  )
  
  p_fc <- plot_forecast_summary(fc_summary, response_var, forecast_date)
  ggplot2::ggsave(
    file.path(out_dir, paste0("ssm_forecast_", tag, ".png")),
    p_fc, width = 10, height = 5, dpi = 160
  )
  
  # Format EFI challenge file and submit
  forecast_file_submit <- write_and_submit_forecast(
    forecast_df   = forecast_out,
    forecast_date = forecast_date,
    model_id      = model_id,
    out_dir       = out_dir,
    task          = task,
    do_submit     = do_submit
  )
  
  message("Local forecast saved to: ", forecast_file_local)
  message("Submission-format forecast saved to: ", forecast_file_submit)
  
  invisible(list(
    fit           = fit,
    prepped       = prepped,
    restart_old   = restart_state,
    restart_new   = restart_new,
    forecast      = forecast_out,
    forecast_full = forecast,
    forecast_sum  = fc_summary
  ))
}

# Submit Forecast
write_and_submit_forecast <- function(forecast_df,
                                      forecast_date,
                                      model_id,
                                      out_dir = "outputs/advanced_ssm",
                                      task = "aquatics",
                                      do_submit = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  forecast_out <- forecast_df %>%
    dplyr::mutate(
      project_id = "neon4cast",
      duration   = "P1D",
      # Format variable name
      variable = dplyr::case_when(
        variable == "chlora_cci" ~ "chla",
        TRUE ~ variable
      ),
      datetime = as.POSIXct(datetime, tz = "UTC"),
      reference_datetime = as.POSIXct(reference_datetime, tz = "UTC")
    ) %>%
    dplyr::select(
      project_id,
      model_id,
      datetime,
      reference_datetime,
      duration,
      site_id,
      family,
      parameter,
      variable,
      prediction
    )
  
  forecast_file <- file.path(
    out_dir,
    paste0(task, "-", as.Date(forecast_date), "-", model_id, ".csv.gz")
  )
  
  readr::write_csv(forecast_out, forecast_file)
  message("Forecast file written to: ", forecast_file)
  
  if (do_submit) {
    neon4cast::submit(
      forecast_file = forecast_file,
      metadata = NULL,
      ask = FALSE
    )
    message("Forecast submitted successfully.")
  } else {
    message("do_submit = FALSE, skipped submission.")
  }
  
  invisible(forecast_file)
}