
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(readr)
})

safe_sd <- function(x, fallback = 1) {
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || is.na(s) || s == 0) fallback else s
}

safe_var <- function(x, fallback = 1e-6) {
  v <- stats::var(x, na.rm = TRUE)
  if (!is.finite(v) || is.na(v) || v <= 0) fallback else v
}

clamp_positive <- function(x, eps = 1e-6) {
  pmax(x, eps)
}


# Prepare data
prepare_dlm_data <- function(targets_processed,
                             site_select  = 1,
                             response_var = "chlora_cci",
                             min_obs      = 100,
                             end_date     = NULL) {
  
  df <- targets_processed$target_wide_daily %>%
    dplyr::filter(site_id == site_select) %>%
    dplyr::mutate(datetime = as.Date(datetime)) %>%
    dplyr::arrange(datetime)
  
  hist_end <- max(df$datetime, na.rm = TRUE)
  
  if (is.null(end_date)) {
    end_date <- hist_end
  } else {
    end_date <- as.Date(end_date)
    if (end_date < hist_end) end_date <- hist_end
  }
  
  df <- df %>%
    tidyr::complete(datetime = seq(min(datetime), end_date, by = "day")) %>%
    dplyr::mutate(site_id = site_select)
  
  if (!response_var %in% names(df)) {
    stop("response_var '", response_var, "' not found in wide data.")
  }
  
  y_raw <- df[[response_var]]
  
  if (sum(!is.na(y_raw)) < min_obs) {
    stop(
      "Too few observed values (", sum(!is.na(y_raw)),
      ") in ", response_var, " — need at least ", min_obs
    )
  }
  
  y_center <- mean(y_raw, na.rm = TRUE)
  y_scale  <- safe_sd(y_raw, fallback = 1)
  
  df <- df %>%
    dplyr::mutate(
      doy  = lubridate::yday(datetime),
      y    = (.data[[response_var]] - y_center) / y_scale,
      sin1 = sin(2 * pi * doy / 365.25),
      cos1 = cos(2 * pi * doy / 365.25),
      sin2 = sin(4 * pi * doy / 365.25),
      cos2 = cos(4 * pi * doy / 365.25)
    )
  
  
  list(
    data         = df,
    response_var = response_var,
    y_center     = y_center,
    y_scale      = y_scale
  )
}

#Estimate fixed model parameters from historical data
estimate_enkf_static_params <- function(prepped,
                                        forecast_date = NULL,
                                        sigma_obs_floor = 0.10,
                                        sigma_proc_floor = 0.05) {
  
  dat <- prepped$data
  
  if (!is.null(forecast_date)) {
    forecast_date <- as.Date(forecast_date)
    dat <- dat %>% dplyr::filter(datetime <= forecast_date)
  }
  
  dat_lag <- dat %>%
    dplyr::mutate(y_lag = dplyr::lag(y)) %>%
    dplyr::filter(
      !is.na(y),
      !is.na(y_lag),
      !is.na(sin1), !is.na(cos1),
      !is.na(sin2), !is.na(cos2)
    )
  
  if (nrow(dat_lag) < 15) {
    beta0   <- 0
    phi     <- 0.8
    beta_s1 <- 0
    beta_c1 <- 0
    beta_s2 <- 0
    beta_c2 <- 0
    sigma_proc <- 0.30
    sigma_obs  <- 0.20
  } else {
    fit <- stats::lm(
      y ~ y_lag + sin1 + cos1 + sin2 + cos2,
      data = dat_lag
    )
    
    cf <- stats::coef(fit)
    cf[is.na(cf)] <- 0
    
    beta0   <- unname(cf["(Intercept)"])
    phi     <- unname(cf["y_lag"])
    beta_s1 <- unname(cf["sin1"])
    beta_c1 <- unname(cf["cos1"])
    beta_s2 <- unname(cf["sin2"])
    beta_c2 <- unname(cf["cos2"])
    
    if (!is.finite(phi)) phi <- 0.8
    
    proc_resid <- stats::residuals(fit)
    sigma_proc <- max(sigma_proc_floor, safe_sd(proc_resid, fallback = 0.30))
    
    # observation noise in scaled space:
    sigma_obs <- max(sigma_obs_floor, min(0.75 * sigma_proc, 0.50))
  }
  
  list(
    beta0      = beta0,
    beta_x     = phi - 1,
    beta_s1    = beta_s1,
    beta_c1    = beta_c1,
    beta_s2    = beta_s2,
    beta_c2    = beta_c2,
    sigma_proc = sigma_proc,
    sigma_obs  = sigma_obs
  )
}

# EnKF one-step forecast
enkf_forecast_step <- function(x_prev_analysis,
                               params,
                               sin1, cos1, sin2, cos2) {
  
  n_members <- length(x_prev_analysis)
  
  mu_t <- params$beta0 +
    (1 + params$beta_x) * x_prev_analysis +
    params$beta_s1 * sin1 +
    params$beta_c1 * cos1 +
    params$beta_s2 * sin2 +
    params$beta_c2 * cos2
  
  x_forecast <- rnorm(
    n_members,
    mean = mu_t,
    sd   = params$sigma_proc
  )
  
  list(
    mu_t       = mu_t,
    x_forecast = x_forecast
  )
}

# EnKF one-step analysis
enkf_analysis_step <- function(x_forecast,
                               y_obs = NA_real_,
                               sigma_obs) {
  
  # if no observation today: analysis = forecast
  if (is.na(y_obs) || !is.finite(y_obs)) {
    return(list(
      x_analysis = x_forecast,
      K          = 0,
      P_f        = safe_var(x_forecast),
      P_a        = safe_var(x_forecast)
    ))
  }
  
  P_f <- safe_var(x_forecast, fallback = 1e-6)
  R   <- sigma_obs^2
  
  K <- P_f / (P_f + R)
  
  # stochastic EnKF
  y_pert <- rnorm(length(x_forecast), mean = y_obs, sd = sigma_obs)
  
  x_analysis <- x_forecast + K * (y_pert - x_forecast)
  
  list(
    x_analysis = x_analysis,
    K          = K,
    P_f        = P_f,
    P_a        = safe_var(x_analysis)
  )
}

# Initialize ensemble at day 1
initialize_enkf_state <- function(prepped,
                                  params,
                                  n_members = 31L) {
  
  dat <- prepped$data
  y_obs <- dat$y[!is.na(dat$y)]
  
  if (length(y_obs) == 0) {
    stop("No observed data available to initialize EnKF.")
  }
  
  y0 <- dat$y[1]
  
  if (is.na(y0)) {
    y0_mean <- mean(y_obs, na.rm = TRUE)
    y0_sd   <- safe_sd(y_obs, fallback = 1)
  } else {
    y0_mean <- y0
    y0_sd   <- max(params$sigma_obs, 0.10)
  }
  
  rnorm(n_members, mean = y0_mean, sd = y0_sd)
}

# Run EnKF over historical period up to forecast_date
run_enkf_assimilation_history <- function(prepped,
                                          forecast_date,
                                          n_members = 31L,
                                          seed = 123) {
  
  set.seed(seed)
  
  dat <- prepped$data %>%
    dplyr::filter(datetime <= as.Date(forecast_date)) %>%
    dplyr::arrange(datetime)
  
  if (nrow(dat) < 2) {
    stop("Need at least 2 days up to forecast_date for EnKF history run.")
  }
  
  params <- estimate_enkf_static_params(
    prepped = list(
      data = dat,
      response_var = prepped$response_var,
      y_center = prepped$y_center,
      y_scale = prepped$y_scale
    ),
    forecast_date = forecast_date
  )
  
  n_t <- nrow(dat)
  
  forecast_ens <- matrix(NA_real_, nrow = n_members, ncol = n_t)
  analysis_ens <- matrix(NA_real_, nrow = n_members, ncol = n_t)
  
  K_vec  <- rep(NA_real_, n_t)
  Pf_vec <- rep(NA_real_, n_t)
  Pa_vec <- rep(NA_real_, n_t)
  
  # t = 1
  x0 <- initialize_enkf_state(prepped = list(data = dat), params = params, n_members = n_members)
  
  # treat day 1 as forecast = init, then optionally assimilate y[1]
  forecast_ens[, 1] <- x0
  
  a1 <- enkf_analysis_step(
    x_forecast = forecast_ens[, 1],
    y_obs      = dat$y[1],
    sigma_obs  = params$sigma_obs
  )
  
  analysis_ens[, 1] <- a1$x_analysis
  K_vec[1]  <- a1$K
  Pf_vec[1] <- a1$P_f
  Pa_vec[1] <- a1$P_a
  
  # t = 2..n_t
  for (t in 2:n_t) {
    fc <- enkf_forecast_step(
      x_prev_analysis = analysis_ens[, t - 1],
      params = params,
      sin1   = dat$sin1[t],
      cos1   = dat$cos1[t],
      sin2   = dat$sin2[t],
      cos2   = dat$cos2[t]
    )
    
    forecast_ens[, t] <- fc$x_forecast
    
    an <- enkf_analysis_step(
      x_forecast = forecast_ens[, t],
      y_obs      = dat$y[t],
      sigma_obs  = params$sigma_obs
    )
    
    analysis_ens[, t] <- an$x_analysis
    K_vec[t]  <- an$K
    Pf_vec[t] <- an$P_f
    Pa_vec[t] <- an$P_a
  }
  
  list(
    data         = dat,
    params       = params,
    forecast_ens = forecast_ens,
    analysis_ens = analysis_ens,
    last_analysis = analysis_ens[, n_t],
    diagnostics = tibble::tibble(
      datetime = dat$datetime,
      y        = dat$y,
      K        = K_vec,
      P_f      = Pf_vec,
      P_a      = Pa_vec
    )
  )
}

# Summarize historical EnKF fit
summarize_dlm_enkf_fit <- function(prepped, enkf_fit) {
  
  dat <- enkf_fit$data
  
  fc_unscaled <- enkf_fit$forecast_ens * prepped$y_scale + prepped$y_center
  an_unscaled <- enkf_fit$analysis_ens * prepped$y_scale + prepped$y_center
  
  fc_ci <- apply(
    fc_unscaled, 2, quantile,
    probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
    na.rm = TRUE
  )
  
  an_ci <- apply(
    an_unscaled, 2, quantile,
    probs = c(0.025, 0.25, 0.50, 0.75, 0.975),
    na.rm = TRUE
  )
  
  tibble::tibble(
    datetime    = dat$datetime,
    site_id     = dat$site_id,
    observation = dat[[prepped$response_var]],
    
    forecast_lower_95 = fc_ci[1, ],
    forecast_lower_50 = fc_ci[2, ],
    forecast_median   = fc_ci[3, ],
    forecast_upper_50 = fc_ci[4, ],
    forecast_upper_95 = fc_ci[5, ],
    
    analysis_lower_95 = an_ci[1, ],
    analysis_lower_50 = an_ci[2, ],
    analysis_median   = an_ci[3, ],
    analysis_upper_50 = an_ci[4, ],
    analysis_upper_95 = an_ci[5, ]
  )
}

# =========================================================
# 8) Extract restart state from EnKF historical analysis
# =========================================================
extract_restart_state_from_enkf <- function(prepped,
                                            enkf_fit,
                                            forecast_date,
                                            n_members = 31L) {
  
  if (n_members != nrow(enkf_fit$analysis_ens)) {
    stop("n_members does not match EnKF ensemble size.")
  }
  
  params <- enkf_fit$params
  x_last <- enkf_fit$last_analysis
  
  tibble::tibble(
    site_id            = as.character(enkf_fit$data$site_id[1]),
    variable           = prepped$response_var,
    parameter          = seq_len(n_members),
    reference_datetime = as.POSIXct(as.Date(forecast_date)),
    
    x_last_scaled = x_last,
    
    beta0   = params$beta0,
    beta_x  = params$beta_x,
    beta_s1 = params$beta_s1,
    beta_c1 = params$beta_c1,
    beta_s2 = params$beta_s2,
    beta_c2 = params$beta_c2,
    
    sigma_proc = params$sigma_proc,
    sigma_obs  = params$sigma_obs,
    
    y_center = prepped$y_center,
    y_scale  = prepped$y_scale
  )
}

# Save / load restart state
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
  if (!dir.exists(out_dir)) return(NULL)
  
  patt <- paste0("^restart_", site_id, "_", variable, "_.*\\.rds$")
  files <- list.files(out_dir, pattern = patt, full.names = TRUE)
  
  if (length(files) == 0) return(NULL)
  
  latest <- files[order(file.info(files)$mtime, decreasing = TRUE)][1]
  message("Loading restart state from: ", latest)
  readRDS(latest)
}

# Update restart with today's observation via EnKF
update_restart_with_today_obs_enkf <- function(restart_state,
                                               prepped,
                                               forecast_date,
                                               seed = 123) {
  set.seed(seed)
  
  forecast_date <- as.Date(forecast_date)
  dat <- prepped$data
  
  today_row <- dat %>% dplyr::filter(datetime == forecast_date)
  
  if (nrow(today_row) != 1) {
    stop("forecast_date not found in prepped$data: ", forecast_date)
  }
  
  # consistency check
  if (!all.equal(restart_state$y_center[1], prepped$y_center, tolerance = 1e-8)) {
    warning("restart_state y_center differs from current prepped y_center.")
  }
  if (!all.equal(restart_state$y_scale[1], prepped$y_scale, tolerance = 1e-8)) {
    warning("restart_state y_scale differs from current prepped y_scale.")
  }
  
  params <- list(
    beta0      = restart_state$beta0[1],
    beta_x     = restart_state$beta_x[1],
    beta_s1    = restart_state$beta_s1[1],
    beta_c1    = restart_state$beta_c1[1],
    beta_s2    = restart_state$beta_s2[1],
    beta_c2    = restart_state$beta_c2[1],
    sigma_proc = restart_state$sigma_proc[1],
    sigma_obs  = restart_state$sigma_obs[1]
  )
  
  prev_date <- as.Date(restart_state$reference_datetime[1])
  
  # if restart already reflects today, do nothing
  if (forecast_date <= prev_date) {
    message("Restart state already at or beyond forecast_date; skipping daily update.")
    return(list(
      restart_state = restart_state,
      analysis_summary = NULL
    ))
  }
  
  # march forward day by day from prev_date + 1 to forecast_date
  dates_to_update <- seq(prev_date + 1, forecast_date, by = "day")
  
  x_a <- restart_state$x_last_scaled
  
  analysis_log <- vector("list", length(dates_to_update))
  
  for (k in seq_along(dates_to_update)) {
    d <- dates_to_update[k]
    row_d <- dat %>% dplyr::filter(datetime == d)
    
    if (nrow(row_d) != 1) {
      stop("Missing row in prepped$data for date: ", d)
    }
    
    fc <- enkf_forecast_step(
      x_prev_analysis = x_a,
      params = params,
      sin1   = row_d$sin1[1],
      cos1   = row_d$cos1[1],
      sin2   = row_d$sin2[1],
      cos2   = row_d$cos2[1]
    )
    
    an <- enkf_analysis_step(
      x_forecast = fc$x_forecast,
      y_obs      = row_d$y[1],
      sigma_obs  = params$sigma_obs
    )
    
    x_f <- fc$x_forecast
    x_a <- an$x_analysis
    
    fc_unscaled <- x_f * prepped$y_scale + prepped$y_center
    an_unscaled <- x_a * prepped$y_scale + prepped$y_center
    
    analysis_log[[k]] <- tibble::tibble(
      datetime    = d,
      observation = row_d[[prepped$response_var]][1],
      
      forecast_lower_95 = quantile(fc_unscaled, 0.025, na.rm = TRUE),
      forecast_lower_50 = quantile(fc_unscaled, 0.25,  na.rm = TRUE),
      forecast_median   = quantile(fc_unscaled, 0.50,  na.rm = TRUE),
      forecast_upper_50 = quantile(fc_unscaled, 0.75,  na.rm = TRUE),
      forecast_upper_95 = quantile(fc_unscaled, 0.975, na.rm = TRUE),
      
      analysis_lower_95 = quantile(an_unscaled, 0.025, na.rm = TRUE),
      analysis_lower_50 = quantile(an_unscaled, 0.25,  na.rm = TRUE),
      analysis_median   = quantile(an_unscaled, 0.50,  na.rm = TRUE),
      analysis_upper_50 = quantile(an_unscaled, 0.75,  na.rm = TRUE),
      analysis_upper_95 = quantile(an_unscaled, 0.975, na.rm = TRUE),
      
      K   = an$K,
      P_f = an$P_f,
      P_a = an$P_a
    )
  }
  
  restart_new <- restart_state %>%
    dplyr::mutate(
      reference_datetime = as.POSIXct(forecast_date),
      x_last_scaled      = x_a
    )
  
  list(
    restart_state    = restart_new,
    analysis_summary = dplyr::bind_rows(analysis_log)
  )
}

# 11) Forecast from saved restart state
#     Keep function name for compatibility
forecast_from_restart_state <- function(restart_state,
                                        forecast_date,
                                        horizon_days = 14L,
                                        model_id = "bu4cast_coastal_ssm",
                                        seed = 123) {
  
  set.seed(seed)
  
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
      
      # propagate latent state
      x_t <<- rnorm(1, mean = mu_t, sd = sig_proc)
      
      # posterior predictive in observation space
      pred_scaled <- rnorm(1, mean = x_t, sd = sig_obs)
      pred <- pred_scaled * y_scale + y_center
      
      # optional non-negativity constraint for chlorophyll-like variables
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
        
        # internal latent state for next restart
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
      title    = paste("Iterative EnKF forecast:", response_var, "—", forecast_date),
      subtitle = "50% and 95% posterior predictive intervals",
      x = "Date", y = response_var
    ) +
    ggplot2::theme_bw()
}

# Submit forecast
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
    if (!requireNamespace("neon4cast", quietly = TRUE)) {
      stop("Package 'neon4cast' is required for submission but is not installed.")
    }
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

# Main daily iterative runner
run_iterative_enkf_daily <- function(targets_processed,
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
  
  forecast_date <- as.Date(forecast_date)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  restart_dir  <- file.path(out_dir, "restart_states")
  analysis_dir <- file.path(out_dir, "analysis_outputs")
  dir.create(restart_dir,  recursive = TRUE, showWarnings = FALSE)
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  
  tag <- paste0("site", site_select, "_", response_var, "_", forecast_date)
  
  site_id_chr <- as.character(site_select)
  
  raw_site_df <- targets_processed$target_wide_daily %>%
    dplyr::filter(site_id == site_select) %>%
    dplyr::mutate(datetime = as.Date(datetime))
  
  last_obs_date <- max(raw_site_df$datetime[!is.na(raw_site_df[[response_var]])], na.rm = TRUE)
  
  prepped <- prepare_dlm_data(
    targets_processed = targets_processed,
    site_select       = site_select,
    response_var      = response_var,
    min_obs           = min_obs,
    end_date          = forecast_date
  )
  
  if (forecast_date > last_obs_date) {
    message(
      "Latest observed ", response_var, " is on ", last_obs_date,
      "; propagating EnKF forward without observations through ", forecast_date, "."
    )
  }
  
  restart_state   <- NULL
  posterior_df    <- NULL
  analysis_today  <- NULL
  enkf_fit        <- NULL
  
  if (!force_refit) {
    restart_state <- load_restart_state(
      site_id  = site_id_chr,
      variable = response_var,
      out_dir  = restart_dir
    )
  }

  if (is.null(restart_state)) {
    message("No restart state found. Running historical EnKF assimilation...")
    
    enkf_fit <- run_enkf_assimilation_history(
      prepped       = prepped,
      forecast_date = forecast_date,
      n_members     = n_members
    )
    
    posterior_df <- summarize_dlm_enkf_fit(prepped, enkf_fit)
    
    restart_state <- extract_restart_state_from_enkf(
      prepped       = prepped,
      enkf_fit      = enkf_fit,
      forecast_date = forecast_date,
      n_members     = n_members
    )
    
    analysis_today <- posterior_df %>%
      dplyr::filter(datetime == forecast_date)
    
  } else {
    message("Restart state found. Updating analysis to forecast_date with EnKF...")
    
    up <- update_restart_with_today_obs_enkf(
      restart_state = restart_state,
      prepped       = prepped,
      forecast_date = forecast_date
    )
    
    restart_state  <- up$restart_state
    analysis_today <- up$analysis_summary
  }
  
  # save current restart after analysis
  save_restart_state(
    restart_state = restart_state,
    out_dir       = restart_dir
  )
  
  # save analysis summary
  analysis_file <- NULL
  if (!is.null(analysis_today) && nrow(analysis_today) > 0) {
    analysis_file <- file.path(
      analysis_dir,
      paste0("analysis_", tag, ".csv")
    )
    readr::write_csv(analysis_today, analysis_file)
    message("Analysis summary written to: ", analysis_file)
  }
  
  # Forecast from today's analysis ensemble
  forecast <- forecast_from_restart_state(
    restart_state = restart_state,
    forecast_date = forecast_date,
    horizon_days  = horizon_days,
    model_id      = model_id
  )
  
  forecast_out <- forecast %>%
    dplyr::select(
      model_id, datetime, reference_datetime,
      site_id, family, parameter, variable, prediction
    )
  
  fc_summary <- summarize_forecast(forecast_out)
  
  # local forecast summary file
  forecast_summary_file <- file.path(
    out_dir,
    paste0("forecast_summary_", tag, ".csv")
  )
  readr::write_csv(fc_summary, forecast_summary_file)
  
  # update restart for tomorrow using forecast trajectory
  restart_new <- update_restart_state_from_forecast(
    restart_state = restart_state,
    forecast_df   = forecast,
    use_last_forecast_day = use_last_forecast_day_for_restart
  )
  
  save_restart_state(
    restart_state = restart_new,
    out_dir       = restart_dir
  )
  
  # EFI-format submission
  forecast_file_submit <- write_and_submit_forecast(
    forecast_df   = forecast_out,
    forecast_date = forecast_date,
    model_id      = model_id,
    out_dir       = out_dir,
    task          = task,
    do_submit     = do_submit
  )
  
  
  invisible(list(
    enkf_fit              = enkf_fit,
    prepped               = prepped,
    restart_after_analysis = restart_state,
    restart_after_forecast = restart_new,
    analysis_today        = analysis_today,
    analysis_file         = analysis_file,
    forecast              = forecast_out,
    forecast_full         = forecast,
    forecast_sum          = fc_summary,
    forecast_summary_file = forecast_summary_file,
    forecast_submit_file  = forecast_file_submit
  ))
}