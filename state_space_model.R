library(tidyverse)
library(lubridate)
library(rjags)
library(coda)
library(dplyr)

source("01_download_data.R")
source("02_process_data.R")
source("03_make_visuals.R")

targets_raw <- download_targets()
targets_processed <- process_targets(targets_raw)

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
  
  if (!response_var %in% names(df))
    stop("response_var '", response_var, "' not found in wide data.")
  
  y_raw <- df[[response_var]]
  
  if (sum(!is.na(y_raw)) < min_obs)
    stop("Too few observed values (", sum(!is.na(y_raw)), ") in ", response_var,
         " — need at least ", min_obs)
  
  
  y_center <- mean(y_raw, na.rm = TRUE)
  y_scale  <- sd(y_raw,   na.rm = TRUE)
  if (!is.finite(y_scale) || y_scale == 0) y_scale <- 1
  
  df <- df %>%
    dplyr::mutate(
      doy  = lubridate::yday(datetime),
      y    = (!!rlang::sym(response_var) - y_center) / y_scale,
      sin1 = sin(2 * pi * doy / 365.25),
      cos1 = cos(2 * pi * doy / 365.25),
      sin2 = sin(4 * pi * doy / 365.25),
      cos2 = cos(4 * pi * doy / 365.25)
    )
  
  message("  Series length  : ", nrow(df), " days")
  message("  Observed y     : ", sum(!is.na(df$y)),
          " (", round(100 * mean(!is.na(df$y)), 1), "%)")
  message("  Missing y      : ", sum(is.na(df$y)),
          " — will be sampled by JAGS")
  
  list(
    data         = df,
    response_var = response_var,
    y_center     = y_center,
    y_scale      = y_scale
  )
}


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
    file       = textConnection(dlm_ssm_string),
    data       = jags_data,
    inits      = inits,
    n.chains   = n_chains,
    quiet      = TRUE
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
    n.iter         = n_iter,
    thin           = thin,
    progress.bar   = "none"
  )
  
  list(model = jm, samples = samp)
}


summarize_dlm_fit <- function(prepped, samples) {
  
  draws  <- as.matrix(samples)
  x_cols <- grep("^x\\[", colnames(draws))
  
  x_draws <- draws[, x_cols, drop = FALSE] * prepped$y_scale + prepped$y_center
  
  ci <- apply(x_draws, 2, quantile,
              probs = c(0.025, 0.25, 0.50, 0.75, 0.975), na.rm = TRUE)
  
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

generate_ssm_forecast <- function(prepped,
                                  samples,
                                  forecast_date,
                                  horizon_days = 35L,
                                  n_members    = 31L) {
  
  draws      <- as.matrix(samples)
  n_draws    <- nrow(draws)
  x_cols     <- grep("^x\\[", colnames(draws))
  n_t        <- length(x_cols)
  
  set.seed(123)
  idx <- sample(n_draws, n_members, replace = FALSE)
  
  x_last_scaled <- draws[idx, x_cols[n_t]]
  
  beta0   <- draws[idx, "beta0"]
  beta_x  <- draws[idx, "beta_x"]
  beta_s1 <- draws[idx, "beta_s1"]
  beta_c1 <- draws[idx, "beta_c1"]
  beta_s2 <- draws[idx, "beta_s2"]
  beta_c2 <- draws[idx, "beta_c2"]
  sig_proc <- draws[idx, "sigma_proc"]
  sig_obs  <- draws[idx, "sigma_obs"]
  
  forecast_dates <- seq(as.Date(forecast_date) + 1,
                        as.Date(forecast_date) + horizon_days,
                        by = "day")
  
  purrr::map_dfr(seq_len(n_members), function(m) {
    
    x_t <- x_last_scaled[m]
    
    purrr::map_dfr(seq_along(forecast_dates), function(i) {
      
      fdate <- forecast_dates[i]
      doy   <- lubridate::yday(fdate)
      
      sin1 <- sin(2 * pi * doy / 365.25)
      cos1 <- cos(2 * pi * doy / 365.25)
      sin2 <- sin(4 * pi * doy / 365.25)
      cos2 <- cos(4 * pi * doy / 365.25)
      
      mu_t <- beta0[m] +
        (1 + beta_x[m]) * x_t +
        beta_s1[m] * sin1 + beta_c1[m] * cos1 +
        beta_s2[m] * sin2 + beta_c2[m] * cos2
      
      x_t <<- rnorm(1, mu_t, sig_proc[m])  
      
      pred_scaled <- rnorm(1, x_t, sig_obs[m])
      pred <- pred_scaled * prepped$y_scale + prepped$y_center
      pred <- max(pred, 0)  
      
      tibble::tibble(
        datetime           = as.POSIXct(fdate),
        reference_datetime = as.POSIXct(forecast_date),
        site_id            = as.character(prepped$data$site_id[1]),
        family             = "ensemble",
        parameter          = as.character(m),
        variable           = prepped$response_var,
        prediction         = pred,
        model_id           = "bu4cast_coastal_ssm"
      )
    })
  })
}

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

run_advanced_ssm <- function(targets_processed,
                             site_select   = 1,
                             response_var  = "chlora_cci",
                             forecast_date = Sys.Date(),
                             horizon_days  = 35L,
                             n_members     = 31L,
                             out_dir       = "outputs/advanced_ssm") {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  tag <- paste0("site", site_select, "_", response_var)
  
  message("\n=== SSM: ", tag, " ===")
  
  prepped <- prepare_dlm_data(
    targets_processed = targets_processed,
    site_select       = site_select,
    response_var      = response_var
  )
  
  fit <- fit_dlm_ssm(prepped)
  
  posterior_df <- summarize_dlm_fit(prepped, fit$samples)
  
  message("Generating ", horizon_days, "-day forecast with ",
          n_members, " ensemble members...")
  forecast <- generate_ssm_forecast(
    prepped       = prepped,
    samples       = fit$samples,
    forecast_date = forecast_date,
    horizon_days  = horizon_days,
    n_members     = n_members
  )
  
  p_fit <- plot_dlm_fit(posterior_df, response_var)
  ggplot2::ggsave(
    file.path(out_dir, paste0("ssm_fit_", tag, ".png")),
    p_fit, width = 12, height = 5.5, dpi = 160
  )
  
  fc_summary <- forecast |>
    dplyr::mutate(datetime = as.Date(datetime)) |>
    dplyr::group_by(datetime) |>
    dplyr::summarise(
      p05 = quantile(prediction, 0.05),
      p25 = quantile(prediction, 0.25),
      p50 = median(prediction),
      p75 = quantile(prediction, 0.75),
      p95 = quantile(prediction, 0.95),
      .groups = "drop"
    )
  
  p_fc <- ggplot2::ggplot(fc_summary, ggplot2::aes(x = datetime)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = p05, ymax = p95),
                         fill = "skyblue", alpha = 0.25) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = p25, ymax = p75),
                         fill = "steelblue", alpha = 0.40) +
    ggplot2::geom_line(ggplot2::aes(y = p50),
                       color = "blue4", linewidth = 0.8) +
    ggplot2::labs(
      title    = paste("SSM forecast:", response_var, "—", forecast_date),
      subtitle = "50% and 95% posterior predictive intervals",
      x = "Date", y = response_var
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    file.path(out_dir, paste0("ssm_forecast_", tag, ".png")),
    p_fc, width = 10, height = 5, dpi = 160
  )
  
  saveRDS(
    list(prepped = prepped, posterior = posterior_df,
         samples = fit$samples, forecast = forecast),
    file.path(out_dir, paste0("ssm_", tag, ".rds"))
  )
  
  message("\nParameter summary:")
  print(summary(fit$samples[, c(
    "beta0", "beta_x", "beta_s1", "beta_c1",
    "beta_s2", "beta_c2", "sigma_obs", "sigma_proc"
  )]))
  
  graphics::plot(fit$samples[, c(
    "beta_x", "sigma_obs", "sigma_proc"
  )])
  
  invisible(list(
    prepped   = prepped,
    fit       = fit,
    posterior = posterior_df,
    forecast  = forecast,
    plots     = list(fit = p_fit, forecast = p_fc)
  ))
}



res <- run_advanced_ssm(
  targets_processed = targets_processed,
  site_select       = 1,
  response_var      = "chlora_cci",
  forecast_date     = Sys.Date(),
  horizon_days      = 35L,
  n_members         = 31L,
  out_dir           = "outputs/advanced_ssm"
)

draws <- as.matrix(res$fit$samples)

param_cols <- c("beta0", "beta_x", "beta_s1", "beta_c1",
                "beta_s2", "beta_c2", "sigma_obs", "sigma_proc")

param_table <- t(apply(draws[, param_cols], 2, function(x) {
  c(Mean  = round(mean(x), 4),
    SD    = round(sd(x),   4),
    `2.5%`  = round(quantile(x, 0.025), 4),
    `50%`   = round(quantile(x, 0.500), 4),
    `97.5%` = round(quantile(x, 0.975), 4))
}))

