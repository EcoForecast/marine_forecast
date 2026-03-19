library(tidyverse)
library(lubridate)
library(rjags)
library(coda)


# state-space model
dlm_ssm_string <- "
model {

  # observation model
  for (t in 1:n) {
    y[t] ~ dnorm(x[t], tau_obs)
  }

  # process model
  x[1] ~ dnorm(mu[1], tau_init)

  for (t in 2:n) {
    mu[t] <- beta0 +
             (1 + beta_x) * x[t-1] +
             beta_z1 * z1[t] +
             beta_z2 * z2[t] +
             beta_s1 * sin1[t] +
             beta_c1 * cos1[t] +
             beta_s2 * sin2[t] +
             beta_c2 * cos2[t]

    x[t] ~ dnorm(mu[t], tau_proc)
  }

  mu[1] <- beta0 +
           beta_s1 * sin1[1] +
           beta_c1 * cos1[1] +
           beta_s2 * sin2[1] +
           beta_c2 * cos2[1]

  # priors
  beta0   ~ dnorm(0, 1)
  beta_x  ~ dnorm(0, 1)
  beta_z1 ~ dnorm(0, 1)
  beta_z2 ~ dnorm(0, 1)
  beta_s1 ~ dnorm(0, 1)
  beta_c1 ~ dnorm(0, 1)
  beta_s2 ~ dnorm(0, 1)
  beta_c2 ~ dnorm(0, 1)

  sigma_obs  ~ dunif(0, 5)
  sigma_proc ~ dunif(0, 5)
  sigma_init ~ dunif(0, 5)

  tau_obs  <- pow(sigma_obs, -2)
  tau_proc <- pow(sigma_proc, -2)
  tau_init <- pow(sigma_init, -2)
}
"

# Prepare one site from target_wide_daily
prepare_dlm_data <- function(targets_processed,
                             site_select = 1,
                             response_var = "chlora_cci",
                             covar1 = "chlora_modis",
                             covar2 = "chlora_buoy",
                             min_obs = 100) {
  
  df <- targets_processed$target_wide_daily %>%
    dplyr::filter(site_id == site_select) %>%
    dplyr::mutate(datetime = as.Date(datetime)) %>%
    dplyr::arrange(datetime) %>%
    tidyr::complete(datetime = seq(min(datetime), max(datetime), by = "day")) %>%
    dplyr::mutate(site_id = site_select)
  
  if (!response_var %in% names(df)) stop("response_var not found.")
  if (!covar1 %in% names(df)) stop("covar1 not found.")
  if (!covar2 %in% names(df)) stop("covar2 not found.")
  
  y_raw <- df[[response_var]]
  z1_raw <- df[[covar1]]
  z2_raw <- df[[covar2]]
  
  if (sum(!is.na(y_raw)) < min_obs) {
    stop("Too few observed values in response series.")
  }
  
  # simple scaling
  y_center <- mean(y_raw, na.rm = TRUE)
  y_scale  <- sd(y_raw, na.rm = TRUE)
  
  z1_center <- mean(z1_raw, na.rm = TRUE)
  z1_scale  <- sd(z1_raw, na.rm = TRUE)
  
  z2_center <- mean(z2_raw, na.rm = TRUE)
  z2_scale  <- sd(z2_raw, na.rm = TRUE)
  
  if (!is.finite(y_scale) || y_scale == 0) y_scale <- 1
  if (!is.finite(z1_scale) || z1_scale == 0) z1_scale <- 1
  if (!is.finite(z2_scale) || z2_scale == 0) z2_scale <- 1
  
  df <- df %>%
    dplyr::mutate(
      doy = yday(datetime),
      y  = ( .data[[response_var]] - y_center ) / y_scale,
      z1 = ( .data[[covar1]]     - z1_center ) / z1_scale,
      z2 = ( .data[[covar2]]     - z2_center ) / z2_scale,
      # fill missing covariates with 0 after scaling = mean imputation
      z1 = ifelse(is.na(z1), 0, z1),
      z2 = ifelse(is.na(z2), 0, z2),
      sin1 = sin(2*pi*doy/365.25),
      cos1 = cos(2*pi*doy/365.25),
      sin2 = sin(4*pi*doy/365.25),
      cos2 = cos(4*pi*doy/365.25)
    )
  
  list(
    data = df,
    response_var = response_var,
    y_center = y_center,
    y_scale = y_scale
  )
}


# Fit model
fit_dlm_ssm <- function(prepped,
                        n_chains = 3,
                        burn_in = 5000,
                        n_iter = 10000,
                        thin = 5) {
  
  dat <- prepped$data
  y_obs <- dat$y[!is.na(dat$y)]
  n <- nrow(dat)
  
  inits <- replicate(n_chains, list(
    x = ifelse(
      is.na(dat$y),
      rnorm(n, mean(y_obs), sd(y_obs)),
      dat$y + rnorm(n, 0, 0.05)
    ),
    beta0 = rnorm(1, 0, 0.2),
    beta_x = rnorm(1, 0, 0.2),
    beta_z1 = rnorm(1, 0, 0.2),
    beta_z2 = rnorm(1, 0, 0.2),
    beta_s1 = rnorm(1, 0, 0.2),
    beta_c1 = rnorm(1, 0, 0.2),
    beta_s2 = rnorm(1, 0, 0.2),
    beta_c2 = rnorm(1, 0, 0.2),
    sigma_obs = runif(1, 0.05, 1),
    sigma_proc = runif(1, 0.05, 1),
    sigma_init = runif(1, 0.05, 1)
  ), simplify = FALSE)
  
  jags_data <- list(
    n = n,
    y = dat$y,
    z1 = dat$z1,
    z2 = dat$z2,
    sin1 = dat$sin1,
    cos1 = dat$cos1,
    sin2 = dat$sin2,
    cos2 = dat$cos2
  )
  
  jm <- jags.model(
    file = textConnection(dlm_ssm_string),
    data = jags_data,
    inits = inits,
    n.chains = n_chains,
    quiet = TRUE
  )
  
  update(jm, n.iter = burn_in, progress.bar = "none")
  
  samp <- coda.samples(
    model = jm,
    variable.names = c(
      "x", "beta0", "beta_x", "beta_z1", "beta_z2",
      "beta_s1", "beta_c1", "beta_s2", "beta_c2",
      "sigma_obs", "sigma_proc", "sigma_init"
    ),
    n.iter = n_iter,
    thin = thin,
    progress.bar = "none"
  )
  
  list(model = jm, samples = samp)
}

# Summarize posterior
summarize_dlm_fit <- function(prepped, samples) {
  draws <- as.matrix(samples)
  x_cols <- grep("^x\\[", colnames(draws))
  x_draws <- draws[, x_cols, drop = FALSE] * prepped$y_scale + prepped$y_center
  
  ci <- apply(x_draws, 2, quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  
  prepped$data %>%
    dplyr::transmute(
      datetime,
      site_id,
      observation = .data[[prepped$response_var]]
    ) %>%
    dplyr::mutate(
      lower_95 = ci[1, ],
      median   = ci[2, ],
      upper_95 = ci[3, ]
    )
}

# Plot Function
plot_dlm_fit <- function(posterior_df, response_var = "chlora_cci") {
  ggplot(posterior_df, aes(x = datetime)) +
    geom_ribbon(aes(ymin = lower_95, ymax = upper_95),
                fill = "skyblue", alpha = 0.35) +
    geom_line(aes(y = median), color = "blue4", linewidth = 0.8) +
    geom_point(aes(y = observation), color = "black", size = 0.9, alpha = 0.7, na.rm = TRUE) +
    labs(
      title = paste("Advanced state-space model:", unique(posterior_df$site_id), "/", response_var),
      subtitle = "95% credible interval and posterior median latent state",
      x = "Date",
      y = response_var
    ) +
    theme_bw()
}

# Runner
run_advanced_ssm <- function(targets_processed,
                             site_select = 1,
                             response_var = "chlora_cci",
                             covar1 = "chlora_modis",
                             covar2 = "chlora_buoy",
                             out_dir = "outputs/advanced_ssm") {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  prepped <- prepare_dlm_data(
    targets_processed = targets_processed,
    site_select = site_select,
    response_var = response_var,
    covar1 = covar1,
    covar2 = covar2
  )
  
  fit <- fit_dlm_ssm(prepped)
  posterior_df <- summarize_dlm_fit(prepped, fit$samples)
  
  p <- plot_dlm_fit(posterior_df, response_var = response_var)
  
  ggsave(
    file.path(out_dir, paste0("advanced_ssm_site", site_select, "_", response_var, ".png")),
    p, width = 12, height = 5.5, dpi = 160
  )
  
  saveRDS(
    list(prepped = prepped, posterior = posterior_df, samples = fit$samples),
    file.path(out_dir, paste0("advanced_ssm_site", site_select, "_", response_var, ".rds"))
  )
  
  print(summary(fit$samples))
  plot(fit$samples[, c("beta_x", "beta_z1", "beta_z2", "sigma_obs", "sigma_proc")])
  
  invisible(list(prepped = prepped, fit = fit, posterior = posterior_df, plot = p))
}

res <- run_advanced_ssm(
  targets_processed = targets_processed,
  site_select = 1,
  response_var = "chlora_cci",
  covar1 = "chlora_modis",
  covar2 = "chlora_buoy",
  out_dir = "outputs/advanced_ssm"
)