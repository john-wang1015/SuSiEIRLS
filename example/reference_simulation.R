## Reference simulation generators for small SuSiEIRLS benchmarks.
##
## These examples use small synthetic data, not locus-scale real genotype
## simulation. The default signal allocation follows the current benchmark
## convention: var(etaZ) = 0.2 and var(etaX) = 0.05.

reference_benchmark_config <- function(ns = 500L,
                                       n_rep = 50L,
                                       families = c("binary_pg", "nb",
                                                    "cox", "tw", "beta"),
                                       p = 10L,
                                       q = 3L,
                                       true_idx = c(2L, 7L),
                                       coverage = 0.95,
                                       min_abs_corr = 0.5,
                                       var_eta_z = 0.2,
                                       var_eta_x = 0.05) {
  list(
    ns = as.integer(ns),
    n_rep = as.integer(n_rep),
    families = families,
    p = as.integer(p),
    q = as.integer(q),
    true_idx = as.integer(true_idx),
    L = length(true_idx),
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    var_eta_z = var_eta_z,
    var_eta_x = var_eta_x,
    oracle = list(
      nb_theta = 2,
      tw_p = 1.5,
      tw_phi = 1.2,
      beta_phi = 8,
      cox_base_hazard = 0.05,
      cox_target_censor = 0.5
    )
  )
}

reference_ar_cov <- function(p, rho = 0.3) {
  toeplitz(rho ^ (0:(p - 1)))
}

reference_scale_to_var <- function(x, target) {
  x <- as.numeric(scale(x, center = TRUE, scale = FALSE))
  vx <- stats::var(x)
  if (!is.finite(vx) || vx <= 0) {
    stop("Cannot scale a zero-variance vector.")
  }
  x * sqrt(target / vx)
}

reference_eta_components <- function(X, Z, cfg) {
  beta <- rep(0, ncol(X))
  beta[cfg$true_idx] <- c(1, -1)[seq_along(cfg$true_idx)]
  alpha <- seq(0.7, -0.4, length.out = ncol(Z))

  eta_z <- reference_scale_to_var(as.numeric(Z %*% alpha), cfg$var_eta_z)
  eta_x <- reference_scale_to_var(as.numeric(X %*% beta), cfg$var_eta_x)
  list(eta_z = eta_z, eta_x = eta_x, eta = eta_z + eta_x)
}

reference_rtweedie <- function(mu, phi, p) {
  if (exists("rTweedie", envir = asNamespace("mgcv"), mode = "function")) {
    return(mgcv::rTweedie(mu = mu, phi = phi, p = p))
  }
  if (requireNamespace("tweedie", quietly = TRUE)) {
    return(tweedie::rtweedie(length(mu), mu = mu, phi = phi, power = p))
  }
  stop("No Tweedie RNG is available. Install tweedie or use mgcv with rTweedie.")
}

reference_make_response <- function(dat, family, cfg) {
  eta <- dat$eta

  if (identical(family, "binary_pg")) {
    y <- stats::rbinom(nrow(dat), 1L, stats::plogis(-0.2 + eta))
    return(list(y = y, trials = NULL, status = NULL, time = NULL))
  }

  if (identical(family, "nb")) {
    y <- stats::rnbinom(
      nrow(dat), mu = exp(log(4) + eta), size = cfg$oracle$nb_theta
    )
    return(list(y = y, trials = NULL, status = NULL, time = NULL))
  }

  if (identical(family, "tw")) {
    y <- reference_rtweedie(
      exp(log(3) + eta), phi = cfg$oracle$tw_phi, p = cfg$oracle$tw_p
    )
    return(list(y = y, trials = NULL, status = NULL, time = NULL))
  }

  if (identical(family, "beta")) {
    mu <- stats::plogis(-0.4 + eta)
    phi <- cfg$oracle$beta_phi
    y <- stats::rbeta(nrow(dat), mu * phi, (1 - mu) * phi)
    y <- pmin(pmax(y, .Machine$double.eps), 1 - .Machine$double.eps)
    return(list(y = y, trials = NULL, status = NULL, time = NULL))
  }

  if (identical(family, "cox")) {
    base <- cfg$oracle$cox_base_hazard
    t_true <- -log(stats::runif(nrow(dat))) / (base * exp(eta))
    censor_rate <- exp(stats::uniroot(
      function(log_rate) {
        mean(1 - exp(-exp(log_rate) * t_true)) -
          cfg$oracle$cox_target_censor
      },
      interval = c(-20, 20)
    )$root)
    c_time <- stats::rexp(nrow(dat), rate = censor_rate)
    time <- pmin(t_true, c_time)
    status <- as.integer(t_true <= c_time)
    y <- survival::Surv(time, status)
    return(list(y = y, trials = NULL, status = status, time = time))
  }

  stop("Unsupported reference family: ", family)
}

simulate_reference_dataset <- function(n, family, seed, cfg) {
  set.seed(seed)

  X <- MASS::mvrnorm(
    n = n, mu = rep(0, cfg$p),
    Sigma = reference_ar_cov(cfg$p, 0.3)
  )
  X <- scale(X)
  colnames(X) <- paste0("X", seq_len(cfg$p))

  Z <- scale(matrix(stats::rnorm(n * cfg$q), n, cfg$q))
  colnames(Z) <- paste0("Z", seq_len(cfg$q))

  eta <- reference_eta_components(X, Z, cfg)
  dat <- data.frame(
    X, Z,
    eta_z = eta$eta_z,
    eta_x = eta$eta_x,
    eta = eta$eta,
    check.names = FALSE
  )
  rsp <- reference_make_response(dat, family, cfg)

  list(
    family = family,
    X = X,
    Z = Z,
    y = rsp$y,
    trials = rsp$trials,
    time = rsp$time,
    status = rsp$status,
    eta_z = eta$eta_z,
    eta_x = eta$eta_x,
    eta = eta$eta,
    true_idx = cfg$true_idx,
    var_eta_z = stats::var(eta$eta_z),
    var_eta_x = stats::var(eta$eta_x),
    var_eta_total = stats::var(eta$eta),
    censor_rate = if (identical(family, "cox")) mean(rsp$status == 0L) else NA_real_
  )
}
