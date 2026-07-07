.conquer_default_bandwidth <- function(n, p, h) {
  h <- as.numeric(h)[1L]
  if (is.finite(h) && h > 0) return(h)
  max(((log(n) + p) / n)^0.4, 0.05)
}

.conquer_validate_args <- function(tau, kernel, h, residual_variance_bounds,
                                   fi_floor_quantile, fi_floor_min) {
  if (!requireNamespace("conquer", quietly = TRUE)) {
    stop("The conquer package is required for SuSiEConquer().")
  }
  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) ||
      tau <= 0 || tau >= 1) {
    stop("tau must be a finite numeric scalar in (0, 1).")
  }
  match.arg(kernel, c("Gaussian", "logistic", "uniform",
                     "parabolic", "triangular"))
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h)) {
    stop("h must be a finite numeric scalar.")
  }
  if (!is.numeric(residual_variance_bounds) ||
      length(residual_variance_bounds) != 2L ||
      !is.finite(residual_variance_bounds[1L]) ||
      residual_variance_bounds[1L] < 0 ||
      residual_variance_bounds[2L] <= residual_variance_bounds[1L]) {
    stop("residual_variance_bounds must be c(lower, upper) with lower >= 0.")
  }
  if (!is.numeric(fi_floor_quantile) || length(fi_floor_quantile) != 1L ||
      !is.finite(fi_floor_quantile) || fi_floor_quantile < 0 ||
      fi_floor_quantile >= 0.5) {
    stop("fi_floor_quantile must be a finite numeric scalar in [0, 0.5).")
  }
  if (!is.numeric(fi_floor_min) || length(fi_floor_min) != 1L ||
      !is.finite(fi_floor_min) || fi_floor_min < 0) {
    stop("fi_floor_min must be a nonnegative finite numeric scalar.")
  }
  invisible(TRUE)
}

.conquer_fit_explicit <- function(y, D, tau, kernel, h, checkSing,
                                  tol, iteMax, stepBounded, stepMax) {
  y <- as.numeric(y)
  n <- length(y)
  if (is.null(D) || ncol(as.matrix(D)) == 0L) {
    h_used <- .conquer_default_bandwidth(n, 0L, h)
    a <- as.numeric(stats::quantile(y, probs = tau, names = FALSE))
    return(list(
      coeff = a,
      ite = 0L,
      residual = y - a,
      bandwidth = h_used,
      tau = tau,
      kernel = kernel,
      n = n,
      p = 0L,
      intercept_only = TRUE
    ))
  }

  D <- as.matrix(D)
  if (nrow(D) != n) stop("nrow(D) must equal length(y).")
  if (ncol(D) >= n) {
    stop("The conquer nuisance design must have fewer columns than rows.")
  }
  if (min(apply(D, 2, stats::sd)) == 0) {
    stop("The conquer nuisance design contains a constant non-intercept column.")
  }

  conquer::conquer(
    X = D, Y = y, tau = tau, kernel = kernel, h = h,
    checkSing = checkSing, tol = tol, iteMax = iteMax,
    stepBounded = stepBounded, stepMax = stepMax, ci = "none"
  )
}

.conquer_predictor_data <- function(Z = NULL, Xextra = NULL, n = NULL) {
  if (is.null(n)) {
    n <- if (!is.null(Z)) nrow(Z) else nrow(Xextra)
  }
  out <- NULL
  if (!is.null(Z) && ncol(as.matrix(Z)) > 0L) out <- as.matrix(Z)
  if (!is.null(Xextra) && ncol(as.matrix(Xextra)) > 0L) {
    Xextra <- as.matrix(Xextra)
    out <- if (is.null(out)) Xextra else cbind(out, Xextra)
  }
  if (is.null(out)) out <- matrix(nrow = n, ncol = 0L)
  out
}

.conquer_kernel_eval <- function(u, kernel, h) {
  kernel <- match.arg(kernel, c("Gaussian", "logistic", "uniform",
                               "parabolic", "triangular"))
  if (kernel == "Gaussian") {
    Phi <- stats::pnorm(u)
    fi <- stats::dnorm(u) / h
  } else if (kernel == "logistic") {
    Phi <- stats::plogis(u)
    fi <- stats::dlogis(u) / h
  } else if (kernel == "uniform") {
    Phi <- ifelse(u <= -1, 0, ifelse(u >= 1, 1, 0.5 + 0.5 * u))
    fi <- ifelse(abs(u) < 1, 0.5 / h, 0)
  } else if (kernel == "parabolic") {
    Phi <- ifelse(u <= -1, 0,
                  ifelse(u >= 1, 1, 0.5 + 0.75 * u - 0.25 * u^3))
    fi <- ifelse(abs(u) < 1, 0.75 * (1 - u^2) / h, 0)
  } else {
    Phi <- numeric(length(u))
    left <- u > -1 & u < 0
    right <- u >= 0 & u < 1
    Phi[u <= -1] <- 0
    Phi[u >= 1] <- 1
    Phi[left] <- 0.5 + u[left] + 0.5 * u[left]^2
    Phi[right] <- 0.5 + u[right] - 0.5 * u[right]^2
    fi <- ifelse(abs(u) < 1, (1 - abs(u)) / h, 0)
  }
  list(Phi = as.numeric(Phi), fi = as.numeric(fi))
}

.conquer_standardize_fi <- function(fi, kernel, fi_floor_quantile = 0,
                                    fi_floor_min = 0,
                                    fi_normalize = TRUE,
                                    fi_normalize_method = c("mean", "median"),
                                    keep_zero_weight = TRUE) {
  fi_normalize_method <- match.arg(fi_normalize_method)
  fi_raw <- as.numeric(fi)
  fi_raw[!is.finite(fi_raw) | fi_raw < 0] <- 0
  compact <- kernel %in% c("uniform", "parabolic", "triangular")
  keep_zero <- isTRUE(keep_zero_weight) && compact
  positive <- fi_raw > 0
  if (!any(positive)) stop("All conquer kernel densities are zero.")

  fi_mod <- fi_raw
  scale_val <- 1
  if (isTRUE(fi_normalize)) {
    if (fi_normalize_method == "mean") {
      scale_val <- mean(fi_mod[positive])
    } else {
      scale_val <- stats::median(fi_mod[positive])
    }
    if (!is.finite(scale_val) || scale_val <= 0) {
      stop("Could not compute a valid fi normalization scale.")
    }
    fi_mod[positive] <- fi_mod[positive] / scale_val
  }

  positive_mod <- fi_mod > 0
  floor_val <- fi_floor_min
  if (fi_floor_quantile > 0) {
    q_floor <- as.numeric(stats::quantile(
      fi_mod[positive_mod], probs = fi_floor_quantile,
      names = FALSE, type = 7
    ))
    if (is.finite(q_floor) && q_floor > floor_val) {
      floor_val <- q_floor
    }
  }
  if (floor_val > 0) {
    fi_mod[positive_mod] <- pmax(fi_mod[positive_mod], floor_val)
    if (!keep_zero) fi_mod[!positive_mod] <- floor_val
  }

  list(
    fi = fi_mod,
    fi_raw = fi_raw,
    fi_floor = floor_val,
    fi_scale = scale_val,
    n_zero_weight = sum(fi_mod <= 0)
  )
}

.conquer_extract_working <- function(fit, y,
                                     fi_floor_quantile = 0,
                                     fi_floor_min = 0,
                                     fi_normalize = TRUE,
                                     fi_normalize_method = c("mean", "median"),
                                     keep_zero_weight = TRUE) {
  y <- as.numeric(y)
  fi_normalize_method <- match.arg(fi_normalize_method)
  h <- as.numeric(fit$bandwidth)[1L]
  if (!is.finite(h) || h <= 0) stop("Invalid conquer bandwidth.")
  tau <- as.numeric(fit$tau)[1L]
  kernel <- fit$kernel
  eta <- y - as.numeric(fit$residual)
  u <- (eta - y) / h
  k <- .conquer_kernel_eval(u, kernel = kernel, h = h)
  der <- k$Phi - tau
  fw <- .conquer_standardize_fi(
    fi = k$fi,
    kernel = kernel,
    fi_floor_quantile = fi_floor_quantile,
    fi_floor_min = fi_floor_min,
    fi_normalize = fi_normalize,
    fi_normalize_method = fi_normalize_method,
    keep_zero_weight = keep_zero_weight
  )
  fi <- fw$fi

  good <- is.finite(eta) & is.finite(der) & is.finite(fi) & fi > 0
  if (!any(good)) stop("All conquer IRLS weights are zero.")
  z <- numeric(length(y))
  z[good] <- eta[good] - der[good] / fi[good]
  W <- numeric(length(y))
  W[good] <- fi[good]

  v <- tau * (1 - tau) / fi[good]
  v0 <- as.numeric(stats::median(v[is.finite(v) & v > 0]))
  if (!is.finite(v0) || v0 <= 0) stop("Could not compute a valid conquer v0.")
  wd <- sum(W^2)
  if (!is.finite(wd) || wd <= 0) stop("All conquer IRLS weights are zero.")

  list(
    pseudo_response = z,
    W_diag = W,
    Phi = k$Phi,
    fi = fi,
    der = der,
    eta = eta,
    h = h,
    v0 = v0,
    n_eff = (sum(W)^2) / wd,
    n_active = sum(W > 0),
    fi_raw = fw$fi_raw,
    fi_floor = fw$fi_floor,
    fi_scale = fw$fi_scale,
    n_zero_weight = fw$n_zero_weight
  )
}

.conquer_greedy_warm_start <- function(X, y, Z, tau, kernel, h, checkSing,
                                       tol, iteMax, stepBounded, stepMax,
                                       L.init = 1) {
  p <- ncol(X)
  k_init <- init_k_from_L(L.init, p)
  selected <- integer(0)
  available <- rep(TRUE, p)
  Xinit <- NULL

  D <- .conquer_predictor_data(Z, Xinit, n = length(y))
  fit <- .conquer_fit_explicit(
    y = y, D = D, tau = tau, kernel = kernel, h = h,
    checkSing = checkSing, tol = tol, iteMax = iteMax,
    stepBounded = stepBounded, stepMax = stepMax
  )

  for (step in seq_len(k_init)) {
    j <- select_by_residual_score(
      X = X, residual = as.numeric(fit$residual), available = available
    )
    if (is.na(j)) break
    selected <- c(selected, j)
    available[j] <- FALSE
    Xinit <- X[, selected, drop = FALSE]
    colnames(Xinit) <- paste0("InitX", seq_along(selected))
    D <- .conquer_predictor_data(Z, Xinit, n = length(y))
    fit <- .conquer_fit_explicit(
      y = y, D = D, tau = tau, kernel = kernel, h = h,
      checkSing = checkSing, tol = tol, iteMax = iteMax,
      stepBounded = stepBounded, stepMax = stepMax
    )
  }

  list(fit = fit, X_fit = Xinit)
}

#' SuSiE from conquer smoothed quantile IRLS sufficient statistics
#'
#' @param X Numeric n by p matrix of main variables.
#' @param Y Numeric response vector.
#' @param Z Optional nuisance covariate matrix.
#' @param tau Quantile level passed to \code{conquer::conquer()}.
#' @param kernel Kernel passed to \code{conquer::conquer()}.
#' @param h Bandwidth passed to \code{conquer::conquer()}; non-positive values
#'   use conquer's default rule.
#' @param checkSing,tol,iteMax,stepBounded,stepMax Main fitting arguments
#'   passed to \code{conquer::conquer()}.
#' @param L,max.iter,min.iter,max.eps,susie.iter Main SuSiE iteration controls.
#' @param residual_variance Initial SuSiE residual variance. If \code{NULL},
#'   uses \code{residual_variance_init_scale * v0}.
#' @param residual_variance_init_scale Multiplier for the initial residual
#'   variance when \code{residual_variance} is \code{NULL}.
#' @param residual_variance_bounds Multipliers for lower and upper SuSiE
#'   residual variance bounds relative to \code{v0}.
#' @param fi_floor_quantile Optional lower quantile used to floor positive
#'   kernel densities after optional normalization. The default \code{0} does
#'   not clip.
#' @param fi_floor_min Optional absolute lower bound for positive kernel
#'   densities after optional normalization. The default \code{0} does not clip.
#' @param fi_normalize If \code{TRUE}, rescale positive kernel densities.
#' @param fi_normalize_method Standardize positive kernel densities by their
#'   mean or median.
#' @param keep_zero_weight For compact kernels, keep exact zero densities as
#'   zero working weights.
#' @param projection_solver Solver for the weighted Schur projection.
#' @param susie_n Sample size passed to \code{susieR::susie_ss()}. Use
#'   \code{"n_eff"} for the weighted effective local size, \code{"n"} for the
#'   original sample size, or \code{"n_active"} for the number of nonzero-weight
#'   observations.
#' @param ... Additional arguments passed to \code{susieR::susie_ss()}.
#' @export
#' @importFrom conquer conquer
SuSiEConquer <- function(X, Y, Z = NULL,
                         tau = 0.5,
                         kernel = c("parabolic", "Gaussian", "logistic",
                                    "uniform", "triangular"),
                         h = 0,
                         checkSing = FALSE,
                         tol = 1e-4,
                         iteMax = 5000,
                         stepBounded = TRUE,
                         stepMax = 100,
                         L = 10,
                         max.iter = 5,
                         min.iter = 2,
                         max.eps = 1e-4,
                         susie.iter = 100,
                         verbose = TRUE,
                         n_threads = 1,
                         coverage = 0.9,
                         estimate_residual_variance = TRUE,
                         residual_variance = NULL,
                         residual_variance_init_scale = 0.5,
                         residual_variance_bounds = c(0.5, 2),
                         fi_floor_quantile = 0,
                         fi_floor_min = 0,
                         fi_normalize = TRUE,
                         fi_normalize_method = c("mean", "median"),
                         keep_zero_weight = TRUE,
                         projection_solver = c("ginv", "ridge"),
                         susie_n = c("n_eff", "n", "n_active"),
                         scaled_prior_variance = 1,
                         L.init = 1,
                         refit_noncs = TRUE,
                         noncs_var = 0.2,
                         suff_block_size = 10000L,
                         ...) {
  kernel <- match.arg(kernel)
  fi_normalize_method <- match.arg(fi_normalize_method)
  projection_solver <- match.arg(projection_solver)
  susie_n <- match.arg(susie_n)
  .conquer_validate_args(
    tau = tau, kernel = kernel, h = h,
    residual_variance_bounds = residual_variance_bounds,
    fi_floor_quantile = fi_floor_quantile,
    fi_floor_min = fi_floor_min
  )

  X <- as.matrix(X)
  Y <- as.numeric(Y)
  n <- length(Y)
  p <- ncol(X)
  if (nrow(X) != n) stop("nrow(X) must equal length(Y).")
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(p))
  suff_block_size <- validate_suff_block_size(suff_block_size)

  if (is.null(Z)) {
    Z <- matrix(nrow = n, ncol = 0L)
    ZI <- matrix(1, nrow = n, ncol = 1L)
    colnames(ZI) <- "Intercept"
  } else {
    if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1L)
    Z <- as.matrix(Z)
    if (nrow(Z) != n) stop("nrow(Z) must equal nrow(X).")
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    ZI <- cbind(Intercept = 1, Z)
  }

  warm <- .conquer_greedy_warm_start(
    X = X, y = Y, Z = Z, tau = tau, kernel = kernel, h = h,
    checkSing = checkSing, tol = tol, iteMax = iteMax,
    stepBounded = stepBounded, stepMax = stepMax, L.init = L.init
  )
  fit_final <- warm$fit
  X_fit <- warm$X_fit

  alpha <- clean_coef(fit_final$coeff[seq_len(ncol(ZI))])
  g <- numeric(0)
  beta <- rep(0, p)
  fitX <- NULL
  XCS <- NULL
  early_no_cs <- FALSE
  last_work <- NULL

  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    alpha_prev <- alpha

    work <- .conquer_extract_working(
      fit_final, Y,
      fi_floor_quantile = fi_floor_quantile,
      fi_floor_min = fi_floor_min,
      fi_normalize = fi_normalize,
      fi_normalize_method = fi_normalize_method,
      keep_zero_weight = keep_zero_weight
    )
    last_work <- work
    v0 <- work$v0
    rv_init <- if (is.null(residual_variance)) {
      residual_variance_init_scale * v0
    } else {
      as.numeric(residual_variance)[1L]
    }
    if (!is.finite(rv_init) || rv_init <= 0) {
      stop("Initial residual variance must be positive.")
    }
    rv_lower <- residual_variance_bounds[1L] * v0
    rv_upper <- residual_variance_bounds[2L] * v0
    susie_n_value <- switch(
      susie_n,
      n_eff = work$n_eff,
      n = n,
      n_active = work$n_active
    )
    susie_n_value <- max(2, as.numeric(susie_n_value))

    S <- which(is.finite(work$W_diag) & work$W_diag > 0 &
                 is.finite(work$pseudo_response))
    if (!length(S)) {
      stop("No nonzero-weight observations in the conquer local set.")
    }

    suff <- weighted_residual_suffstats(
      X = X[S, , drop = FALSE],
      y = work$pseudo_response[S],
      ZI = ZI[S, , drop = FALSE],
      weights = work$W_diag[S],
      n_threads = n_threads,
      block_size = suff_block_size,
      projection_solver = projection_solver
    )

    fitX <- susieR::susie_ss(
      XtX = suff$XtX, Xty = suff$Xty, yty = suff$yty,
      n = susie_n_value, L = L,
      scaled_prior_variance = scaled_prior_variance,
      estimate_residual_variance = estimate_residual_variance,
      residual_variance = rv_init,
      residual_variance_lowerbound = rv_lower,
      residual_variance_upperbound = rv_upper,
      max_iter = susie.iter,
      estimate_prior_method = "EM",
      coverage = coverage, ...
    )
    rm(suff)

    beta <- clean_coef(stats::coef(fitX)[-1L])
    CSdt <- summary(fitX)$vars
    cs_indices <- sort(unique(CSdt$cs[CSdt$cs > 0]))

    if (!length(cs_indices)) {
      if (iter <= min.iter) {
        noncs_res <- build_no_cs_noncs_refit_term(X, fitX)
        if (is.null(noncs_res)) {
          early_no_cs <- TRUE
          if (verbose) {
            cat("No credible set detected; returning current no-CS fit.\n")
          }
          break
        }
        XCS <- matrix(noncs_res, ncol = 1L)
        colnames(XCS) <- "Main_CS_noncs"
        XCS_refit <- XCS
      } else {
        early_no_cs <- TRUE
        if (verbose) {
          cat("No credible set detected; returning current no-CS fit.\n")
        }
        break
      }
    } else {
      Alpha_filtered <- fitX$alpha * 0
      for (i in cs_indices) {
        vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
        Alpha_filtered[i, vars_in_cs_i] <- fitX$alpha[i, vars_in_cs_i]
      }
      Alpha_filtered <- Alpha_filtered * sign(fitX$mu)
      XCS <- CppMatrix::matrixMultiply(X, as.matrix(Alpha_filtered),
                                       transB = TRUE)
      XCS <- XCS[, cs_indices, drop = FALSE]
      if (is.null(dim(XCS))) XCS <- matrix(XCS, ncol = 1L)
      colnames(XCS) <- paste0("Main_CS", cs_indices)
      XCS_refit <- XCS

      if (isTRUE(refit_noncs)) {
        noncs_term <- build_noncs_refit_term(
          X = X, fitX = fitX, CSdt = CSdt, cs_indices = cs_indices,
          XCS = XCS, noncs_var = noncs_var
        )
        if (!is.null(noncs_term)) {
          XCS_refit <- cbind(XCS_refit, Main_CS_noncs = noncs_term)
        }
      }
    }

    D <- .conquer_predictor_data(Z, XCS_refit, n = n)
    fit_final <- .conquer_fit_explicit(
      y = Y, D = D, tau = tau, kernel = kernel, h = h,
      checkSing = checkSing, tol = tol, iteMax = iteMax,
      stepBounded = stepBounded, stepMax = stepMax
    )
    X_fit <- XCS_refit

    alpha <- clean_coef(fit_final$coeff[seq_len(ncol(ZI))])
    err <- max(sqrt(mean((beta - beta_prev)^2)),
               sqrt(mean((alpha - alpha_prev)^2)))
    g[iter] <- err

    if (verbose) {
      cat(sprintf(
        "Iteration %d: err = %.3e, n_susie = %.1f, n_eff = %.1f, n_active = %d, h = %.4g, v0 = %.4g, sigma2 = %.4g\n",
        iter, err, susie_n_value, work$n_eff, work$n_active, work$h, v0,
        fitX$sigma2
      ))
    }

    if (err < max.eps && iter > min.iter) {
      if (verbose) cat("Converged!\n")
      break
    }
  }

  MainIndex <- if (is.null(fitX)) NULL else {
    Identifying_MainEffect(fitX, colnames(X))
  }
  if (!is.null(XCS)) {
    D <- .conquer_predictor_data(Z, XCS, n = n)
    fit_final <- .conquer_fit_explicit(
      y = Y, D = D, tau = tau, kernel = kernel, h = h,
      checkSing = checkSing, tol = tol, iteMax = iteMax,
      stepBounded = stepBounded, stepMax = stepMax
    )
  }

  last_err <- if (length(g)) utils::tail(g, 1L) else Inf
  out <- list(
    iter = if (exists("iter")) iter else 0L,
    error = g,
    converged = early_no_cs ||
      (exists("iter") && iter < max.iter && last_err < max.eps),
    fitX = fitX,
    fitJoint = fit_final,
    main_index = MainIndex,
    JointCoef = fit_final$coeff,
    n_eff = if (!is.null(last_work)) last_work$n_eff else NA_real_,
    n_active = if (!is.null(last_work)) last_work$n_active else NA_integer_,
    susie_n = susie_n,
    susie_n_value = if (exists("susie_n_value")) susie_n_value else NA_real_,
    tau = tau,
    kernel = kernel,
    bandwidth = if (!is.null(fit_final$bandwidth)) fit_final$bandwidth else NA_real_,
    v0 = if (!is.null(last_work)) last_work$v0 else NA_real_,
    fi_floor = if (!is.null(last_work)) last_work$fi_floor else NA_real_,
    fi_scale = if (!is.null(last_work)) last_work$fi_scale else NA_real_,
    n_zero_weight = if (!is.null(last_work)) last_work$n_zero_weight else NA_integer_,
    residual_variance = if (!is.null(fitX)) fitX$sigma2 else NA_real_,
    residual_variance_bounds = if (!is.null(last_work)) {
      last_work$v0 * residual_variance_bounds
    } else {
      c(NA_real_, NA_real_)
    },
    Phi = if (!is.null(last_work)) last_work$Phi else NULL,
    fi = if (!is.null(last_work)) last_work$fi else NULL,
    fi_raw = if (!is.null(last_work)) last_work$fi_raw else NULL,
    W_diag = if (!is.null(last_work)) last_work$W_diag else NULL,
    pseudo_response = if (!is.null(last_work)) last_work$pseudo_response else NULL,
    early_no_cs = early_no_cs
  )
  class(out) <- c("SuSiEConquer", "list")
  out
}
