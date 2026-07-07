.qgam_validate_args <- function(qu, lsig, qgam_control, qgam_argGam,
                                qgam_discrete) {
  if (!requireNamespace("qgam", quietly = TRUE)) {
    stop("The qgam package is required for family = \"qgam\".")
  }
  if (!is.numeric(qu) || length(qu) != 1L || !is.finite(qu) ||
      qu <= 0 || qu >= 1) {
    stop("qu must be a finite numeric scalar in (0, 1).")
  }
  if (!is.null(lsig)) {
    stop("lsig must be NULL. Run_QGAM estimates calibrated lsig at every qgam refit.")
  }
  if (!is.list(qgam_control)) stop("qgam_control must be a list.")
  if (!is.null(qgam_argGam) && !is.list(qgam_argGam)) {
    stop("qgam_argGam must be NULL or a list.")
  }
  if (!is.logical(qgam_discrete) || length(qgam_discrete) != 1L ||
      is.na(qgam_discrete)) {
    stop("qgam_discrete must be TRUE or FALSE.")
  }
  invisible(TRUE)
}

.qgam_fit_explicit <- function(response, rhs, data, qu = 0.5, lsig = NULL,
                               qgam_control = list(), qgam_argGam = NULL,
                               qgam_discrete = FALSE) {
  .qgam_validate_args(
    qu = qu, lsig = lsig, qgam_control = qgam_control,
    qgam_argGam = qgam_argGam, qgam_discrete = qgam_discrete
  )
  qgam::qgam(
    form = .mgcv_explicit_formula(response, rhs),
    data = data,
    qu = qu,
    discrete = qgam_discrete,
    lsig = NULL,
    control = qgam_control,
    argGam = qgam_argGam
  )
}

.qgam_prepare_response <- function(y) {
  if (is.matrix(y) || inherits(y, "Surv")) {
    stop("qgam IRLS currently supports a numeric response vector.")
  }
  y <- as.numeric(y)
  list(data = data.frame(y = y), response = "y", n = length(y))
}

.qgam_fit_init <- function(X, response_info, Z, selected, qu, lsig,
                           qgam_control, qgam_argGam, qgam_discrete) {
  Xinit <- NULL
  if (length(selected) > 0L) {
    Xinit <- X[, selected, drop = FALSE]
    colnames(Xinit) <- paste0("InitX", seq_along(selected))
  }
  pred <- .mgcv_predictor_data(Z, Xinit, n = response_info$n)
  dat <- cbind(response_info$data, pred)
  .qgam_fit_explicit(
    response = response_info$response,
    rhs = colnames(pred),
    data = dat,
    qu = qu,
    lsig = lsig,
    qgam_control = qgam_control,
    qgam_argGam = qgam_argGam,
    qgam_discrete = qgam_discrete
  )
}

.qgam_greedy_warm_start <- function(X, response_info, Z, qu = 0.5,
                                    lsig = NULL, qgam_control = list(),
                                    qgam_argGam = NULL,
                                    qgam_discrete = FALSE, L.init = 1,
                                    init_cor_method = NULL) {
  p <- ncol(X)
  k_init <- init_k_from_L(L.init, p)
  selected <- integer(0)
  available <- rep(TRUE, p)
  fit <- .qgam_fit_init(
    X = X, response_info = response_info, Z = Z,
    selected = selected, qu = qu, lsig = lsig,
    qgam_control = qgam_control, qgam_argGam = qgam_argGam,
    qgam_discrete = qgam_discrete
  )

  for (step in seq_len(k_init)) {
    r <- stats::residuals(fit, type = "response")
    j <- select_by_residual_score(X = X, residual = r, available = available)
    if (is.na(j)) break
    selected <- c(selected, j)
    available[j] <- FALSE
    fit <- .qgam_fit_init(
      X = X, response_info = response_info, Z = Z,
      selected = selected, qu = qu, lsig = lsig,
      qgam_control = qgam_control, qgam_argGam = qgam_argGam,
      qgam_discrete = qgam_discrete
    )
  }

  fit
}

.qgam_lsig <- function(fit) {
  if (!is.null(fit$lsig) && length(fit$lsig) &&
      is.finite(as.numeric(fit$lsig[1L]))) {
    return(as.numeric(fit$lsig[1L]))
  }
  theta <- tryCatch(
    log(fit$family$getTheta(TRUE)),
    error = function(e) fit$lsig
  )
  if (is.null(theta) || !length(theta) || !is.finite(theta[1L])) {
    stop("Could not extract theta (lsig) from qgam object.")
  }
  as.numeric(theta[1L])
}

.qgam_extract_working <- function(fit, weight_cutoff = 0.005) {
  eta <- as.numeric(fit$linear.predictors)
  y_work <- as.numeric(fit$y)
  theta <- .qgam_lsig(fit)
  wt <- if (is.null(fit$prior.weights)) {
    rep(1, length(y_work))
  } else {
    as.numeric(fit$prior.weights)
  }

  Dval <- fit$family$Dd(y_work, eta, theta, wt)
  pseudo_response <- eta - Dval$Dmu / Dval$Dmu2
  W_diag <- as.numeric(Dval$Dmu2) / 2
  phi0 <- exp(theta)

  bad <- !is.finite(pseudo_response) | !is.finite(W_diag) | W_diag <= 0
  if (mean(bad) > 0.9) stop("Too many invalid qgam working observations.")
  if (any(bad)) {
    W_diag[bad] <- 0
    pseudo_response[bad] <- 0
  }
  W_diag <- robust_weight(W_diag, cutoff = weight_cutoff)
  weight_denom <- sum(W_diag^2)
  if (!is.finite(weight_denom) || weight_denom <= 0) {
    stop("All qgam working weights are zero.")
  }
  # qgam's Dmu2 already includes 1 / exp(lsig); SuSiE gets the dispersion-free
  # weight and uses exp(lsig) as the fixed residual variance.
  W_diag <- W_diag * phi0

  list(
    pseudo_response = pseudo_response,
    W_diag = W_diag,
    phi0 = phi0,
    n_eff = (sum(W_diag)^2) / weight_denom,
    lsig = theta
  )
}

#' Quantile-GAM IRLS-SuSiE path
#' @inheritParams SuSiE_IRLS
#' @param y Numeric response vector.
#' @param estimate_residual_variance Ignored in the qgam path. SuSiE residual
#'   variance is fixed to the current \code{exp(lsig)}.
#' @param residual_variance Ignored in the qgam path. The fixed SuSiE residual
#'   variance is extracted from the current calibrated qgam fit.
#' @param residual_variance_lowerbound Ignored in the qgam path.
#' @param residual_variance_upperbound Ignored in the qgam path.
#' @export
Run_QGAM <- function(X, y, Z = NULL, weight_cutoff = 0.005,
                     qu = 0.5, lsig = NULL,
                     qgam_control = list(),
                     qgam_argGam = NULL,
                     qgam_discrete = FALSE,
                     L, max.iter, min.iter, max.eps, susie.iter,
                     verbose = TRUE, n_threads = 1, coverage = 0.9,
                     estimate_residual_variance = FALSE,
                     residual_variance = NULL, scaled_prior_variance = 1,
                     residual_variance_lowerbound = 0.1,
                     residual_variance_upperbound = 1,
                     L.init = 1,
                     init_cor_method = NULL,
                     refit_noncs = TRUE,
                     noncs_var = 0.2,
                     suff_block_size = 10000L, ...) {

  n <- NROW(y)
  p <- ncol(X)
  suff_block_size <- validate_suff_block_size(suff_block_size)
  .qgam_validate_args(
    qu = qu, lsig = lsig, qgam_control = qgam_control,
    qgam_argGam = qgam_argGam, qgam_discrete = qgam_discrete
  )
  if (!isTRUE(verbose) && is.null(qgam_control$progress)) {
    qgam_control$progress <- FALSE
  }
  response_info <- .qgam_prepare_response(y)
  if (response_info$n != nrow(X)) stop("Length(y) must equal nrow(X).")

  if (is.null(Z)) {
    Z <- matrix(nrow = n, ncol = 0)
    ZI <- matrix(1, nrow = n, ncol = 1)
    colnames(ZI) <- "Intercept"
  } else {
    if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1)
    if (nrow(Z) != n) stop("nrow(Z) must equal nrow(X).")
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
    ZI <- cbind(Intercept = 1, Z)
  }

  fit_final <- .qgam_greedy_warm_start(
    X = X, response_info = response_info, Z = Z,
    qu = qu, lsig = lsig, qgam_control = qgam_control,
    qgam_argGam = qgam_argGam, qgam_discrete = qgam_discrete,
    L.init = L.init, init_cor_method = init_cor_method
  )

  alpha <- clean_coef(stats::coef(fit_final)[seq_len(ncol(ZI))])
  g <- numeric(0)
  beta <- rep(0, p)
  beta_prev <- beta
  alpha_prev <- alpha * 0
  fitX <- NULL
  XCS <- NULL
  early_no_cs <- FALSE

  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    alpha_prev <- alpha

    work <- .qgam_extract_working(fit_final, weight_cutoff = weight_cutoff)
    suff <- weighted_residual_suffstats(
      X = X,
      y = work$pseudo_response,
      ZI = ZI,
      weights = work$W_diag,
      n_threads = n_threads,
      block_size = suff_block_size
    )

    fitX <- susieR::susie_ss(
      XtX = suff$XtX, Xty = suff$Xty, yty = suff$yty,
      n = max(n / 2, work$n_eff), L = L,
      scaled_prior_variance = scaled_prior_variance,
      estimate_residual_variance = FALSE,
      residual_variance = work$phi0,
      max_iter = susie.iter,
      estimate_prior_method = "EM",
      coverage = coverage, ...
    )
    rm(suff)

    beta <- clean_coef(stats::coef(fitX)[-1])
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
        XCS <- matrix(noncs_res, ncol = 1)
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
      if (is.null(dim(XCS))) XCS <- matrix(XCS, ncol = 1)
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

    pred <- .mgcv_predictor_data(Z, XCS_refit)
    Data <- cbind(response_info$data, pred)
    fit_final <- .qgam_fit_explicit(
      response = response_info$response,
      rhs = colnames(pred),
      data = Data,
      qu = qu,
      lsig = NULL,
      qgam_control = qgam_control,
      qgam_argGam = qgam_argGam,
      qgam_discrete = qgam_discrete
    )

    alpha <- clean_coef(stats::coef(fit_final)[seq_len(ncol(ZI))])
    err <- max(sqrt(mean((beta - beta_prev)^2)),
               sqrt(mean((alpha - alpha_prev)^2)))
    g[iter] <- err

    if (verbose) {
      cat(sprintf("Iteration %d: err = %.3e, n_eff = %.1f, lsig = %.4f\n",
                  iter, err, work$n_eff, work$lsig))
    }

    if (err < max.eps && iter > min.iter) {
      if (verbose) cat("Converged!\n")
      break
    }
  }

  MainIndex <- if (is.null(fitX)) NULL else Identifying_MainEffect(fitX, colnames(X))
  if (!is.null(XCS)) {
    pred <- .mgcv_predictor_data(Z, XCS)
    Data <- cbind(response_info$data, pred)
    fit_final <- .qgam_fit_explicit(
      response = response_info$response,
      rhs = colnames(pred),
      data = Data,
      qu = qu,
      lsig = NULL,
      qgam_control = qgam_control,
      qgam_argGam = qgam_argGam,
      qgam_discrete = qgam_discrete
    )
  }

  G <- tryCatch(summary(fit_final)$p.table, error = function(e) NULL)
  if (!is.null(G)) MainIndex <- safe_add_p(MainIndex, G)
  fit_final <- clean_model_environment(fit_final)

  if (verbose && length(g)) {
    plot(g, type = "o", col = "black", pch = 16,
         xlab = "Iteration",
         ylab = "Max Parameter Change",
         main = "Convergence Trace (qgam)")
    for (i in seq_along(g)) {
      graphics::text(x = i, y = g[i],
                     labels = formatC(g[i], format = "e", digits = 1),
                     pos = 3, cex = 0.7, col = "red")
    }
  }

  last_err <- if (length(g)) utils::tail(g, 1) else Inf
  list(
    iter = if (exists("iter")) iter else 0,
    error = g,
    converged = early_no_cs || (exists("iter") && iter < max.iter && last_err < max.eps),
    fitX = fitX,
    fitJoint = fit_final,
    main_index = MainIndex,
    JointCoef = G,
    n_eff = if (exists("work")) work$n_eff else NA_real_,
    qu = qu,
    lsig = .qgam_lsig(fit_final),
    residual_variance = exp(.qgam_lsig(fit_final)),
    early_no_cs = early_no_cs
  )
}
