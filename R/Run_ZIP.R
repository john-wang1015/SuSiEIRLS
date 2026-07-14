.zip_is_family <- function(family) {
  inherits(family, "family") &&
    grepl("zero inflated poisson", tolower(paste(family$family, collapse = " ")),
          fixed = TRUE)
}

.zip_validate_family <- function(family) {
  if (!.zip_is_family(family)) {
    stop("family must be mgcv::ziP() for the ZIP path.")
  }
  if (!is.function(family$Dd) || !is.function(family$getTheta)) {
    stop("mgcv::ziP() family must provide Dd() and getTheta().")
  }
  if (!identical(family$link, "identity")) {
    stop("Only the identity ziP link is supported by mgcv::ziP().")
  }
  invisible(TRUE)
}

.zip_prepare_response <- function(y) {
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop("ziP y must contain finite counts.")
  if (any(y < 0)) stop("ziP y must be non-negative.")
  if (any(abs(y - round(y)) > sqrt(.Machine$double.eps))) {
    stop("ziP y must contain integer counts.")
  }
  if (min(y) == 0 && max(y) == 1) {
    stop("Using ziP for binary data makes no sense.")
  }
  list(data = data.frame(y = y), response = "y", n = length(y))
}

.zip_theta <- function(fit_or_family, transformed = FALSE) {
  family <- if (inherits(fit_or_family, "family")) fit_or_family else fit_or_family$family
  if (!is.function(family$getTheta)) return(NULL)
  theta <- family$getTheta(transformed)
  if (is.null(theta) || any(!is.finite(theta))) return(NULL)
  as.numeric(theta)
}

.zip_extract_working <- function(fit, weight_cutoff = 0.005,
                                 zip_info = c("expected", "observed")) {
  zip_info <- match.arg(zip_info)
  family <- fit$family
  .zip_validate_family(family)

  gamma <- as.numeric(fit$linear.predictors)
  y_work <- as.numeric(fit$y)
  theta <- .zip_theta(fit, transformed = FALSE)
  if (is.null(theta)) stop("Could not extract raw ziP theta from the fitted model.")

  w0 <- fit$prior.weights
  if (is.null(w0)) w0 <- rep(1, length(gamma))
  dd <- family$Dd(y_work, gamma, theta, wt = as.numeric(w0), level = 0)

  score <- -0.5 * as.numeric(dd$Dmu)
  h_obs <- 0.5 * as.numeric(dd$Dmu2)
  h_exp <- 0.5 * as.numeric(dd$EDmu2)
  h <- if (identical(zip_info, "expected")) h_exp else h_obs

  bad <- !is.finite(gamma) | !is.finite(score) | !is.finite(h) | h <= 0
  if (mean(bad) > 0.9) stop("Too many invalid ziP working observations.")
  h[bad] <- NA_real_
  W_diag <- robust_weight(h, cutoff = weight_cutoff)

  pseudo_response <- numeric(length(gamma))
  good <- is.finite(gamma) & is.finite(score) & is.finite(W_diag) & W_diag > 0
  if (mean(!good) > 0.9) stop("All ziP working weights are zero.")
  pseudo_response[good] <- gamma[good] + score[good] / W_diag[good]
  W_diag[!good] <- 0

  weight_denom <- sum(W_diag^2)
  if (!is.finite(weight_denom) || weight_denom <= 0) {
    stop("All ziP working weights are zero.")
  }

  list(
    pseudo_response = pseudo_response,
    W_diag = W_diag,
    n_eff = (sum(W_diag)^2) / weight_denom,
    theta_raw = theta,
    theta = .zip_theta(fit, transformed = TRUE),
    min_weight = min(W_diag[W_diag > 0]),
    med_weight = stats::median(W_diag[W_diag > 0]),
    max_weight = max(W_diag)
  )
}

#' Zero-inflated Poisson IRLS-SuSiE path
#' @inheritParams SuSiE_IRLS
#' @param family An `mgcv::ziP()` family object.
#' @param zip_info Curvature used for the local quadratic. `"expected"`
#'   uses the Fisher information returned by `mgcv::ziP()`; `"observed"`
#'   uses the observed Hessian.
#' @importFrom mgcv gam bam ziP
#' @keywords internal
#' @noRd
Run_ZIP <- function(X, y, Z = NULL, weight_cutoff = 0.005,
                    family = mgcv::ziP(),
                    mgcv_model = NULL,
                    zip_info = c("expected", "observed"),
                    L, max.iter, min.iter, max.eps, susie.iter,
                    verbose = TRUE, n_threads = 1, coverage = 0.9,
                    estimate_residual_variance = TRUE,
                    residual_variance = 0.5, prior_variance = 1,
                    estimate_prior_variance = TRUE,
                    residual_variance_lowerbound = 0.1,
                    residual_variance_upperbound = 1,
                    L.init = 1,
                    init_cor_method = NULL,
                    refit_noncs = TRUE,
                    noncs_var = 0.2,
                    noncs_max_abs_cor = 0.9,
                    suff_block_size = 10000L, ...) {

  run_start <- proc.time()[["elapsed"]]
  n <- NROW(y)
  p <- ncol(X)
  estimate_prior_variance <- .validate_estimate_prior_variance(
    estimate_prior_variance
  )
  prior_variance <- .validate_prior_variance(prior_variance)
  zip_info <- match.arg(zip_info)
  suff_block_size <- validate_suff_block_size(suff_block_size)
  .zip_validate_family(family)
  response_info <- .zip_prepare_response(y)
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

  fit_final <- .mgcv_greedy_warm_start(
    X = X, response_info = response_info, Z = Z, family = family,
    L.init = L.init, init_cor_method = init_cor_method,
    mgcv_model = mgcv_model
  )

  alpha <- clean_coef(stats::coef(fit_final)[seq_len(ncol(ZI))])
  g <- numeric(0)
  beta <- rep(0, p)
  beta_prev <- beta
  alpha_prev <- alpha * 0
  fitX <- NULL
  XCS <- NULL
  V_main <- numeric(0)
  early_no_cs <- FALSE
  work <- NULL

  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    alpha_prev <- alpha

    work <- .zip_extract_working(
      fit_final, weight_cutoff = weight_cutoff, zip_info = zip_info
    )
    suff <- weighted_residual_suffstats(
      X = X,
      y = work$pseudo_response,
      ZI = ZI,
      weights = work$W_diag,
      n_threads = n_threads,
      block_size = suff_block_size
    )

    updateV <- if (iter <= min.iter) 2 else 3
    V_main[iter] <- updateV
    fitX <- susieR::susie_ss(
      XtX = suff$XtX, Xty = suff$Xty, yty = suff$yty,
      n = max(n / 2, work$n_eff), L = L,
      scaled_prior_variance = updateV,
      estimate_prior_variance = iter > min.iter &&
        isTRUE(estimate_prior_variance),
      estimate_residual_variance = estimate_residual_variance,
      residual_variance = residual_variance,
      residual_variance_lowerbound = residual_variance_lowerbound,
      residual_variance_upperbound = residual_variance_upperbound,
      max_iter = susie.iter,
      estimate_prior_method = "optim",
      coverage = coverage, ...
    )
    rm(suff)

    beta <- clean_coef(stats::coef(fitX)[-1L])
    CSdt <- summary(fitX)$vars
    cs_indices <- sort(unique(CSdt$cs[CSdt$cs > 0]))

    if (!length(cs_indices)) {
      if (iter <= min.iter) {
        noncs_res <- build_no_cs_noncs_refit_term(
          X, fitX, cor_design = Z,
          noncs_max_abs_cor = noncs_max_abs_cor
        )
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
      XCS <- CppMatrix::matrixMultiply(X, as.matrix(Alpha_filtered), transB = TRUE)
      XCS <- XCS[, cs_indices, drop = FALSE]
      if (is.null(dim(XCS))) XCS <- matrix(XCS, ncol = 1)
      colnames(XCS) <- paste0("Main_CS", cs_indices)
      XCS_refit <- XCS

      if (isTRUE(refit_noncs)) {
        noncs_term <- build_noncs_refit_term(
          X = X, fitX = fitX, CSdt = CSdt, cs_indices = cs_indices,
          XCS = XCS, noncs_var = noncs_var,
          noncs_max_abs_cor = noncs_max_abs_cor, cor_design = Z
        )
        if (!is.null(noncs_term)) {
          XCS_refit <- cbind(XCS_refit, Main_CS_noncs = noncs_term)
        }
      }
    }

    pred <- .mgcv_predictor_data(Z, XCS_refit)
    Data <- cbind(response_info$data, pred)
    fit_final <- .mgcv_fit_explicit(
      response_info$response, colnames(pred), Data,
      family, mgcv_model = mgcv_model
    )


    alpha <- clean_coef(stats::coef(fit_final)[seq_len(ncol(ZI))])
    err <- max(sqrt(mean((beta - beta_prev)^2)),
               sqrt(mean((alpha - alpha_prev)^2)))
    g[iter] <- err

    if (verbose) {
      cat(sprintf(
        "Iteration %d: err = %.3e, n_eff = %.1f, med_w = %.3g\n",
        iter, err, work$n_eff, work$med_weight
      ))
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
    fit_final <- .mgcv_fit_explicit(
      response_info$response, colnames(pred), Data,
      family, mgcv_model = mgcv_model
    )
  }

  G <- tryCatch(summary(fit_final)$p.table, error = function(e) NULL)
  if (!is.null(G)) MainIndex <- safe_add_p(MainIndex, G)
  fit_final <- clean_model_environment(fit_final)

  if (verbose && length(g)) {
    plot(g, type = "o", col = "black", pch = 16,
         xlab = "Iteration",
         ylab = "Max Parameter Change",
         main = "ZIP Convergence Trace")
    for (i in seq_along(g)) {
      graphics::text(x = i, y = g[i],
                     labels = formatC(g[i], format = "e", digits = 1),
                     pos = 3, cex = 0.7, col = "red")
    }
  }

  list(
    diagnostics = make_diagnostics(
      if (exists("iter")) iter else 0L, g, run_start
    ),
    fitX = fitX,
    fitJoint = fit_final,
    main_index = MainIndex,
    JointCoef = G,
    n_eff = if (!is.null(work)) work$n_eff else NA_real_,
    theta = .zip_theta(fit_final, transformed = TRUE),
    theta_raw = .zip_theta(fit_final, transformed = FALSE),
    zip_info = zip_info,
    prior_variance_main = V_main
  )
}
