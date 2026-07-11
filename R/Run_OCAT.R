.ocat_is_family <- function(family) {
  inherits(family, "family") &&
    grepl("ordered categorical", tolower(paste(family$family, collapse = " ")),
          fixed = TRUE)
}

.ocat_family_ncat <- function(family) {
  if (!.ocat_is_family(family)) return(NULL)
  if (!is.function(family$getTheta)) return(NULL)
  alpha <- family$getTheta(TRUE)
  if (is.null(alpha) || !length(alpha)) return(NULL)
  length(alpha) + 1L
}

.ocat_prepare_response <- function(y, family = NULL) {
  ncat <- .ocat_family_ncat(family)

  if (is.factor(y)) {
    if (is.null(ncat)) ncat <- nlevels(y)
    if (nlevels(y) != ncat) {
      stop("Number of y levels does not match the ocat family categories.")
    }
    y_ord <- ordered(y, levels = levels(y))
    y_int <- as.integer(y_ord)
  } else {
    y_num <- as.numeric(y)
    if (any(!is.finite(y_num))) stop("Ordinal y must contain finite class labels.")
    if (any(abs(y_num - round(y_num)) > sqrt(.Machine$double.eps))) {
      stop("Ordinal y must be integer class labels or an ordered factor.")
    }
    y_int <- as.integer(y_num)
    if (is.null(ncat)) ncat <- max(y_int)
    if (any(y_int < 1L) || any(y_int > ncat)) {
      stop("Ordinal y class labels must be in 1, ..., R.")
    }
    y_ord <- ordered(y_int, levels = seq_len(ncat))
  }

  if (ncat < 2L) stop("Ordinal y must have at least two categories.")
  list(y = y_ord, y_int = y_int, ncat = ncat)
}

.ocat_predictor_data <- function(Z = NULL, Xextra = NULL, n = NULL) {
  if (is.null(n)) {
    n <- if (!is.null(Z)) nrow(Z) else nrow(Xextra)
  }
  out <- data.frame(row.names = seq_len(n))
  if (!is.null(Z) && ncol(Z) > 0L) {
    Zdf <- as.data.frame(Z)
    colnames(Zdf) <- paste0("Z", seq_len(ncol(Z)))
    out <- cbind(out, Zdf)
  }
  if (!is.null(Xextra) && ncol(as.matrix(Xextra)) > 0L) {
    Xdf <- as.data.frame(Xextra)
    out <- cbind(out, Xdf)
  }
  out
}

.ocat_backtick <- function(x) {
  paste0("`", gsub("`", "``", x, fixed = TRUE), "`")
}

.ocat_formula <- function(rhs) {
  if (!length(rhs)) return(stats::as.formula("y ~ 1"))
  stats::as.formula(paste("y ~", paste(.ocat_backtick(rhs), collapse = " + ")))
}

.ocat_fit_explicit <- function(y, pred, alpha_start = NULL,
                               clm_link = "logit") {
  dat <- cbind(data.frame(y = y), pred)
  rhs <- colnames(pred)
  start <- NULL
  if (!is.null(alpha_start)) {
    start <- c(alpha_start, rep(0, length(rhs)))
  }
  ordinal::clm(.ocat_formula(rhs), data = dat, link = clm_link,
               threshold = "flexible", start = start, model = TRUE)
}

.ocat_linear_predictor <- function(fit, pred, n) {
  if (is.null(fit$beta) || !length(fit$beta)) return(rep(0, n))
  beta <- as.numeric(fit$beta)
  names(beta) <- names(fit$beta)
  Xp <- as.matrix(pred[, names(beta), drop = FALSE])
  as.numeric(CppMatrix::matrixVectorMultiply(Xp, beta))
}

.ocat_nuisance_coef <- function(fit, q) {
  zb <- if (q > 0L && length(fit$beta) >= q) {
    clean_coef(fit$beta[seq_len(q)])
  } else {
    numeric(0)
  }
  c(clean_coef(fit$alpha), zb)
}

.ocat_coef_table <- function(fit) {
  G <- summary(fit)$coefficients
  if (is.null(G) || is.null(dim(G))) return(NULL)
  G
}

.ocat_link_parts <- function(t, clm_link) {
  if (identical(clm_link, "logit")) {
    Fv <- stats::plogis(t)
    fv <- stats::dlogis(t)
    fp <- fv * (1 - 2 * Fv)
  } else if (identical(clm_link, "probit")) {
    Fv <- stats::pnorm(t)
    fv <- stats::dnorm(t)
    fp <- -t * fv
  } else if (identical(clm_link, "cauchit")) {
    Fv <- stats::pcauchy(t)
    fv <- stats::dcauchy(t)
    fp <- -2 * t / (pi * (1 + t^2)^2)
  } else if (identical(clm_link, "cloglog")) {
    et <- exp(pmin(t, 35))
    Fv <- 1 - exp(-et)
    fv <- exp(pmin(t - et, 700))
    fp <- fv * (1 - et)
  } else if (identical(clm_link, "loglog")) {
    ent <- exp(pmin(-t, 35))
    Fv <- exp(-ent)
    fv <- exp(pmin(-t - ent, 700))
    fp <- fv * (ent - 1)
  } else {
    stop("Unsupported clm_link. Use logit, probit, cauchit, cloglog, or loglog.")
  }
  ii <- !is.finite(t)
  if (any(ii)) {
    Fv[ii & t < 0] <- 0
    Fv[ii & t > 0] <- 1
    fv[ii] <- 0
    fp[ii] <- 0
  }
  Fv <- pmin(pmax(Fv, 0), 1)
  fv[!is.finite(fv)] <- 0
  fp[!is.finite(fp)] <- 0
  list(F = Fv, f = fv, fp = fp)
}

.ocat_solve_with_ridge <- function(A, B, ridge = 1e-6) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  rr <- if (is.finite(ridge) && ridge > 0) {
    ridge * c(1, 100, 10000, 1e6)
  } else {
    c(1e-8, 1e-6, 1e-4, 1e-2)
  }
  for (r in rr) {
    Ar <- A
    diag(Ar) <- diag(Ar) + r
    if (qr(Ar, tol = sqrt(.Machine$double.eps))$rank == ncol(Ar)) {
      return(CppMatrix::matrixSolve(Ar, B))
    }
  }
  Ar <- A
  diag(Ar) <- diag(Ar) + utils::tail(rr, 1)
  MASS::ginv(Ar) %*% B
}

.ocat_prob_parts <- function(y_int, eta, alpha, clm_link = "logit",
                             eps = 1e-12) {
  n <- length(y_int)
  K <- length(alpha)
  cuts <- c(-Inf, alpha, Inf)

  lo <- cuts[y_int]
  hi <- cuts[y_int + 1L]
  tl <- lo - eta
  tu <- hi - eta

  lp <- .ocat_link_parts(tl, clm_link = clm_link)
  up <- .ocat_link_parts(tu, clm_link = clm_link)
  Fl <- lp$F
  Fu <- up$F
  fl <- lp$f
  fu <- up$f
  fpl <- lp$fp
  fpu <- up$fp

  pr <- pmax(Fu - Fl, eps)
  A <- fl - fu
  B <- fpu - fpl
  u_eta <- A / pr
  h_eta <- A^2 / pr^2 - B / pr

  D <- matrix(0, n, K)
  C <- matrix(0, n, K)
  E <- matrix(0, n, K)

  ii <- which(y_int <= K)
  if (length(ii)) {
    jj <- y_int[ii]
    D[cbind(ii, jj)] <- fu[ii]
    C[cbind(ii, jj)] <- -fpu[ii]
    E[cbind(ii, jj)] <- fpu[ii]
  }

  ii <- which(y_int > 1L)
  if (length(ii)) {
    jj <- y_int[ii] - 1L
    D[cbind(ii, jj)] <- D[cbind(ii, jj)] - fl[ii]
    C[cbind(ii, jj)] <- C[cbind(ii, jj)] + fpl[ii]
    E[cbind(ii, jj)] <- E[cbind(ii, jj)] - fpl[ii]
  }

  h_eta_th <- sweep(D, 1L, A / pr^2, "*") -
    sweep(C, 1L, 1 / pr, "*")

  list(
    u_eta = u_eta,
    h_eta = h_eta,
    h_eta_th = h_eta_th,
    D = D,
    E = E,
    pr = pr
  )
}

ocat_suffstats <- function(X, y_int, eta, Z, alpha,
                           clm_link = "logit",
                           n_threads = 1, ridge = 1e-6,
                           block_size = 10000L) {
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  n <- nrow(X)
  K <- length(alpha)
  q <- ncol(Z)

  pp <- .ocat_prob_parts(y_int = y_int, eta = eta, alpha = alpha,
                         clm_link = clm_link)
  h <- pmax(as.numeric(pp$h_eta), 1e-8)
  sw <- sqrt(h)

  Xh <- X * sw
  XtX <- SuSiE4I::blockwise_crossprod(Xh, n_threads = n_threads,
                                      block_size = block_size)
  Eh <- matrix(eta * sw, ncol = 1)
  XtE <- CppMatrix::matrixMultiply(Xh, Eh, transA = TRUE)
  XtU <- as.numeric(CppMatrix::matrixMultiply(
    X, matrix(pp$u_eta, ncol = 1), transA = TRUE
  ))
  XtT <- CppMatrix::matrixMultiply(X, pp$h_eta_th, transA = TRUE)

  EtT <- CppMatrix::matrixMultiply(matrix(eta, ncol = 1), pp$h_eta_th,
                                   transA = TRUE)
  EtU <- sum(eta * pp$u_eta)

  if (q > 0L) {
    Zh <- Z * sw
    XtZ <- CppMatrix::matrixMultiply(Xh, Zh, transA = TRUE)
    EtZ <- CppMatrix::matrixMultiply(Eh, Zh, transA = TRUE)
    ZtZ <- CppMatrix::matrixMultiply(Zh, Zh, transA = TRUE)
    ZtT <- CppMatrix::matrixMultiply(Z, pp$h_eta_th, transA = TRUE)
    ZtU <- as.numeric(CppMatrix::matrixMultiply(
      Z, matrix(pp$u_eta, ncol = 1), transA = TRUE
    ))
  }

  TtT <- matrix(0, K, K)
  for (i in seq_len(n)) {
    di <- pp$D[i, ]
    TtT <- TtT + tcrossprod(di) / pp$pr[i]^2
    TtT <- TtT - diag(pp$E[i, ] / pp$pr[i], K)
  }
  TtU <- colSums(pp$D / pp$pr)

  if (q > 0L) {
    HNN <- rbind(
      cbind(ZtZ, ZtT),
      cbind(t(ZtT), TtT)
    )
    HXN <- cbind(XtZ, XtT)
    HNE <- rbind(t(EtZ), t(EtT))
    UN <- c(ZtU, TtU)
  } else {
    HNN <- TtT
    HXN <- XtT
    HNE <- t(EtT)
    UN <- TtU
  }

  HNN_inv_HNX <- .ocat_solve_with_ridge(HNN, t(HXN), ridge = ridge)
  HNN_inv_HNE <- .ocat_solve_with_ridge(HNN, HNE, ridge = ridge)
  HNN_inv_UN <- .ocat_solve_with_ridge(HNN, matrix(UN, ncol = 1), ridge = ridge)

  XtX <- XtX - CppMatrix::matrixMultiply(HXN, HNN_inv_HNX)
  XtE <- as.numeric(XtE - CppMatrix::matrixMultiply(HXN, HNN_inv_HNE))
  XtU <- XtU - as.numeric(CppMatrix::matrixMultiply(HXN, HNN_inv_UN))
  Xty <- XtE + XtU

  XtX <- (XtX + t(XtX)) / 2
  diag(XtX) <- diag(XtX) + ridge

  dimnames(XtX) <- list(colnames(X), colnames(X))
  names(Xty) <- colnames(X)

  list(
    XtX = XtX,
    Xty = Xty,
    yty = n - 1,
    n_eff = n,
    min_pr = min(pp$pr),
    min_h = min(h),
    med_h = stats::median(h),
    max_h = max(h),
    EtU = EtU
  )
}

#' Ordered-categorical score IRLS-SuSiE path
#' @inheritParams SuSiE_IRLS
#' @param Z An n by q matrix or vector of covariates. If NULL, only the
#'   ordered-categorical threshold nuisance parameters are projected out.
#' @param family An ordered-categorical family object, typically
#'   \code{mgcv::ocat(R = )}, or the string \code{"clm"} / \code{"ocat"} /
#'   \code{"ordinal"}.
#' @param clm_link Link used by \code{ordinal::clm()} and the local
#'   cumulative-link score statistics. Supported values are \code{"logit"},
#'   \code{"probit"}, \code{"cauchit"}, \code{"cloglog"}, and
#'   \code{"loglog"}.
#' @param ridge Diagonal ridge added to the nuisance and projected information
#'   matrices for numerical stability.
#' @export
Run_OCAT <- function(X, y, Z = NULL,
                     family = NULL,
                     clm_link = c("logit", "probit", "cauchit",
                                  "cloglog", "loglog"),
                     L, max.iter, min.iter, max.eps, susie.iter,
                     verbose = TRUE, n_threads = 1, coverage = 0.9,
                     estimate_residual_variance = TRUE,
                     residual_variance = 0.5, scaled_prior_variance = 1,
                     residual_variance_lowerbound = 0.1,
                     residual_variance_upperbound = 1,
                     ridge = 1e-6,
                     L.init = 1,
                     init_cor_method = NULL,
                     refit_noncs = TRUE,
                     noncs_var = 0.2,
                     suff_block_size = 10000L, ...) {

  n <- NROW(y)
  p <- ncol(X)
  clm_link <- match.arg(clm_link)
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(p))
  suff_block_size <- validate_suff_block_size(suff_block_size)
  y_info <- .ocat_prepare_response(y, family = family)
  y <- y_info$y
  y_int <- y_info$y_int

  if (is.null(Z)) {
    Z <- matrix(nrow = n, ncol = 0)
  } else {
    if (is.null(dim(Z))) Z <- matrix(Z, ncol = 1)
    if (nrow(Z) != n) stop("nrow(Z) must equal nrow(X).")
    colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  }
  q <- ncol(Z)

  pred <- .ocat_predictor_data(Z = Z, n = n)
  fit_final <- .ocat_fit_explicit(y, pred, clm_link = clm_link)
  eta <- .ocat_linear_predictor(fit_final, pred, n = n)

  alpha <- .ocat_nuisance_coef(fit_final, q = q)
  g <- numeric(0)
  beta <- rep(0, p)
  beta_prev <- beta
  alpha_prev <- alpha * 0
  early_no_cs <- FALSE
  fitX <- NULL
  XCS <- NULL
  stat <- NULL

  for (iter in seq_len(max.iter)) {
    beta_prev <- beta
    alpha_prev <- alpha

    stat <- ocat_suffstats(
      X = X, y_int = y_int, eta = eta, Z = Z, alpha = fit_final$alpha,
      clm_link = clm_link,
      n_threads = n_threads, ridge = ridge, block_size = suff_block_size
    )

    fitX <- susieR::susie_ss(
      XtX = stat$XtX, Xty = stat$Xty, yty = stat$yty,
      n = n, L = L,
      scaled_prior_variance = scaled_prior_variance,
      estimate_residual_variance = estimate_residual_variance,
      residual_variance = residual_variance,
      residual_variance_lowerbound = residual_variance_lowerbound,
      residual_variance_upperbound = residual_variance_upperbound,
      max_iter = susie.iter,
      estimate_prior_method = "optim",
      coverage = coverage, ...
    )

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
          XCS = XCS, noncs_var = noncs_var
        )
        if (!is.null(noncs_term)) {
          XCS_refit <- cbind(XCS_refit, Main_CS_noncs = noncs_term)
        }
      }
    }

    pred_refit <- .ocat_predictor_data(Z = Z, Xextra = XCS_refit, n = n)
    fit_final <- .ocat_fit_explicit(
      y, pred_refit, alpha_start = fit_final$alpha, clm_link = clm_link
    )
    eta <- .ocat_linear_predictor(fit_final, pred_refit, n = n)

    alpha <- .ocat_nuisance_coef(fit_final, q = q)
    err <- max(sqrt(mean((beta - beta_prev)^2)),
               sqrt(mean((alpha - alpha_prev)^2)))
    g[iter] <- err

    if (verbose) {
      cat(sprintf(
        "Iteration %d: err = %.3e, sigma2 = %.3g, min_pr = %.3g, med_h = %.3g\n",
        iter, err, fitX$sigma2, stat$min_pr, stat$med_h
      ))
    }

    if (err < max.eps && iter > min.iter) {
      if (verbose) cat("Converged!\n")
      break
    }
  }

  MainIndex <- if (is.null(fitX)) NULL else Identifying_MainEffect(fitX, colnames(X))
  if (early_no_cs) {
    G <- .ocat_coef_table(fit_final)
    MainIndex <- safe_add_p(MainIndex, G)
    fit_final <- clean_model_environment(fit_final)
    return(list(
      iter = iter,
      error = g,
      converged = FALSE,
      fitX = fitX,
      fitJoint = fit_final,
      main_index = MainIndex,
      JointCoef = G,
      clm_link = clm_link,
      n_eff = if (!is.null(stat)) stat$n_eff else NA_real_
    ))
  }

  if (!is.null(XCS)) {
    pred <- .ocat_predictor_data(Z = Z, Xextra = XCS, n = n)
    fit_final <- .ocat_fit_explicit(
      y, pred, alpha_start = fit_final$alpha, clm_link = clm_link
    )
  }

  G <- .ocat_coef_table(fit_final)
  MainIndex <- safe_add_p(MainIndex, G)
  fit_final <- clean_model_environment(fit_final)

  if (verbose && length(g)) {
    plot(g, type = "o", col = "black", pch = 16,
         xlab = "Iteration",
         ylab = "Max Parameter Change",
         main = sprintf("Convergence Trace (CLM, %s)", clm_link))
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
    converged = exists("iter") && iter < max.iter && last_err < max.eps,
    fitX = fitX,
    fitJoint = fit_final,
    main_index = MainIndex,
    JointCoef = G,
    clm_link = clm_link,
    n_eff = if (!is.null(stat)) stat$n_eff else NA_real_,
    score_diagnostics = if (!is.null(stat)) {
      stat[c("min_pr", "min_h", "med_h", "max_h")]
    } else {
      NULL
    }
  )
}

#' @rdname Run_OCAT
#' @export
Run_CLM <- Run_OCAT
