Run_Cox <- function(X, y, status, Z = NULL,
                    L, max.iter, min.iter, max.eps, susie.iter, pip.thres = 5e-3,
                    verbose = TRUE, n_threads = 1, coverage = 0.5,
                    estimate_residual_variance = FALSE, scaled_prior_variance = 0.5,
                    residual_variance = 1, ridge = 1e-6, ...) {

  n = length(y)
  p = ncol(X)

  # ============================================
  # Handle Z edge cases
  # ============================================
  if (is.null(Z)) {
    Z = matrix(nrow = n, ncol = 0)
    ZI = matrix(1, nrow = n, ncol = 1)
    colnames(ZI) = "Intercept"

  } else {
    if (is.null(dim(Z))) {
      Z = matrix(Z, ncol = 1)
    }

    if (is.null(colnames(Z))) {
      colnames(Z) = paste0("Z", seq_len(ncol(Z)))
    }

    ZI = cbind(1, Z)
    colnames(ZI)[1] = "Intercept"
  }

  # Survival outcome
  surv_y = survival::Surv(y, status)

  # ============================================
  # Initial Cox fit with covariates only
  # ============================================
  if (ncol(Z) == 0) {
    fit_final = survival::coxph(surv_y ~ 1, ties = "breslow")
  } else {
    fit_final = survival::coxph(surv_y ~ Z, ties = "breslow")
  }
  alpha = coef(fit_final)

  # Initialize tracking variables
  g = c()
  beta = rep(0, p)
  beta_prev = beta
  alpha_prev = alpha * 0

  # ============================================
  # Main iteration loop
  # ============================================
  for (iter in 1:max.iter) {
    beta_prev = beta
    alpha_prev = alpha

    ## ===== Cox score-based sufficient statistics with binary-style projection =====

    # Current linear predictor
    eta = fit_final$linear.predictors

    # Same projection logic as Run_Binary:
    # ZI always contains intercept; ZI is used only for projection.
    q = ncol(ZI)

    XZE = cbind(X, eta, ZI)

    ss  = cox_suffstat(X = XZE, eta = eta, time = y,
                       status = as.integer(status), n_threads = n_threads)
    a     = as.numeric(ss$a)
    B     = as.matrix(ss$B)
    XZEty = as.numeric(ss$Xty)        # 这里抠 X 块 = X'M
    n_eff = ss$d

    XZEa = XZE * sqrt(a)
    A    = blockwise_crossprod(XZEa, n_threads = n_threads)
    BtB  = blockwise_crossprod(B,    n_threads = n_threads)
    XZEtXZE = A - BtB
    XZEtXZE = (XZEtXZE + t(XZEtXZE)) / 2
    idxX = seq_len(p)
    idxE = p + 1L
    idxZ = p + 1L + seq_len(q)

    # 信息块
    XtX = XZEtXZE[idxX, idxX, drop = FALSE]
    XtE = XZEtXZE[idxX, idxE, drop = FALSE]   # X'W eta
    XtZ = XZEtXZE[idxX, idxZ, drop = FALSE]
    EtZ = XZEtXZE[idxE, idxZ, drop = FALSE]
    ZtZ = XZEtXZE[idxZ, idxZ, drop = FALSE]
    ZtX = XZEtXZE[idxZ, idxX, drop = FALSE]
    ZtE = XZEtXZE[idxZ, idxE, drop = FALSE]

    XtM = XZEty[idxX]                          # X'M，天然垂直于 Z

    # 对 (X, eta) 块整体投影消 Z
    diag(ZtZ) = diag(ZtZ) + ridge
    Zinv_ZtX = solve(ZtZ, ZtX)
    Zinv_ZtE = solve(ZtZ, ZtE)

    XtX_proj = XtX - matrixMultiply(XtZ, Zinv_ZtX)
    XtE_proj = as.vector(XtE - matrixVectorMultiply(XtZ, Zinv_ZtE))

    # 合成
    Xty = XtE_proj + XtM
    XtX = XtX_proj

    XtX = (XtX + t(XtX)) / 2
    diag(XtX) = diag(XtX) + ridge

    dXtX = diag(XtX)
    zhat = Xty / sqrt(dXtX)
    R = cov2cor(XtX)
    R = R / 2 + t(R) / 2

    # Run SuSiE-RSS on the score statistics.
    fitX <- susieR::susie_rss(
      z = zhat, R = R, n = n_eff, L = L,
      residual_variance = 1,
      estimate_residual_variance = estimate_residual_variance,
      max_iter = susie.iter,
      estimate_prior_method = "EM",
      coverage = coverage, ...
    )

    beta = coef(fitX)[-1]
    beta.cs = group.pip.filter(
      pip.summary = summary(fitX)$var,
      xQTL.cred.thres = coverage,
      xQTL.pip.thres = pip.thres
    )
    pip.alive = beta.cs$ind.keep
    beta[-pip.alive] = 0

    # Extract credible sets using summary information
    CSdt <- summary(fitX)$vars
    cs_indices <- unique(CSdt$cs[CSdt$cs > 0])
    cs_indices = sort(cs_indices)

    if (length(cs_indices) == 0) {
      warning("No credible set detected at iteration ", iter)
      break
    }

    Alpha_filtered <- fitX$alpha * 0
    for (i in cs_indices) {
      vars_in_cs_i <- CSdt$variable[CSdt$cs == i]
      Alpha_filtered[i, vars_in_cs_i] <- fitX$alpha[i, vars_in_cs_i]
    }

    # Align within-CS SNP directions while preserving PIP weights.
    Alpha_filtered <- Alpha_filtered * sign(fitX$mu)
    XCS <- matrixMultiply(X, t(as.matrix(Alpha_filtered)))
    XCS <- XCS[, cs_indices, drop = FALSE]

    if (is.null(dim(XCS))) {
      XCS <- matrix(XCS, ncol = 1)
    }

    colnames(XCS) <- paste0("Main_CS", cs_indices)
    XCS <- as.matrix(XCS)

    # ============================================
    # Refit Cox with selected credible sets
    # ============================================
    if (ncol(Z) == 0) {
      Data = data.frame(XCS)
    } else {
      Data = cbind(Z, XCS)
      Data = as.data.frame(Data)
    }

    fit_final = survival::coxph(surv_y ~ ., data = Data, ties = "breslow")

    # Extract covariate coefficients only
    if (ncol(Z) == 0) {
      alpha = numeric(0)
    } else {
      alpha = coef(fit_final)[seq_len(ncol(Z))]
    }

    # Check convergence
    err = max(
      sqrt(mean((beta - beta_prev)^2)),
      if (length(alpha)) sqrt(mean((alpha - alpha_prev)^2)) else 0
    )
    g[iter] = err

    if (verbose) {
      cat(sprintf("Iteration %d: err = %.3e, events = %d\n", iter, err, n_eff))
      cat("iter", iter, "cs:", cs_indices, "eta sd:", sd(eta),
          "beta nonzero:", which(beta != 0), "\n")
    }

    if (err < max.eps && iter > min.iter) {
      if (verbose) cat("Converged!\n")
      break
    }
  }

  # ============================================
  # Post-processing
  # ============================================
  MainIndex = Identifying_MainEffect(fitX, colnames(X))
  G = summary(fit_final)$coefficients[,-2]
  MainIndex <- safe_add_p(MainIndex, G)

  if (verbose) {
    plot(g, type = "o", col = "black", pch = 16,
         xlab = "Iteration",
         ylab = "Max Parameter Change",
         main = "Convergence Trace (Cox PH, Breslow)")
    for (i in seq_along(g)) {
      text(x = i, y = g[i],
           labels = formatC(g[i], format = "e", digits = 1),
           pos = 3, cex = 0.7, col = "red")
    }
  }

  AA = list(
    iter = iter,
    error = g,
    converged = (iter < max.iter && err < max.eps),
    fitX = fitX,
    fitJoint = fit_final,
    main_index = MainIndex,
    JointCoef = G
  )

  return(AA)
}
