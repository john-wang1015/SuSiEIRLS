#' SuSiE_IRLS: SuSiE for non-Gaussian outcomes via IRLS + SuSiE
#'
#' Iterates between IRLS-based pseudo-response construction and SuSiE fine-mapping
#' on covariate-adjusted data. At each outer iteration the algorithm:
#' 1) fits a GLM to obtain pseudo responses and weights, 2) projects out covariates,
#' 3) runs SuSiE on X, 4) refits a GLM with selected variables, and checks convergence.
#'
#' For survival outcomes, pass a `survival::Surv` object as `y`. The Cox path
#' builds score-based sufficient statistics and runs SuSiE-SS without IRLS
#' weights.
#'
#' @param X An n by p numeric matrix of predictors.
#' @param y Response vector, or a `survival::Surv` object for Cox PH.
#' @param Z An n by q matrix or vector of covariates. If NULL, only an intercept is used.
#' @param family A GLM or mgcv family object, such as
#'   `binomial(link = "logit")` or `mgcv::nb(theta = NULL)`. The strings
#'   `"negbin"`, `"zip"`, and `"ocat"` select their corresponding paths.
#'   Cumulative-link models use `"clm_logit"`, `"clm_probit"`,
#'   `"clm_cloglog"`, `"clm_loglog"`, or `"clm_cauchit"`. An
#'   `mgcv::ocat(R = )` family object is also accepted. This argument is ignored
#'   when `y` is a `Surv` object.
#' @param mgcv_model Either `NULL`, `"gam"`, or `"bam"` for ordinary GLM and
#'   ZIP refits. `NULL` uses `gam` when `n < 50000` and `bam` otherwise.
#' @param L Number of single effects in SuSiE. Default 10.
#' @param L.init Number of SNPs used in the initial low-dimensional warm start.
#'   Default 1.
#' @param max.iter Maximum outer iterations. Default 15.
#' @param min.iter Minimum iterations before convergence checks and prior
#'   variance updates. Default 3.
#' @param max.eps Convergence threshold on max parameter change. Default 1e-5.
#' @param susie.iter Maximum iterations for each SuSiE fit. Default 30.
#' @param verbose Logical flag for progress printing. Default TRUE.
#' @param n_threads Integer number of threads for internal parallel blocks. Default 4.
#' @param coverage Credible set coverage level in SuSiE. Default 0.9.
#' @param estimate_residual_variance Logical for SuSiE residual variance estimation. Default TRUE.
#' @param residual_variance Initial or fixed residual variance. Default 0.5.
#' @param residual_variance_lowerbound Lower bound when estimating residual variance.
#'   Default 0.1.
#' @param residual_variance_upperbound Upper bound when estimating residual variance.
#'   If NULL, defaults to 1.
#' @param prior_variance Prior variance supplied to SuSiE after the warm-up
#'   iterations. The warm-up value is fixed at 2. Default 1.
#' @param estimate_prior_variance Logical. If TRUE, estimate the prior variance
#'   with SuSiE's `"optim"` method after `min.iter`; otherwise keep
#'   `prior_variance` fixed after `min.iter`.
#' @param weight_cutoff Quantile in (0, 0.05) to clip extreme IRLS weights. Default 0.005.
#' @param ridge Diagonal ridge added to the Cox information matrix for positive
#'   definiteness. Used only in the Cox path. Default 1e-6.
#' @param zip_theta Optional raw two-parameter `theta` vector passed to
#'   `mgcv::ziP(theta = )` when `family = "zip"`. If NULL, mgcv estimates it.
#' @param zip_b Non-negative `b` parameter passed to `mgcv::ziP()` when
#'   `family = "zip"`. Default 0.
#' @param zip_info Curvature used by the `mgcv::ziP()` local quadratic.
#'   `"expected"` uses Fisher information; `"observed"` uses the
#'   observed Hessian.
#' @param init_cor_method Deprecated and ignored. The greedy warm start now
#'   ranks variables by `abs(crossprod(X, residual))`.
#' @param refit_noncs Logical. If TRUE, add a one-dimensional non-CS residual
#'   summary variable to the refit model when the current credible-set summary
#'   leaves enough posterior mean variation outside the CS terms. This variable
#'   is used only to improve the next linear predictor estimate and is not
#'   reported as a credible set. Default TRUE.
#' @param noncs_var Minimum non-CS variance fraction required to add the
#'   non-CS residual summary variable. For example, `noncs_var = 0.2`
#'   adds it when the CS summary explains less than 80% of the posterior mean
#'   linear predictor variance. Default 0.2.
#' @param noncs_max_abs_cor Maximum allowed absolute correlation between a
#'   non-CS term and applicable existing refit covariates, including `Z` and
#'   credible-set summary variables. Correlations equal to or above the
#'   threshold are rejected. The effective threshold is capped at 0.9;
#'   smaller user-supplied values remain effective. Default 0.9.
#' @param scale_data Logical. If TRUE, standardize `X` with
#'   `SuSiE4I::large_scale()` and center and scale non-binary columns of `Z`;
#'   binary columns of `Z` remain on their original scale. If FALSE, `X` and
#'   non-binary columns of `Z` are assumed to already
#'   be standardized. Default TRUE.
#' @param suff_block_size Row block size for weighted sufficient-statistic
#'   crossproducts. Larger values can be faster for small-to-moderate p when
#'   memory is sufficient. Default 10000.
#' @param ... Additional arguments passed to the SuSiE fitting routine.
#'
#' @return A list containing the main-effect SuSiE fit, final joint model,
#'   discovery and coefficient summaries, convergence information, and the
#'   prior variance used at each outer iteration. `diagnostics` contains the
#'   number of outer iterations, final convergence eps, and runtime in seconds.
#'
#' @importFrom stats var lm coef glm
#' @importFrom susieR susie_ss coef.susie
#' @importFrom survival coxph Surv
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixCor
#' @importFrom mgcv gam bam nb tw betar scat ziP
#' @importFrom ordinal clm
#' @export
SuSiE_IRLS <- function(X, Z = NULL, y,
                       family = binomial(link = "logit"),
                       mgcv_model = NULL,
                       n_threads = 4, L = 10, coverage = 0.9,
                       estimate_residual_variance = TRUE, residual_variance = 0.5,
                       residual_variance_lowerbound = 0.1,
                       residual_variance_upperbound = NULL,
                       prior_variance = 1,
                       estimate_prior_variance = TRUE,
                       max.iter = 15, max.eps = 1e-5, min.iter = 3,
                       weight_cutoff = 0.005,
                       susie.iter = 30,
                       ridge = 1e-6,
                       zip_theta = NULL,
                       zip_b = 0,
                       zip_info = c("expected", "observed"),
                       L.init = 1,
                       init_cor_method = NULL,
                       refit_noncs = TRUE,
                       noncs_var = 0.2,
                       noncs_max_abs_cor = 0.9,
                       scale_data = TRUE,
                       suff_block_size = 10000L,
                       verbose = TRUE, ...) {

  # ---- basic checks ----
  if (is.null(X)) stop("X must not be NULL.")
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.")
  if (ncol(X) == 0) stop("X has zero columns.")
  if (is.null(colnames(X))) colnames(X) <- paste0("X", seq_len(ncol(X)))

  n <- nrow(X)

  if (is.null(y)) stop("y must be provided.")
  if (NROW(y) != n) stop("Length(y) must equal nrow(X).")

  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (!is.numeric(Z)) stop("Z must be numeric.")
    if (nrow(Z) != n) stop("nrow(Z) must equal nrow(X).")
    if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(ncol(Z)))
  }

  if (!is.logical(scale_data) || length(scale_data) != 1L || is.na(scale_data)) {
    stop("scale_data must be TRUE or FALSE.")
  }

  if (!is.numeric(weight_cutoff) || length(weight_cutoff) != 1L || !is.finite(weight_cutoff)) {
    stop("weight_cutoff must be a finite numeric scalar.")
  }
  if (weight_cutoff <= 0) weight_cutoff <- 1e-6
  if (weight_cutoff >= 0.05) weight_cutoff <- 0.049
  noncs_max_abs_cor <- validate_noncs_max_abs_cor(noncs_max_abs_cor)

  # ---- optional standardization ----
  if (isTRUE(scale_data)) {
    x_dimnames <- dimnames(X)
    X <- SuSiE4I:::large_scale(X, n_threads = n_threads, center = TRUE, scale = TRUE)
    X <- as.matrix(X)
    dimnames(X) <- x_dimnames

    if (!is.null(Z)) {
      z_dimnames <- dimnames(Z)
      Z <- apply(Z, 2L, function(x) {
        vals <- unique(stats::na.omit(x))
        if (length(vals) <= 2L) return(as.numeric(x))
        as.numeric(scale(x))
      })
      Z <- as.matrix(Z)
      dimnames(Z) <- z_dimnames
    }
  } else if (isTRUE(verbose)) {
    cat("scale_data = FALSE: assuming X is column-centered and standardized; non-binary columns of Z should also be standardized.\n")
  }

  # ---- helpers ----
  family_string <- if (is.character(family) && length(family) == 1L) {
    tolower(family)
  } else {
    NULL
  }
  is_negbin_flag <- !is.null(family_string) &&
    family_string %in% c("negbin", "nb", "negative.binomial")
  is_zip_flag <- (!is.null(family_string) &&
    family_string %in% c("zip", "zero.inflated.poisson",
                         "zero_inflated_poisson", "zero inflated poisson",
                         "zero-inflated-poisson", "zeroinflatedpoisson")) ||
    .zip_is_family(family)
  clm_links <- c("logit", "probit", "cloglog", "loglog", "cauchit")
  clm_link <- NULL
  if (!is.null(family_string) && startsWith(family_string, "clm_")) {
    clm_link <- sub("^clm_", "", family_string)
    if (!clm_link %in% clm_links) {
      stop("Unsupported CLM family. Use clm_logit, clm_probit, clm_cloglog, clm_loglog, or clm_cauchit.")
    }
  }
  if (identical(family_string, "clm")) {
    stop("family = 'clm' is incomplete. Specify clm_logit, clm_probit, clm_cloglog, clm_loglog, or clm_cauchit.")
  }
  is_clm_flag <- !is.null(clm_link)
  is_ocat_flag <- (!is.null(family_string) &&
    family_string %in% c("ocat", "ordinal", "ordered", "ordered.categorical",
                         "ordered_categorical")) ||
    .ocat_is_family(family)
  # Cox is identified by a Surv-typed response; family is then ignored.
  is_cox_flag <- inherits(y, "Surv")
  estimate_prior_variance <- .validate_estimate_prior_variance(
    estimate_prior_variance
  )
  prior_variance <- .validate_prior_variance(prior_variance)
  zip_info <- match.arg(zip_info)
  rv_upper_default <- if (is.null(residual_variance_upperbound)) 1 else residual_variance_upperbound
  suff_block_size <- validate_suff_block_size(suff_block_size)

  # ---- dispatch ----
  if (is_cox_flag) {
    # Surv carries time and status; split and hand off to Run_Cox.
    surv_time   <- as.numeric(y[, 1])
    surv_status <- as.integer(y[, 2])
    return(
      Run_Cox(
        X = X, y = surv_time, status = surv_status, Z = Z,
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie.iter = susie.iter,
        verbose = verbose,
        n_threads = n_threads,
        coverage = coverage,
        prior_variance = prior_variance,
        estimate_prior_variance = estimate_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        ridge = ridge,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  if (is_negbin_flag) {
    return(
      Run_GLM(
        X = X, y = y, Z = Z,
        family = mgcv::nb(theta = NULL),
        mgcv_model = mgcv_model,
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie.iter = susie.iter,
        verbose = verbose,
        n_threads = n_threads,
        coverage = coverage,
        weight_cutoff = weight_cutoff,
        prior_variance = prior_variance,
        estimate_prior_variance = estimate_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  if (is_zip_flag) {
    if (!is.null(zip_theta)) {
      if (!is.numeric(zip_theta) || length(zip_theta) != 2L ||
          any(!is.finite(zip_theta))) {
        stop("zip_theta must be NULL or a finite numeric vector of length 2.")
      }
    }
    if (!is.numeric(zip_b) || length(zip_b) != 1L ||
        !is.finite(zip_b) || zip_b < 0) {
      stop("zip_b must be a non-negative finite numeric scalar.")
    }
    zip_family <- if (is.character(family)) {
      mgcv::ziP(theta = zip_theta, b = zip_b)
    } else {
      family
    }
    return(
      Run_ZIP(
        X = X, y = y, Z = Z,
        family = zip_family,
        mgcv_model = mgcv_model,
        zip_info = zip_info,
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie.iter = susie.iter,
        verbose = verbose,
        n_threads = n_threads,
        coverage = coverage,
        weight_cutoff = weight_cutoff,
        prior_variance = prior_variance,
        estimate_prior_variance = estimate_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  if (is_clm_flag) {
    return(
      Run_CLM(
        X = X, y = y, Z = Z, family = NULL, clm_link = clm_link,
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie.iter = susie.iter,
        verbose = verbose,
        n_threads = n_threads,
        coverage = coverage,
        prior_variance = prior_variance,
        estimate_prior_variance = estimate_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        ridge = ridge,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  if (is_ocat_flag) {
    return(
      Run_OCAT(
        X = X, y = y, Z = Z, family = family, clm_link = "logit",
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie.iter = susie.iter,
        verbose = verbose,
        n_threads = n_threads,
        coverage = coverage,
        prior_variance = prior_variance,
        estimate_prior_variance = estimate_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        ridge = ridge,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        noncs_max_abs_cor = noncs_max_abs_cor,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  if (is.character(family)) {
    stop("Unsupported family string. Use negbin, zip, ocat, a clm_<link> family, or a GLM/mgcv family object.")
  }

  # general GLM (e.g., Poisson, Gaussian with non-default link, etc.)
  return(
    Run_GLM(
      X = X, y = y, Z = Z, family = family,
      mgcv_model = mgcv_model,
      L = L, max.iter = max.iter, min.iter = min.iter, max.eps = max.eps,
      susie.iter = susie.iter, verbose = verbose, n_threads = n_threads,
      coverage = coverage, weight_cutoff = weight_cutoff,
      prior_variance = prior_variance,
      estimate_prior_variance = estimate_prior_variance,
      estimate_residual_variance = estimate_residual_variance,
      residual_variance = residual_variance,
      residual_variance_lowerbound = residual_variance_lowerbound,
      residual_variance_upperbound = rv_upper_default,
      L.init = L.init,
      init_cor_method = init_cor_method,
      refit_noncs = refit_noncs,
      noncs_var = noncs_var,
      noncs_max_abs_cor = noncs_max_abs_cor,
      suff_block_size = suff_block_size,
      ...
    )
  )
}
