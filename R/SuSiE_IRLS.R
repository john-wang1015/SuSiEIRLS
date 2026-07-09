#' SuSiE_IRLS: SuSiE for non-Gaussian outcomes via IRLS + SuSiE
#'
#' Iterates between IRLS-based pseudo-response construction and SuSiE fine-mapping
#' on covariate-adjusted data. At each outer iteration the algorithm:
#' 1) fits a GLM to obtain pseudo responses and weights, 2) projects out covariates,
#' 3) runs SuSiE on X, 4) refits a GLM with selected variables, and checks convergence.
#'
#' For survival outcomes, pass a \code{survival::Surv} object as \code{y}; the
#' algorithm then dispatches to a Cox proportional-hazards fine-mapping routine
#' (\code{Run_Cox}) that builds score-based sufficient statistics and runs
#' SuSiE-SS. No IRLS weights are used in the Cox path.
#'
#' @param X An n by p numeric matrix of predictors.
#' @param y Response vector, or a \code{survival::Surv} object for Cox PH.
#' @param Z An n by q matrix or vector of covariates. If NULL, only an intercept is used.
#' @param family A GLM or mgcv family object (e.g., \code{binomial(link="logit")}
#'   or \code{mgcv::nb(theta = NULL)}), the string \code{"negbin"} for
#'   \code{mgcv::nb(theta = NULL)}, the string \code{"zip"} for
#'   \code{mgcv::ziP()}, or \code{"clm"} /
#'   \code{"ocat"} / \code{mgcv::ocat(R = )} for ordered-categorical
#'   cumulative-link
#'   outcomes. Ignored when \code{y} is a
#'   \code{Surv} object (Cox PH is used).
#' @param logit_method Method for binomial-logit outcomes. \code{"pg"} uses
#'   \code{Run_Binary}; \code{"glm"} uses the general GLM IRLS path.
#' @param clm_link Link for ordered-categorical \code{family = "clm"}
#'   outcomes. Default \code{"logit"}. Other supported values are
#'   \code{"probit"}, \code{"cauchit"}, \code{"cloglog"}, and \code{"loglog"}.
#' @param L Number of single effects in SuSiE. Default 10.
#' @param L.init Number of SNPs used in the initial low-dimensional warm start.
#'   Default 1.
#' @param max.iter Maximum outer iterations. Default 15.
#' @param min.iter Minimum iterations before convergence check. Default 3.
#' @param max.eps Convergence threshold on max parameter change. Default 1e-5.
#' @param susie.iter Maximum iterations for each SuSiE fit. Default 300.
#' @param verbose Logical flag for progress printing. Default TRUE.
#' @param n_threads Integer number of threads for internal parallel blocks. Default 4.
#' @param coverage Credible set coverage level in SuSiE. Default 0.9.
#' @param estimate_residual_variance Logical for SuSiE residual variance estimation. Default TRUE.
#' @param residual_variance Initial or fixed residual variance. Default 0.5.
#' @param residual_variance_lowerbound Lower bound when estimating residual variance.
#'   Default 0.1; the PG-logit final correction uses 0.5 when this argument is not supplied.
#' @param residual_variance_upperbound Upper bound when estimating residual variance.
#'   If NULL, defaults to 1.
#' @param scaled_prior_variance Prior variance for SuSiE single effects. Default 1.
#' @param weight_cutoff Quantile in (0, 0.05) to clip extreme IRLS weights. Default 0.005.
#' @param ridge Diagonal ridge added to the Cox information matrix for positive
#'   definiteness. Used only in the Cox path. Default 1e-6.
#' @param zip_theta Optional raw two-parameter \code{theta} vector passed to
#'   \code{mgcv::ziP(theta = )} when \code{family = "zip"}. If NULL, theta is
#'   estimated by \code{mgcv}.
#' @param zip_b Non-negative \code{b} parameter passed to \code{mgcv::ziP()} when
#'   \code{family = "zip"}. Default 0.
#' @param zip_info Curvature used by the \code{mgcv::ziP()} local quadratic.
#'   \code{"expected"} uses Fisher information; \code{"observed"} uses the
#'   observed Hessian.
#' @param init_cor_method Deprecated and ignored. The greedy warm start now
#'   ranks variables by `abs(crossprod(X, residual))`.
#' @param refit_noncs Logical. If TRUE, add a one-dimensional non-CS residual
#'   summary variable to the refit model when the current credible-set summary
#'   leaves enough posterior mean variation outside the CS terms. This variable
#'   is used only to improve the next linear predictor estimate and is not
#'   reported as a credible set. Default TRUE.
#' @param noncs_var Minimum non-CS variance fraction required to add the
#'   non-CS residual summary variable. For example, \code{noncs_var = 0.2}
#'   adds it when the CS summary explains less than 80% of the posterior mean
#'   linear predictor variance. Default 0.2.
#' @param scale_data Logical. If TRUE, standardize \code{X} with
#'   \code{SuSiE4I::large_scale()} and center/scale non-binary columns of
#'   \code{Z}; binary columns of \code{Z} are left on their original scale. If
#'   FALSE, \code{X} and non-binary columns of \code{Z} are assumed to already
#'   be standardized. Default TRUE.
#' @param suff_block_size Row block size for weighted sufficient-statistic
#'   crossproducts. Larger values can be faster for small-to-moderate p when
#'   memory is sufficient. Default 10000.
#' @param ... Additional arguments passed to the SuSiE fitting routine.
#'
#' @return A list with elements:
#'   \item{iter}{Number of iterations completed}
#'   \item{error}{Max parameter change per iteration}
#'   \item{converged}{Logical}
#'   \item{fitX}{SuSiE results for main effects}
#'   \item{fitJoint}{Final model fit}
#'   \item{main_index}{Summary table of identified effects with PIPs and p-values}
#'   \item{JointCoef}{Coefficient table from the final joint model}
#'
#' @importFrom stats var lm coef glm
#' @importFrom susieR susie_ss coef.susie
#' @importFrom survival coxph Surv
#' @importFrom CppMatrix matrixMultiply matrixVectorMultiply matrixCor
#' @importFrom mgcv gam bam nb tw betar scat ziP
#' @importFrom ordinal clm
#' @importFrom SuSiE4I blockwise_crossprod large_scale
#' @export
SuSiE_IRLS <- function(X, Z = NULL, y = NULL,
                       family = binomial(link = "logit"),
                       n_threads = 4, L = 10, coverage = 0.9,
                       estimate_residual_variance = TRUE, residual_variance = 0.5,
                       residual_variance_lowerbound = 0.1,
                       residual_variance_upperbound = NULL,
                       scaled_prior_variance = 1,
                       max.iter = 15, max.eps = 1e-5, min.iter = 4,
                       weight_cutoff = 0.005,
                       susie.iter = 30,
                       ridge = 1e-6,
                       zip_theta = NULL,
                       zip_b = 0,
                       zip_info = c("expected", "observed"),
                       logit_method = c("pg", "glm"),
                       clm_link = c("logit", "probit", "cauchit",
                                    "cloglog", "loglog"),
                       L.init = 1,
                       init_cor_method = NULL,
                       refit_noncs = TRUE,
                       noncs_var = 0.2,
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

  # ---- optional standardization ----
  if (isTRUE(scale_data)) {
    x_dimnames <- dimnames(X)
    X <- SuSiE4I::large_scale(X, n_threads = n_threads, center = TRUE, scale = TRUE)
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
  is_logit_binomial <- function(fam) {
    inherits(fam, "family") && identical(fam$family, "binomial") && identical(fam$link, "logit")
  }
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
  is_clm_flag <- !is.null(family_string) &&
    family_string %in% c("clm", "cumulative.link", "cumulative_link",
                         "cumulative")
  is_ocat_flag <- (!is.null(family_string) &&
    family_string %in% c("ocat", "ordinal", "ordered", "ordered.categorical",
                         "ordered_categorical")) ||
    .ocat_is_family(family)
  # Cox is identified by a Surv-typed response; family is then ignored.
  is_cox_flag <- inherits(y, "Surv")
  logit_method <- match.arg(logit_method)
  clm_link <- match.arg(clm_link)
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
        scaled_prior_variance = scaled_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        ridge = ridge,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
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
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie.iter = susie.iter,
        verbose = verbose,
        n_threads = n_threads,
        coverage = coverage,
        weight_cutoff = weight_cutoff,
        scaled_prior_variance = scaled_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
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
        scaled_prior_variance = scaled_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  if (is_clm_flag || is_ocat_flag) {
    return(
      Run_OCAT(
        X = X, y = y, Z = Z, family = family, clm_link = clm_link,
        L = L,
        max.iter = max.iter,
        min.iter = min.iter,
        max.eps = max.eps,
        susie.iter = susie.iter,
        verbose = verbose,
        n_threads = n_threads,
        coverage = coverage,
        scaled_prior_variance = scaled_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = residual_variance_lowerbound,
        residual_variance_upperbound = rv_upper_default,
        ridge = ridge,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  if (is.character(family)) {
    stop("Unsupported family string. Use \"negbin\", \"zip\", \"clm\", \"ocat\", or a GLM/mgcv family object.")
  }

  if (is_logit_binomial(family) && identical(logit_method, "pg")) {
    pg_logit_rv_lower <- if (missing(residual_variance_lowerbound)) 0.5 else residual_variance_lowerbound
    pg_logit_rv_upper <- if (is.null(residual_variance_upperbound)) 1 else residual_variance_upperbound
    return(
      Run_Binary(
        X = X, y = y, Z = Z, family = family,
        L = L, max.iter = max.iter, min.iter = min.iter, max.eps = max.eps,
        susie.iter = susie.iter, verbose = verbose, n_threads = n_threads,
        coverage = coverage, weight_cutoff = weight_cutoff,
        scaled_prior_variance = scaled_prior_variance,
        estimate_residual_variance = estimate_residual_variance,
        residual_variance = residual_variance,
        residual_variance_lowerbound = pg_logit_rv_lower,
        residual_variance_upperbound = pg_logit_rv_upper,
        L.init = L.init,
        init_cor_method = init_cor_method,
        refit_noncs = refit_noncs,
        noncs_var = noncs_var,
        suff_block_size = suff_block_size,
        ...
      )
    )
  }

  # general GLM (e.g., Poisson, Gaussian with non-default link, etc.)
  return(
    Run_GLM(
      X = X, y = y, Z = Z, family = family,
      L = L, max.iter = max.iter, min.iter = min.iter, max.eps = max.eps,
      susie.iter = susie.iter, verbose = verbose, n_threads = n_threads,
      coverage = coverage, weight_cutoff = weight_cutoff,
      scaled_prior_variance = scaled_prior_variance,
      estimate_residual_variance = estimate_residual_variance,
      residual_variance = residual_variance,
      residual_variance_lowerbound = residual_variance_lowerbound,
      residual_variance_upperbound = rv_upper_default,
      L.init = L.init,
      init_cor_method = init_cor_method,
      refit_noncs = refit_noncs,
      noncs_var = noncs_var,
      suff_block_size = suff_block_size,
      ...
    )
  )
}
