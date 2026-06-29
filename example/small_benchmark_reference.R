## Small reference benchmark:
##   binary(PG), negative binomial, Cox, Tweedie, and beta outcomes.
##
## The comparison is CS-based and uses the standard benchmark helper in
## example/evaluation.R:
##   coverage = 0.95
##   min_abs_corr = 0.5
##
## By default this runs n = 500 and 50 replicates per family. Override with:
##   BENCH_NS=300,500
##   BENCH_NREP=10
##   BENCH_FAMILIES=binary_pg,nb,cox,tw,beta
##   BENCH_METHODS=irls,ibss
##   BENCH_OUT_PREFIX=tmp/my_benchmark

suppressPackageStartupMessages({
  library(MASS)
  library(mgcv)
  library(survival)
  library(susieR)
  library(logisticsusie)
  library(SuSiEIRLS)
})

source("example/evaluation.R")
source("example/reference_simulation.R")

split_env <- function(name, default) {
  val <- Sys.getenv(name, "")
  if (!nzchar(val)) return(default)
  strsplit(val, ",", fixed = TRUE)[[1]]
}

cfg <- reference_benchmark_config(
  ns = as.integer(split_env("BENCH_NS", "500")),
  n_rep = as.integer(Sys.getenv("BENCH_NREP", "50")),
  families = split_env("BENCH_FAMILIES",
                       c("binary_pg", "nb", "cox", "tw", "beta")),
  coverage = as.numeric(Sys.getenv("BENCH_COVERAGE", "0.95")),
  min_abs_corr = as.numeric(Sys.getenv("BENCH_MIN_ABS_CORR", "0.5")),
  var_eta_z = as.numeric(Sys.getenv("BENCH_VAR_ETA_Z", "0.2")),
  var_eta_x = as.numeric(Sys.getenv("BENCH_VAR_ETA_X", "0.05"))
)

methods <- split_env("BENCH_METHODS", c("irls", "ibss"))
out_prefix <- Sys.getenv(
  "BENCH_OUT_PREFIX",
  "tmp/reference_small_benchmark_minabs05"
)

bt <- function(x) paste0("`", gsub("`", "``", x, fixed = TRUE), "`")

explicit_formula <- function(response, rhs, offset_name = NULL) {
  terms <- c(bt(rhs), if (!is.null(offset_name)) paste0("offset(", offset_name, ")"))
  if (!length(terms)) return(stats::as.formula(paste(response, "~ 1")))
  stats::as.formula(paste(response, "~", paste(terms, collapse = " + ")))
}

theta_of <- function(fit) {
  if (!is.function(fit$family$getTheta)) return(NULL)
  theta <- tryCatch(fit$family$getTheta(TRUE), error = function(e) NULL)
  if (is.null(theta) || any(!is.finite(theta))) return(NULL)
  as.numeric(theta)
}

zfit_family <- function(family) {
  if (identical(family, "binary_pg")) return(stats::binomial())
  if (identical(family, "nb")) return(mgcv::nb(theta = NULL))
  if (identical(family, "tw")) return(mgcv::tw(theta = NULL))
  if (identical(family, "beta")) return(mgcv::betar(link = "logit", theta = NULL))
  stop("Unsupported Z-fit family: ", family)
}

fixed_family_from_theta <- function(family, theta) {
  if (identical(family, "binary_pg")) return(stats::binomial())
  if (is.null(theta)) return(zfit_family(family))
  if (identical(family, "nb")) return(do.call(mgcv::nb, list(theta = theta)))
  if (identical(family, "tw")) return(do.call(mgcv::tw, list(theta = theta)))
  if (identical(family, "beta")) {
    return(do.call(mgcv::betar, list(link = "logit", theta = theta)))
  }
  stop("Unsupported fixed family: ", family)
}

irls_family <- function(family) {
  if (identical(family, "binary_pg")) return(stats::binomial(link = "logit"))
  if (identical(family, "nb")) return("negbin")
  if (identical(family, "cox")) return(stats::binomial(link = "logit"))
  if (identical(family, "tw")) return(mgcv::tw())
  if (identical(family, "beta")) return(mgcv::betar(link = "logit"))
  stop("Unsupported IRLS family: ", family)
}

response_frame <- function(dat_obj, family) {
  if (identical(family, "cox")) {
    return(list(
      data = data.frame(time = dat_obj$time, status = dat_obj$status),
      response = "survival::Surv(time, status)"
    ))
  }
  list(data = data.frame(y = dat_obj$y), response = "y")
}

fit_reference_z_model <- function(dat_obj, family) {
  Zdf <- as.data.frame(dat_obj$Z)
  colnames(Zdf) <- colnames(dat_obj$Z)

  if (identical(family, "cox")) {
    dat <- data.frame(time = dat_obj$time, status = dat_obj$status, Zdf)
    fit <- survival::coxph(
      explicit_formula("survival::Surv(time, status)", colnames(Zdf)),
      data = dat, ties = "breslow"
    )
    return(list(
      eta = as.numeric(stats::predict(fit, type = "lp")),
      theta = NULL,
      family = NULL,
      fit = fit
    ))
  }

  rsp <- response_frame(dat_obj, family)
  dat <- cbind(rsp$data, Zdf)
  fit <- mgcv::bam(
    explicit_formula(rsp$response, colnames(Zdf)),
    data = dat, family = zfit_family(family), method = "fREML"
  )
  theta <- theta_of(fit)
  list(
    eta = as.numeric(fit$linear.predictors),
    theta = theta,
    family = fixed_family_from_theta(family, theta),
    fit = fit
  )
}

ser_empty <- function(v0) {
  list(mu = 0, var = v0, lbf = -1e30, prior_variance = v0, intercept = 0)
}

make_mgcv_ser_fun <- function(family, eta_z, family_object) {
  fam <- family_object
  etaZ_hat <- as.numeric(eta_z)

  function(x, y, e, prior_variance, estimate_intercept = 0, ...) {
    v0 <- prior_variance
    off <- as.numeric(e)
    dat <- data.frame(y = y, x = x, etaZ_hat = etaZ_hat, off = off)

    fit <- tryCatch(
      mgcv::bam(y ~ x + etaZ_hat + offset(off),
                data = dat, family = fam, method = "fREML"),
      error = function(e) e
    )
    if (inherits(fit, "error")) return(ser_empty(v0))

    co <- tryCatch(summary(fit)$p.table, error = function(e) NULL)
    if (is.null(co) || !("x" %in% rownames(co))) return(ser_empty(v0))
    bhat <- co["x", "Estimate"]
    s <- co["x", "Std. Error"]
    if (!is.finite(bhat) || !is.finite(s) || s <= 0) return(ser_empty(v0))

    z <- bhat / s
    lbf_wake <- 0.5 * log(s^2 / (v0 + s^2)) +
      0.5 * z^2 * v0 / (v0 + s^2)

    fit0 <- tryCatch(
      mgcv::bam(y ~ etaZ_hat + offset(off),
                data = dat, family = fam, method = "fREML"),
      error = function(e) NULL
    )
    lrt <- if (is.null(fit0)) 0 else as.numeric(2 * (logLik(fit) - logLik(fit0)))
    lbf <- lbf_wake - 0.5 * z^2 + 0.5 * lrt

    v1 <- 1 / (1 / v0 + 1 / s^2)
    mu1 <- v1 * bhat / s^2
    list(mu = mu1, var = v1, lbf = lbf,
         prior_variance = mu1^2 + v1, intercept = 0)
  }
}

make_cox_ser_fun <- function(eta_z) {
  etaZ_hat <- as.numeric(eta_z)

  function(x, y, e, prior_variance, estimate_intercept = 0, ...) {
    v0 <- prior_variance
    off <- as.numeric(e)
    dat <- data.frame(
      time = as.numeric(y[, 1]),
      status = as.integer(y[, 2]),
      x = x,
      etaZ_hat = etaZ_hat,
      off = off
    )

    fit <- tryCatch(
      survival::coxph(
        survival::Surv(time, status) ~ x + etaZ_hat + offset(off),
        data = dat, ties = "breslow"
      ),
      error = function(e) e
    )
    if (inherits(fit, "error")) return(ser_empty(v0))

    co <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    if (is.null(co) || !("x" %in% rownames(co))) return(ser_empty(v0))
    bhat <- co["x", "coef"]
    s <- co["x", "se(coef)"]
    if (!is.finite(bhat) || !is.finite(s) || s <= 0) return(ser_empty(v0))

    z <- bhat / s
    lbf_wake <- 0.5 * log(s^2 / (v0 + s^2)) +
      0.5 * z^2 * v0 / (v0 + s^2)

    fit0 <- tryCatch(
      survival::coxph(
        survival::Surv(time, status) ~ etaZ_hat + offset(off),
        data = dat, ties = "breslow"
      ),
      error = function(e) NULL
    )
    lrt <- if (is.null(fit0)) 0 else as.numeric(2 * (logLik(fit) - logLik(fit0)))
    lbf <- lbf_wake - 0.5 * z^2 + 0.5 * lrt

    v1 <- 1 / (1 / v0 + 1 / s^2)
    mu1 <- v1 * bhat / s^2
    list(mu = mu1, var = v1, lbf = lbf,
         prior_variance = mu1^2 + v1, intercept = 0)
  }
}

fit_irls_reference <- function(dat_obj, family, cfg) {
  fit <- benchmark_quiet_eval(SuSiEIRLS::SuSiE_IRLS(
    X = dat_obj$X, Z = dat_obj$Z, y = dat_obj$y,
    family = irls_family(family),
    logit_method = if (identical(family, "binary_pg")) "pg" else "glm",
    L = cfg$L, L.init = 1L,
    max.iter = 5L, min.iter = 2L, max.eps = 1e-4,
    susie.iter = 100L,
    coverage = cfg$coverage,
    min_abs_corr = cfg$min_abs_corr,
    n_threads = 1L,
    verbose = FALSE,
    estimate_residual_variance = TRUE,
    residual_variance = 0.5,
    residual_variance_lowerbound = 0.1,
    residual_variance_upperbound = 1
  ))

  benchmark_eval_cs(fit$main_index, dat_obj$true_idx, ncol(dat_obj$X))
}

fit_ibss_reference <- function(dat_obj, family, cfg) {
  zfit <- fit_reference_z_model(dat_obj, family)
  uni_fun <- if (identical(family, "cox")) {
    make_cox_ser_fun(zfit$eta)
  } else {
    make_mgcv_ser_fun(family, zfit$eta, zfit$family)
  }

  fit <- benchmark_quiet_eval(logisticsusie::ibss_from_ser(
    X = dat_obj$X, y = dat_obj$y, L = cfg$L,
    tol = 1e-4, maxit = 100L, num_cores = 1L,
    ser_function = logisticsusie::ser_from_univariate(uni_fun)
  ))

  main_index <- benchmark_ibss_x_main_index(
    fit, X_aug = dat_obj$X, p = ncol(dat_obj$X),
    coverage = cfg$coverage, min_abs_corr = cfg$min_abs_corr
  )
  benchmark_eval_cs(main_index, dat_obj$true_idx, ncol(dat_obj$X))
}

run_one_method <- function(method, dat_obj, family, cfg) {
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch(
    if (identical(method, "irls")) {
      fit_irls_reference(dat_obj, family, cfg)
    } else {
      fit_ibss_reference(dat_obj, family, cfg)
    },
    error = function(e) c(
      power = NA_real_, fdr = NA_real_, n_cs = NA_real_,
      mean_cs_size = NA_real_, lbf_cs1 = NA_real_,
      log10bf_cs1 = NA_real_, error = conditionMessage(e)
    )
  )
  elapsed <- proc.time()[["elapsed"]] - t0

  data.frame(
    method = method,
    power = as.numeric(out[["power"]]),
    fdr = as.numeric(out[["fdr"]]),
    n_cs = as.numeric(out[["n_cs"]]),
    mean_cs_size = as.numeric(out[["mean_cs_size"]]),
    lbf_cs1 = as.numeric(out[["lbf_cs1"]]),
    log10bf_cs1 = as.numeric(out[["log10bf_cs1"]]),
    time_sec = elapsed,
    error = if ("error" %in% names(out)) as.character(out[["error"]]) else NA_character_,
    stringsAsFactors = FALSE
  )
}

rows <- list()
k <- 0L
for (n in cfg$ns) {
  for (family in cfg$families) {
    for (rep in seq_len(cfg$n_rep)) {
      dat_obj <- simulate_reference_dataset(
        n = n, family = family,
        seed = 100000 + n * 100 + rep + match(family, cfg$families) * 10000,
        cfg = cfg
      )

      for (method in methods) {
        k <- k + 1L
        res <- run_one_method(method, dat_obj, family, cfg)
        res$n <- n
        res$family <- family
        res$rep <- rep
        res$coverage <- cfg$coverage
        res$min_abs_corr <- cfg$min_abs_corr
        res$var_eta_z <- dat_obj$var_eta_z
        res$var_eta_x <- dat_obj$var_eta_x
        res$var_eta_total <- dat_obj$var_eta_total
        res$censor_rate <- dat_obj$censor_rate
        rows[[k]] <- res
      }

      if (rep %% 10L == 0L) {
        cat("done n", n, "family", family, "rep", rep, "\n")
      }
    }
  }
}

per_run <- do.call(rbind, rows)
summary_df <- benchmark_summary(per_run)

dir.create(dirname(out_prefix), showWarnings = FALSE, recursive = TRUE)
per_run_file <- paste0(out_prefix, "_per_run.csv")
summary_file <- paste0(out_prefix, "_summary.csv")
write.csv(per_run, per_run_file, row.names = FALSE)
write.csv(summary_df, summary_file, row.names = FALSE)

cat("\nWrote:\n")
cat("  ", per_run_file, "\n", sep = "")
cat("  ", summary_file, "\n\n", sep = "")
print(summary_df)
