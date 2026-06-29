## Shared evaluation helpers for SuSiEIRLS benchmarks.
##
## Standard benchmark rule:
##   - CS-based power/FDR
##   - coverage = 0.95 unless overridden
##   - min_abs_corr = 0.5, matching susieR::susie_get_cs() default purity
##
## Older exploratory scripts sometimes used min_abs_cor = 0.1. That is not the
## standard benchmark rule used here.

benchmark_quiet_eval <- function(expr) {
  capture.output(value <- eval.parent(substitute(expr)))
  value
}

benchmark_empty_main_index <- function() {
  data.frame(
    Index = integer(0),
    CS = character(0),
    PIP = numeric(0),
    lbf = numeric(0),
    log10BF = numeric(0),
    Pvalue = numeric(0)
  )
}

benchmark_normalize_main_index <- function(main_index, p) {
  if (is.null(main_index) || !nrow(main_index)) {
    return(benchmark_empty_main_index())
  }

  idx <- suppressWarnings(as.integer(main_index$Index))
  keep <- is.finite(idx) & idx >= 1L & idx <= p
  if (!any(keep)) return(benchmark_empty_main_index())

  out <- main_index[keep, , drop = FALSE]
  out$Index <- as.integer(out$Index)
  if (!("CS" %in% names(out))) out$CS <- "Main_CS1"
  if (!("PIP" %in% names(out))) out$PIP <- NA_real_
  if (!("lbf" %in% names(out))) out$lbf <- NA_real_
  if (!("log10BF" %in% names(out))) out$log10BF <- out$lbf / log(10)
  if (!("Pvalue" %in% names(out))) out$Pvalue <- NA_real_
  row.names(out) <- NULL
  out
}

benchmark_cs1_lbf <- function(main_index, p) {
  main_index <- benchmark_normalize_main_index(main_index, p)
  if (!nrow(main_index)) {
    return(c(lbf_cs1 = NA_real_, log10bf_cs1 = NA_real_))
  }

  cs_names <- unique(as.character(main_index$CS))
  cs_num <- suppressWarnings(as.integer(gsub("[^0-9]", "", cs_names)))
  first_cs <- cs_names[
    order(ifelse(is.finite(cs_num), cs_num, seq_along(cs_names)))[1L]
  ]
  ii <- which(main_index$CS == first_cs)[1L]

  lbf <- suppressWarnings(as.numeric(main_index$lbf[ii]))
  log10bf <- suppressWarnings(as.numeric(main_index$log10BF[ii]))
  if (!is.finite(log10bf) && is.finite(lbf)) log10bf <- lbf / log(10)
  if (!is.finite(lbf) && is.finite(log10bf)) lbf <- log10bf * log(10)
  c(lbf_cs1 = lbf, log10bf_cs1 = log10bf)
}

benchmark_eval_cs <- function(main_index, true_idx, p) {
  main_index <- benchmark_normalize_main_index(main_index, p)
  bf <- benchmark_cs1_lbf(main_index, p)

  if (!nrow(main_index)) {
    return(c(
      power = 0,
      fdr = 0,
      n_cs = 0,
      mean_cs_size = NA_real_,
      lbf_cs1 = bf[["lbf_cs1"]],
      log10bf_cs1 = bf[["log10bf_cs1"]]
    ))
  }

  cs_split <- split(main_index$Index, main_index$CS)
  detected <- sort(unique(main_index$Index))
  cs_has_signal <- vapply(
    cs_split,
    function(ii) any(ii %in% true_idx),
    logical(1)
  )

  c(
    power = mean(true_idx %in% detected),
    fdr = mean(!cs_has_signal),
    n_cs = length(cs_split),
    mean_cs_size = mean(vapply(cs_split, length, integer(1))),
    lbf_cs1 = bf[["lbf_cs1"]],
    log10bf_cs1 = bf[["log10bf_cs1"]]
  )
}

benchmark_ibss_x_main_index <- function(fit, X_aug, p, coverage = 0.95,
                                        min_abs_corr = 0.5) {
  if (is.null(fit)) return(benchmark_empty_main_index())

  fit_susie <- fit
  class(fit_susie) <- c("susie", "list")
  cs_out <- tryCatch(
    susieR::susie_get_cs(
      fit_susie, X = X_aug, coverage = coverage,
      min_abs_corr = min_abs_corr
    ),
    error = function(e) NULL
  )
  if (is.null(cs_out) || !length(cs_out$cs)) {
    return(benchmark_empty_main_index())
  }

  effect_index <- cs_out$cs_index
  if (is.null(effect_index)) {
    effect_index <- suppressWarnings(
      as.integer(gsub("[^0-9]", "", names(cs_out$cs)))
    )
  }
  if (length(effect_index) != length(cs_out$cs) ||
      any(!is.finite(effect_index))) {
    effect_index <- seq_along(cs_out$cs)
  }

  cs_list <- lapply(cs_out$cs, function(ii) ii[ii <= p])
  keep <- vapply(cs_list, length, integer(1)) > 0
  cs_list <- cs_list[keep]
  effect_index <- effect_index[keep]
  if (!length(cs_list)) return(benchmark_empty_main_index())

  idx <- unlist(cs_list)
  cs_id <- rep(seq_along(cs_list), times = vapply(cs_list, length, integer(1)))

  lbf_src <- if (!is.null(fit$lbf_ser)) fit$lbf_ser else fit$lbf
  lbf <- rep(NA_real_, length(cs_list))
  if (!is.null(lbf_src) && length(lbf_src) >= max(effect_index)) {
    lbf <- suppressWarnings(as.numeric(lbf_src[effect_index]))
  }

  out <- data.frame(
    Index = idx,
    CS = paste0("Main_CS", cs_id),
    stringsAsFactors = FALSE
  )
  out$PIP <- fit$pip[out$Index]
  out$lbf <- lbf[cs_id]
  out$log10BF <- out$lbf / log(10)
  out$Pvalue <- NA_real_
  out
}

benchmark_mean_finite <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ok <- is.finite(x)
  if (!any(ok)) return(NA_real_)
  mean(x[ok])
}

benchmark_summary <- function(per_run) {
  rows <- lapply(
    split(per_run, list(per_run$n, per_run$family, per_run$method), drop = TRUE),
    function(x) {
      data.frame(
        n = x$n[1],
        family = x$family[1],
        method = x$method[1],
        n_run = nrow(x),
        n_failed = sum(!is.na(x$error)),
        mean_power = benchmark_mean_finite(x$power),
        mean_fdr = benchmark_mean_finite(x$fdr),
        mean_n_cs = benchmark_mean_finite(x$n_cs),
        mean_cs_size = benchmark_mean_finite(x$mean_cs_size),
        mean_lbf_cs1 = benchmark_mean_finite(x$lbf_cs1),
        mean_log10bf_cs1 = benchmark_mean_finite(x$log10bf_cs1),
        mean_time_sec = benchmark_mean_finite(x$time_sec),
        var_eta_z = benchmark_mean_finite(x$var_eta_z),
        var_eta_x = benchmark_mean_finite(x$var_eta_x),
        var_eta_total = benchmark_mean_finite(x$var_eta_total),
        stringsAsFactors = FALSE
      )
    }
  )
  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out[order(out$n, out$family, out$method), , drop = FALSE]
}
