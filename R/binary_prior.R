.binary_response <- function(y, family) {
  if (!inherits(family, "family") ||
      !identical(family$family, "binomial")) {
    return(NULL)
  }
  if (is.matrix(y)) {
    if (ncol(y) != 2L || any(!is.finite(y)) || any(y < 0) ||
        any(rowSums(y) != 1)) return(NULL)
    return(as.numeric(y[, 1L]))
  }
  y <- as.numeric(y)
  if (any(!is.finite(y)) || any(!y %in% c(0, 1))) return(NULL)
  y
}

.binary_response_link <- function(y, family) {
  if (is.null(.binary_response(y, family))) return(NULL)
  family$link
}

.binary_prior_weights <- function(prior_weights, p) {
  if (is.null(prior_weights)) return(rep(1 / p, p))
  if (!is.numeric(prior_weights) || length(prior_weights) != p ||
      any(!is.finite(prior_weights)) || any(prior_weights < 0) ||
      sum(prior_weights) <= 0) {
    stop("prior_weights must be non-negative, finite, and have length ncol(X).")
  }
  prior_weights / sum(prior_weights)
}

.validate_estimate_prior_variance <- function(x) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("estimate_prior_variance must be TRUE or FALSE.")
  }
  x
}

.binary_loglik_derivatives <- function(b, x, y, eta, link) {
  t <- eta + x * b
  if (identical(link, "logit")) {
    p <- stats::plogis(t)
    q <- p * (1 - p)
    softplus <- pmax(t, 0) + log1p(exp(-abs(t)))
    return(c(
      l0 = sum(y * t - softplus),
      l1 = sum(x * (y - p)),
      l2 = -sum(x^2 * q),
      l3 = -sum(x^3 * q * (1 - 2 * p)),
      l4 = -sum(x^4 * q * (1 - 6 * q))
    ))
  }

  s <- 2 * y - 1
  z <- s * t
  log_Phi <- stats::pnorm(z, log.p = TRUE)
  r <- exp(stats::dnorm(z, log = TRUE) - log_Phi)
  r1 <- -r * (z + r)
  A <- (z + r) * (z + 2 * r) - 1
  g3 <- r * A
  g4 <- r * (-(z + r) * A + 2 * z + 3 * r -
               r * (z + r) * (3 * z + 4 * r))
  tail <- z < -10
  if (any(tail)) {
    a <- -z[tail]
    r1[tail] <- -1 + a^-2 - 6 * a^-4 + 50 * a^-6 -
      518 * a^-8 + 6354 * a^-10
    g3[tail] <- 2 * a^-3 - 24 * a^-5 + 300 * a^-7 -
      4144 * a^-9 + 63540 * a^-11
    g4[tail] <- 6 * a^-4 - 120 * a^-6 + 2100 * a^-8 -
      37296 * a^-10 + 698940 * a^-12
  }
  sx <- s * x
  c(
    l0 = sum(log_Phi),
    l1 = sum(sx * r),
    l2 = sum(sx^2 * r1),
    l3 = sum(sx^3 * g3),
    l4 = sum(sx^4 * g4)
  )
}

.binary_laplace4_lbf <- function(V, x, y, eta, link, null_loglik = NULL) {
  if (!is.finite(V) || V < 0) return(-Inf)
  if (V == 0) return(0)
  if (is.null(null_loglik)) {
    null_loglik <- .binary_loglik_derivatives(0, x, y, eta, link)["l0"]
  }

  b <- 0
  converged <- FALSE
  for (i in seq_len(50L)) {
    d <- .binary_loglik_derivatives(b, x, y, eta, link)
    H <- -d["l2"] + 1 / V
    score <- d["l1"] - b / V
    if (!is.finite(H) || H <= 0 || !is.finite(score)) return(-Inf)
    step <- score / H
    if (!is.finite(step)) return(-Inf)
    b_new <- b + step
    h_now <- d["l0"] - b^2 / (2 * V)
    for (k in seq_len(25L)) {
      h_new <- .binary_loglik_derivatives(b_new, x, y, eta, link)["l0"] -
        b_new^2 / (2 * V)
      if (is.finite(h_new) && h_new >= h_now) break
      step <- step / 2
      b_new <- b + step
    }
    if (!is.finite(h_new) || h_new < h_now) return(-Inf)
    if (abs(step) <= 1e-10 * (1 + abs(b_new))) {
      b <- b_new
      converged <- TRUE
      break
    }
    b <- b_new
  }
  if (!converged) return(-Inf)

  d <- .binary_loglik_derivatives(b, x, y, eta, link)
  H <- -d["l2"] + 1 / V
  if (!is.finite(H) || H <= 0 || any(!is.finite(d))) return(-Inf)
  correction <- d["l4"] / (8 * H^2) + 5 * d["l3"]^2 / (24 * H^3)
  if (!is.finite(correction) || correction <= -1) return(-Inf)
  h <- d["l0"] - b^2 / (2 * V) - 0.5 * log(2 * pi * V)
  as.numeric(h + 0.5 * log(2 * pi / H) + log1p(correction) - null_loglik)
}

.binary_laplace4_objective <- function(V, X, y, eta, link,
                                       prior_weights = NULL) {
  p <- ncol(X)
  pi_j <- .binary_prior_weights(prior_weights, p)
  null_loglik <- .binary_loglik_derivatives(
    0, X[, 1L], y, eta, link
  )["l0"]
  lbf <- vapply(seq_len(p), function(j) {
    .binary_laplace4_lbf(
      V, X[, j], y, eta, link, null_loglik = null_loglik
    )
  }, numeric(1))
  keep <- pi_j > 0 & is.finite(lbf)
  if (!any(keep)) return(-Inf)
  z <- log(pi_j[keep]) + lbf[keep]
  zmax <- max(z)
  zmax + log(sum(exp(z - zmax)))
}

.estimate_binary_prior_variance <- function(X, y, eta, link,
                                            initial_variance = 1,
                                            prior_weights = NULL) {
  if (!link %in% c("logit", "probit")) return(initial_variance)
  if (!is.numeric(initial_variance) || length(initial_variance) != 1L ||
      !is.finite(initial_variance) || initial_variance < 0) {
    stop("scaled_prior_variance must be a non-negative finite scalar.")
  }
  if (initial_variance == 0) return(0)

  log_bounds <- c(-30, 15)
  objective <- function(log_V) {
    value <- .binary_laplace4_objective(
      exp(log_V), X, y, eta, link, prior_weights
    )
    if (is.finite(value)) value else -.Machine$double.xmax
  }
  opt <- stats::optim(
    par = min(max(log(initial_variance), log_bounds[1]), log_bounds[2]),
    fn = function(log_V) -objective(log_V), method = "Brent",
    lower = log_bounds[1], upper = log_bounds[2]
  )
  opt$objective <- objective(opt$par)
  initial_value <- objective(log(initial_variance))
  if (initial_value > opt$objective) {
    opt$par <- log(initial_variance)
    opt$objective <- initial_value
  }
  if (!is.finite(opt$objective) || opt$objective <= 0) return(0)
  exp(opt$par)
}

.binary_prior_for_fit <- function(X, y, eta, family, estimate_prior_variance,
                                  scaled_prior_variance,
                                  prior_weights = NULL) {
  link <- .binary_response_link(y, family)
  if (is.null(link) || !estimate_prior_variance ||
      !link %in% c("logit", "probit")) {
    return(scaled_prior_variance)
  }
  .estimate_binary_prior_variance(
    X = X, y = .binary_response(y, family), eta = as.numeric(eta), link = link,
    initial_variance = scaled_prior_variance,
    prior_weights = prior_weights
  )
}
