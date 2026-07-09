require(SuSiEIRLS)

set.seed(101)
n <- 250L
p <- 20L
q <- 3L
true_idx <- c(4L, 11L)

X <- matrix(stats::rnorm(n * p), n, p)
Z <- matrix(stats::rnorm(n * q), n, q)
X <- as.matrix(scale(X))
Z <- as.matrix(scale(Z))
colnames(X) <- paste0("X", seq_len(p))
colnames(Z) <- paste0("Z", seq_len(q))

beta <- rep(0, p)
beta[true_idx] <- c(0.45, -0.45)
alpha <- c(0.35, -0.25, 0.15)
eta <- as.numeric(X %*% beta + Z %*% alpha)
a <- stats::uniroot(
  function(a0) mean(stats::plogis(a0 + eta)) - 0.35,
  interval = c(-20, 20)
)$root
y <- stats::rbinom(n, size = 1L, prob = stats::plogis(a + eta))

fit <- SuSiEIRLS::SuSiE_IRLS(
  X = X, Z = Z, y = y,
  family = stats::binomial(link = "logit"),
  logit_method = "pg",
  L = 2L, L.init = 1L,
  max.iter = 3L, min.iter = 1L, max.eps = 1e-3,
  susie.iter = 80L,
  n_threads = 1L,
  verbose = FALSE,
  scale_data = FALSE
)

pip <- fit$fitX$pip[seq_len(p)]
if (length(pip) != p || any(!is.finite(pip))) stop("binary smoke test failed")
print(data.frame(family = "binary", iter = fit$iter,
                 n_cs = nrow(fit$main_index),
                 top_signal_pip = max(pip[true_idx])))
