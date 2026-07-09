require(SuSiEIRLS)
require(mgcv)

set.seed(103)
n <- 300L
p <- 20L
q <- 3L
true_idx <- c(6L, 15L)

X <- matrix(stats::rnorm(n * p), n, p)
Z <- matrix(stats::rnorm(n * q), n, q)
X <- as.matrix(scale(X))
Z <- as.matrix(scale(Z))
colnames(X) <- paste0("X", seq_len(p))
colnames(Z) <- paste0("Z", seq_len(q))

beta <- rep(0, p)
beta[true_idx] <- c(0.30, -0.30)
alpha <- c(0.25, -0.15, 0.10)
eta <- as.numeric(X %*% beta + Z %*% alpha)
a <- stats::uniroot(
  function(a0) mean(exp(a0 + eta)) - 2.5,
  interval = c(-20, 20)
)$root
mu <- exp(a + eta)
y <- stats::rnbinom(n, size = 2.5, mu = mu)

fit <- SuSiEIRLS::SuSiE_IRLS(
  X = X, Z = Z, y = y,
  family = mgcv::nb(theta = NULL, link = "log"),
  L = 2L, L.init = 1L,
  max.iter = 3L, min.iter = 1L, max.eps = 1e-3,
  susie.iter = 80L,
  n_threads = 1L,
  verbose = FALSE,
  scale_data = FALSE
)

pip <- fit$fitX$pip[seq_len(p)]
if (length(pip) != p || any(!is.finite(pip))) stop("NB smoke test failed")
print(data.frame(family = "nb", iter = fit$iter,
                 n_cs = nrow(fit$main_index),
                 top_signal_pip = max(pip[true_idx]),
                 theta_hat = fit$theta))
