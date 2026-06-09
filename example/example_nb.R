library(CppMatrix)
library(susieR)
library(devtools)
library(MASS)
library(ncvreg)
library(glmnet)
library(logisticsusie)
options(bitmapType = "cairo")
load_all("~/SuSiEGLM/")

# ----------------------------
# 工具函数：协方差、索引映射
# ----------------------------
ARcov <- function(p, rho){
  s <- rho ^ (0:(p-1))
  toeplitz(s)
}
CScov <- function(p, rho){
  diag(p)*(1-rho) + matrix(rho, p, p)
}

selected_from_main_index <- function(main_index, p, x_prefix = "SNP"){
  if (is.null(main_index) || nrow(main_index) == 0) return(integer(0))
  if ("Index" %in% names(main_index)) {
    idx <- as.integer(main_index$Index)
    idx[!is.na(idx) & idx >= 1 & idx <= p]
  } else if ("Variable" %in% names(main_index)) {
    v <- as.character(main_index$Variable)
    as.integer(gsub(paste0("^", x_prefix), "", v))
  } else integer(0)
}

susie_to_main_index <- function(fit_susie, X, p, coverage = 0.95, min_abs_cor = 0.1) {
  class(fit_susie) <- c("susie", "list")
  cs_out <- susie_get_cs(fit_susie, X = X, coverage = coverage, min_abs_cor = min_abs_cor)
  cs_list <- lapply(cs_out$cs, function(idx) idx[idx <= p])
  cs_list <- cs_list[sapply(cs_list, length) > 0]
  if (length(cs_list) == 0) return(data.frame())
  main_index <- data.frame(
    Index = unlist(cs_list),
    Variable = colnames(X)[unlist(cs_list)],
    CS = rep(paste0("Main_CS", 1:length(cs_list)),
             times = sapply(cs_list, length)),
    stringsAsFactors = FALSE
  )
  main_index$PIP <- fit_susie$pip[main_index$Index]
  main_index <- main_index[order(main_index$CS, main_index$Index), ]
  row.names(main_index) <- NULL
  return(main_index)
}

# ----------------------------
# 单变量 SER（Poisson）给 IBSS 用
# Wakefield 近似 + LRT 修正，带 offset(e)
# ----------------------------
poisson_uni_fun <- function(x, y, e, prior_variance,
                            estimate_intercept = 0, ...) {
  v0  <- prior_variance
  fit <- MASS::glm.nb(y ~ x + offset(e))
  co  <- summary(fit)$coefficients
  bhat <- if ("x" %in% rownames(co)) co["x","Estimate"] else co[2, "Estimate"]
  s    <- if ("x" %in% rownames(co)) co["x","Std. Error"] else co[2, "Std. Error"]
  z    <- bhat / s

  # Wakefield approximate BF + LRT tweak（与你 logistic 版本对应）
  lbf_wake <- 0.5 * log(s^2 / (v0 + s^2)) + 0.5 * z^2 * v0 / (v0 + s^2)
  fit0 <- MASS::glm.nb(y ~ 1 + offset(e))
  lrt  <- as.numeric(2 * (logLik(fit) - logLik(fit0)))
  lbf  <- lbf_wake - 0.5 * z^2 + 0.5 * lrt

  v1 <- 1 / (1 / v0 + 1 / s^2)
  mu1 <- v1 * bhat / s^2
  list(mu = mu1, var = v1, lbf = lbf,
       prior_variance = mu1^2 + v1,
       intercept = 0)
}

# ----------------------------
# 你的评估/其它函数
# ----------------------------
source("~/SuSiEGLM/example/main_tptn.R")
source("~/SuSiEGLM/example/otherfunction.R")

n  <- 500
p  <- 10
nb_size <- 5  # 负二项 'size' (θ); 你可调整过度离散度

Main_TP <- Main_TN <- TIME <- matrix(0, 200, 4)

for(iter in 1:100){
  # --------------------------------------
  # 生成 X, Z, beta, alpha, eta
  # --------------------------------------
  R  <- ARcov(p, 0.5)
  X  <- mvrnorm(n = n, mu = runif(p, 0, 1), R)
  X  <- scale(X, center = TRUE, scale = FALSE) # （可选）去均值，避免截距偏移
  Z  <- mvrnorm(n = n, mu = runif(5, 0, 1), diag(5))
  colnames(Z) <- paste0("UKBB", 1:5)

  alpha0 <- rnorm(5, 0, 1/sqrt(5))
  beta0  <- rep(0, p)
  ind    <- sort(sample(p, 3))
  beta0[ind] <- 0.5

  eta <- matrixVectorMultiply(Z, alpha0) + matrixVectorMultiply(X, beta0)
  mu  <- exp(eta/2+1)

  # --------------------------------------
  # 负二项模拟（NB2）：y ~ NB(size=nb_size, mu=mu)
  # Var(Y)=mu + mu^2/nb_size
  # 你可以自行“调 eta 的 scale”来控制均值/方差
  # --------------------------------------
  y <- rnbinom(n = n, size = nb_size, mu = mu)

  # --------------------------------------
  # 1) SuSiE_IRLS（若支持 NB）; 否则切到 Poisson 作参考
  # --------------------------------------
  t1 <- Sys.time()
  fit_SuSiEGLM08 <- tryCatch({
    SuSiE_IRLS(
      X = X, Z = Z, y = y, L = 5, n_threads = 4,estimate_residual_variance = F,
      max.iter = 15, verbose = TRUE, susie.iter = 500,
      family = "negbin"
    )
  }, error = function(e) {
    message("❌ SuSiEGLM08 failed at iteration ", iter, ": ", e$message)
    NULL
  })
  t2 <- Sys.time()
  t_SuSiEGLM08 <- difftime(t2, t1, units = "secs")

  # --------------------------------------
  # 2) glmnet：只能用 Poisson（glmnet 无 quasipoisson/NB）
  #    注意：我们把 Z 并到设计矩阵，不用 offset
  # --------------------------------------
  t1 <- Sys.time()
  fit_lasso_cv <- cv.glmnet(
    x = cbind(X, Z), y = y,
    family = quasipoisson()
  )
  t2 <- Sys.time()
  t_lasso_cv <- difftime(t2, t1, units = "secs")

  # --------------------------------------
  # 3) ncvreg：按你的要求都用 Poisson
  # --------------------------------------
  t1 <- Sys.time()
  fit_mcp_cv <- cv.ncvreg(X = cbind(X, Z), y = y, family = "poisson")
  t2 <- Sys.time()
  t_mcp_cv <- difftime(t2, t1, units = "secs")

  # --------------------------------------
  # 4) IBSS（官方 SER）：换成 Poisson 的单变量 SER
  # # --------------------------------------
  t1 <- Sys.time()
  fit_ibss <- ibss_from_ser(
    X = cbind(X, Z), y = y, L = 10,
    tol = 1e-4, maxit = 100, num_cores = 8,
    ser_function = ser_from_univariate(poisson_uni_fun)
  )
  out_ibss <- susie_to_main_index(
    fit_ibss, X = cbind(X, Z), p = p,
    coverage = 0.95, min_abs_cor = 0.1
  )
  t2 <- Sys.time()
  t_ibss <- as.numeric(difftime(t2, t1, units = "secs"))

  # --------------------------------------
  # 指标计算
  # --------------------------------------
  if (is.null(out_ibss) || nrow(out_ibss) == 0) {
    tptn_ibss <- list(tp = 0, tn = 1)
  } else {
    tptn_ibss <- main_tptn(true_main_index = ind, main_index = out_ibss)
  }
  if (is.null(fit_SuSiEGLM08)) {
    tptn_SuSiEGLM08 <- list(tp = 0, tn = 1)
  } else {
    tptn_SuSiEGLM08 <- main_tptn(true_main_index = ind, main_index = fit_SuSiEGLM08$main_index)
  }

  tptn_lasso_cv  <- tptn_evaulate(true_main_index = ind,
                                  hat_beta = fit_lasso_cv$glmnet.fit$beta[1:p, which(fit_lasso_cv$lambda == fit_lasso_cv$lambda.1se)])
  tptn_mcp_cv    <- tptn_evaulate(true_main_index = ind,
                                  hat_beta = fit_mcp_cv$fit$beta[2:(p+1), which.min(fit_mcp_cv$cve)])

  TIME[iter,]    <- c(t_SuSiEGLM08, t_lasso_cv, t_mcp_cv, t_ibss)
  Main_TP[iter,] <- c(tptn_SuSiEGLM08$tp, tptn_lasso_cv$tp, tptn_mcp_cv$tp, tptn_ibss$tp)
  Main_TN[iter,] <- c(tptn_SuSiEGLM08$tn, tptn_lasso_cv$tn, tptn_mcp_cv$tn, tptn_ibss$tn)

  if(iter %% 10 == 0) print(iter)
}
colMeans(Main_TP)
colMeans(Main_TN)
colMeans(TIME)
