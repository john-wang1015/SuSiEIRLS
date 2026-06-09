library(CppMatrix)
library(susieR)
library(devtools)
library(MASS)
library(ncvreg)
library(glmnet)
library(logisticsusie)
load_all("~/SuSiEGLM/")
ARcov=function(p,rho){
s=c(1:p)
for(i in 1:p){
s[i]=rho^(i-1)
}
return(toeplitz(s))
}
CScov=function(p,rho){
diag(p)*(1-rho)+matrix(rho,p,p)
}
logit_uni_fun <- function(x, y, e, prior_variance,
                          estimate_intercept = 0, ...) {
  v0  <- prior_variance
  fit <- glm(y ~ x + offset(e), family = binomial())
  co  <- summary(fit)$coefficients
  bhat <- if ("x" %in% rownames(co)) co["x", "Estimate"] else co[2, "Estimate"]
  s    <- if ("x" %in% rownames(co)) co["x", "Std. Error"] else co[2, "Std. Error"]
  z    <- bhat / s
  lbf_wake <- 0.5 * log(s^2 / (v0 + s^2)) + 0.5 * z^2 * v0 / (v0 + s^2)
  fit0 <- glm(y ~ 1 + offset(e), family = binomial())
  lrt  <- as.numeric(2 * (logLik(fit) - logLik(fit0)))
  lbf  <- lbf_wake - 0.5 * z^2 + 0.5 * lrt
  v1   <- 1 / (1 / v0 + 1 / s^2)
  mu1  <- v1 * bhat / s^2
  list(mu = mu1, var = v1, lbf = lbf,
       prior_variance = mu1^2 + v1,
       intercept = 0)
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

source("~/SuSiEGLM/example/main_tptn.R")
source("~/SuSiEGLM/example/otherfunction.R")
n=1000
p=10
Main_TP=Main_TN=TIME=matrix(0,100,6)

for(iter in 1:100){
###############################################################################
R=ARcov(10,0.5)
X=mvrnorm(n=n,mu=runif(p,0,1),R)
Xmean=colMeans(X)
Z=mvrnorm(n=n,mu=runif(5,0,1),diag(5))
colnames(Z)=paste0("UKBB",1:5)
alpha0=c(1:5)*0
alpha0[1:5]=rnorm(5,0,1/sqrt(5))
beta0=rep(0,p)
ind=sample(p,3)
ind=sort(ind)
beta0[ind]=0.5
eta=matrixVectorMultiply(Z,alpha0)+matrixVectorMultiply(X,beta0)
y=rbinom(n=n,size=1,prob=1/(1+exp(-eta)))

t1=Sys.time()
fit_SuSiEGLM08 <- tryCatch(
{
SuSiE_IRLS(
X = X, Z = Z, y = y, L = 5, n_threads = 4,
max.iter = 15, verbose = TRUE, susie.iter = 1000, family = binomial(link="logit")
)
},
error = function(e) {
message("❌ SuSiEGLM08 failed at iteration ", iter, ": ", e$message)
NULL
}
)
t2 <- Sys.time()
t_SuSiEGLM08 <- difftime(t2, t1, units = "secs")

t1=Sys.time()
fit_lasso_cv=cv.glmnet(x=cbind(X,Z),y=y,family=binomial(link="logit"))
t2=Sys.time()
t_lasso_cv=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit_mcp_cv=cv.ncvreg(X=cbind(X,Z),y=y,family="binomial")
t2=Sys.time()
t_mcp_cv=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit_lasso_bic=ncvreg(X=cbind(X,Z),y=y,family="binomial",penalty="lasso")
t2=Sys.time()
t_lasso_bic=difftime(t2, t1, units = "secs")

t1=Sys.time()
fit_mcp_bic=ncvreg(X=cbind(X,Z),y=y,family="binomial",penalty="MCP")
t2=Sys.time()
t_mcp_bic=difftime(t2, t1, units = "secs")

## 4) 官方 IBSS (logisticsusie) ------------------------------------
t1 <- Sys.time()
fit_ibss <- ibss_from_ser(
  X = cbind(X, Z), y = y, L = 10,
  tol = 1e-4, maxit = 100,num_cores = 8,
  ser_function = ser_from_univariate(logit_uni_fun)
)
out_ibss <- susie_to_main_index(
  fit_ibss, X = cbind(X, Z), p = p,
  coverage = 0.95, min_abs_cor = 0.1
)
t2 <- Sys.time()
t_ibss <- as.numeric(difftime(t2, t1, units = "secs"))

if (is.null(out_ibss) || nrow(out_ibss) == 0) {
  tptn_ibss <- list(tp = 0, tn = 1)
} else {
  tptn_ibss <- main_tptn(
    true_main_index = ind,
    main_index = out_ibss
  )
}
if (is.null(fit_SuSiEGLM08)) {
  tptn_SuSiEGLM08 <- list(tp = 0, tn = 1)
} else {
  tptn_SuSiEGLM08 <- main_tptn(
    true_main_index = ind,
    main_index = fit_SuSiEGLM08$main_index
  )
}

tptn_lasso_cv=tptn_evaulate(true_main_index=ind,hat_beta=fit_lasso_cv$glmnet.fit$beta[1:p,which(fit_lasso_cv$lambda==fit_lasso_cv$lambda.1se)])
tptn_lasso_bic=tptn_evaulate(true_main_index=ind,hat_beta=fit_lasso_bic$beta[2:(p+1),which.min(BIC(fit_lasso_bic))])
tptn_mcp_cv=tptn_evaulate(true_main_index=ind,hat_beta=fit_mcp_cv$fit$beta[2:(p+1),which.min(fit_mcp_cv$cve)])
tptn_mcp_bic=tptn_evaulate(true_main_index=ind,hat_beta=fit_mcp_bic$beta[2:(p+1),which.min(BIC(fit_mcp_bic))])

TIME[iter,]=c(t_SuSiEGLM08,t_lasso_cv,t_mcp_cv,t_lasso_bic,t_mcp_bic,t_ibss)
Main_TP[iter,]=c(tptn_SuSiEGLM08$tp,tptn_lasso_cv$tp,tptn_mcp_cv$tp,tptn_lasso_bic$tp,tptn_mcp_bic$tp,tptn_ibss$tp)
Main_TN[iter,]=c(tptn_SuSiEGLM08$tn,tptn_lasso_cv$tn,tptn_mcp_cv$tn,tptn_lasso_bic$tn,tptn_mcp_bic$tn,tptn_ibss$tn)

if(iter%%10==0) print(iter)

}

print( colMeans(Main_TP))
print( colMeans(Main_TN))
print( colMeans(TIME))
