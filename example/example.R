library(ncvreg)
library(devtools)
library(logisticsusie)
coef.binsusie <- function(fit, data) {
  b <- colSums(get_alpha(fit) * get_mu(fit))
  return(b)
}
document()
data("Heart")
Z=Heart$X[,c("famhist","age")]
X=Heart$X[,setdiff(colnames(Heart$X),c("famhist","age"))]
X[,"tobacco"]=log(X[,"tobacco"]+1)
X[,"alcohol"]=log(X[,"alcohol"]+1)
y=Heart$y
fit1=ncvreg(X=Heart$X,y=Heart$y,family="binomial")
fit3=SuSiE_IRLS(X=X,Z=Z,y=y,L=5,family="binomial",coverage=0.7)

data("Prostate")
Z=Prostate$X[,c("age","svi")]
X=Prostate$X[,setdiff(colnames(Prostate$X),c("age","svi"))]
y=exp(Prostate$y)
fit1=ncvreg(X=Prostate$X,y=exp(Prostate$y),family="poisson")
fit3=SuSiE_IRLS(X=X,Z=Z,y=y,L=5,family=Gamma(),coverage=0.95,estimate_residual_variance=F)
fit1$beta[,which.min(BIC(fit1))]
fit3$main_index
fit3$fitX$sigma2
summary(fit3$fitJoint)

