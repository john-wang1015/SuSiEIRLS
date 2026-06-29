# SuSiEIRLS

`SuSiEIRLS` fits SuSiE-type fine-mapping models for non-Gaussian outcomes by
constructing approximate Gaussian sufficient statistics and passing them to
`susieR::susie_ss()`. The package separates three cases: general GLM/extended
GLM outcomes, binary outcomes with a Polya-Gamma working path, and Cox
proportional-hazards outcomes.

## GLM And Extended GLM Branch

For a GLM with linear predictor

$$
\eta = Z\alpha + X\beta,
$$

the algorithm uses the usual IRLS expansion around the current estimate. At
iteration \(t\), with mean \(\mu^{(t)}\), derivative
\(d\mu / d\eta\), and variance function \(V(\mu)\), the working response is

$$
z^{(t)}
  =
  \eta^{(t)}
  +
  \frac{y - \mu^{(t)}}{d\mu / d\eta},
$$

with working weights

$$
w_i^{(t)}
  =
  \left\{
    V(\mu_i^{(t)})
    \left(\frac{d\eta_i}{d\mu_i}\right)^2
  \right\}^{-1}.
$$

After projecting out the covariates \(Z\), the package constructs weighted
sufficient statistics \((X^\top W X, X^\top W z, z^\top W z)\) and runs SuSiE
on this local Gaussian approximation. If the IRLS expansion were taken at the
true linear predictor, the residual variance in this working Gaussian problem
would be one.

In finite samples, especially in early outer iterations, the current
\(\eta^{(t)}\) is not exact. We therefore allow the SuSiE residual variance to
move within a conservative interval, by default

$$
\sigma^2 \in [0.1, 1].
$$

The upper bound keeps an incidental overestimate of the working residual
variance from unnecessarily lowering power, while the lower bound prevents the
working regression from becoming too aggressive. Negative binomial, Tweedie,
beta, and other supported `mgcv` families follow this same working-IRLS logic;
additional distributional parameters are estimated by the `mgcv` fit and then
carried into the subsequent working updates.

## Binary Branch

For binary logistic outcomes, the package provides a Polya-Gamma-enhanced path.
Using the identity

$$
p(y_i \mid \eta_i)
  \propto
  \exp\{\kappa_i \eta_i\}
  \int
  \exp\left(-\frac{\omega_i \eta_i^2}{2}\right)
  p(\omega_i)\,d\omega_i,
  \qquad
  \kappa_i = y_i - \frac12,
$$

the binary likelihood is represented through a conditionally Gaussian working
problem. This path is used to improve signal detection, particularly when the
ordinary IRLS approximation is weak.

The final reported evidence is then recomputed through the standard IRLS
working likelihood. This final correction places the log Bayes factors on the
ordinary IRLS scale, so the PG step acts as a power-enhancing search path rather
than as a different reporting scale. If that final IRLS correction is unstable,
the package returns the PG fit and marks the returned main-effect summary with
`status = "PG IRLS"`.

## Cox Branch

For Cox proportional-hazards outcomes,

$$
\lambda(t \mid X, Z)
  =
  \lambda_0(t)
  \exp\{Z\alpha + X\beta\},
$$

there is no observed Gaussian response \(Y\) analogous to the GLM working
response. The package therefore constructs score-based sufficient statistics
from the Cox partial likelihood. In this branch, \(X^\top X\) is replaced by the
projected information matrix and \(X^\top y\) by the projected Cox score.

Because no literal working response is constructed, the residual scale is not
the same theoretical unit-variance object used in the GLM and extended-GLM IRLS
branches. The implementation therefore leaves a small degree of freedom for
SuSiE to estimate the residual variance, again using conservative bounds. This
keeps the Cox score approximation flexible without treating it as if it came
from an exact Gaussian response model.
