# psmqte

`psmqte` estimates treatment effects by propensity-score matching in R. It is designed for empirical workflows where users want quantile treatment effects as well as the more familiar average treatment effects from the same matching design.

The package currently supports:

- quantile treatment effects (QTE),
- quantile treatment effects on the treated (QTT),
- average treatment effects (ATE), and
- average treatment effects on the treated (ATT).

Matching is nearest-neighbor matching with replacement on an estimated or user-supplied propensity score. The default propensity score is estimated by logistic regression.

## Installation

Install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("Li-Lu-SoB/psmqte")
```

## Quick Start

The package includes a self-contained synthetic demo data set. It is useful for
checking installation, learning the function interface, and reproducing examples
without downloading external data.

```r
library(psmqte)

dat <- psmqte_example_data()

fit <- psmqte(
  Y ~ X1 + X2 + X3 + age + size,
  treatment = "D",
  data = dat,
  taus = seq(0.1, 0.9, by = 0.1),
  M = c(1, 2, 4, 8),
  se = FALSE
)

head(summary(fit))
plot(fit, estimand = "qte")
plot(fit, estimand = "qtt")
```

Set `se = TRUE` to compute feasible influence-function standard errors and
pointwise confidence intervals:

```r
fit_se <- psmqte(
  Y ~ X1 + X2 + X3 + age + size,
  treatment = "D",
  data = dat,
  taus = c(0.25, 0.5, 0.75),
  M = c(1, 4),
  se = TRUE
)

summary(fit_se)
```

## Average Effects

Use `psmate()` to estimate ATE and ATT with the same style of formula interface:

```r
avg_fit <- psmate(
  Y ~ X1 + X2 + X3 + age + size,
  treatment = "D",
  data = dat,
  M = c(1, 2, 4, 8)
)

summary(avg_fit)
```

## Using Your Own Data

Your data should contain:

- one numeric outcome variable;
- one binary treatment variable coded as `0/1`, `TRUE/FALSE`, or
  `control/treated`;
- observed covariates used to estimate the propensity score.

Example:

```r
fit <- psmqte(
  outcome ~ x1 + x2 + x3,
  treatment = "treated",
  data = my_data,
  taus = seq(0.1, 0.9, by = 0.1),
  M = c(1, 2, 4),
  se = TRUE
)

summary(fit)
```

To use a different propensity-score specification from the outcome formula,
provide `ps_formula`:

```r
fit <- psmqte(
  outcome ~ x1 + x2,
  treatment = "treated",
  data = my_data,
  ps_formula = treated ~ x1 + x2 + x3 + x4,
  M = c(1, 2, 4)
)
```

To match on an externally estimated propensity score, provide `ps_score`:

```r
my_data$ps <- fitted(glm(treated ~ x1 + x2 + x3, data = my_data,
                         family = binomial()))

fit <- psmqte(
  outcome ~ x1 + x2 + x3,
  treatment = "treated",
  data = my_data,
  ps_score = my_data$ps,
  M = c(1, 2, 4)
)
```

## Diagnostics

The fitted object stores propensity-score overlap and match-reuse diagnostics:

```r
fit$diagnostics$overlap
fit$diagnostics$reuse
```

For covariate balance, use `psm_balance()`:

```r
psm_balance(dat, treatment = "D", covariates = c("X1", "X2", "X3", "age"))
```

The matching matrix and reuse counts are also available directly:

```r
match_obj <- psm_match(dat$p_true, dat$D, M = 4)
head(match_obj$match_matrix)
head(match_obj$K_all)
```

## Example Data

`psmqte_example_data()` loads the package's built-in demo data:

```r
demo <- psmqte_example_data()
```

The demo data are synthetic and generated specifically for this package. They
include observed variables (`Y`, `D`, `X1`, `X2`, `X3`, `age`, `size`) and
simulation-only variables (`Y0`, `Y1`, `p_true`) that are useful for teaching
and checking examples.

Fresh simulated data can also be generated:

```r
sim <- psmqte_example_data("simulated", n = 2000, seed = 123)
```

An optional real-data example is available through the `causaldata` package:

```r
install.packages("causaldata")
nhefs <- psmqte_example_data("nhefs")
```

## Main Functions

- `psmqte()` estimates QTE and QTT curves.
- `psmate()` estimates ATE and ATT.
- `psm_match()` constructs nearest-neighbor matches and reuse counts.
- `psm_balance()` computes standardized mean differences.
- `m_grid_diverging()` creates a diverging-M grid.
- `psmqte_example_data()` loads example data.
- `simulate_psmqte_data()` generates synthetic data.

## Inference

For QTE and QTT, `psmqte()` can compute feasible influence-function standard
errors using local same-arm nearest-neighbor estimates of conditional CDF terms.
When the propensity score is estimated internally by logistic regression, the
default uses a first-step score adjustment.

For QTT, the estimated-propensity-score adjustment uses a practical
finite-difference approximation to the target-drift term. For ATE and ATT,
`psmate()` reports large-sample matching standard errors.

## References

Abadie, A. and Imbens, G. W. (2016). Matching on the estimated propensity
score. *Econometrica*, 84(2), 781-807.

He, Y. and Han, F. (2024). On propensity score matching with a diverging number
of matches. *Biometrika*, 111(4), 1421-1428.
