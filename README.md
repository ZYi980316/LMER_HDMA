# LMER_HDMA

High-dimensional mediation analysis with a random intercept using **lme4** for `X → M` (alpha path) and **hdi** (lasso/ridge projection) for `M → Y` (beta path). Supports **Gaussian** and **binomial** outcomes.

> ⚠️ Interpretation note: when `family = "binomial"`, effects are on the log-odds scale. The product `alpha*beta` and `% total effect = (alpha*beta) / gamma` are **link-scale** decompositions and are not guaranteed to equal natural indirect/direct effects without additional identification assumptions.

## Features
- Optional **SIS** pre-screening to reduce mediators from `p` to `d` (default: `d = 2n/log n` for Gaussian; `d = n/(2 log n)` for binomial).
- **HDI** (`hdi::lasso.proj` or `ridge.proj`) for p-values of `M → Y`.
- **Random-intercept** `lmer/glmer` for `X → M` and `X → Y`.
- Combine p-values via `max` (default) or **Fisher's method**.
- Optional **parallel** alpha-fitting via `foreach`/`doParallel`.

## Installation
```r
# install.packages(c("devtools"))
# devtools::install_github("yourname/lmer_hdma")
