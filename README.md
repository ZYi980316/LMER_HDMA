# LMER_HDMA

High-dimensional mediation analysis with a random intercept using **lme4** for `X → M` (alpha path) and **hdi** (lasso/ridge projection) for `M → Y` (beta path). Supports **Gaussian** and **binomial** outcomes.

> ⚠️ Interpretation note: when `family = "binomial"`, effects are on the log-odds scale. The product `alpha*beta` and `% total effect = (alpha*beta) / gamma` are **link-scale** decompositions and are not guaranteed to equal natural indirect/direct effects without additional identification assumptions.

## Features
- Optional **SIS** pre-screening to reduce mediators from `p` to `d` (default: `d = 2n/log n` for Gaussian; `d = n/(2 log n)` for binomial).
- **HDI** (`hdi::lasso.proj` or `ridge.proj`) for p-values of `M → Y`.
- **Random-intercept** `lmer/glmer` for `X → M` and `X → Y`.
- Combine p-values via `max` (default) or **Fisher's method**.
- Optional **parallel** alpha-fitting via `foreach`/`doParallel`.

## Uasage

library(lme4); library(hdi)

set.seed(1)  
n <- 200; p <- 500; G <- 20  
group <- sample(letters[1:G], n, replace = TRUE)  
X <- rnorm(n)  

# mediators with a few true signals  
M <- matrix(rnorm(n*p), n, p)  
colnames(M) <- paste0("M", 1:p)  
M[,1] <- 0.5*X + rnorm(n)  
M[,2] <- -0.6*X + rnorm(n)  

# outcome (Gaussian)    
beta <- c(0.8, -0.9, rep(0, p-2))  
Y_lin <- 0.7*X + M %*% beta + rnorm(n, sd = 1)  
Y <- as.numeric(Y_lin)  

fit <- lmer_hdma(
  X = X, Y = Y, M = M, group = group,
  family = "gaussian", method = "lasso",
  topN = NULL, parallel = FALSE, verbose = TRUE, seed = 123
)  

str(fit$results)  
head(fit$results[order(fit$results$p_joint), ], 10)  

