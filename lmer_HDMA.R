#' High-dimensional mediation with random intercept (lmer + HDI)
#' @param X Numeric vector (exposure), length = n
#' @param Y Numeric (gaussian) or binary 0/1 (binomial) outcome, length = n
#' @param M Matrix/data.frame of mediators (n x p), numeric columns
#' @param group Grouping factor/id for random intercept (length = n)
#' @param COV.XM Covariates for M ~ X (data.frame/matrix, n rows) (optional)
#' @param COV.MY Covariates for Y ~ M + X (data.frame/matrix, n rows) (optional). Default: same as COV.XM
#' @param family "gaussian" or "binomial" (for Y model / HDI link)
#' @param method "lasso" or "ridge" (HDI engine)
#' @param topN Optional integer; if NULL, use SIS size d = ceiling(2*n/log(n)) for gaussian or ceiling(n/(2*log(n))) for binomial
#' @param parallel Logical; fit alpha models in parallel (default FALSE)
#' @param ncore Integer; number of cores when parallel = TRUE
#' @param alpha Numeric; p-value threshold for HDI screening (default 0.05)
#' @param combine_p "max" (default) or "fisher" to combine p_alpha & p_beta
#' @param standardize Logical; standardize M columns before HDI (default TRUE)
#' @param seed Optional integer for reproducibility (default NULL)
#' @param verbose Logical
#' @return A list with: results data.frame and metadata
lmer_hdma <- function(
    X, Y, M, group,
    COV.XM = NULL, COV.MY = COV.XM,
    family = c("gaussian", "binomial"),
    method = c("lasso", "ridge"),
    topN = NULL, parallel = FALSE, ncore = 1,
    alpha = 0.05, combine_p = c("max","fisher"),
    standardize = TRUE, seed = NULL, verbose = FALSE
) {
  # ---- deps ----
  pkgs <- c("hdi","glmnet","lme4","lmerTest","foreach","doParallel")
  for (pkg in pkgs) if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required but not installed.")
  }
  library(lme4)
  library(lmerTest)
  family  <- match.arg(family)
  method  <- match.arg(method)
  combine_p <- match.arg(combine_p)
  
  # ---- coerce & checks ----
  M <- as.matrix(M)
  storage.mode(M) <- "numeric"
  n <- nrow(M); p <- ncol(M)
  if (is.null(colnames(M))) colnames(M) <- paste0("M", seq_len(p))
  stopifnot(length(X) == n, length(Y) == n, length(group) == n)
  if (!is.null(COV.XM) && nrow(COV.XM) != n) stop("'COV.XM' must have n rows.")
  if (!is.null(COV.MY) && nrow(COV.MY) != n) stop("'COV.MY' must have n rows.")
  
  # remove rows with NA
  aux <- list(X=X, Y=Y, group=group)
  if (!is.null(COV.XM)) aux$COV.XM <- COV.XM
  if (!is.null(COV.MY)) aux$COV.MY <- COV.MY
  keep <- complete.cases(cbind(X, Y, group, M,
                               if (is.null(COV.XM)) NULL else COV.XM,
                               if (is.null(COV.MY)) NULL else COV.MY))
  if (sum(!keep) > 0) {
    if (verbose) message("Removed ", sum(!keep), " rows with missing values.")
    X <- X[keep]; Y <- Y[keep]; group <- group[keep]
    M <- M[keep, , drop = FALSE]
    if (!is.null(COV.XM)) COV.XM <- COV.XM[keep, , drop = FALSE]
    if (!is.null(COV.MY)) COV.MY <- COV.MY[keep, , drop = FALSE]
    n <- nrow(M)
  }
  
  # standardize M (optional)
  if (standardize) {
    M <- scale(M)
  }
  
  # sane names
  colnames(M) <- make.names(colnames(M), unique = TRUE)
  M_names <- colnames(M)
  if (!is.null(COV.XM)) colnames(COV.XM) <- make.names(colnames(COV.XM), unique = TRUE)
  if (!is.null(COV.MY)) colnames(COV.MY) <- make.names(colnames(COV.MY), unique = TRUE)
  
  # group as factor
  group <- factor(group)
  
  # default SIS size d
  d <- if (is.null(topN)) {
    if (family == "binomial") ceiling(n / (2 * log(n))) else ceiling(2 * n / log(n))
  } else {
    min(topN, p)
  }
  
  # ---- SIS (keep top d mediators if p>d) ----
  idx_all <- seq_len(p)
  if (p > d) {
    # combine |cor(M, X)| and |cor(M, Y)| ranks
    y_num <- if (family == "binomial") as.numeric(Y) else Y
    r1 <- suppressWarnings(abs(stats::cor(M, X, use = "pairwise.complete.obs")))
    r2 <- suppressWarnings(abs(stats::cor(M, y_num, use = "pairwise.complete.obs")))
    r1[is.na(r1)] <- 0; r2[is.na(r2)] <- 0
    rank1 <- rank(-as.numeric(r1), ties.method = "min")
    rank2 <- rank(-as.numeric(r2), ties.method = "min")
    score <- pmin(rank1, rank2)
    keep_sis <- order(score)[seq_len(d)]
    M_sis <- M[, keep_sis, drop = FALSE]
    M_sis_names <- colnames(M_sis)
    if (verbose) message("SIS reduced mediators from ", p, " to ", ncol(M_sis), ".")
  } else {
    M_sis <- M
    M_sis_names <- M_names
  }
  
  # ---- HDI step: Y ~ M + X (+ COV.MY) ----
  XM <- cbind(M_sis, X = X)
  if (!is.null(COV.MY)) {
    # expand covariates for MY model
    COV.MY.mm <- model.matrix(~ ., data = as.data.frame(COV.MY))[, -1, drop = FALSE]
    XM_COV <- cbind(XM, COV.MY.mm)
  } else {
    XM_COV <- XM
  }
  
  if (!is.null(seed)) set.seed(seed)
  hdi_fit <- switch(method,
                    lasso = hdi::lasso.proj(XM_COV, Y, family = family),
                    ridge = hdi::ridge.proj(XM_COV, Y, family = family))
  
  P_hdi <- hdi_fit$pval
  bhat  <- hdi_fit$bhat
  # first columns correspond to mediators after cbind(M_sis, X, cov)
  p_M <- ncol(M_sis)
  beta_p <- as.numeric(P_hdi[seq_len(p_M)])
  beta_hat <- as.numeric(bhat[seq_len(p_M)])
  
  sel <- which(beta_p <= alpha)
  if (length(sel) == 0) {
    warning("No mediators passed HDI threshold (p <= ", alpha, ").")
    return(list(
      results = NULL,
      meta = list(n = n, p = p, p_after_SIS = p_M, method = method, family = family,
                  alpha = alpha, combine_p = combine_p)
    ))
  }
  if (verbose) message("HDI selected ", length(sel), " mediators.")
  
  selected_M  <- M_sis[, sel, drop = FALSE]
  selected_nm <- M_sis_names[sel]
  
  # ---- alpha step: M_i ~ X + COV.XM + (1|group) ----
  get_alpha <- function(m_i) {
    df <- data.frame(Mi = m_i, X = X, group = group)
    if (!is.null(COV.XM)) df <- cbind(df, COV.XM)
    # build formula safely
    rhs <- c("X", if (!is.null(COV.XM)) colnames(COV.XM) else NULL)
    fml <- as.formula(paste("Mi ~", paste(rhs, collapse = " + "), "+ (1 | group)"))
    fit <- tryCatch(lmerTest::lmer(fml, data = df), error = function(e) e)
    if (inherits(fit, "error")) return(c(alpha = NA, p = NA))
    cf <- summary(fit)$coefficients
    if (!"X" %in% rownames(cf)) return(c(alpha = NA, p = NA))
    c(alpha = unname(cf["X","Estimate"]),
      p     = unname(cf["X","Pr(>|t|)"]))
  }
  
  if (parallel && ncore > 1) {
    doParallel::registerDoParallel(ncore)
    alpha_mat <- foreach::foreach(i = seq_len(ncol(selected_M)), .combine = rbind) %dopar% {
      get_alpha(selected_M[, i])
    }
    doParallel::stopImplicitCluster()
  } else {
    alpha_mat <- t(apply(selected_M, 2, get_alpha))
  }
  
  alpha_hat <- as.numeric(alpha_mat[, "alpha"])
  alpha_p   <- as.numeric(alpha_mat[, "p"])
  
  # ---- combine p-values ----
  if (combine_p == "max") {
    p_joint <- pmax(alpha_p, beta_p[sel], na.rm = TRUE)
  } else {
    # Fisher's method (df = 2*k with k=2)
    X2 <- -2 * (log(alpha_p) + log(beta_p[sel]))
    p_joint <- stats::pchisq(X2, df = 4, lower.tail = FALSE)
  }
  
  # ---- total effect gamma: Y ~ X + COV.MY + (1|group) ----
  Y_cov <- data.frame(Y = Y, X = X, group = group)
  if (!is.null(COV.MY)) Y_cov <- cbind(Y_cov, COV.MY)
  gamma_fit <- tryCatch({
    if (family == "binomial") {
      summary(lme4::glmer(Y ~ . + (1 | group), data = Y_cov, family = stats::binomial()))
    } else {
      summary(lme4::lmer(Y ~ . + (1 | group), data = Y_cov))
    }
  }, error = function(e) stop("Total-effect model failed: ", e$message))
  
  gamma_est <- tryCatch(unname(gamma_fit$coefficients["X","Estimate"]),
                        error = function(e) NA_real_)
  pct_total <- if (!is.na(gamma_est) && abs(gamma_est) > 1e-8) {
    (alpha_hat * beta_hat[sel]) / gamma_est * 100
  } else {
    rep(NA_real_, length(sel))
  }
  
  results <- data.frame(
    mediator = selected_nm,
    alpha = alpha_hat,
    beta  = beta_hat[sel],
    gamma = gamma_est,
    ab    = alpha_hat * beta_hat[sel],
    pct_total = pct_total,
    p_alpha = alpha_p,
    p_beta  = beta_p[sel],
    p_joint = p_joint,
    stringsAsFactors = FALSE
  )
  
  if (verbose) message("Done! (", Sys.time(), ")")
  list(
    results = results,
    meta = list(n = n, p = p, p_after_SIS = p_M, method = method, family = family,
                alpha = alpha, combine_p = combine_p, standardize = standardize)
  )
}
