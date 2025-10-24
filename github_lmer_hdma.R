####### Load three internal functions before using hdma function to do mediation analysis.
####### Three internal functions are established by YinanZheng (2016) and posted in the website of 
####### https://github.com/YinanZheng/HIMA/blob/master/R/utils.R.  
####### There is need to download the script of utils.R from the website and then load into your workplace of R. 
lmer_hdma <- function(X, Y, M, group,
                      COV.XM = NULL, COV.MY = COV.XM,
                      family = c("gaussian", "binomial"),
                      method = c("lasso", "ridge"),
                      topN = NULL, parallel = FALSE, ncore = 1, verbose = FALSE) {
  ####################################################################################################################################
  #########################################                   Function body                ###########################################
  ####################################################################################################################################
  ####### INPUT
  ####### X : Independent variable that is a vector
  ####### Y : Dependent variable that is a vector and can be either continuous or binary variable
  ####### M : High-dimensional mediators that can be either data.frame or matrix. Rows represent samples, columns represent variables
  ####### COV.XM : a data.frame or matrix of covariates dataset for testing the association X ~ M. Default = NULL. 
  #######          If the covariates contain mixed types, please make sure all categorical variables are properly transformed into factor
  #######          type.
  ####### COV.MY : a data.frame or matrix of covariates dataset for testing the association Y ~ M. Using covariates should be careful.
  #######          If the cavariables are not specified, the covariates for Y ~ M are the same with that of M ~ X.
  ####### family : either 'gaussian' or 'binomial', relying on the type of outcome (Y). See hdi package.
  ####### method : either "lasso" or "ridge" to estimate the effect of M -> Y.
  ####### topN : an integer can be used to set the number of top markers by the method of sure independent screening. Default = NULL.
  #######        If topN is NULL, it will be either ceiling(n/log(n)) if family = 'gaussian', or ceiling(n/(2*log(n))) if family = 
  #######	      'binomial', where n is the sample size. If the sample size is greater than topN (pre-specified or calculated), all
  #######        markers will be harbored in the test.
  ####### verbose : logical. Default = FALSE.
  ####################################################################################################################################
  ####### checking the necessary packages 
  pkgs <- c("hdi", "MASS", "doParallel", "foreach", "iterators", "glmnet", "lme4","lmerTest")
  invisible(lapply(pkgs, require, character.only = TRUE))
  
  family <- match.arg(family)
  method <- match.arg(method)
  
  n <- nrow(M)
  p <- ncol(M)
  if (is.null(colnames(M))) colnames(M) <- paste0("M", seq_len(p))
  M_names <- colnames(M)
  
  d <- if (is.null(topN)) {
    if (family == "binomial") ceiling(n / (2 * log(n))) else ceiling(2 * n / log(n))
  } else {
    min(topN, p)
  }
  
  XM <- cbind(M, X)
  
  if (!is.null(COV.MY)) {
    COV.MY <- model.matrix(~ ., data = data.frame(COV.MY))[, -1, drop = FALSE]
    if (verbose) message("Adjusting for covariates: ", paste(colnames(COV.MY), collapse = ", "))
    XM_COV <- cbind(XM, COV.MY)
  } else {
    XM_COV <- XM
  }
  #  High-dimensional Inference (HDI) 
  set.seed(1029)
  fit <- switch(method,
                lasso = lasso.proj(XM_COV, Y, family = family),
                ridge = ridge.proj(XM_COV, Y, family = family))
  
  P_hdi <- fit$pval
  selected_index <- which(M_names %in% names(P_hdi[P_hdi <= 0.05]))
  if (length(selected_index) == 0) {
    warning("No mediators passed the significance threshold (p <= 0.05)")
    return(NULL)
  }
  if (verbose) message("Significant mediators selected: ", paste(M_names[selected_index], collapse = ", "))
  
  selected_M <- M[, selected_index, drop = FALSE]
  
  # Estimate alpha using lmer (X -> M with random effect)
  alpha_hat <- numeric(length(selected_index))
  alpha_p <- numeric(length(selected_index))
  
  for (i in seq_along(selected_index)) {
    mediator_i <- selected_M[, i]
    df_alpha <- data.frame(M = mediator_i, X = X, group = group)
    if (!is.null(COV.XM)) df_alpha <- cbind(df_alpha, COV.XM)
    
    # Construct the formula:M ~ X + covariates + (1 | group)
    fixed_part <- paste("X", if (!is.null(COV.XM)) paste0("+", paste(colnames(COV.XM), collapse = " + ")) else "")
    formula_alpha <- as.formula(paste("M ~", fixed_part, "+ (1 | group)"))
    
    # Fit the lmer model
    fit_alpha <- tryCatch({
      lmer(formula_alpha, data = df_alpha)
    }, error = function(e) {
      warning(paste0("Failed to fit alpha model for mediator ", M_names[selected_index[i]], ": ", e$message))
      return(NULL)
    })
    
    if (!is.null(fit_alpha)) {
      coef_summary <- summary(fit_alpha)$coefficients
      if ("X" %in% rownames(coef_summary)) {
        alpha_hat[i] <- coef_summary["X", "Estimate"]
        alpha_p[i] <- coef_summary["X", "Pr(>|t|)"]
      } else {
        alpha_hat[i] <- NA
        alpha_p[i] <- NA
      }
    } else {
      alpha_hat[i] <- NA
      alpha_p[i] <- NA
    }
  }
  beta_hat <- fit$bhat[selected_index]
  beta_p <- P_hdi[selected_index]
  ab_est <- alpha_hat * beta_hat
  p_joint <- apply(rbind(beta_p, alpha_p), 2, max)
  
  if (!is.null(COV.MY)) {
    Y_cov <- data.frame(Y = Y, X = X, COV.MY)
  } else {
    Y_cov <- data.frame(Y = Y, X = X)
  }
  df_lmm <- data.frame(Y_cov, group = group)
  df_lmm$group <- as.factor(df_lmm$group)
  
  gamma_model <- tryCatch({
    if (family == "binomial") {
      summary(glmer(Y ~ . + (1 | group), data = df_lmm, family = binomial()))
    } else {
      summary(lmer(Y ~ . + (1 | group), data = df_lmm))
    }
  }, error = function(e) {
    stop("Model fitting for total effect failed: ", e$message)
  })
  
  gamma_est <- gamma_model$coefficients["X", "Estimate"]
  
  results <- data.frame(
    mediator = M_names[selected_index],
    alpha = alpha_hat,
    beta = beta_hat,
    gamma = gamma_est,
    `alpha*beta` = ab_est,
    `%total effect` = ab_est / gamma_est * 100,
    `P.value` = p_joint,
    check.names = FALSE
  )
  
  message("Done! (", Sys.time(), ")")
  return(results)
}
