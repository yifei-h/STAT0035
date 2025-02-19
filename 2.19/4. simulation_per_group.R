# Simulation for each group
simulation_per_group <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05) {
  library(car)
  source("1. MC_sample.R")
  
  count_group0 <- 0
  count_group1 <- 0
  coefficients_group0 <- list()
  coefficients_group1 <- list()
  
  for (i in 1:mc_iter) {
    sample_result <- MC_sample(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq)
    fit <- sample_result$fit
    coefs <- sample_result$coefs
    coefficients_group0[[i]] <- coefs
    coefficients_group1[[i]] <- coefs
    
    p_value_group0 <- tryCatch({
      as.numeric(coefs["A", "Pr(>|t|)"])
    }, error = function(e) NA)
    
    # Use tryCatch for linearHypothesis to prevent failure
    p_value_group1 <- tryCatch({
      as.numeric(linearHypothesis(fit, c("A + A:Z = 0"))$Pr[2])  
    }, error = function(e) NA) 
    
    if (!is.na(p_value_group0) && !is.null(p_value_group0) && p_value_group0 < alpha) {
      count_group0 <- count_group0 + 1
    }
    if (!is.na(p_value_group1) && !is.null(p_value_group1) && p_value_group1 < alpha) {
      count_group1 <- count_group1 + 1
    }
  }
  assurance_group0 <- count_group0 / mc_iter
  assurance_group1 <- count_group1 / mc_iter
  
  return(list(
    assurance_group0 = assurance_group0,
    assurance_group1 = assurance_group1
  ))
}
