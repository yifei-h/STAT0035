# Simulation
simulation <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05) {
  source("1. MC_sample.R")
  # Monte Carlo simulation
  count <- 0
  coefficients <- list()
  
  for (i in 1:mc_iter) {
    results <- MC_sample(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq)
    coefs <- results$coefs
    coefficients[[i]] <- coefs
    
    p_value <- tryCatch({
      if ("A" %in% rownames(coefs)) {
        coefs["A", "Pr(>|t|)"]
      } else {
        NA
      }
    }, error = function(e) NA) 
    
    if (!is.na(p_value) && p_value < alpha) {
      count <- count + 1
    }
  }
  
  assurance <- count / mc_iter
  return(list(
    coefficients = coefficients,
    assurance = assurance))
}