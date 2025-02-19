# Simulation for the whole population
simulation_whole <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05) {
  source("1. MC_sample.R")
  count <- 0
  for (i in 1:mc_iter) {
    results <- MC_sample(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq)
    coefs <- results$coefs
    p_value <- tryCatch({
      if ("A" %in% rownames(coefs)) {
        coefs["A", "Pr(>|t|)"]
      } else {
        NA
      }
    }, error = function(e) NA) 
    if (p_value < alpha) {
      count <- count + 1
    }
  }
  assurance <- count / mc_iter
  return(list(assurance_population = assurance))
}
