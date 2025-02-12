# simulation_per_group: simulation for each group
simulation_per_group <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05) {
  count_group0 <- 0
  count_group1 <- 0
  coefficients_group0 <- list()
  coefficients_group1 <- list()
  
  for (i in 1:mc_iter) {
    coefs <- MC_sample(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq)
    coefficients_group0[[i]] <- coefs
    coefficients_group1[[i]] <- coefs
    p_value_group0 <- coefs["A", "Pr(>|t|)"]
    p_value_group1 <- coefs["A", "Pr(>|t|)"] #figureout the true p-value
    
    if (p_value_group0 < alpha) {
      count_group0 <- count_group0 + 1
    }
    if (p_value_group1 < alpha) {
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

results_group <- simulation_per_group(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter)
print(results_group)