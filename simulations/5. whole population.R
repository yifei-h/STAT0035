# simulation_whole: simulation for the whole population
simulation_whole <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05) {
  count <- 0
  for (i in 1:mc_iter) {
    coefs <- MC_sample(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq)
    p_value <- coefs["A", "Pr(>|t|)"]
    if (p_value < alpha) {
      count <- count + 1
    }
  }
  assurance <- count / mc_iter
  return(list(assurance_population = assurance))
}

result_whole <- simulation_whole(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter)
print(result_whole)