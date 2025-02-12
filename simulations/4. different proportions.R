# simulation_varying_proportions: simulation for different proportions of the two groups
simulation_varying_proportions <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05, proportions = seq(0, 1, by = 0.1)) {
  results <- data.frame(Proportion_Group1 = numeric(), Assurance = numeric())
  for (prop in proportions) {
    n_group1 <- round(n * prop)  # Number of individuals in Group 1
    n_group0 <- n - n_group1
    count <- 0
    for (i in 1:mc_iter) {
      coefs <- MC_sample(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq)
      p_value <- coefs["A", "Pr(>|t|)"] 
      if (p_value < alpha) {
        count <- count + 1
      }
    }
    assurance <- count / mc_iter
    results <- rbind(results, data.frame(Proportion_Group1 = prop, Assurance = assurance))
  }
  return(results)
}

results_proportion <- simulation_varying_proportions(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter)
print(results_proportion)

# Visualisation
library(ggplot2)
ggplot(results_proportion, aes(x = Proportion_Group1, y = Assurance)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Assurance vs. Proportion of Group 1",
       x = "Proportion of Group 1",
       y = "Assurance")