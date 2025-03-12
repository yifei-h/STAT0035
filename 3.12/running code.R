r_files <- c("1. MC_sample.R", "2. simulation.R", "3. sample_size.R", "4. simulation_per_group.R", "5. simulation_varying_proportions.R")

for (file in r_files) {
  source(file)
}

set.seed(123)
n <- 50
p <- 2
beta <- c(1.5, -0.5)
tau0 <- 1.5
tau1 <- 1.0
gamma <- 1.0
sigsq <- 1.5
mc_iter <- 1000

# Run simulation
simulation_results <- simulation(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq, mc_iter = mc_iter)
assurance <- simulation_results$assurance
cat("Assurance:", assurance, "\n")

# Analyse coefficients
coefficients <- do.call(rbind, lapply(simulation_results$coefficients, function(res) res[, "Estimate"]))
mean_coefficients <- colMeans(coefficients, na.rm = TRUE)
print(mean_coefficients)

target <- 0.8
max_n <- 2000

# Run sample_size
result <- sample_size(target, p, beta, tau0, tau1, gamma, sigsq, mc_iter, max_n)
print(result)

# Run simulation_per_group
results_group <- simulation_per_group(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter)
print(results_group)

# Run simulation_varying_proportions
results_proportion <- simulation_varying_proportions(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter)
print(results_proportion)

# Visualisation
library(ggplot2)
ggplot(results_proportion, aes(x = Proportion_Group1, y = Assurance_Group1)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Assurance vs. Proportion of Group 1",
       x = "Proportion of Group 1",
       y = "Assurance")
