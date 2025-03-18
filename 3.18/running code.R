r_files <- c("1. MC_sample.R", "2. simulation.R", "3. sample_size.R", "4. simulation_per_group.R", "5. simulation_varying_proportions.R")

for (file in r_files) {
  source(file)
}

set.seed(123)
n <- 50
p <- 2
beta <- c(1.5, -0.5)
tau0 <- 2.0
tau1 <- -2.0
gamma <- 2.0 #assurance decreases as it increases
sigsq <- 3.0 #assurance decreases as it increases
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

#Create data frames
n_values <- seq(50, 500, by = 50)
tau0_values <- seq(-3, 3, length.out = 10)
tau1_values <- seq(-3, 3, length.out = 10)
prop_values <- seq(0, 1, by = 0.1)
gamma_values <- seq(0, 3, length.out = 5)
sigsq_values <- seq(0.5, 3, length.out = 5)

assurance_values <- numeric(length(n_values))

assurance_tau <- expand.grid(tau0 = tau0_values, tau1 = tau1_values)
assurance_tau$assurance <- NA

assurance_prop <- data.frame(Proportion_Group1 = prop_values, Assurance = NA)

assurance_gamma_sigsq <- expand.grid(gamma = gamma_values, sigsq = sigsq_values)
assurance_gamma_sigsq$assurance <- NA

# Run simulations for different n
for (i in seq_along(n_values)) {
  result <- simulation(n = n_values[i], p = p, beta = beta, 
                       tau0 = tau0_values[i], tau1 = tau1_values[i], 
                       gamma = gamma, sigsq = sigsq, mc_iter = mc_iter)
  assurance_values[i] <- result$assurance
}
df_n <- data.frame(Sample_Size = n_values, Assurance = assurance_values)

# Run simulations for tau0 and tau1
for (i in 1:nrow(assurance_tau)) {
  sim_result <- simulation(n = 50, p = 2, beta = c(1.5, -0.5),
                           tau0 = assurance_tau$tau0[i], tau1 = assurance_tau$tau1[i],
                           gamma = 1, sigsq = 1.5, mc_iter = 1000)
  assurance_tau$assurance[i] <- sim_result$assurance
}

# Run simulations for proportions
results_proportion <- simulation_varying_proportions(n = 50, p = 2, beta = c(1.5, -0.5),
                                                     tau0 = 2.0, tau1 = -2.0,
                                                     gamma = 1, sigsq = 1.5, mc_iter = 1000)
assurance_prop$Assurance <- results_proportion$Assurance_Group1

# Run simulations for gamma and sigsq
for (i in 1:nrow(assurance_gamma_sigsq)) {
  sim_result <- simulation(n = 50, p = 2, beta = c(1.5, -0.5),
                           tau0 = 2.0, tau1 = -2.0,
                           gamma = assurance_gamma_sigsq$gamma[i], sigsq = assurance_gamma_sigsq$sigsq[i],
                           mc_iter = 1000)
  assurance_gamma_sigsq$assurance[i] <- sim_result$assurance
}

# Plot Assurance vs Sample Size
ggplot(df_n, aes(x = Sample_Size, y = Assurance)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  theme_minimal() +
  labs(title = "Assurance vs. Sample Size (n)",
       x = "Sample Size (n)",
       y = "Assurance")
ggsave(filename = "assurance vs sample size.png", plot = last_plot(), width = 8, height = 6, dpi = 300)

# Plot assurance for tau0 and tau1
ggplot(assurance_tau, aes(x = tau0, y = tau1, fill = assurance)) +
  geom_tile() +
  geom_text(aes(label = round(assurance, 2)), color = "white", size = 4) + 
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Assurance vs Tau0 and Tau1", x = "Tau0", y = "Tau1")
ggsave(filename = "assurance vs tau.png", plot = last_plot(), width = 8, height = 6, dpi = 300)

# Plot assurance for different proportions of Groups
ggplot(assurance_prop, aes(x = Proportion_Group1, y = Assurance)) +
  geom_line(color = "blue", size = 1, aes(linetype = "Group 1")) + 
  geom_point(color = "blue", size = 2) +
  geom_line(aes(x = 1 - Proportion_Group1, y = Assurance, linetype = "Group 0"), 
            color = "red", size = 1) + 
  geom_point(aes(x = 1 - Proportion_Group1, y = Assurance), 
             color = "red", size = 2) +
  theme_minimal() +
  labs(title = "Assurance vs Proportion of Groups", 
       x = "Proportion of Group 1", 
       y = "Assurance",
       linetype = "Group") +
  scale_linetype_manual(values = c("Group 1" = "solid", "Group 0" = "dashed"))
ggsave(filename = "assurance vs proportion of groups.png", plot = last_plot(), width = 8, height = 6, dpi = 300)

# Plot assurance for gamma and sigsq
ggplot(assurance_gamma_sigsq, aes(x = gamma, y = sigsq, fill = assurance)) +
  geom_tile() +
  geom_text(aes(label = round(assurance, 2)), color = "white", size = 4) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  labs(title = "Assurance vs Gamma and Sigma^2", x = "Gamma", y = "Sigma^2")
ggsave(filename = "assurance vs gamma and sigma^2.png", plot = last_plot(), width = 8, height = 6, dpi = 300)
