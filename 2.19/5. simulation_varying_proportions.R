# Simulation for different proportions of the two groups
simulation_varying_proportions <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05, proportions = seq(0, 1, by = 0.1)) {
  source("1. MC_sample.R")
  results <- data.frame(Proportion_Group1 = numeric(), 
                        Assurance_Group0 = numeric(),
                        Assurance_Group1 = numeric())
  for (prop in proportions) {
    n_group1 <- floor(n * prop)  # Number of individuals in Group 1
    n_group0 <- n - n_group1
    count_group0 <- 0
    count_group1 <- 0
    for (i in 1:mc_iter) {
      # Generate data separately for each group
      Xn_group0 <- matrix(rnorm(n_group0 * p), nrow = n_group0, ncol = p)
      Xn_group1 <- matrix(rnorm(n_group1 * p), nrow = n_group1, ncol = p)
      A_group0 <- rbinom(n_group0, 1, 0.5)
      A_group1 <- rbinom(n_group1, 1, 0.5)
      Z_group0 <- rep(0, n_group0)
      Z_group1 <- rep(1, n_group1)
      
      # Combine data for the whole population
      Xn <- rbind(Xn_group0, Xn_group1)
      A <- c(A_group0, A_group1)
      Z <- c(Z_group0, Z_group1)
      
      epsilon <- rnorm(n, mean = 0, sd = sqrt(sigsq))
      y <- Xn %*% beta + A * tau0 + A * Z * tau1 + Z * gamma + epsilon
      fit <- lm(y ~ Xn + A + A:Z + Z)
      
      results_fit <- MC_sample(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq)
      coefs <- results_fit$coefs
      
      p_value_group0 <- tryCatch({
        if ("A" %in% rownames(coefs)) {
          coefs["A", "Pr(>|t|)"]  # Group 0 treatment effect significance
        } else {
          NA
        }
      }, error = function(e) NA)
      
      p_value_group1 <- tryCatch({
        if ("A:Z" %in% rownames(coefs)) {
          coefs["A:Z", "Pr(>|t|)"]  # Interaction term significance
        } else {
          NA
        }
      }, error = function(e) NA)
      
      if (!is.na(p_value_group0) && p_value_group0 < alpha) {
        count_group0 <- count_group0 + 1
      }
      if (!is.na(p_value_group1) && p_value_group1 < alpha) {
        count_group1 <- count_group1 + 1
      }
    }
    
    assurance_group0 <- count_group0 / mc_iter
    assurance_group1 <- count_group1 / mc_iter
    new_row <- data.frame(
      Proportion_Group1 = as.numeric(prop), 
      Assurance_Group0 = as.numeric(assurance_group0),
      Assurance_Group1 = as.numeric(assurance_group1))
    results <- rbind(results, new_row)
  }
  return(results)
}
