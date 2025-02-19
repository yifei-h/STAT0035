# Hepler function to simulate one dataset
MC_sample <- function(n, p, beta, tau0, tau1, gamma, sigsq){
  #' n: number of individuals
  #' p: number of covariates
  #' beta: effect coefficients
  #' tau0: treatment effect for Group 0
  #' tau1: additional treatment effect for Group 1  
  #' gamma: baseline for Group 1
  #' sigqs: variance of noise
  
  Xn <- matrix(rnorm(n * p), nrow = n, ncol = p) # Design matrix
  A <- rbinom(n, 1, 0.5) # Treatment indicator (binary)
  Z <- rbinom(n, 1, 0.5) # Group indicator (binary)
  epsilon <- rnorm(n, mean = 0, sd = sqrt(sigsq)) # Noise
  
  # Simulate y
  y <- Xn %*% beta + A * tau0 + A * Z * tau1 + Z * gamma + epsilon
  fit <- lm(y ~ Xn + A + A * Z + Z)
  return(list(fit = fit, coefs = summary(fit)$coefficients))
}


# Simulation
simulation <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05) {
  
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

set.seed(123)
n <- 50
p <- 2
beta <- c(1.5, -0.5)
tau0 <- 1.5
tau1 <- 1.0
gamma <- 1.0
sigsq <- 1.5
mc_iter <- 1000

simulation_results <- simulation(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq, mc_iter = mc_iter)
assurance <- simulation_results$assurance
cat("Assurance:", assurance, "\n")

# Analyse coefficients
coefficients <- do.call(rbind, lapply(simulation_results$coefficients, function(res) res[, "Estimate"]))
mean_coefficients <- colMeans(coefficients, na.rm = TRUE)
print(mean_coefficients)



# Find the sample size for a given assurance
sample_size <- function(target, p, beta, tau0, tau1, gamma, sigsq, mc_iter, max_n, alpha = 0.05) {
  n <- 20
  while (n <= max_n) {
    print(n)
    my_sim <- simulation(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq, mc_iter = mc_iter, alpha = alpha)
    assur <- my_sim$assurance
    if (assur >= target) {
      return(list(
        Sample_Size = n,
        Assurance = assurance))
    }
    n <- n + 10
  }
  return("No sample size within the range achieves the target assurance.")
}

target <- 0.8
max_n <- 2000

result <- sample_size(target, p, beta, tau0, tau1, gamma, sigsq, mc_iter, max_n)
print(result)



# Simulation for each group
simulation_per_group <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05) {
  library(car)
  
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

results_group <- simulation_per_group(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter)
print(results_group)



# Simulation for different proportions of the two groups
simulation_varying_proportions <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05, proportions = seq(0, 1, by = 0.1)) {
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



# Simulation for the whole population
simulation_whole <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05) {
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

result_whole <- simulation_whole(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter)
print(result_whole)

