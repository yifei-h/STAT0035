# Define the simulation model, return the assurance and model coefficients

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
  #fit <- lm(y ~ Xn + A + A * Z + Z)
  fit <- lm(y ~ Xn + A + Z)
  
  return(summary(fit)$coefficients)
}


# simulation: run the simulation model, return the assurance and coefficients
simulation <- function(n, p, beta, tau0, tau1, gamma, sigsq, mc_iter, alpha = 0.05) {
  
  # Monte Carlo simulation
  count <- 0
  coefficients <- list()
  
  for (i in 1:mc_iter) {
    coefs <- MC_sample(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq)
    coefficients[[i]] <- coefs
    #print(coefs)
    p_value <- coefs["A", "Pr(>|t|)"]
    if (p_value < alpha) {
      count <- count + 1
    }
  }
  
  # Calculate assurance
  assurance <- count / mc_iter
  # Return results
  return(list(
    coefficients = coefficients,
    assurance = assurance))
}

set.seed(123)
n <- 50
p <- 2
beta <- c(1.5, -0.5)
tau0 <- 2.0
tau1 <- -1.0
gamma <- 1.0
sigsq <- 1.5
mc_iter <- 1000

simulation_results <- simulation(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq, mc_iter = mc_iter)

# Return assurance
assurance <- simulation_results$assurance
cat("Assurance:", assurance, "\n")

# Analyse coefficients
coefficients <- do.call(rbind, lapply(simulation_results$coefficients, function(res) res[, "Estimate"]))
mean_coefficients <- colMeans(coefficients, na.rm = TRUE)
print(mean_coefficients)