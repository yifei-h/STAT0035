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
