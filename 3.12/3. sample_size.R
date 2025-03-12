# Find the sample size for a given assurance
sample_size <- function(target, p, beta, tau0, tau1, gamma, sigsq, mc_iter, max_n, alpha = 0.05) {
  source("2. simulation.R")
  n <- 20
  error_count <- 0
  while (n <= max_n) {
    print(n)
    my_sim <- simulation(n = n, p = p, beta = beta, tau0 = tau0, tau1 = tau1, gamma = gamma, sigsq = sigsq, mc_iter = mc_iter, alpha = alpha)
    p_values <- sapply(my_sim$coefficients, function(coefs) {
      tryCatch({
        coefs["A1", "Pr(>|t|)"]
      }, error = function(e) {
        error_count <<- error_count + 1  # Track errors
        NA
      })
    })
    
    count <- sum(!is.na(p_values) & p_values < alpha)
    assurance <- count / mc_iter
    
    if (assurance >= target) {
      return(list(
        Sample_Size = n,
        Assurance = assurance,
        Errors = error_count # Return number of errors encountered
      ))
    }
    n <- n + 10
  }
  return("No sample size within the range achieves the target assurance.")
}
