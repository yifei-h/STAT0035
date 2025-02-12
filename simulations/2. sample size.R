# sample size: find the sample size for a given assurance
sample_size <- function(target, p, beta, tau0, tau1, gamma, sigsq, mc_iter, max_n, alpha = 0.05) {
  n <- 10
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
