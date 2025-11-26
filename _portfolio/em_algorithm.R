# Hardy-Weinberg Equilibrium initial genotype frequencies
hardy_weinberg_freq <- function(q) {
  c((1 - q)^2, 2*q*(1 - q), q^2)
}
# E-Step
e_step <- function(data, params, alpha) {
  mu <- params$mu
  sigma <- params$sigma
  
  weights <- sapply(1:3, function(k) {
    alpha[k] * dnorm(data, mean = mu[k], sd = sqrt(sigma))
  })
  
  weights <- weights / rowSums(weights)
  return(weights)
}

# M-Step
m_step <- function(data, weights) {
  Nk <- colSums(weights)
  alpha_new <- Nk / sum(Nk)
  mu_new <- colSums(weights * data) / Nk
  
  # common variance
  sigma_new <- sum(sapply(1:3, function(k) sum(weights[,k] * (data - mu_new[k])^2)) / sum(Nk))
  
  return(list(alpha = alpha_new, mu = mu_new, sigma = sigma_new))
}

# Compute log-likelihood
log_likelihood <- function(data, params, alpha) {
  mu <- params$mu
  sigma <- params$sigma
  
  ll <- sum(log(rowSums(sapply(1:3, function(k) {
    alpha[k] * dnorm(data, mean = mu[k], sd = sqrt(sigma))
  }))))
  
  return(ll)
}

# Estimate MAF-q from alpha
estimate_maf <- function(alpha){
  return((alpha[2] + 2*alpha[3]) / 2)
}

# EM algorithm main function
em_algorithm <- function(data, maf_init = 0.4, tol = 1e-4, max_iter = 1000) {
  alpha <- hardy_weinberg_freq(maf_init)
  mu    <- quantile(data, probs = c(0.25, 0.5, 0.75))
  sigma <- var(data)
  
  params <- list(mu = mu, sigma = sigma)
  loglik_old <- -Inf
  loglik_hist <- numeric(max_iter)  
  
  for (iter in 1:max_iter) {
    weights <- e_step(data, params, alpha)
    updated <- m_step(data, weights)
    params$mu <- updated$mu
    params$sigma <- updated$sigma
    alpha <- updated$alpha
    
    loglik_new <- log_likelihood(data, params, alpha)
    loglik_hist[iter] <- loglik_new
    
    if (abs(loglik_new - loglik_old) < tol) {
      loglik_hist <- loglik_hist[1:iter] 
      break
    }
    loglik_old <- loglik_new
  }
  
  maf_estimate <- estimate_maf(alpha)
  
  return(list(
    mu_estimates = params$mu,
    sigma_estimate = params$sigma,
    alpha_estimates = alpha,
    maf_estimate = maf_estimate,
    loglik = loglik_hist,
    iterations = iter
  ))
}

