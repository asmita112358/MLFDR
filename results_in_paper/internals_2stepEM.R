##composite alt debugging
#source("~/Library/CloudStorage/OneDrive-JohnsHopkins/MLFDR_Plos/codes/MLFDR/2step_EM_revised.R", echo = FALSE)
#source("~/Library/CloudStorage/OneDrive-JohnsHopkins/MLFDR_Plos/codes/MLFDR/MLFDR.R", echo = FALSE)
library(NMOF)

LL.data.1 = function(lambda, mu, var, x)
{
  k = length(lambda)
  m = length(x)
  t = matrix(nrow = m, ncol = k)
  
  # for(i in 1:m)
  # {
  #   for(j in 1:k)
  #   {
  #     t[i,j] = lambda[j]*dnorm(x[i], mu[j], sqrt(var[i,j]))
  #   }
  # }
  
  for(j in 1:k){
    t[,j] = lambda[j] * dnorm(x, mu[j], sqrt(var[,j]))
  }
  return(sum(log(rowSums(t))))
}


LL.complete.1 = function(kappa, lambda, mu, var_coeff, x, z)
{
  k = length(lambda)
  m = length(var_coeff)
  t <- matrix(nrow = m, ncol = k)
  kappa = c(0, kappa)
  # for(i in 1:m)
  # {
  #   for(j in 1:k)
  #   {
  #     var_mat[i,j] = var_coeff[i] + kappa[j]
  #     t[i,j] = max(dnorm(x[i], mu[j], sqrt(var_mat[i,j])), 9e-324)
  #   }
  # }
  
  for(j in 1:k){
    var_j <- var_coeff + kappa[j]
    t[,j] <- dnorm(x, mu[j], sqrt(var_j))
  }
  t[t==0] <- 9e-324
  return(-sum(z*log(t)))
}


# objective_function <- function(kappa_values, lambda, mu, var_coeff, x, z) {
#   return(-LL.complete.1(kappa_values, lambda, mu, var_coeff, x, z))  # Negate for maximization
# }

EM_fun.1 <- function(coeff, var_coeff, k, epsilon = 1e-02, maxit = 10000, lambda.init = NULL, mu.init = NULL, method = "multicore", mc.cores = detectCores()-1){
  mc_settings <- list(
    mc.cores = detectCores() - 1,        # Use 4 cores
    mc.set.seed = TRUE,  # Set seed for reproducibility
    mc.preschedule = FALSE # Do not preschedule tasks
  )
  if(is.null(lambda.init)){
    lambda.init = runif(k)
    lambda.init = lambda.init/sum(lambda.init)
  }
  lambda = lambda.init
  if(k != length(lambda.init)) message("length of lambda.init is different from k, k is assigned as length(lambda.init)")
  k = length(lambda)
  m = length(coeff)
  
  ##Initialize mean vector
  if(is.null(mu.init)){
    probs = seq(0.1,0.98,length.out = k-1)
    mu.init = c(0, quantile(coeff, probs))
  }
  mu =  mu.init 
  if(k!= length(mu))message("length of mu.init is different from k")
  kappa = c(0, rep(1,k-1))
  var_mat <- matrix(nrow = m, ncol = k)
  
  for(j in 1:k){
    var_mat[,j] <- var_coeff + kappa[j]
  }
  
  diff = 2
  iter = 0
  
  ll <- LL.data.1(lambda, mu, var_mat, coeff)
  w = matrix(nrow = m, ncol = k)
  
  
  while(diff > epsilon & iter < maxit){
    
    z = matrix(nrow = m, ncol = k)
    for(j in 1:k){
      z[,j] <- lambda[j] * dnorm(coeff, mu[j], sqrt(var_mat[,j]))
    }
    z = z/rowSums(z)
    
    #Update probabilities of each cluster
    lambda.new <- colMeans(z)
    
    #Update mu
    mu.new = c()
    for(j in 1:k){
      w[,j] = z[,j]/(var_coeff + kappa[j])
      mu.new[j] = sum(coeff*w[,j])/sum(w[,j])
    }
    
    mu.new[1] = 0
    #Update variances
    lower_bounds <- c(rep(0.01, k - 1))  # First kappa is fixed at 0
    upper_bounds <- c(rep(10, k - 1))  # Adjust upper limit as needed
    
   grid_results <- tryCatch(
        gridSearch(
          fun = LL.complete.1,
          lambda = lambda.new,
          mu = mu.new,
          var_coeff = var_coeff,
          x = coeff,
          z = z,
          lower = lower_bounds,
          upper = upper_bounds,
          method = "multicore",
          mc.control = list(mc.silent = TRUE, mc.cores = mc.cores)
          ),
        error = function(e) NULL
      )
      
      kappa.new <- if (is.null(grid_results)) {
        kappa
      } else {
        c(0, grid_results$minlevels)
      }
    
       
    #Update variance matrix
    
    var_mat.new <- matrix(nrow = m, ncol = k)
    for(j in 1:k){
      var_mat.new[,j] = var_coeff + kappa.new[j]
    }
    
    #Reassign all parameters
    lambda <- lambda.new
    mu <- mu.new
    kappa <- kappa.new
    var_mat <- var_mat.new
    newobsloglik <- LL.data.1(lambda, mu, var_mat, coeff)
    
    diff <- newobsloglik - ll
    ll <- newobsloglik
    iter <- iter + 1
    if(diff < 0){
      cat("WARNING! log-likelihood has decreased!", "\n")
    }
  }
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  cat("number of iterations=", iter, "\n")
  a = list(coeff = coeff, lambda = lambda, mu = mu, var_mat = var_mat, 
           loglik = newobsloglik, posterior = z)
  #class(a) = "mixEM"
  a
}


pi.est.comp <- function(alpha, beta, mu, theta, var_mat.alpha, var_mat.beta){
  
  m = length(alpha)
  d1 <- length(mu) - 1
  d2 <- length(theta) - 1
  k <- (1 + d1) * (1 + d2)
  
  indices = expand.grid(v = 0:d2, u = 0:d1)[2:1]
  
  z = matrix(nrow = m, ncol = k)
  pi.init = rep(1,k)
  pi.new = runif(k)
  pi.new = pi.new/sum(pi.new)
  j = 1
  while(sum((pi.init - pi.new)^2) > 1e-6){
    pi.init <- pi.new
    for(u in 1:(d1+1)){
      for(v in 1:(d2+1)){
        #print(c(u,v,j))
        z[,j] <- pi.init[j] * dnorm(alpha, mu[u], sqrt(var_mat.alpha[,u])) * dnorm(beta, theta[v], sqrt(var_mat.beta[,v]))
        j = j+1
      }
    }
    j = 1
    z <- z/rowSums(z)
    pi.new <- apply(z, 2, mean)
  }
  return(data.frame(indices, pi.new))
  
}

