#' Internal function for MLFDR, mediation analysis with localFDR
#' @importFrom stats dnorm optimize quantile runif
#' @importFrom NMOF gridSearch
#' @importFrom parallel detectCores mclapply
#' @noRd
LL.data = function(coeff_mat, mu, sigma, lambda)
{
  k = length(mu)
  m = nrow(coeff_mat)
  t = matrix(nrow = m, ncol = k)

  for(j in 1:k){
    t[,j] = lambda[j]* dnorm(coeff_mat[,1], mu[[j]][1], sqrt(sigma[j,,1,1])) * dnorm(coeff_mat[,2], mu[[j]][2],sqrt(sigma[j,,2,2]))
  }
  return(sum(log(rowSums(t))))
}

LL.complete.v2 <- function(kappa, psi, var_alpha, var_beta, coeff_mat, mu.new, z){
  k = 4
  m = length(var_alpha)
  t = matrix(nrow = m, ncol = 4)
  #Define sigma

  sigma <- array(0,dim = c(k, m, 2,2))
  for(i in 1:3){
    sigma[i, ,1,1] <- var_alpha
    sigma[i, ,2,2] <- var_beta
  }
  sigma[2, ,1,1] <- var_alpha + kappa
  sigma[3, ,2,2] <- var_beta + psi
  sigma[4, ,1,1] <- var_alpha + kappa
  sigma[4, ,2,2] <- var_beta + psi

  for(i in 1:k){
    temp1= dnorm(coeff_mat[,1], mu.new[[i]][1], sqrt(sigma[i,,1,1]))*dnorm(coeff_mat[,2], mu.new[[i]][2], sqrt(sigma[i,,2,2]))
    t[,i] = pmax(temp1, rep(9e-321,m))
  }
  return(sum(z*log(t)))

}


EM_fun <- function(coeff_mat, k = 4, var_alpha, var_beta ,lambda.init = NULL, kappa.init = NULL, psi.init = NULL,
                   kappa_int = NULL, psi_int = NULL,
                   epsilon = 1e-02, maxit = 10000, verbose = FALSE)
{
  if(is.null(lambda.init)){
   lambda = c(0.7, 0.1, 0.1, 0.1)
  }else{
    lambda = lambda.init
  }
  coeff_mat <- as.matrix(coeff_mat)
  m <- nrow(coeff_mat)
  p <- ncol(coeff_mat)
  mu.init = quantile(coeff_mat[,1], 0.99)
  theta.init = quantile(coeff_mat[,2], 0.99)
  if(is.null(kappa.init)){
    kappa = 1
  }else{
    kappa = kappa.init
  }

  if(is.null(psi.init)){
    psi = 1
    }else{
      psi = psi.init
    }
  if(is.null(kappa_int)){
    kappa_int = c(min(0.1, min(var_alpha)), max(10, max(var_alpha)))
  }
  if(is.null(psi_int)){
    psi_int = c(min(0.1, min(var_beta)), max(10, max(var_beta)))
  }

  mu = list(c(0,0), c(mu.init, 0), c(0, theta.init), c(mu.init, theta.init))
  sigma <- array(0,dim = c(k, m, 2,2))

  for(i in 1:3){
    sigma[i, ,1,1] <- var_alpha
    sigma[i, ,2,2] <- var_beta
  }

  sigma[2, ,1,1] <- var_alpha + kappa
  sigma[3, ,2,2] <- var_beta + psi
  sigma[4, ,1,1] <- var_alpha + kappa
  sigma[4, ,2,2] <- var_beta + psi



  diff <- 1
  iter <- 0


  ll <- LL.data(coeff_mat, mu, sigma, lambda)
  restarts <- 0
  while (diff > epsilon & iter < maxit) {

    ##Compute Q
    z = matrix(nrow = m, ncol = k)
    for(j in 1:k){
      z[,j] = lambda[j]*dnorm(coeff_mat[,1], mu[[j]][1], sqrt(sigma[j,,1,1]))*dnorm(coeff_mat[,2], mu[[j]][2],sqrt(sigma[j,,2,2]))
    }
    z = z/rowSums(z)

    lambda.new <- apply(z, 2, mean)
    w = (z[,2] + z[,4])/(var_alpha + kappa)
    v = (z[,3] + z[,4])/(var_beta + psi)


    m.new = sum(coeff_mat[,1]*w)/sum(w)
    theta.new = sum(coeff_mat[,2]*v)/sum(v)
    mu.new <- list(c(0,0), c(m.new, 0), c(0, theta.new), c(m.new, theta.new))

    #Update kappa and psi

    kappa.new = optimize(LL.complete.v2, interval = kappa_int, psi = psi, var_alpha = var_alpha, var_beta = var_beta, coeff_mat = coeff_mat, mu.new = mu.new,z = z, maximum = TRUE)$maximum
    psi.new = optimize(LL.complete.v2, interval = psi_int, kappa = kappa.new, var_alpha = var_alpha, var_beta = var_beta, coeff_mat = coeff_mat, mu.new = mu.new,z = z, maximum = TRUE)$maximum

    ##Update sigma

    sigma.new <- array(0,dim = c(k, m, 2,2))
    for(i in 1:k)
    {
      for(j in 1:m)
      {
        sigma.new[i,j,1,1] = var_alpha[j] + ifelse(i == 2||i == 4, 1, 0)*kappa.new
        sigma.new[i,j,2,2] = var_beta[j] + ifelse(i == 3||i == 4, 1, 0)*psi.new
      }

    }


    #Reassign all parameters

    lambda <- lambda.new
    mu <- mu.new
    kappa <- kappa.new
    psi <- psi.new
    sigma <- sigma.new
    newobsloglik <- LL.data(coeff_mat, mu, sigma, lambda)
    if(verbose){
      cat("iteration=", iter, "loglik=", newobsloglik, "\n")
    }
    if(newobsloglik == -Inf){
      diff = 5
    }else{
      diff = newobsloglik - ll
      ll <- newobsloglik
    }
    iter <- iter +1


  }

  if (iter == maxit) {
    warning("WARNING! NOT CONVERGENT!", "\n")
  }

  if(verbose) cat("number of iterations=", iter, "\n")
  a = list(coeff_mat = coeff_mat, lambda = lambda, mu = mu, sigma = sigma,
           loglik = newobsloglik, posterior = z, all.loglik = ll, restarts = restarts)
  #class(a) = "mixEM"
  a
}
LL.data.1 = function(lambda, mu, var, x)
{
  k = length(lambda)
  m = length(x)
  t = matrix(nrow = m, ncol = k)

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

  for(j in 1:k){
    var_j <- var_coeff + kappa[j]
    t[,j] <- dnorm(x, mu[j], sqrt(var_j))
  }
  t[t==0] <- 9e-324
  return(-sum(z*log(t)))
}


EM_fun.1 <- function(coeff, var_coeff, k, epsilon = 1e-02, maxit = 10000, lambda.init = NULL, mu.init = NULL, method = "multicore", mc.cores = detectCores()-1, verbose = FALSE){
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
      NMOF::gridSearch(
        fun = LL.complete.1,
        lambda = lambda.new,
        mu = mu.new,
        var_coeff = var_coeff,
        x = coeff,
        z = z,
        lower = lower_bounds,
        upper = upper_bounds,
        printDetail = verbose,
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
      warning("WARNING! log-likelihood has decreased!", "\n")
    }
  }
  if (iter == maxit) {
    warning("WARNING! NOT CONVERGENT!", "\n")
  }
  if(verbose) cat("number of iterations=", iter, "\n")
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

