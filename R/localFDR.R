#' localFDR: A function that fits a Gaussian Mixture model to the mediation coefficients and returns the localFDR estimated from the mixture model
#'
#' @param alpha a vector of estimated alpha coefficients from the first equation of mediation analysis
#' @param beta a vector of estimated beta coefficients from the second equation of mediation analysis
#' @param var_alpha a vector of estimated variances for alpha coefficients
#' @param var_beta a vector of estimated variances for beta coefficients
#' @param lambda.init initial values of the proportion of mixture, must sum to 1
#' @param kappa.init inital value of kappa, the variance of the prior of alpha under the alternative
#' @param psi.init inital value of psi, the variance of the prior of beta under the alternative
#' @param psi_int gridSearch interval for psi, to be used in optimize function
#' @param kappa_int gridSearch interval for kappa, to be used in optimize function
#' @param twostep logical, whether to use two-step MLFDR
#' @param k number of mixture components, default is 4. Used in one-step MLFDR.
#' @param d1 number of non-null components for alpha in two-step MLFDR.
#' @param d2 number of non-null components for beta in two-step MLFDR.
#' @param eps stopping criteria for EM algorithm
#' @param verbose logical, whether to print the log-likelihood at each iteration. Default is TRUE.
#' @param method either "unicore" or "multicore", specifies whether parallelization should be used when estimating variances in two step EM. Is not used unless twostep = TRUE
#' @returns A list consisting of two elements - lfdr, a vector of local false discovery rates and pi, the estimated mixture proportions
#' @export
#'
#' @examples
#' n = 100
#' m = 1000
#' pi = c(0.4, 0.1, 0.3, 0.2)
#' X = rbinom(n, 1, 0.1)
#' M = matrix(nrow = m, ncol = n)
#' Y = matrix(nrow = m, ncol = n)
#' gamma = sample(1:4, m, replace = TRUE, prob = pi)
#' alpha = vector()
#' beta = vector()

#' vec1 = rnorm(m, 0.05, 1)
#' vec2 = rnorm(m, -0.5, 2)

#' alpha = ((gamma==2) + (gamma == 4))*vec1
#' beta = ((gamma ==3) + (gamma == 4))*vec2



#' alpha_hat = vector()
#' beta_hat = vector()
#' var_alpha = c()
#' var_beta = c()
#' p1 = vector()
#' p2 = vector()
#' for(i in 1:m)
#' {
#'   M[i,] = alpha[i]*X + rnorm(n)
#'   Y[i,] = beta[i]*M[i,]  + rnorm(1,0.5)*X + rnorm(n)
#'   obj1 = lm(M[i,] ~  X )
#'   obj2 = lm(Y[i,] ~  M[i,] + X)

#'   table1 = coef(summary(obj1))
#'   table2 = coef(summary(obj2))


#'   alpha_hat[i] = table1["X",1]
#'   beta_hat[i] = table2["M[i, ]",1]
#'   p1[i] = table1["X",4]
#'   p2[i] = table2["M[i, ]",4]
#'   var_alpha[i] = table1["X",2]^2
#'   var_beta[i] = table2["M[i, ]",2]^2
#' }

#' lfdr <- localFDR(alpha_hat, beta_hat, var_alpha, var_beta, twostep = FALSE)

localFDR <- function(alpha, beta, var_alpha, var_beta, lambda.init = NULL,
                     kappa.init = 1, psi.init = 1, psi_int = NULL, kappa_int = NULL,
                     twostep = FALSE,
                     k = 4, d1 = NULL, d2 = NULL, eps = 1e-02,
                     verbose = TRUE, method = "unicore"){
  if(!twostep){
    x = cbind(alpha, beta)
    m = nrow(x)
    fit = EM_fun(x, k = 4, var_alpha, var_beta,lambda.init = lambda.init,
                 kappa.init = kappa.init, psi.init = psi.init,
                 kappa_int = kappa_int, psi_int = psi_int,
                 epsilon = eps, verbose = verbose)
    pi = fit$lambda
    mu = fit$mu
    k = length(mu)
    sigma = fit$sigma
    lfdr = vector()
    t = matrix(nrow = m, ncol = k)
    for(j in 1:k){
      t[,j] = pi[j] * dnorm(x[,1], mu[[j]][1], sqrt(sigma[j,,1,1])) * dnorm(x[,2], mu[[j]][2], sqrt(sigma[j,,2,2]))
    }
    lfdr <- (t[,1] + t[,2] + t[,3])/rowSums(t)
    pi_em <- expand.grid(H_alpha = c(0,1), H_beta = c(0,1))
    pi_em = data.frame(pi_em, prob = pi)
  }else{
    if(is.null(d1)) {stop("Please specify d1 for two-step MLFDR")}
    if(is.null(d2)) {stop("Please specify d2 for two-step MLFDR")}
    m = length(alpha)
    fit_alpha = EM_fun.1(alpha, var_alpha, k = d1 + 1, epsilon = eps, maxit = 10000, method = method, verbose = verbose)
    fit_beta = EM_fun.1(beta, var_beta, k = d2 + 1, epsilon = eps, maxit = 10000, method = method, verbose = verbose)
    p_em = fit_alpha$lambda
    q_em = fit_beta$lambda
    mu_em = fit_alpha$mu
    theta_em = fit_beta$mu
    var_mat.alpha = fit_alpha$var_mat
    var_mat.beta = fit_beta$var_mat
    pi_em = pi.est.comp(alpha, beta, mu_em, theta_em, var_mat.alpha, var_mat.beta)
    pi_est <- pi_em$pi.new
    lfdr_em = c()
    z_em = matrix(nrow = m, ncol = (1+d1)*(1 + d2))

    j = 0
    for(u in 1:(d1+1)){
      for(v in 1:(d2 + 1)){
        j = j+1
        z_em[,j] = 10* pi_est[j] * dnorm(alpha, mu_em[u], sqrt(var_mat.alpha[,u])) * dnorm(beta, theta_em[v], sqrt(var_mat.beta[,v]))
      }
    }
    h0 <- (pi_em$u)*(pi_em$v) == 0
    lfdr <- rowSums(z_em[,h0])/rowSums(z_em)
  }
  return(list(lfdr = lfdr, pi = pi_em))
}
