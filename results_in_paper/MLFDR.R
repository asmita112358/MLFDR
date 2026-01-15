source("~/Library/CloudStorage/OneDrive-JohnsHopkins/MLFDR_Plos/codes/MLFDR/internals.R", echo = FALSE)
#source("~/Library/CloudStorage/OneDrive-JohnsHopkins/MLFDR_Plos/codes/MLFDR/2step_EM_revised.R", echo = FALSE)
source("~/Library/CloudStorage/OneDrive-JohnsHopkins/MLFDR_Plos/codes/MLFDR/internals_2stepEM.R", echo = FALSE)
localFDR <- function(alpha, beta, var_alpha, var_beta, lambda.init = c(0.7, 0.1, 0.1, 0.1),
                  kappa.init = 1, psi.init = 1, psi_int = NULL, kappa_int = NULL,
                  twostep = FALSE, 
                  k = 4, d1 = NULL, d2 = NULL, eps = 1e-02,
                  verbose = TRUE){
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
   
  }else{
    if(is.null(d1)) {stop("Please specify d1 for two-step MLFDR")}
    if(is.null(d2)) {stop("Please specify d2 for two-step MLFDR")}
    m = length(alpha)
    fit_alpha = EM_fun.1(alpha, var_alpha, k = d1 + 1, epsilon = eps, maxit = 10000)
    fit_beta = EM_fun.1(beta, var_beta, k = d2 + 1, epsilon = eps, maxit = 10000)
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
  return(lfdr = lfdr)
  }
MLFDR <- function(lfdr, size = size){
  m = length(lfdr)
  st.lfdr<-sort(lfdr)
  k=1
  while(k<m && ((1/k)*sum(st.lfdr[1:k])) <= size){
    k=k+1
  }
  k<-k-1
  lfdrk<-st.lfdr[k]
  reject<- lfdr<=lfdrk
  return(reject)
}
