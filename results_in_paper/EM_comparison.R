source(
  "~/Library/CloudStorage/OneDrive-JohnsHopkins/MLFDR_revision2/codes/MDACT_code/funcs.R",
  echo = FALSE
)
#source("~/Library/CloudStorage/OneDrive-JohnsHopkins/MLFDR_revision2/codes/MLFDR.R", echo = FALSE)
#libraries
library(HDMT)
library(locfdr)
library(qvalue)
library(emdbook)
library(parallel)
library(MLFDR)




sim.size = function(tau, p, q, size = 0.05)
{
  X = rbinom(n, 1, 0.1)
  Z = rnorm(n, 0, sd = 1)
  M = matrix(nrow = m, ncol = n)
  Y = matrix(nrow = m, ncol = n)
  #Y = vector()
  #gamma1 = sample(1:2, m, replace = T, prob = p)
  #gamma2 = sample(1:2, m, replace = T, prob = q)
  gamma = sample(1:4, m, replace = T, prob = c(0.4, 0.2, 0.2, 0.2))
  #del = rnorm(m, 0, 0.5)
  alpha = vector()
  beta = vector()
  tn = vector()
  tp = vector()
  vec1 = rnorm(m, 0.05*tau, kap)
  vec2 = rnorm(m, -0.5*tau, psi)
  
  #alpha = 0 + (gamma1 == 2)*vec1
  #beta = 0 + (gamma2 == 2)*vec2
  alpha = (gamma == 2 | gamma == 4) * vec1
  beta = (gamma == 3 | gamma == 4) * vec2
  
  tn = alpha*beta == 0
  tp = alpha*beta != 0
 
  for (i in 1:m)
  {
    
    #w/o confounders
    M[i, ] = alpha[i] * X + rnorm(n)
    Y[i, ] = beta[i] * M[i, ]  + rnorm(1, 0.5) * X + rnorm(n)
    
    # with confounder
    #M[i,] = alpha[i]*X + runif(1,0,0.5)*Z + rnorm(n)
    #Y[i,] = beta[i]*M[i,]  + runif(1,0,0.5)*Z + rnorm(1,0.5)*X + rnorm(n)
    
    
  }
  
  
  ##Estimate coefficients
  alpha_hat = vector()
  beta_hat = vector()
  var_alpha = c()
  var_beta = c()
  p1 = vector()
  p2 = vector()
  for (i in 1:m) {
    #w/o confounders
    obj1 = lm(M[i, ] ~  X)
    obj2 = lm(Y[i, ] ~  M[i, ] + X)
    
    #with confounders
    # obj1 = lm(M[i,] ~ -1 + X + Z)
    # obj2 = lm(Y[i,] ~ -1 + M[i,] + Z + X)
    
    table1 = coef(summary(obj1))
    table2 = coef(summary(obj2))
    
    
    alpha_hat[i] = table1["X", 1]
    beta_hat[i] = table2["M[i, ]", 1]
    p1[i] = table1["X", 4]
    p2[i] = table2["M[i, ]", 4]
    var_alpha[i] = table1["X", 2]^2
    var_beta[i] = table2["M[i, ]", 2]^2
  }
  
  p1[p1 == 0] <- 1e-17
  p2[p2 == 0] <- 1e-17
  input_pvalues = cbind(p1, p2)
  pmax = apply(input_pvalues, 1, max)
  
  ##null estimation
  
  nullprop = null_estimation(input_pvalues)
  
  #hdmt
  res_hdmt <- tryCatch({
    fdr_hdmt = HDMT::fdr_est(
      nullprop$alpha00,
      nullprop$alpha01,
      nullprop$alpha10,
      nullprop$alpha1,
      nullprop$alpha2,
      input_pvalues
    )
    threshhold = max(pmax[fdr_hdmt <= size])
    rej_hdmt = pmax <= threshhold
    fdr_hdmt = sum(rej_hdmt * tn) / max(1, sum(rej_hdmt))
    pow_hdmt = sum(rej_hdmt * tp) / sum(tp)
    list(fdr_hdmt = fdr_hdmt, pow_hdmt = pow_hdmt)
  }, error = function(e) {
    fdr_hdmt = NA
    pow_hdmt = NA
    list(fdr_hdmt = fdr_hdmt, pow_hdmt = pow_hdmt)
  })
  
  fdr_hdmt <- res_hdmt$fdr_hdmt
  pow_hdmt <- res_hdmt$pow_hdmt
  
  
  #mdact
  res <- tryCatch({
    p.mdact <- balancing_DACT_control_DR_adjust(p1,
                                                p2,
                                                nullprop,
                                                significance_upper = 0.05,
                                                control.method = "FDR")
    rej_mdact <- rep(FALSE, m)
    rej_mdact[p.mdact] <- TRUE
    fdr_mdact <- sum(rej_mdact * tn) / max(1, sum(rej_mdact))
    pow_mdact <- sum(rej_mdact * tp) / sum(tp)
    list(fdr_mdact = fdr_mdact, pow_mdact = pow_mdact)
  }, error = function(e) {
    list(fdr_mdact = NA, pow_mdact = NA)
  })
  
  fdr_mdact <- res$fdr_mdact
  pow_mdact <- res$pow_mdact
  
  #MLFDR
  
  lfdr <- localFDR(alpha_hat, beta_hat, var_alpha, var_beta, twostep = FALSE)
  
  rej_MLFDR <- MLFDR(lfdr, size = size)
  fdr_mlfdr <- sum(rej_MLFDR * tn) / max(1, sum(rej_MLFDR))
  pow_mlfdr <- sum(rej_MLFDR * tp) / sum(tp)
  
  #MLFDR2
  
  lfdr2 <- localFDR(alpha_hat, beta_hat, var_alpha, var_beta, twostep = TRUE, d1 = 1, d2 = 1)
  rej_MLFDR2 <- MLFDR(lfdr2, size = size)
  fdr_mlfdr2 <- sum(rej_MLFDR2 * tn) / max(1, sum(rej_MLFDR2))
  pow_mlfdr2 <- sum(rej_MLFDR2 * tp) / sum(tp)
  
  
  out = c(
    fdr_mlfdr2,
    fdr_mlfdr,
    fdr_hdmt,
    fdr_mdact,
    pow_mlfdr2,
    pow_mlfdr,
    pow_hdmt,
    pow_mdact
  )
  names(out) <- c(
    "fdr_mlfdr2",
    "fdr_mlfdr",
    "fdr_hdmt",
    "fdr_mdact",
    "pow_mlfdr2",
    "pow_mlfdr",
    "pow_hdmt",
    "pow_mdact"
  )
  return(out)
}



kap = 1
psi = 2
tau = seq(0.1, 1.9, by = 0.2)

##simulation
m = 1000
n = 300
n.sim = 250
#pi = c(0.4,0.28, 0.3, 0.02)
sim.res = matrix(nrow = length(tau), ncol = 8)
p = c(0.6, 0.4)
q = c(0.7, 0.3)
#pi = c(0.88, 0.05, 0.05, 0.02)
sim.se = matrix(nrow = length(tau), ncol = 8)
colnames(sim.res) <- c("fdr_mlfdr2",
                       "fdr_mlfdr",
                       "fdr_hdmt",
                       "fdr_mdact",
                       "pow_mlfdr2",
                       "pow_mlfdr",
                       "pow_hdmt",
                       "pow_mdact")
colnames(sim.se) <- c("fdr_mlfdr2",
                      "fdr_mlfdr",
                      "fdr_hdmt",
                      "fdr_mdact",
                      "pow_mlfdr2",
                      "pow_mlfdr",
                      "pow_hdmt",
                      "pow_mdact")
for (k1 in 1:length(tau)) {
  # Parallelize the inner loop with mclapply
  temp <- mclapply(
    1:n.sim,
    function(s) {
      sim.size(tau[k1], p, q)
    },
    mc.cores = detectCores() - 1  # Use available cores minus one
  ) |>
    do.call(what = rbind)
  
  sim.res[k1,]  <- colMeans(temp, na.rm = TRUE)
  sim.se[k1,] <- apply(temp, 2, sd, na.rm = TRUE)
  print(k1)
}
write.csv(sim.res, paste0("results/mean_m",m,"n",n,"EM_comp.csv"))
write.csv(sim.se, paste0("results/se_m",m,"n",n,"EM_comp.csv"))
sim.res



beepr::beep(4)
