##SVA
##DO NOT RUN LOCALLY WITHOUT CHANGING DETECTCORES
library(sva)
library(limma)
library(bladderbatch)
library(MLFDR)
library(HDMT)
library(parallel)
#setwd('~/mlfdr')
#source('helper_mlfdr.R')
#source("~/Library/CloudStorage/OneDrive-JohnsHopkins/MLFDR_revision2/codes/MDACT_code/funcs.R", echo = FALSE)

# data(bladderdata)
# pheno = pData(bladderEset)
# edata = exprs(bladderEset)
# mod = model.matrix(~as.factor(cancer), data=pheno)
# mod0 = model.matrix(~1,data=pheno)
# n.sv = num.sv(edata,mod,method="leek")
# n.sv
# svobj = sva(edata,mod,mod0,n.sv=n.sv)
# pValues = f.pvalue(edata,mod,mod0)
# qValues = p.adjust(pValues,method="BH")
# modSv = cbind(mod,svobj$sv)
# mod0Sv = cbind(mod0,svobj$sv)
# pValuesSv = f.pvalue(edata,modSv,mod0Sv)
# qValuesSv = p.adjust(pValuesSv,method="BH")


uni.mclapply <- function(X, FUN, ...,
                         mc.cores = max(parallel::detectCores()-2, 1),
                         force_windows = FALSE) {
  os_type <- if (force_windows) "windows" else .Platform$OS.type

  if (os_type == "unix") {
    out <- tryCatch(
      parallel::mclapply(X, FUN, ..., mc.cores = mc.cores),
      error = function(e) {
        message("mclapply failed, falling back to lapply: ", e$message)
        lapply(X, FUN, ...)
      }
    )
  } else {
    out <- tryCatch({
      future::plan(future::multisession, workers = mc.cores)
      future.apply::future_lapply(X, FUN, ..., future.seed = TRUE)
    }, error = function(e) {
      message("future_lapply failed, falling back to lapply: ", e$message)
      lapply(X, FUN, ...)
    })
  }

  return(out)
}

sim.size <- function(tau,  pi, size = 0.05)
{
X = rbinom(n, 1, 0.1)
M = matrix(nrow = m, ncol = n)
Y = matrix(nrow = m, ncol = n)
Z = rnorm(n, 2, sd = 1.2) #observed confounder
Z1 <- rnorm(n) #unobserved confounder 1
Z2 <- runif(n) #unobserved confounder 2
cli <- data.frame(X, Z)
theta = rnorm(m, 2, 0.25)
gamma = sample(1:4, m, replace = T, prob = pi)
alpha = vector()
beta = vector()

vec1 = rnorm(m, 0.05*tau, 1)
vec2 = rnorm(m, -0.5*tau, 2)

alpha = ((gamma==2) + (gamma == 4))*vec1
beta = ((gamma ==3) + (gamma == 4))*vec2

tn = alpha*beta == 0
tp = alpha*beta != 0


for(i in 1:m)
{
  M[i,] = alpha[i]*X + runif(1)*Z + 0.4*Z1 - 0.5*Z2 + rnorm(n)
  
}
indices = sample(1:m, 20)
theta_ind = theta[indices]
M_ind = scale(as.vector(theta_ind %*% M[indices,]))
for(i in 1:m)
{
  Y[i,] = beta[i]*M[i,]  + rnorm(1,0.5)*X +runif(1)*Z - 0.5*Z1 + rnorm(1,0.5)*M_ind  + rnorm(n)  
}


mod = model.matrix(~ X + Z, data = cli)
mod0 = model.matrix(~ Z, data = cli)
n.sv1 = num.sv(M,mod,method="be")
svobj = sva(M, mod, mod0, n.sv = n.sv1)
modSv = cbind(mod,svobj$sv)

  
n.sv2 = num.sv(Y,mod,method="be")
svobj2 = sva(Y, mod, mod0, n.sv = n.sv2)
modSv2 = cbind(mod, svobj2$sv)

alpha_hat = vector()
beta_hat = vector()
var_alpha = c()
var_beta = c()
p1 = vector()
p2 = vector()


for(i in 1:m){
  obj1 = lm(M[i,] ~  -1 + modSv)
  #mod = model.matrix(~ M[i,] + X + Z)
  #mod0 = model.matrix(~ X + Z)
  #n.sv2[i] = num.sv(Y,mod,method="be")
  #svobj2 = sva(Y, mod, mod0, n.sv = n.sv2[i])
  #modSv2 = cbind(mod, svobj2$sv)
  obj2 <- lm(Y[i,] ~ -1 + M[i,] + modSv2)
  table1 = coef(summary(obj1))
  table2 = coef(summary(obj2))
  alpha_hat[i] = table1["modSvX",1]
  beta_hat[i] = table2["M[i, ]",1]
  p1[i] = table1["modSvX",4]
  p2[i] = table2["M[i, ]",4]
  var_alpha[i] = table1["modSvX",2]^2
  var_beta[i] = table2["M[i, ]",2]^2
}

#Change the above to mclapply
# results <- uni.mclapply(1:m, function(i) {
#   obj1 = lm(M[i,] ~  -1 + modSv)
#   mod = model.matrix(~ M[i,] + X + Z)
#   mod0 = model.matrix(~ X + Z)
#   n.sv2 = num.sv(Y,mod,method="be")
#   svobj2 = sva(Y, mod, mod0, n.sv = n.sv2)
#   modSv2 = cbind(mod, svobj2$sv)
#   obj2 <- lm(Y[i,] ~ -1 + modSv2)
#   table1 = coef(summary(obj1))
#   table2 = coef(summary(obj2))
#   list(
#     alpha_hat = table1["modSvX",1],
#     beta_hat = table2["modSv2M[i, ]",1],
#     p1 = table1["modSvX",4],
#     p2 = table2["modSv2M[i, ]",4],
#     var_alpha = table1["modSvX",2]^2,
#     var_beta = table2["modSv2M[i, ]",2]^2)
# }, mc.cores = 12, force_windows = TRUE)
# 
# # Extract results from the list
# for (i in 1:m) {
#   alpha_hat[i] <- results[[i]]$alpha_hat
#   beta_hat[i] <- results[[i]]$beta_hat
#   p1[i] <- results[[i]]$p1
#   p2[i] <- results[[i]]$p2
#   var_alpha[i] <- results[[i]]$var_alpha
#   var_beta[i] <- results[[i]]$var_beta
# }

  
p1[p1==0] <- 1e-17 ; p2[p2==0] <- 1e-17
input_pvalues = cbind(p1, p2)
pmax = apply(input_pvalues, 1, max)

##null estimation

nullprop = null_estimation(input_pvalues)

#hdmt
res_hdmt <- tryCatch({
  fdr_hdmt = HDMT::fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,
                           nullprop$alpha1,nullprop$alpha2,input_pvalues)
  threshhold = max(pmax[fdr_hdmt<= size])
  rej_hdmt = pmax <= threshhold
  fdr_hdmt = sum(rej_hdmt*tn)/max(1,sum(rej_hdmt))
  pow_hdmt = sum(rej_hdmt*tp)/sum(tp)
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
  p.mdact <- balancing_DACT_control_DR_adjust(p1, p2, nullprop, significance_upper = size, control.method = "FDR")
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

lfdr <- localFDR(alpha_hat, beta_hat, var_alpha, var_beta, twostep = FALSE)
rej_MLFDR <- MLFDR(lfdr$lfdr, size = size)
fdr_MLFDR <- sum(rej_MLFDR*tn)/max(1,sum(rej_MLFDR))
pow_MLFDR <- sum(rej_MLFDR*tp)/sum(tp)

out <- c(fdr_MLFDR,  fdr_hdmt,fdr_mdact, pow_MLFDR, pow_hdmt, pow_mdact)
names(out) <- c("FDR_MLFDR", "FDR_HDMT", "FDR_MDACT", "Power_MLFDR", "Power_HDMT", "Power_MDACT")
return(out)
}
#kap = 1
#psi = 4
tau = seq(0.1, 1.9, by = 0.2)


##Sparse alternative
m = 1000
n = 300
n.sim = 100
#pi = c(0.4,0.28, 0.3, 0.02)
sim.res = matrix(nrow = length(tau), ncol = 6)
pi = c(0.88, 0.05, 0.05, 0.02)
sim.se = matrix(nrow = length(tau), ncol = 6)
colnames(sim.res) <- c("fdr_mlfdr", "fdr_hdmt", "fdr_mdact", "pow_mlfdr", "pow_hdmt", "pow_mdact")
colnames(sim.se) <- c("fdr_mlfdr", "fdr_hdmt", "fdr_mdact", "pow_mlfdr", "pow_hdmt", "pow_mdact")
for (k1 in 1:length(tau)) {
  temp <- uni.mclapply(1:n.sim, 
                   function(x) sim.size(tau[k1], pi, size = 0.05), 
                   mc.cores = detectCores()-1, force_windows = TRUE)|>
    do.call(what = rbind)
  
  sim.res[k1,]  <- colMeans(temp, na.rm = TRUE)
  sim.se[k1,] <- apply(temp, 2, sd, na.rm = TRUE)
  print(k1)
  print(sim.res[k1,])
  #beepr::beep()
}
write.csv(sim.res, paste0("results/mean_m",m,"n",n,"sparse_pl.csv"))
write.csv(sim.se, paste0("results/se_m",m,"n",n,"sparse_pl.csv"))
sim.res
beepr::beep(4)


system.time({u = sim.size(tau=0.1, pi = c(0.88,0.05, 0.05, 0.02), size = 0.05)})


###Dense alternative
m = 1000
n = 300
n.sim = 100
#pi = c(0.4,0.28, 0.3, 0.02)
sim.res = matrix(nrow = length(tau), ncol = 6)
pi = c(0.4,0.2, 0.2, 0.2)
sim.se = matrix(nrow = length(tau), ncol = 6)
colnames(sim.res) <- c("fdr_mlfdr", "fdr_hdmt", "fdr_mdact", "pow_mlfdr", "pow_hdmt", "pow_mdact")
colnames(sim.se) <- c("fdr_mlfdr", "fdr_hdmt", "fdr_mdact", "pow_mlfdr", "pow_hdmt", "pow_mdact")
for (k1 in 1:length(tau)) {
  temp <- uni.mclapply(1:n.sim, 
                   function(x) sim.size(tau[k1], pi, size = 0.05), 
                   mc.cores = detectCores() - 1, force_windows = TRUE)|>
    do.call(what = rbind)
  
  sim.res[k1,]  <- colMeans(temp, na.rm = TRUE)
  sim.se[k1,] <- apply(temp, 2, sd, na.rm = TRUE)
  print(k1)
  print(sim.res[k1])
  beepr::beep()
}
write.csv(sim.res, paste0("results/mean_m",m,"n",n,"dense_pl.csv"))
write.csv(sim.se, paste0("results/se_m",m,"n",n,"dense_pl.csv"))
sim.res
