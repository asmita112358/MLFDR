source("funcs.R")
options(warn = -1)

Parallelsim_mu_fdr <- function(ws,pras,sim.num,control.method,sig_level,num_runs){
  library(parallel)
  
  SizePowerSim = function(ws,pras,sim.num,control.method,sig_level){
    source("funcs.R")
    #10 01 00 11
    combinations <- as.matrix(expand.grid(rep(list(c(0, 1)),2)))[c(2,3,1,4),]
    # Calculate the number of units to be assigned to each set based on probabilities
    sim.num_per_set <- round(ws * sim.num)
    sets <- vector("list", length = length(ws))
    current_index <- 1
    for (i in 1:length(ws)) {
      if (sim.num_per_set[i] > 0) {
        sets[[i]] <- current_index:(current_index + sim.num_per_set[i] - 1)
        current_index <- current_index + sim.num_per_set[i]
      } else {
        sets[[i]] <- integer(0)  # Empty set
      }
    }
    
    true_set <- unlist(sets[4])
    
    
    Z.M = c()
    Z.Y = c()
    p.M = c()
    p.Y = c()
    
    ## directly generate z scores
    # m <- 6
    # delta = 1
    # pras <- c(m/sqrt(1+delta^2),m/sqrt(1+delta^2)*delta)
    # 
    mu_M <- sample(c(1,-1),length(c(unlist(sets[1]),unlist(sets[4]))),replace = T)*rnorm(length(c(unlist(sets[1]),unlist(sets[4]))),pras[1],1)
    mu_Y <- sample(c(1,-1),length(c(unlist(sets[2]),unlist(sets[4]))),replace = T)*rnorm(length(c(unlist(sets[2]),unlist(sets[4]))),pras[2],1)
    Z.M <- rnorm(sim.num,0,1); Z.M[c(unlist(sets[1]),unlist(sets[4]))] <- rnorm(length(c(unlist(sets[1]),unlist(sets[4]))),mu_M,1);
    Z.Y <- rnorm(sim.num,0,1); Z.Y[c(unlist(sets[2]),unlist(sets[4]))] <- rnorm(length(c(unlist(sets[2]),unlist(sets[4]))),mu_Y,1);
    p.M <- 2* (1 - pnorm(abs(Z.M))); p.Y <- 2* (1 - pnorm(abs(Z.Y)))
    p.M[p.M==0] <- 1e-17 ; p.Y[p.Y==0] <- 1e-17
    
    estws <- null_estimation(cbind(p.M,p.Y))
    
    if(control.method == "FDR"){
      result_tab <- matrix(NA,4,2)
    }
    if(control.method == "FWER"){
      result_tab <- matrix(NA,4,1)
    }
    
    p.maxp <- pmax(p.M,p.Y)
    p.MTComp = CompTestER(Z.M,Z.Y)$pp
    if(control.method == "FDR"){
      result_tab[1,] <- index_fpn_fdr_power_func(which(p.adjust(p.maxp,method = "fdr")< sig_level),true_set)
      result_tab[2,] <- index_fpn_fdr_power_func(which(p.adjust(p.MTComp,method = "fdr")< sig_level),true_set)
    }
    if(control.method == "FWER"){
      result_tab[1,] <- index_fpn_fdr_power_func(which(p.adjust(p.maxp,method = "bonferroni")< sig_level),true_set)
      result_tab[2,] <- index_fpn_fdr_power_func(which(p.adjust(p.MTComp,method = "bonferroni")< sig_level),true_set)
    }
    
    p.mdact <- balancing_DACT_control_DR_adjust(p.M,p.Y,estws,significance_upper = sig_level,control.method=control.method)
    result_tab[3,] <- index_fpn_fdr_power_func(p.mdact,true_set)
    
    p.hdmt <- HDMT_control_DR_adjust(p.M,p.Y,estws, significance_upper = sig_level,control.method=control.method)
    result_tab[4,] <- index_fpn_fdr_power_func(p.hdmt,true_set)
    
    rownames(result_tab) <- c("MaxP","MT-comp","MDACT","HDMT")
    
    return(list(result_tab=result_tab))
  }
  
  sim_results <- mclapply(1:num_runs, function(i) {
    SizePowerSim(ws, pras, sim.num, control.method, sig_level)
  }, mc.cores = 20)
  
  sim_results_table <- list()
  
  for (i in 1:num_runs) {
    sim_results_table[[i]] <- sim_results[[i]]$result_tab
  }
  
  average_tab <- as.matrix(apply(array(sapply(sim_results_table, function(x) x), 
                                       dim = c(4, dim(sim_results[[i]]$result_tab)[2], num_runs)), 1:dim(sim_results[[i]]$result_tab)[2], mean))
  rownames(average_tab) <- c("MaxP","MT-comp","MDACT","HDMT")
  if(control.method == "FWER"){
    colnames(average_tab) <- c("FDP")
  }else{
    if(ws[4]==0){
      colnames(average_tab) <- c("FDP","FDP")
    }else{
      colnames(average_tab) <- c("FDP","power")
    }
  }
  
  res <- list(result_tab=average_tab)
  return(res)
}

### control the FDR is the same as control FWER under \pi_11=0

m <- 4
delta <- 1

pras <- c(m/sqrt(1+delta^2),m/sqrt(1+delta^2)*delta)##mean of the z-scores distribution
res_fdr <- Parallelsim_mu_fdr(c(0.2,0.2,0.6,0),pras=pras,sim.num=4.8e5,control.method = "FDR", sig_level=0.1,num_runs=500)
res_fdr[["result_tab"]]


pras <- c(m/sqrt(1+delta^2),m/sqrt(1+delta^2)*delta)##mean of the z-scores distribution
res_fdr <- Parallelsim_mu_fdr(c(0.2,0.2,0.6,0),pras=pras,sim.num=4.8e5,control.method = "FWER", sig_level=0.1,num_runs=1000)
res_fdr[["result_tab"]]

