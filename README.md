
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MLFDR

## Overview

The package MLFDR implements a local-FDR based solution to the high
dimensional mediation analysis problem. For
$H_{0i} = \alpha_i\beta_i = 0$, $\alpha_i$ denotes the coefficient of
the mediator-exposure relationship, and $\beta_i$, that of the
mediator-outcome relationship. MLFDR takes the estimates and standard
errors of $\alpha_i$ and $\beta_i$ and fits a Gaussian Mixture model
accounting for the composite null hypothesis, and estimates a localFDR
based on this Gaussian Mixture model. An adaptive procedure determines
the cutoff for the local FDR.

## Installation

You can install the development version of MLFDR from
[GitHub](https://github.com/asmita112358/MLFDR) with:

``` r
devtools::install_github("asmita112358/MLFDR")
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo asmita112358/MLFDR@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>      checking for file ‘/private/var/folders/mk/vjmxbpp92m750qctwwm_qn440000gn/T/RtmpiSHtim/remotes4944626e5e69/asmita112358-MLFDR-37db225/DESCRIPTION’ ...  ✔  checking for file ‘/private/var/folders/mk/vjmxbpp92m750qctwwm_qn440000gn/T/RtmpiSHtim/remotes4944626e5e69/asmita112358-MLFDR-37db225/DESCRIPTION’
#>   ─  preparing ‘MLFDR’:
#>      checking DESCRIPTION meta-information ...  ✔  checking DESCRIPTION meta-information
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>    Omitted ‘LazyData’ from DESCRIPTION
#>   ─  building ‘MLFDR_0.1.0.tar.gz’
#>      
#> 
```

## Example

This is a basic example which demonstrates the implementation of MLFDR.
For more in-depth examples, please refer to the results_in_paper as well
as (<https://arxiv.org/abs/2402.13933>).

``` r
library(MLFDR)
n = 150
m = 1000
pi = c(0.4, 0.1, 0.3, 0.2)
X = rbinom(n, 1, 0.1)
M = matrix(nrow = m, ncol = n)
Y = matrix(nrow = m, ncol = n)
gamma = sample(1:4, m, replace = T, prob = pi)
alpha = vector()
beta = vector()

vec1 = rnorm(m, 0.05, 1)
vec2 = rnorm(m, -0.5, 2)

alpha = ((gamma==2) + (gamma == 4))*vec1
beta = ((gamma ==3) + (gamma == 4))*vec2

tn = (alpha*beta == 0)
tp = (alpha*beta!= 0)


alpha_hat = vector()
beta_hat = vector()
var_alpha = c()
var_beta = c()
p1 = vector()
p2 = vector()
for(i in 1:m)
{
  M[i,] = alpha[i]*X + rnorm(n)
  Y[i,] = beta[i]*M[i,]  + rnorm(1,0.5)*X + rnorm(n)  
  obj1 = lm(M[i,] ~  X )
  obj2 = lm(Y[i,] ~  M[i,] + X)
  
  table1 = coef(summary(obj1))
  table2 = coef(summary(obj2))
  
  
  alpha_hat[i] = table1["X",1]
  beta_hat[i] = table2["M[i, ]",1]
  p1[i] = table1["X",4]
  p2[i] = table2["M[i, ]",4]
  var_alpha[i] = table1["X",2]^2
  var_beta[i] = table2["M[i, ]",2]^2
}

lfdr <- localFDR(alpha_hat, beta_hat, var_alpha, var_beta, twostep = FALSE)
#> iteration= 0 loglik= -2019.293 
#> iteration= 1 loglik= -1988.861 
#> iteration= 2 loglik= -1984.51 
#> iteration= 3 loglik= -1983.162 
#> iteration= 4 loglik= -1982.668 
#> iteration= 5 loglik= -1982.477 
#> iteration= 6 loglik= -1982.401 
#> iteration= 7 loglik= -1982.37 
#> iteration= 8 loglik= -1982.358 
#> iteration= 9 loglik= -1982.353 
#> number of iterations= 10
print(lfdr$pi)
#>   H_alpha H_beta       prob
#> 1       0      0 0.36500291
#> 2       1      0 0.09294027
#> 3       0      1 0.29634692
#> 4       1      1 0.24570990
rej = MLFDR(lfdr$lfdr, size = 0.05)
fdr <- sum(rej*tn)/max(1,sum(rej))
power <- sum(rej*tp)/sum(tp)

print(paste("FDR:", round(fdr,4)))
#> [1] "FDR: 0.0464"
print(paste("Power:", round(power,4)))
#> [1] "Power: 0.6128"
```
