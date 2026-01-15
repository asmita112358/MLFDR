MDACT_pvalues = function(p.M,p.Y,estws){
  sim.num <- length(p.M)
  
  pi.01.est = max(1e-3, estws$alpha01)
  pi.10.est = max(1e-3, estws$alpha10)
  pi.00.est = max(1e-3, estws$alpha00)
  c <- pi.01.est + pi.10.est + pi.00.est
  pi.11.est <- 1-c
  
  pi.01.est.sd = pi.01.est/c
  pi.10.est.sd = pi.10.est/c
  pi.00.est.sd = pi.00.est/c
  
  Ts <- pi.01.est*p.M + pi.10.est*p.Y + pi.00.est*pmax(p.M,p.Y)^2
  
  F00 <- function(t){
    F_0int <- function(p.y){
      c0 <- pi.00.est * p.y^2 + pi.01.est * p.y
      c1 <- t- pi.10.est * p.y
      c2 <- (t - pi.00.est * p.y^2 - pi.10.est * p.y)/pi.01.est
      m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.y))))/(2*pi.00.est)
      s1 <- pmin(pmax(c2,0), 1)
      s2 <- pmin(pmax(m2,0), 1)
      s <- rep(0, length(p.y))
      s[c0>=c1] <- s1[c0>=c1]
      s[c0<c1] <- s2[c0<c1]
      s
    }
    
    x <- tryCatch(integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps), error = function(e) e)
    if(!inherits(x, "error")){
      res <- integrate(F_0int, 0, 1, rel.tol = .Machine$double.eps)
    }else{
      res <- integrate(F_0int, 0, 1)
    }
    res$value
  }
  
  F01 <- function(t){
    c0 <- pi.00.est * p.Y^2 + pi.01.est * p.Y
    c1 <- t- pi.10.est * p.Y
    c2 <- (t - pi.00.est * p.Y^2 - pi.10.est * p.Y)/pi.01.est
    m2 <- (-pi.01.est+sqrt(pmax(0,pi.01.est^2 + 4*pi.00.est*(t-pi.10.est*p.Y))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)
    
    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.Y)
  }
  
  F10 <- function(t){
    c0 <- pi.00.est * p.M^2 + pi.10.est * p.M
    c1 <- t- pi.01.est * p.M
    c2 <- (t - pi.00.est * p.M^2 - pi.01.est * p.M)/pi.10.est
    m2 <- (-pi.10.est+sqrt(pmax(0,pi.10.est^2 + 4*pi.00.est*(t-pi.01.est*p.M))))/(2*pi.00.est)
    s1 <- pmin(pmax(c2,0), 1)
    s2 <- pmin(pmax(m2,0), 1)
    
    (sum(s1[c0>=c1]) + sum(s2[c0<c1])) / length(p.M)
  }
  
  Fnull <- function(t) {
    (1/(1-pi.11.est)*(pi.10.est/(pi.10.est+pi.11.est)*F10(t) + 
                        pi.01.est/(pi.01.est+pi.11.est)*F01(t) + 
                        (pi.00.est - pi.10.est*(pi.00.est + pi.01.est)/(pi.10.est + pi.11.est) - 
                           pi.01.est*(pi.00.est + pi.10.est)/(pi.01.est + pi.11.est))*F00(t)))
  }
  
  pv <- vapply(Ts, Fnull, numeric(1))
  
  return(pv)
}
