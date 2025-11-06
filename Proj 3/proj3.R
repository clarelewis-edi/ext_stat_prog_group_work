#### ---- Q1 Computing Xtilde, X and S -----------------------------------------
library(splines)
getwd()
#setwd("C:\\Users\\Luke Egan\\OneDrive\\Desktop\\Extended Statistical Programming\\Practical 3")
dat <- read.table("engcov.txt", header = T, stringsAsFactors = T)

dat
spline_func <- function(dat, K){
  # Probability function for days from infection until death
  d <- 1:K
  edur <- 3.151
  sdur <- 0.469
  
  pd <- dlnorm(d, edur, sdur)
  pd <- pd/sum(pd)
  
  # Xtilde
  ks <- seq(min(dat$julian) - 30, max(dat$julian), length = K - 2)
  ks_diff <- ks[2] - ks[1]
  
  lower_ks <- ks[1] - (3:1)*ks_diff
  upper_ks <- ks[K-2] + (1:3)*ks_diff
  
  knots <- c(lower_ks, ks, upper_ks)
  
  Xtilde <- splineDesign(knots, (min(dat$julian)-30):max(dat$julian) , outer.ok = T)
  
  
  # X
  n <- nrow(dat)
  X <- matrix(0, n, K)
  for(i in 1:n){
    
    j_upper <- min(29 + i, 80)
    
    for(j in 1:j_upper){
      if((30 + i -j) <= n){
        X[i,] <- X[i,] + (Xtilde[30 + i - j,] * pd[j])
      }
    }
  }
  
  # S
  S <- crossprod(diff(diag(K), diff = 2))
  
  # Return all important information as a list
  spline_mats <- list(Xtilde = Xtilde, X = X, S = S, pd = pd, knots = knots)
  return(spline_mats)
  
}

#### ---- Q2 ---- ####
K <- 80
list1 <- spline_func(dat, K)
#gamma01 <- rep(log(1), K)
gamma02 <- matrix(log(mean(dat$nhs) / K), K)


y <- dat$nhs
n <- nrow(dat)
X <- list1$X
S <- list1$S
lambda <- 5e-5

beta0 <- exp(gamma02)
mu0 <- list1$X %*% beta0
mu0

# Dropped factorial as independent of parameter of interest. As said in pg 2
log_lik <- function(beta){
  mu0 <- X %*% beta
  log_lik <- sum( y*log(mu0) - mu0)
  return(log_lik)
}
log_lik(beta0)

penalty <- function(beta){
  penalty <- (lambda * t(beta) %*% (S %*% beta)) / 2
  return(penalty)
}
penalty(beta0)

d_log_lik <- function(beta){
  mu <- X %*% beta
  d_log_lik <- diag(as.vector(y / mu - 1)) %*% X %*% diag(beta)
  return(d_log_lik)
}
d_log_lik(beta0)

d_penalty <- function(beta){
  d_penalty <- as.vector(diag(beta) %*% (S %*% beta))
  return(as.vector(d_penalty))
}
d_penalty(beta0)

optim_func <- function(beta){
  return(-log_lik(beta) + penalty(beta))
  
}

optim_grad <- function(beta){
  return(-1 * d_log_lik(beta) + d_penalty(beta))
  
}

penalty(beta0)


optim(par = beta0, fn = optim_func, gr = optim_grad, method = 'BFGS')
