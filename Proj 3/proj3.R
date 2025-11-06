#### ---- Q1 Computing Xtilde, X and S -----------------------------------------
library(splines)
getwd()
setwd("C:\\Users\\Luke Egan\\OneDrive\\Desktop\\Extended Statistical Programming\\Practical 3")
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
gamma01 <- rep(log(1), K)
gamma02 <- rep(log(mean(dat$deaths) / K), K)


y <- dat$nhs
n <- nrow(dat)
X <- list1$X
S <- list1$S
lambda <- 5e-5

beta0 <- exp(gamma02)
mu0 <- list1$X %*% beta0
mu0
class(d_penalty)

# Dropped factorial as independent of parameter of interest. As said in pg 2
log_lik <- function(beta0){
  mu0 <- X %*% beta0
  log_lik <- y * log(mu0) - mu0 
  return(log_lik)
}

penalty <- function(beta0){
  penalty <- (lambda * t(beta0) %*% (S %*% beta0)) / 2
  return(as.vector(penalty))
  
  
}
d_log_lik <- function(beta0){
  mu0 <- X %*% beta0
  d_log_lik <- diag(as.vector((y / (mu0 - 1))))
  return(d_log_lik)
}


d_penalty <- function(beta0){
  d_penalty <- as.vector(diag(beta0) %*% (S %*% beta0))
  return(as.vector(d_penalty))
  
}

optim_func <- function(beta0){
  return(-log_lik(beta0) + penalty(beta0))
  
}

optim_grad <- function(beta0){
  return(-1 * d_log_lik(beta0) + d_penalty(beta0))
  
}

penalty(beta0)


optim(beta0, fn = optim_func(beta0), gr <- optim_grad(beta0))