#### ---- Q1 Computing Xtilde, X and S -----------------------------------------
library(splines)
getwd()
setwd("C:\\Users\\Luke Egan\\OneDrive\\Desktop\\Extended Statistical Programming\\Practical 3")
dat <- read.table("engcov.txt", header = T, stringsAsFactors = T)


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
