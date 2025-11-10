#### ---- Q1 Computing Xtilde, X and S -----------------------------------------
#set.seed(3)
library(splines)
library(ggplot2)
#setwd("C:/Users/Grace Sheahan/ext_stat_prog_group_work/Proj 3")
dat <- read.table("engcov.txt", header = T, stringsAsFactors = T)


spline_func <- function(dat, K){
  # Probability function for days from infection until deah
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
  
  Xtilde <- splineDesign(knots, (min(dat$julian)-30):max(dat$julian))
  
  
  n<- nrow(dat)
  X <- matrix(0, n, K)
  # 
  # X
  # n <- nrow(dat)
  # X <- matrix(0, n, K)
  for(i in 1:n){
    
    first_row <- max(1, i - 50)
    last_row <- min(29 + i, 80 + first_row - 1)
    
    X[i,] <- colSums(Xtilde[first_row:last_row, ] * pd[(last_row - first_row + 1):1] )
  }
  # 
  #   j_upper <- min(29 + i, 80)
  # 
  #   for(j in 1:j_upper){
  #     if((30 + i -j) <= n){
  #       X[i,] <- X[i,] + (Xtilde[30 + i - j,] * pd[j])
  #     }
  #   }
  
  # S
  S <- crossprod(diff(diag(K), diff = 2))
  
  # Return all important information as a list
  spline_mats <- list(Xtilde = Xtilde, X = X, S = S, pd = pd, knots = knots)
  return(spline_mats)
  
}

splines <- spline_func(dat, K = 80)




#### ---- Q(2) ---- ####

K <- 80
#gamma0 <- rep(log(mean(dat$deaths) / K), K)
gamma0 <- rep(0,K)
beta0 <- exp(gamma0)

y <- dat$nhs
n <- nrow(dat)
X <- splines$X
S <- splines$S
lambda0 <- 5e-5
#mu <- X %*% beta

# Dropped factorial as independent of parameter of interest. As said in pg 2

optim_func <- function(y, gamma, X, S, lambda, weights = rep(1, nrow(X)) ){
  beta <- exp(gamma)
  mu <- X %*% beta
  
  log_lik <- sum((y * log(mu) - mu - lgamma(y + 1)) * weights)
  penalty <- 0.5 * (lambda * t(beta) %*% (S %*% beta))
  
  val <- -log_lik + penalty
  return(as.numeric(val))
  
}

optim_grad <- function(y, gamma, X, S, lambda, weights = rep(1, nrow(X))){
  beta <- exp(gamma)
  mu <- X %*% beta
  
  d_log_lik <- apply((diag((as.vector(y/mu - 1))) %*% X %*% diag(beta))*weights, 2, sum)
  d_penalty <- lambda * as.vector(diag(beta) %*% (S %*% beta))
  
  val <- -d_log_lik + d_penalty
  return(val)
}

##### Finite Differencing

eps <- 5e-7
#gamma0 <- rep(log(mean(dat$deaths) / K), K) Already defined
# just calling grad function
grad <- optim_grad(y, gamma0, X, S, lambda0) # only used in finite diff


holding <- numeric(80)
for(i in 1:length(gamma0)){
  gamma1 <- gamma0
  gamma1[i] <- gamma0[i] + eps
  optim_func0 <- optim_func(y, gamma0, X, S, lambda = lambda0)
  optim_func1 <- optim_func(y, gamma1, X, S, lambda = lambda0)
  holding[i] <- (optim_func1 - optim_func0)/eps
  print(holding[i]-grad[i])
  
}


#### Q3 #####
# rename initial opt_vals
optim_vals_1 <- optim(gamma0, optim_func, optim_grad, y = y, X = X, S = S, lambda = lambda0, method = "BFGS")

gamma_hat <- optim_vals_1$par
beta_hat <- exp(gamma_hat)
Xtilde <- splines$Xtilde
mu <- X %*% beta_hat
t <- (min(dat$julian)-30):max(dat$julian)
f <- Xtilde %*% beta_hat

ggplot() +
  geom_line(
    aes(x = t, y = f, col = 'Infection Curve f(t)'),
    linewidth = .75
  ) +
  
  geom_line(
    aes(x = dat$julian, y = mu, col = 'Fitted Deaths μ'),
  ) +
  
  geom_point(
    aes(x = dat$julian, y = dat$nhs, col = 'Observed Deaths'),
    size = 2,
    alpha = .4
  ) +
  
  scale_color_manual(values = c('Infection Curve f(t)' ='dodgerblue',
                                'Fitted Deaths μ' = 'gray11',
                                'Observed Deaths' = 'firebrick')
  ) +
  
  labs(
    x = "Day of Year (2020)",
    y = "Count",
    title = "Daily COVID Deaths and Estimated Infection Curve",
    col = ""
  ) +
  
  theme_light(
    base_size = 13
  ) +
  
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) 


#### Q4 ####
min_BIC <- function(gamma, X, S, y, lambda_vals){
  BIC_vals <- c()
  
  for(i in seq_along(lambda_vals)){
    
    optim_vals <- optim(gamma, optim_func, optim_grad, y = y, X = X, S = S, lambda = lambda_vals[i], method = "BFGS")
    # is it bad practice to redefine an input variable
    optim_gamma <- optim_vals$par
    optim_beta <- exp(optim_gamma)
    mu <- X %*% optim_beta
    
    W <- diag(as.vector(y/(mu^2)))
    H0 <- t(X) %*% (W %*% X)
    H_lambda <- H0 + lambda_vals[i]*S
    EDF <- sum(diag(solve(H_lambda, H0)))
    
    log_lik <- sum(y * log(mu) - mu)
    
    n <- nrow(X)
    
    BIC_vals[i] <- -2*log_lik + log(n)*EDF
    
  }
  ii_min_BIC <- which.min(BIC_vals)
  opt_lambda <- lambda_vals[ii_min_BIC]
  return(opt_lambda)
}

lambda_vals <- exp(seq(-13, -7, length = 50))

lambda_hat <- min_BIC(gamma0, X, S, y, lambda_vals)
optim_vals_2 <- optim(gamma0, optim_func, optim_grad, y = y, X = X, S = S, lambda = lambda_hat, method = "BFGS")

# updated values based on new optimisation
gamma_hat <- optim_vals_2$par
beta_hat <- exp(gamma_hat)
Xtilde <- splines$Xtilde
mu <- X %*% beta_hat
t <- (min(dat$julian)-30):max(dat$julian)
f <- Xtilde %*% beta_hat

# install.packages("matrixcalc")
# library(matrixcalc)
# is.symmetric.matrix(H_lambda)
# is.positive.definite(H_lambda)
# Probably should use QR-decomp

# Q5
#set.seed(3)
n <- nrow(dat)
nb <- 200
f_b <- matrix(0, nb, length(f))

for (i in 1:nb){
  wb <- tabulate(sample(n,replace=TRUE),n) ## non-para bootstrap weights
  optim_vals_b <- optim(gamma0, optim_func, optim_grad, y = y, X = X, S = S, lambda = lambda_hat, weights = wb, method = "BFGS")
  gamma_hat_b <- optim_vals_b$par
  beta_hat_b <- exp(gamma_hat_b)
  f_b[i,] <- Xtilde %*% beta_hat_b
}

f_b_limits <- apply(f_b, 2, quantile, probs = c(0.025, 0.95))


ggplot() +
  geom_line(
    aes(x = t, y = f, col = 'Infection Curve f(t)'),
    linewidth = .75
  ) +
  
  geom_ribbon(
    aes(x = t, ymin = f_b_limits[1,], ymax = f_b_limits[2,],
        fill = 'Confidence Interval'), alpha = 0.5
  ) +
  
  geom_line(
    aes(x = dat$julian, y = mu, col = 'Fitted Deaths μ'),
  ) +
  
  geom_point(
    aes(x = dat$julian, y = dat$nhs, col = 'Observed Deaths'),
    size = 2,
    alpha = .4
  ) +
  
  scale_color_manual(values = c('Infection Curve f(t)' ='dodgerblue',
                                'Fitted Deaths μ' = 'gray11',
                                'Observed Deaths' = 'firebrick')
  ) +
  
  scale_fill_manual(values = 'gray70'
  ) +
  
  labs(
    x = "Day of Year (2020)",
    y = "Count",
    title = "Daily COVID Deaths and Estimated Infection Curve",
    col = "",
    fill = ""
  ) +

  theme_light(
    base_size = 13
  ) +
  
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) 
