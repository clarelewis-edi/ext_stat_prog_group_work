# ------------------------------------------------------------------------------ 
# Clare Lewis (s2879721), Grace Sheahan (s2898645), Luke Egan (s2837709)
#
# Clare: 
# Grace: 
# Luke: 
#
# We all feel that we equally contributed to this project, primarily through 
# in-person collaboration along with some independent coding, roughly completing 
# 1/3 of the work each
#
# Our github repository can be found at;
# https://github.com/clarelewis-edi/ext_stat_prog_group_work
#
# ---- Introduction ------------------------------------------------------------
# 
# This project works with the provided data set 'engcov.txt' on deaths from Covid-19
# in English hospitals against the day of year and aims to use this data to infer
# a model of the new infections per day which resulted in these deaths.

# The 'engcov.txt' file contains records on 150 days of deaths (150 rows of data)
# and contains 5 columns; "date", "deaths", "julian", "gov", "nhs".
# Date - Date for which the death rates were collected
# Julian - The day of the year corresponding to the date
# Deaths/NHS - These rows contain equivalent 

#### ---- Q1 Computing Xtilde, X and S -----------------------------------------
#set.seed(3)
start <- Sys.time()
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
  
  for(i in 1:n){
    
    lower <- max(1, i - 50)
    upper <- min(29 + i, 80 + lower - 1)
    
    X[i,] <- colSums(Xtilde[lower:upper, ] * pd[(upper - lower + 1):1])
  }
  
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
  
  d_log_lik <- apply(((as.vector(y/mu - 1)) * t(beta * t(X)) * weights), 2, sum)
  d_penalty <- lambda * as.vector(diag(beta) %*% (S %*% beta))
  
  val <- -d_log_lik + d_penalty
  return(val)
}

##### Finite Differencing

eps <- 5e-7

# just calling grad function
grad <- optim_grad(y, gamma0, X, S, lambda0) # only used in finite diff

est_grad <- numeric(80)
for(i in 1:length(gamma0)){
  gamma1 <- gamma0
  gamma1[i] <- gamma0[i] + eps
  optim_func0 <- optim_func(y, gamma0, X, S, lambda = lambda0)
  optim_func1 <- optim_func(y, gamma1, X, S, lambda = lambda0)
  est_grad[i] <- (optim_func1 - optim_func0)/eps
  print(est_grad[i]-grad[i])
}

#### Q3 #####
# use our optimisation and gradient functions to get the optimal values of gamma
optim_vals_1 <- optim(gamma0, optim_func, optim_grad, y = y, X = X, S = S, lambda = lambda0, method = "BFGS")

# assign variables based on optimisation output
gamma_hat <- optim_vals_1$par
beta_hat <- exp(gamma_hat)
Xtilde <- splines$Xtilde
mu <- X %*% beta_hat
t <- (min(dat$julian)-30):max(dat$julian)
f <- Xtilde %*% beta_hat


ggplot() +
  
  # Plot the fitted deaths (mu)
  geom_line(
    aes(x = dat$julian, y = mu, col = 'Fitted Deaths μ'),
  ) +
  
  # Plot the actual observed deaths from our dataset
  geom_point(
    aes(x = dat$julian, y = dat$nhs, col = 'Observed Deaths'),
    size = 2,
    alpha = .4
  ) +
  
  # Plot the infection curve f(t)
  geom_line(
    aes(x = t, y = f, col = 'Infection Curve f(t)'),
    linewidth = .75
  ) +
  
  # In order to show the legend with the different plots, we put the colour
  # inside the aes() with the name we want to print and then manually assigning
  # the colours
  scale_color_manual(values = c('Infection Curve f(t)' ='royalblue4',
                                'Fitted Deaths μ' = 'gray11',
                                'Observed Deaths' = 'firebrick')
  ) +
  
  labs(
    x = "Day of Year (2020)",
    y = "Count",
    title = "Daily COVID Deaths and Estimated Infection Curve",
    col = "" # we don't need a "colour" label so assigning this as empty
  )

#### Q4 ####
min_BIC <- function(gamma, X, S, y, lambda_vals){
  BIC_val <- 1000 
  for(i in seq_along(lambda_vals)){
    
    optim_vals <- optim(gamma, optim_func, optim_grad, y = y, X = X, S = S, lambda = lambda_vals[i], method = "BFGS")
    # is it bad practice to redefine an input variable
    optim_gamma <- optim_vals$par
    optim_beta <- exp(optim_gamma)
    mu <- X %*% optim_beta
    
    W <- diag(as.vector(y/(mu^2)))
    H0 <- t(X) %*% (W %*% X)
    H_lambda <- H0 + lambda_vals[i]*S
    cholesk <- chol(H_lambda)
    solved <- backsolve(cholesk,forwardsolve(t(cholesk),H0))
    EDF <- sum(diag(solved))
    
    log_lik <- sum(y * log(mu) - mu)
    
    n <- nrow(X)
    
    BIC_val <- min(BIC_val, -2*log_lik + log(n)*EDF)
    
    if (BIC_val == (-2*log_lik + log(n)*EDF)){
      opt_lambda <- lambda_vals[i]
    }
  }
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

# Q5
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
  
  # We want the confidence intervals to be furthest back on the plot so we plot
  # those first
  # confidence intervals determined from bootstrapping
  geom_ribbon(
    aes(x = t, ymin = f_b_limits[1,], ymax = f_b_limits[2,],
        fill = 'Confidence Interval'), alpha = 0.5
  ) +
  
  # Plot our fitted deaths (mu) next
  geom_line(
    aes(x = dat$julian, y = mu, col = 'Fitted Deaths μ'),
  ) +
  
  # Add our actual observed deaths on top of that
  geom_point(
    aes(x = dat$julian, y = dat$nhs, col = 'Observed Deaths'),
    size = 2,
    alpha = .4
  ) +
  
  # Finally add the infection curve f(t) updated with our optimal values 
  geom_line(
    aes(x = t, y = f, col = 'Infection Curve f(t)'),
    linewidth = .75
  ) +
  
  # same as previous graph, want to create a legend that includes each plot element
  # so manually assigning colours here so r will create it
  scale_color_manual(values = c('Infection Curve f(t)' ='royalblue4',
                                'Fitted Deaths μ' = 'gray11',
                                'Observed Deaths' = 'firebrick')
  ) +
  
  # doing the same thing as the manual colour but for the fill of the confidence interva
  scale_fill_manual(values = 'lightskyblue'
  ) +
  
  labs(
    x = "Day of Year (2020)",
    y = "Count",
    title = "Daily COVID Deaths and Estimated Infection Curve",
    # in this case we want the colour and fill titles both to be removed
    col = "",
    fill = ""
  ) +
  
  # change the theme to make the plots easier to read
  theme_light(
    base_size = 13
  ) +
  
  # Make title bold and move the legend to the top of the graph
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) 
Sys.time() - start