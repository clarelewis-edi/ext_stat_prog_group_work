# -----------------------------------------------------------------------------# 
# Clare Lewis (s2879721), Grace Sheahan (s2898645), Luke Egan (s2837709)
#
# Clare: Collaborated on all sections, taking a lead on questions 2 and 6.
# Grace: Collaborated on all sections, taking a lead on questions 3 and 5
# Luke: Collaborated on all sections, taking a lead on questions 1 and 4
#
# We all feel that we equally contributed to this project, almost entirely through 
# in-person collaboration along with some independent coding, roughly completing 
# 1/3 of the work each
#
# Our github repository can be found at;
# https://github.com/clarelewis-edi/ext_stat_prog_group_work
#
# ---- Introduction ------------------------------------------------------------
# 
# This project works with the provided data set 'engcov.txt' on deaths from 
# Covid-19 in English hospitals against the day of year and aims to use this data
# to infer a model of the new infections per day which resulted in these deaths.

# The 'engcov.txt' file contains records on 150 days of deaths (150 rows of data)
# and contains 5 columns; "date", "deaths", "julian", "gov", "nhs".
# Date - Date for which the death rates were collected
# Julian - The day of the year corresponding to the date
# Deaths/NHS - These rows contain equivalent data on the number of deaths from 
# Covid reported by the NHS

# ---- Computing Xtilde, X and S -----------------------------------------------
# Import necessary libraries
library(splines)
library(ggplot2)

# Read in the dataset
dat <- read.table("engcov.txt", header = T, stringsAsFactors = T)

# SPLINE_FUNC
# Constructs the spline-based matrices that will be used to fit the deconvolution
# model to the Covid death data
#
# Inputs:
#     - dat: Dataframe containing a column "julian" which contains the days for
#            observed deaths
#
#     - K: Total number of spline basis functions
#
# Outputs:
#     List containing:
#     - Xtilde: Spline basis matrix over the infection times
#     
#     - X: Model matrix for deaths
#
#     - S: Penalty matrix with smoothing parameter for smoother model fit

spline_func <- function(dat, K){
  
  # ---- Probability function for days from infection until death                
  
  # It is given that there is data to suggest that if d is the interval from
  # infection to death then log(d) ~ N(3.151, 0.469^2)
  d <- 1:80
  edur <- 3.151
  sdur <- 0.469
  
  # Calculate model probability of each fatal disease duration for each day
  pd <- dlnorm(d, edur, sdur)
  # Normalise probabilities
  pd <- pd/sum(pd)
  
  # ---- Xtilde

  # Sequence for inner K-2 knots covering the interval over which f(t) 
  # is evaluated. This extends to 30 days before the first observed death (over
  # the infection period)
  ks <- seq(min(dat$julian) - 30, max(dat$julian), length = K - 2)
  # Interval size between knots
  ks_diff <- ks[2] - ks[1]
  
  # Add boundary knots to the start and end of the sequence. Accounts for boundary
  # splines, the tails of which could effect the edge cases
  lower_ks <- ks[1] - (3:1)*ks_diff
  upper_ks <- ks[K-2] + (1:3)*ks_diff
  
  # Complete sequence of knots 
  knots <- c(lower_ks, ks, upper_ks)
  
  # Matrix of spline coefficients over period of infections 
  Xtilde <- splineDesign(knots, (min(dat$julian)-30):max(dat$julian))
  
  # ---- X

  # Initialise X
  n <- nrow(dat)
  X <- matrix(0, n, K)
  
  
  for(i in 1:n){
    # Define the range for the plausible infection period for deaths on day i
    lower <- max(1, i - 50)
    upper <- min(29 + i, 80 + lower - 1)
    
    # Builds the ith row of the model matrix for the deaths based on the formula given
    X[i,] <- colSums(Xtilde[lower:upper, ] * pd[(upper - lower + 1):1])
  }
  
  # ---- S
  
  # Compute the penalty matrix
  S <- crossprod(diff(diag(K), diff = 2))
  
  # Return all important information as a list
  spline_mats <- list(Xtilde = Xtilde, X = X, S = S)
  return(spline_mats)
}

K <- 80
splines <- spline_func(dat, K)

# PEN_NLL
# Calculates the negative penalised log likelihood function
#
# Inputs:
#     - y: number of deaths per day
#     - X: model matrix for deaths
#     - S: penalty matrix
#     - gamma: log of the coefficients of the infection model
#     - lambda: smoothing parameter
#     - weights: weighting of the rows of dat (default to all ones except in
#                the case of bootstrapping)
#
# Outputs:
#     - val: vector of the negative penalised log likelihood values

pen_nll <- function(y, gamma, X, S, lambda, weights = rep(1, nrow(X)) ){
  # To ensure that beta is positive, it is defined as the exponential of gamma
  beta <- exp(gamma)
  mu <- X %*% beta
  
  # log likelihood and penalty calculated based on formulas provided
  log_lik <- sum((y * log(mu) - mu - lgamma(y + 1)) * weights)
  penalty <- 0.5 * (lambda * crossprod(beta, S %*% beta))

  val <- -log_lik + penalty
  return(as.numeric(val))
}

# GRAD_PEN_NLL
# Calculate the gradient function based on the derivatives of the log likelihood
# and penalty
#
# Inputs: Same as pen_nll
#
# Outputs:
#     - val: vector of negative derivative of the penalised log likelihood

grad_pen_nll <- function(y, gamma, X, S, lambda, weights = rep(1, nrow(X))){
  beta <- exp(gamma)
  mu <- X %*% beta
  
  # derivatives of log likelihood and penalty based on formulas provided
  d_log_lik <- apply(((as.vector(y/mu - 1)) * t(beta * t(X)) * weights), 2, sum)
  d_penalty <- lambda * as.vector(diag(beta) %*% (S %*% beta))
  
  val <- -d_log_lik + d_penalty
  return(val)
}

# ---- Finite Differencing -----------------------------------------------------

# Initial estimate of gamma values
gamma0 <- rep(0,K)

# Defining values that will be used as inputs in functions
y <- dat$nhs # Number of deaths per day
n <- nrow(dat) # Number of days, over which data is collected
X <- splines$X # Model matrix for deaths
S <- splines$S # Penalty matrix
lambda0 <- 5e-5 # Initial estimate of smoothing parameter

eps <- 5e-7 # Finite differencing interval

# Calling the gradient function plugging in the above as inputs
grad <- grad_pen_nll(y, gamma0, X, S, lambda0)

# Run pen_nll to find negative log likelihood corresponding to gamma0
pen_nll0 <- pen_nll(y, gamma0, X, S, lambda = lambda0)

est_grad <- numeric(80)
for(i in 1:length(gamma0)){
  gamma1 <- gamma0
  gamma1[i] <- gamma0[i] + eps # Increase gamma0[i] by eps
  
  # Run pen_nll to find negative log likelihood corresponding to gamma1
  pen_nll1 <- pen_nll(y, gamma1, X, S, lambda = lambda0)
  est_grad[i] <- (pen_nll1 - pen_nll0)/eps # Approximate the gradient
  # Print difference between estimated and actual gradient 
  # print(est_grad[i]-grad[i]) # (Should ≈ 0)
}

# ---- Initial fit of model ----------------------------------------------------

# Use optimisation and gradient functions to find the optimal values of gamma
optim_vals_1 <- optim(gamma0, pen_nll, grad_pen_nll, y = y, X = X, S = S, 
                      lambda = lambda0, method = "BFGS")

# Assign variables based on optimisation output
gamma_hat <- optim_vals_1$par
beta_hat <- exp(gamma_hat)

Xtilde <- splines$Xtilde # Extract Xtilde from the output of the splines function
mu <- X %*% beta_hat # Fitted deaths
t <- (min(dat$julian)-30):max(dat$julian) # Infection period
f <- Xtilde %*% beta_hat # Modelled infection curve

ggplot() +

  # Plot the fitted deaths (mu)
  geom_line(
    aes(x = dat$julian, y = mu, col = 'Fitted Deaths μ'),
  ) +

  # Plot the actual observed deaths from the dataset
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
    col = "" # We don't need a "colour" label so assigning this as empty
  )


# ---- Optimising Lambda and Gamma Parameters ----------------------------------

# MIN_BIC
# Determines the optimal values based on minimising the BIC criterion
#
# Inputs:
#     - gamma, X, S, y: As above in pen_nll and grad_pen_nll functions
#     - log_lambda_vals: Sequence of log of smoothing parameters over which to 
#                   search for the optimal lambda and corresponding gammas

# Outputs:
#     - hat_params: list containing the optimal lambda value and corresponding
#                 gamma values

min_BIC <- function(gamma, X, S, y, log_lambda_vals){
  BIC_val <- 1000 # Defining an arbitrarily large starting value for BIC for which
  # the BIC can be improved
  # The log of lambda values are used as inputs to ensure positive lambda values
  # The lambda values are found by exponentiating these log values
  lambda_vals <- exp(log_lambda_vals)
  
  for(i in seq_along(lambda_vals)){
    # Run optim to find the optimal gamma values corresponding to the lambda value
    optim_vals <- optim(gamma, pen_nll, grad_pen_nll, y = y, X = X, S = S, 
                        lambda = lambda_vals[i], method = "BFGS")
    optim_gamma <- optim_vals$par 
    optim_beta <- exp(optim_gamma)
    
    mu <- X %*% optim_beta # Calculate mu
    
    # Defining the inputs H0 and W, and the resulting Hessian in order to
    # calculate the BIC
    W <- diag(as.vector(y/(mu^2)))
    H0 <- t(X) %*% (W %*% X)
    H_lambda <- H0 + lambda_vals[i]*S
    
    # Use cholesky decomposition to solve (H_lambda inverse * H0) and calculate
    # the effective degrees of freedom (EDF)
    H_lambda_decomp <- chol(H_lambda)
    trace_input <- backsolve(H_lambda_decomp,forwardsolve(t(H_lambda_decomp),H0))
    EDF <- sum(diag(trace_input))
    
    log_lik <- sum(y * log(mu) - mu) # Calculate the log likelihood
    
    n <- nrow(X)
    
    # Define initial BIC value
    if (i == 1){
      BIC_val <-  -2*log_lik + log(n)*EDF
      }
    # Update the BIC value based on comparison with current minimum value
    else { 
      BIC_val <- min(BIC_val, -2*log_lik + log(n)*EDF)
    }
    # Save the corresponding optimal gamma values if the BIC value has been updated
    if (BIC_val == (-2*log_lik + log(n)*EDF)){  
     opt_lambda <- lambda_vals[i]
     opt_gamma <- optim_gamma
    }
  }
  # Create list of optimal parameters to return
  hat_params <- list(lambda_hat = opt_lambda, gamma_hat = opt_gamma)
  return(hat_params)
}

# Define log lambda values over which to grid search
log_lambda_vals <- seq(-13, -7, length = 50)

# Run and store output of min_BIC, and extract optimal lambda and gamma
hat_params <- min_BIC(gamma_hat, X, S, y, log_lambda_vals)
lambda_hat <- hat_params$lambda_hat
gamma_hat_updated <- hat_params$gamma_hat
beta_hat_updated <- exp(gamma_hat_updated) # Calculate optimal beta

# Update mu and the infection curve with optimal beta
mu_updated <- X %*% beta_hat_updated
f_updated <- Xtilde %*% beta_hat_updated

# ---- Bootstrapping------------------------------------------------------------

# Define the number of bootstraps to complete and a holding matrix for the outputs
n <- nrow(dat)
nb <- 200
f_b <- matrix(0, nb, length(f_updated))

for (i in 1:nb){
  wb <- tabulate(sample(n,replace=TRUE),n) ## Non-parametric bootstrap weights
  
  # Optimise gamma (and hence beta) for the bootstrap sample
  optim_vals_b <- optim(gamma0, pen_nll, grad_pen_nll, y = y, X = X, S = S, 
                        lambda = lambda_hat, weights = wb, method = "BFGS")
  gamma_hat_b <- optim_vals_b$par
  beta_hat_b <- exp(gamma_hat_b)
  # Store the bootstrapped infection model values in the holding matrix
  f_b[i,] <- Xtilde %*% beta_hat_b
}

# Find the confidence limits of the bootstrapped infection curves
f_b_limits <- apply(f_b, 2, quantile, probs = c(0.025, 0.975))

# ---- Final Plot --------------------------------------------------------------

ggplot() +
  
  # We want the confidence limits to be furthest back on the plot so we plot
  # those first confidence limits determined from bootstrapping
  geom_ribbon(
    aes(x = t, ymin = f_b_limits[1,], ymax = f_b_limits[2,],
        fill = 'Confidence Limits'), alpha = 0.5
  ) +
  
  # Plot our fitted deaths (mu) next
  geom_line(
    aes(x = dat$julian, y = mu_updated, col = 'Fitted Deaths μ'),
  ) +
  
  # Add our actual observed deaths on top of that
  geom_point(
    aes(x = dat$julian, y = dat$nhs, col = 'Observed Deaths'),
    size = 2,
    alpha = .4
  ) +
  
  # Finally add the infection curve f(t) updated with our optimal values 
  geom_line(
    aes(x = t, y = f_updated, col = 'Infection Curve f(t)'),
    linewidth = .75
  ) +
  
  # As for the previous graph, want to create a legend that includes each plot 
  # element so manually assigning colours here so r will create it
  scale_color_manual(values = c('Infection Curve f(t)' ='royalblue4',
                                'Fitted Deaths μ' = 'gray11',
                                'Observed Deaths' = 'firebrick')
  ) +
  
  # As for the manual colour but for the fill of the confidence limits
  scale_fill_manual(values = 'lightskyblue'
  ) +
  
  labs(
    x = "Day of Year (2020)",
    y = "Count",
    title = "Daily COVID Deaths and Estimated Infection Curve",
    # Removing the colour and fill titles by leaving blank
    col = "",
    fill = ""
  ) +
  
  # Change the theme to enhance plots readability
  theme_light(
    base_size = 13
  ) +
  
  # Make title bold and move the legend to the top of the graph
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  ) 