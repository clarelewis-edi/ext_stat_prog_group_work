# ------------------------------------------------------------------------------ 
# Clare Lewis (s2879721), Grace Sheahan (s2898645), Luke Egan (s2837709)

# Clare: Completed plotting steps and household function, collaborated on
# other sections
# Grace: Completed original iteration of nseir function and refined plot
# function, collaborated on other sections
# Luke: Completed get net function and refined nseir function, collaborated on
# other sections

# We all feel that we equally contributed to this project, through a mixture of
# independent coding and in-person collaboration, roughly completing 1/3 of the
# work each

#---- Introduction -------------------------------------------------------------
# 
# This project aims to build a SEIR stochastic model of the spread of an epidemic
# within a population with a rudimentary social structure.
# This is done through expansion of the basic SEIR model to account for
# households and set social networks, through both of which infection is more
# likely to spread than through random population mixing.
# 
# In the below code functions are built to model the social structure for a
# population - assigning every individual a household and social network.
# 
# A function is built to simulate the spread of an epidemic through this population
# for a given interval - taking into account the epidemiological factors such as
# infectiousness and recovery probability as well as the societal connections
# created in our model.
# 
# A plot function is then built to provide clear visualisation of the progression
# of the epidemic.
# 
# These functions are then used to model, visualise and compare four scenarios in
# order to better understand the effects of applying sociability and social
# structure to such models. (Our findings on these are present below, at the end 
# of the code).

## BUILDING THE MODEL ##

# ---- Social Structures ----

# Households #

# BUILD H
# Function to allocate households to individuals in the population
#
# Inputs:
# - pop_size: Total population size
# - hmax: Maximum number of people in a household, the default for this value is 5
#
# Output:
# - households: A vector of pop_size length that has an integer value corresponding
# to each individual that identifies what household they are a member of

build.h <- function(pop_size, hmax=5) {
  # the outer sample function avoids ordering the population indices by household
  households <- sample(rep(1:pop_size, sample(1:hmax, pop_size, replace = TRUE))[1:pop_size])
  return(households)
}

# Network #

# GET NET
# Function that creates the non-household network structure for each person
# in the population
#
# Inputs:
# - beta: Vector of sociability parameters of the population
# - h: Vector assigning each person to a household
# - nc: Average number of contacts per person
#
# Output:
# - net_list: A list where each element is a vector containing the indices of
#               individuals that person i has contact with


get.net <- function(beta, h, nc = 15){
  n <- length(beta)
  
  # Only concerned with the probabilities in the upper triangle of the
  # network probability matrix
  # Network probability matrix is symmetric
  # Pr(link i->j) = Pr(link j->i)
  # Also want to avoid links being created twice
  
  # Indices for upper triangle of network probability matrix
  upper_triangle_i <- which(upper.tri(diag(n)), arr.ind = T)
  
  # Vectorised computation of network probabilities
  prob_vec <- (nc/((mean(beta)^2)*(n-1)))*
    (beta[upper_triangle_i[,1]]*beta[upper_triangle_i[,2]])
  
  # The probability of creating a link with people in the same household should be 0
  same_hh <- h[upper_triangle_i[,1]] == h[upper_triangle_i[,2]]
  prob_vec[same_hh] <- 0
  
  # Using the network probabilities between individuals
  # Perform Bernoulli trial for each possible contact to see if link is created
  link_vec <- rbinom(length(prob_vec), size = 1, prob = prob_vec)
  
  # If a link is made i->j, must account for j->i
  indivduals_i <- upper_triangle_i[link_vec == 1,1]
  indivduals_j <- upper_triangle_i[link_vec == 1,2]
  # Matrix with two columns which holds all links created i->j and j->i
  all_links <- rbind(cbind(indivduals_i, indivduals_j), cbind(indivduals_j, indivduals_i))
  
  # Create list of all contacts for each individual
  # Start with individuals with at least one contact
  existing_link_list <- split(all_links[,2], all_links[,1])
  net_list <- as.list(rep(0, n))
  net_list[as.integer(names(existing_link_list))] <- existing_link_list
  # Now add in individuals with no contacts 
  no_link <- setdiff(1:n, as.integer(names(existing_link_list)))
  net_list[no_link] <- list(integer(0))
  
  # Return full list of network links
  return(net_list)

}

# ---- NSEIR Model ---- 

# NSEIR FUNCTION
# Function that takes in population and epidemiological factors and models
# the progression of the epidemic through the population for a given duration.
#
# Inputs:
# - beta: Vector of sociability parameters of the population
# - h: Vector assigning each person to a household
# - alink: List defining the regular contacts of each person (as returned by get.net)
# - alpha: Probability vector of (α_h, α_c, α_r)
#     - α_h: Probability of infecting a member of the same household
#     - α_c: Probability of infecting a member of their regular contacts network
#     - α_r: Constant factor of the probability of infection spreading randomly
#             between two people, irrespective of other relations
# - delta: Probability of moving from the infected to recovered class
# - gamma: Probability of moving from the exposed to infected class
# - nc: Average number of contacts per person
# - nt: Number of days (duration) for which to simulate the model
# - pinf: Proportion of the population that are initially infected
#
# Output:
# - SEIRt: A list  of vectors S, E, I, R, and t, which give the number of people
#          in each state by day, and the corresponding days.

nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, 
                  nc = 15, nt = 100, pinf = .005) {
  
  n <- length(beta) # The population size
  t <- 1:nt # Vector of the days of the simulation 
  
  x <- rep(0, n) # Initialising a state vector with all susceptible
  x[sample(c(1:n), round(pinf*n))] <- 1 # Moving pinf proportion to infected state
  
  # Defining constant factor of probability of random infection between two people
  daily_constant <- (alpha[3] * nc)/((mean(beta)^2)*(n-1)) 
  
  # Defining vectors that will contain the population in the states on each day
  S <- E <- I <- R <- rep(0, nt)
  S[1] <- sum(x == 0) # Initial susceptible population is all non-infected people
  I[1] <- sum(x == 2) # Initial infected population
  
  for (i in 2:nt) { # Looping over days to update and store the state for each day
    
    # Uniform random deviates
    # Two needed to ensure that random infection is distinct from relational infection
    u <- runif(n) 
    v <- runif(n)
    
    prev_infectious <- which(x == 2) # Vector of indices of infected people
    # This assumes that moving to state E occurs the day after contact with state I
    
    x[x == 2 & u < delta] <- 3 # Moves population in I to R with probability delta
    x[x == 1 & u < gamma] <- 2 # Moves population in E to I with probability gamma
    
    # Susceptible people can only become exposed if there were people infected
    # the day prior
    if (length(prev_infectious) > 0) {
      
      # Household
      # Count the number of infected people per household
      infected_hh_member_count <- tabulate(h[prev_infectious], nbins = max(h)) 
      
      # If m in a house are infected, the probability of a susceptible person in
      # the house becoming exposed is '1 - (1 - α_h)^m'
      # Define probability of infections by each household
      alpha_by_hh <- 1 - (1 - alpha[1])^(infected_hh_member_count)  
      # Vector of length n which holds the probability of infection by one's household
      alpha_hh <- alpha_by_hh[h] 
      # Moves population in S to E by their household probability
      x[x == 0 & u < alpha_hh] <- 1 
      
      # Network
      # Vector of those in the network of any infected person
      infected_network_ids <- unlist(lapply(alink[prev_infectious], as.numeric)) 
      # Only proceeds if the infected have people in their network(s)
      if (length(infected_network_ids) > 0) { 
        # Counts the number of infected people in each persons network
        infected_by_network_count <- tabulate(infected_network_ids, nbins = n) 
        # Defines each person's probability of infection by their network
        alpha_net <- 1 - (1 - alpha[2])^(infected_by_network_count) 
        # Moves population in S to E by their network probability
        x[x == 0 & u < alpha_net] <- 1 
      }
      
      # Random mixing
      alpha_daily <- daily_constant * outer(beta[prev_infectious], beta, FUN = "*")
      # Person i's probability of infection from random mixing is 
      #'1 - Pr(not catching through random mixing from anyone)'
      # That is, Pr(not catching) = product of Pr(not catching from each j) 
      
      # Defines each person's probability of infection by random mixing (based
      # on their specific alpha_dailys with those in prev_infectious)
      alpha_rm <- 1 - apply(1 - alpha_daily, 2, prod) 
      # Moves population in S to E by their random mixing probability
      x[x == 0 & v < alpha_rm] <- 1 
    }
    
    # Define the state vector elements as the number of people in the state by day
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
  }
  
  # The number of people in each state by day and the days
  SEIRt <- list(S = S, E = E, I = I, R = R, t = t) 
  return(SEIRt)
}

# ---- Model Visualisation ----


# PLOT OUTPUT:
# Function sets up parameters and produces the plot for the nseir model
# 
# Inputs:
# - beta, h, alink, alpha, delta, gamma, nc, nt, and pinf are all inputs for
# the nseir function and are defined above that function
# - title: Title to be displayed above the plot
#
# Output: 
# - A plot showing the model results

plot.output <- function(beta,h,alink,alpha,delta=.2,gamma=.4,nc=15,nt = 100,
                        pinf = .005, title="NSEIR Model") {
  # Define epi as a run of the nseir model
  epi <- nseir(beta,h,alink,alpha,delta,gamma,nc,nt,pinf)
  # First, plot the susceptible population (gray) and set up some parameters
  # for the plots
  plot(epi$S,type="l",ylim=c(0,max(epi$S)),
       # Title from function input, manual axis labels, and adjust line width
       main=title,xlab="",ylab="Population",col="gray11", lwd = 2)
  points(epi$E,type="l",col="dodgerblue", lwd = 2) # Plot exposed population (blue)
  points(epi$I,type="l",col="firebrick", lwd = 2) # Plot infected population (red)
  points(epi$R,type="l",col="green4", lwd = 2) # Plot recovered population (green)
  # Add grid lines to improve readability
  grid(nx = NULL, ny = NULL,
       lty = 2,
       col = "gray",
       lwd = .5) 
  legend(x="right",
         legend = c("Susceptible", "Exposed", "Infected", "Recovered"),
         bty='n',
         lty=1,
         lwd = 2, seg.len = 1, cex = .75, # symbol and size of legend
         col = c("gray11","dodgerblue","firebrick", "green4"),
         )
  # manually add x axis label to have the ability to adjust the spacing between
  # the axis and the label
  title(xlab="Day", line=1.75)
}

# ---- Running the Model  ----

# Defining Parameters #
pop_size <- 10000
h <- build.h(pop_size)
beta <- runif(pop_size) # Sociability vector
alink <- get.net(runif(pop_size),h,nc = 15) # Building social networks
beta_constant <- rep(sum(beta)/length(beta), length(beta)) # Defining constant beta
# Set plot window up for multiple plots, adjust margins, rotate axes labels to
# be horizontal
par(mfrow=c(2,2), mar=c(4,4,2,1), las = 1) 

# Run 1: Default parameters
plot.output(beta,h,alink,alpha=c(.1,.01,.01), title = "Full Model") 
# Run 2: Random mixing only
plot.output(beta,h,alink,alpha=c(0,0,.04), title = "Only Random Mixing")
# Run 3: Set a constant beta
plot.output(beta = beta_constant,h,alink,alpha=c(.1,.01,.01),
            title = "Full Model with Constant Beta") 
# Run 4: Constant beta and random mixing
plot.output(beta = beta_constant,h,alink,alpha=c(0,0,.04),
            title = "Constant Beta and Random Mixing") 


