# ------------------------------------------------------------------------------ 
# Clare Lewis (s2879721), Grace Sheahan (s2898645), Luke Egan (s2837709)

# Clare: 
# Grace: 
# Luke: 

################################################################################ Comments should not be longer than this line of ##
# ------------------------------------------------------------------------------

# Setting working directory to folder for Assessment 1
# setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work/Proj 2")
# setwd("C:/Users/Grace Sheahan/ext_stat_prog_group_work/Proj 2")
# setwd("C:/Users/Luke Egan/Desktop/Extended Statistical Programming/ext_stat_prog_group_work/Proj 2")

#---- Introduction -------------------------------------------------------------
# 
# This project aims to build a SEIR stochastic model of the spread of an epidemic
# within a population with a rudimentary social structure.
# This is done through expansion of the basic SEIR model to account for
# households and set social networks, through both of which infection is more
# likely to spread than through random population mixing.
# 
# In the below code functions are build to model the social structure for a
# population - assigning every individual a household and social network.
# 
# A function is build to simulate the spread of an epidemic through this population
# for a given interval - taking into account the epidemiological factors such as
# infectiousness and recovery probability as well as the societal connections
# created in our model).
# 
# A plot function is then built to provide clear visualisation of the progression
# of the epidemic.
# 
# These functions are then used to model, visualise and compare four scenarios in
# order to better understand the effects of applying socialbility and social
# structure to such models. (Our findings on these are present below, at the end 
# of the code).

## BUILDING THE MODEL ##

# ---- Social Structures ----
# Households #

hmax <- 5
n <- 10000

h <- sample(rep(1:n, sample(1:hmax, n, replace = TRUE))[1:n])

# Network #

# GET NET
# 
# The inputs for this function: beta, h, nc
# - beta:
# - h: 
# - nc:
# The output for this function: net_list_i
# - net_list_i


get.net <- function(beta, h, nc = 15){
  
  # Create a holding matrix to hold the nxn possible contacts
  holding_matrix = matrix(0, nrow = length(beta), ncol = length(beta))
  # Matrix of probabilities 
  # Each element is the probability of a link being created between person i and person j
  # Symmetric Matrix
  prob_matrix = (nc/((mean(beta)^2)*(length(beta)-1)))*(beta%*%t(beta))
  
  # The probability of creating a link in this way with people in the same household should be 0
  # Could remove after rbinom() but might effect the average number of non-household contacts
  family_link = outer(h, h, FUN = "==")
  prob_matrix[family_link] = 0
  
  # As probability matrix is symmetric, only need to consider the probabilities in upper triangle of the probability matrix
  # If a connection is made between person i and person j, then theres a link between person j and person i
  upper_tri_probs = prob_matrix[upper.tri(prob_matrix)]
  # Given the probabilities of contacts between individuals, use Bernoulli trial for each possible contact to see if link is created
  link_vec = rbinom(length(upper_tri_probs), size = 1, prob = upper_tri_probs)
  
  # If a link is made between person iand person j, reflect to lower triangle of matrix as there is thereby a link between person j and person i
  holding_matrix[upper.tri(holding_matrix)] = link_vec
  holding_matrix = holding_matrix + t(holding_matrix)
  # No link created with oneself
  diag(holding_matrix) = 0
  
  # Remove any links created between family members
  # family_link = outer(h, h, FUN = "==")
  # holding_matrix[family_link] = 0
  
  
  # Change to logical argument for the lapply function
  holding_matrix = holding_matrix == 1
  
  # Return the list of link indices
  net_list = as.list(as.data.frame(holding_matrix))
  net_list_i = lapply(net_list, which)
  return(net_list_i)
  
}

# ---- NSEIR Model ---- 

# NSEIR FUNCTION
# A function that takes in population and epidemiological factors and models
# the progression of the epidemic through the population for a given duration.
#
# The inputs for this function: beta, h, alink, alpha, delta, gamma, nc, nt, pinf
# - beta: Vector of socialbility parameters of the population
# -h: Vector assigning each person to a household
# -alink: List defining the regular contacts of each person (as returned by get.net)
# -alpha: Probability vector of (α_h, α_c, α_r)
#     - α_h: Probability of infecting a member of the same household
#     - α_c: Probability of infecting a member of their regular contacts network
#     - α_r: Constant factor of the probability of infection spreading randomly
#             between two people, irrespective of other relations
# -delta: Probability of moving from the infected to recovered class
# -gamma: Probability of moving from the exposed to infected class
# -nc: Average number of contacts per person
# -nt: Number of days (duration) for which to simulate the model
# -pinf: Proportion of the population that are initially infected
# The output for this function: SEIRt
# - SEIRt: A list  of vectors S, E, I, R, and t, which give the number of people
#              in each state by day, and the corresponding days.

nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, 
                  nc = 15, nt = 100, pinf = .005) {
  
  n = length(beta) # The population size
  t = 1:nt # Vector of the days of the simulation 
  
  x <- rep(0, n) # Initialising a state vector with all suseptible
  x[sample(c(1:n), round(pinf*n))] <- 1 # Moving pinf proportion to infected state
  
  # Basing initial infection based on pinf being prob, should be proportion, above correct.
  #initial_state <- c(0,2) # To start members of pop are only either suseptible or infected
  #initial_prob <- c(1-pinf, pinf)
  #x <- sample(initial_state, replace = TRUE, size = n, prob = initial_prob) # Randomly assigns initial states based on pinf
  
  daily_constant <- (alpha[3] * nc)/((mean(beta)^2)*(n-1)) # Defining constant factor of probability of random infection between two people
  
  S <- E <- I <- R <- rep(0, nt)  # Defining vectors that will contain the pop in the states on each day
  S[1] <- sum(x == 0) # Initial suseptible population is all non-infected people
  I[1] <- sum(x == 2) # Initial infected population
  
  for (i in 2:nt) { # Looping over days to update and store the state for each day
    
    u <- runif(n) # Uniform random deviates
    v <- runif(n) # Two needed to ensure that random infection is distinct from relational infection
    
    prev_infectious <- which(x == 2) # Vector of indices of infected people
    # This assumes that moving to state E occurs the day after contact with state I
    
    x[x == 2 & u < delta] <- 3 # Moves population in I to R with probability delta
    x[x == 1 & u < gamma] <- 2 # Moves population in E to I with probability delta
    
    # Suseptible people can only become exposed if there were people infected
    if (length(prev_infectious) > 0) {
      
      # Household
      infected_hh_member_count <- tabulate(h[prev_infectious], nbins = max(h)) # Counts the number of infected people per household
      # If m in a house are infected, the probabilty of a suseptible person in the house becoming exposed is '1 - (1 - α_h)^m'
      alpha_by_hh <- 1 - (1 - alpha[1])^(infected_hh_member_count)  # Defines probability of infections by each household
      alpha_hh <- alpha_by_hh[h] # Vector of probability of infection for every person by their household
      x[x == 0 & u < alpha_hh] <- 1 # Moves population in S to E by their household probability
      
      # Network
      infected_network_ids <- unlist(lapply(alink[prev_infectious], as.numeric)) # Vector of those in the network of any infected person
      # infected_network_ids <- infected_network_ids[!is.na(infected_network_ids)]
      if (length(infected_network_ids) > 0) { # Only proceeds if the infected have people in their network(s)
        infected_by_network_count <- tabulate(infected_network_ids, nbins = n) # Counts the number of infected people in each persons network
        alpha_net <- 1 - (1 - alpha[2])^(infected_by_network_count) # Defines each person's probability of infection by their network
        x[x == 0 & u < alpha_net] <- 1 # Moves population in S to E by their network probability
      }
      
      # Random
      alpha_daily <- daily_constant * outer(beta[prev_infectious], beta, FUN = "*")
      # Person i's probability of infection from random mixing is '1 - (probability of not catching through random mixing)'
      # That is, probability of not catching = (probability of not catching from j) multiplied for each j in prev_infectious
      alpha_rm <- 1 - apply(1 - alpha_daily, 2, prod) # Defines each person's probability of infection by random mixing (based on their specific alpha_daily's with those in prev_infectious)
      x[x == 0 & v < alpha_rm] <- 1 # Moves population in S to E by their random mixing probability
    }
    
    # Defining each state vector input corresponding to the day, as the number of people in that state
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
  }
  
  SEIRt <- list(S = S, E = E, I = I, R = R, t = t) # The number of people in each state by day, and the corresponding days
  SEIRt # Function returns this list
}

# ---- Model Visualisation ----


# GET NET
# 
# The inputs for this function: beta,h,alink,alpha,title
# - beta:
# - h:
# - alink:
# - alpha: 
# - title:
# The output for this function: 
# - 

plot.output <- function(beta, h, alink, alpha, title) {
  
  epi <- nseir(beta,h,alink,alpha,delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005)
  plot(epi$S,ylim=c(0,max(epi$S)),main=title,xlab="",ylab="N",las=1) ## S black
  points(epi$E,col=4);points(epi$I,col=2);points(epi$R,col=3) ## E (blue) and I (red)
  legend(x="right",legend = c("Suceptible", "Exposed", "Infected", "Recovered"),
         bty='n',
         pch=1, cex = .75,
         col = c("black","blue","red", "green"),
         text.col = c("black","blue","red", "green"))
  title(xlab="Day", line=1.75)
}


################## 
# 5: Setting beta to a vector of U(0,1) random variables, use the model to compare 4 scenarios and plot them next to each other (suitably labelled).
# First: full model with default parameters
# Second: what happens when you remove the household and regular network structure, while keeping the average initial number of infectious contacts per day the same for each person
# by setting ah=ac=0 and ar=0.04
# Third: consider the full model but with the beta vector set to simply contain the average of the previous beta vector for every element
# Fourth: combine the previous two scenarios (constant beta, random mixing)
# Comment on the apparent effect of the household and network structure relative to random mixing

# ---- Running the Model  ----

# Defining Parameters #
# ~~ Should we define h as a function and the set n and hmax down here and run everything to build the structures and models here?

n <- 10000 # Size of population to be modeled
beta <- runif(n) # Sociability vector
alink <- get.net(runif(n),h,nc = 15) # Building social networks
beta_constant <- rep(sum(beta)/length(beta), n) # Defining constant beta (eg. if social distancing enacted)

par(mfcol=c(2,2), mar=c(4,4,2,1)) # Set plot window up for multiple plots

plot.output(beta,h,alink,alpha=c(.1,.01,.01),"Full Model") # Run 1: default params
plot.output(beta,h,alink,alpha=c(0,0,.04),"Random Mixing") # Run 2: remove household and regular network structure, while keeping the average initial number of infectious contacts per day the same for each person by setting ah=ac=0 and ar=0.04
plot.output(beta = beta_constant,h,alink,alpha=c(.1,.01,.01),"Constant Beta") # Run 3: consider the full model but with the beta vector set to simply contain the average of the previous beta vector for every element
plot.output(beta = beta_constant,h,alink,alpha=c(0,0,.04),"Combined") # Run 4: combine the previous two scenarios (constant beta, random mixing)
