# ------------------------------------------------------------------------------ 
# Clare Lewis (s2879721), Grace Sheahan (s2898645), Luke Egan (s2837709)

# Clare: 
# Grace: 
# Luke: 

################################################################################ Comments should not be longer than this line of ##
# ------------------------------------------------------------------------------

# Setting working directory to folder for Assessment 1
# setwd("/Users/clarelewis/Documents/github/ext_stat_prog_group_work/Proj 2")
# setwd("C:/Users/Grace Sheahan/Documents/Extended Statistical Programming/ext_stat_prog_group_work/Proj 2")
# setwd("C:/Users/Luke Egan/Desktop/Extended Statistical Programming/ext_stat_prog_group_work/Proj 2")

#-------------------------------------------------------------------------------
### Background provided info
#-------------------------------------------------------------------------------

################## 
# There are n people in each state S, E, I, R.
# If in state I they have a daily probability of delta of moving to state R.
# If in state E they have a daily probability of gamma of moving to state I.
# If in state S the probability of moving to E is the result of infection by someone in state I.
#   This depends on two characteristics of each person:
#     who they share a household with
#     who is in their network of regular contacts
#   These are modelled as follows:
################## 

# ################## 
# # model household characteristic by dividing the n people into households of sizes between 1 and hmax, assume a uniform distribution of household sizes
# ################## 
# 
# 
# ################## 
# # contact network model: sociability parameter: Bi where person is i, network is created randomly by assigning a link between person i and j
# # with probability: (nc*Bi*Bj)/(Bmean^2*(n-1)) where nc is average number of contacts per person an Bmean is mean social parameter. People in the same household are excluded from such contacts
# ################## 


# There are several ways for person i in the state I to infect a person j in state S
# ##################
# # Way 1: if j is a household member, there is a dialy probability ah of i infecting j
# ##################
# 
# 
# ##################
# # Way 2: if j is in i's regular network of contacts, then there is a daily probability ac of i infecting j
# ##################
# 
# 
# ##################
# # Way 3: irrespective of household or regular network relations there is a daily probability: (ar*nc*Bi*Bj)/(Bmean^2*(n-1)) of i infecting j (random mixing)
# ##################

# # note: the use of the same nc in the random mixing and regular network expressions is a simplification

#-------------------------------------------------------------------------------
################## 
# Assignment Goal: impliment the model and provide illustration of its use to investigate the role of household and network structure on epidemic dynamics.
# code work with any pop size: n up to at least 10,000, but we should test and develop with n = 1000
################## 

#-------------------------------------------------------------------------------
#### Task steps
#-------------------------------------------------------------------------------

################## 
# 1: Write code to produce vector h (should be length n) of integers indicating which household each person belongs to
# Household sizes should be uniformly distributed between 1 and hmax. For example if h1, h56 and h907 are all 13 and these are the only
# 13s in h then people 1, 56, and 907 all live in the same 3 person household
# h can be created with one line of code by careful use of rep and sample. Use hmax = 5 by default
################## 

set.seed(17)

hmax <- 5
n <- 10000

h <- sample(rep(1:n, sample(1:hmax, n, replace = TRUE))[1:n])

################## 
# 2: Write function get.net(beta,nc=15) where beta is the n vector of Bi value for each person. The function should return a list,
# the ith element of which is a vector of indices of the regular (non-household) contacts of a person i. Must be careful to implement
# the model properly so that you do not CREATE any links twice (links need to be RECORDED twice, however)
################## 

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


################## 
# 3: Write the function:
# nseir(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005)
# where beta and h are as above, alink = list defining the regular contacts of each person as returned by get.net,
# nt = number of days to simulate, and pinf = proportion of initial population to randomly start in the I state
# The function should implement the model given above and return a list with elements S,E,I,R and t = total population in each class
# each day and the day, respectively
# Note: while looping through individuals in the I state is probably inevitable, the code can still be made to run in a few seconds for n = 10000,
# especially with careful use of expressions like x[ind1][ind2] <- y in places (assignment to a subvector of a subvector)
################## 

# nseir <- function(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005){
#   ## SEIR stochastic simulation model (S,E,I,R are states 0,1,2,3)
#   # b, h, alink as defined above
#   # alpha = vector of transmition probabilities by relationship
#   # delta = probability of moving from I to R
#   # gamma = probability of moving from E to I
#   # nc = average number of contacts per person
#   # nt = number of days to simulate
#   # proportion of the initial population to randomly start in the I state.
#   
#   t = 1:nt
#   
#   initial_state <- c(0,2) # To start members of pop are only either suseptible or infected
#   initial_prob <- c(1-pinf, pinf)
#   x <- sample(initial_state, replace = TRUE, size = n, prob = initial_prob) # Randomly assigns initial states based on pinf
#   
#   beta_bar <- sum(beta)/n
#   daily_constant <- alpha[3]*nc/(beta_bar^2*(n-1))
#   
#   S <- E <- I <- R <- rep(0, nt)  # Defining vectors that will contain the pop in the states on each day
#   S[1] <- sum(x == 0) # Initial suseptible population is all non-infected people
#   I[1] <- sum(x == 2) # Initial infected population
#   
#   for (i in 2:nt) {
#     u <- runif(n)
#     prev_infectious <- which(x == 2)
#     
#     x[x == 2 & u < delta] <- 3
#     x[x == 1 & u < gamma] <- 2
#     
#     e_if_s <- rep(FALSE, n)
#     for (j in prev_infectious){
#       # v = runif(n)
#       
#       xs <- x
#       
#       household <- rep(0, n)
#       household[which(h == h[j])] <- 1
#       xh <- x
#       xh[runif(n) < alpha[1]*household] <- 1
#       xhs <- which(xh == 1)
#       
#       network <- rep(0, n)
#       network[alink[[j]]] <- 1
#       xn <- x
#       xn[runif(n) < alpha[2]*network] <- 1
#       xns <- which(xn == 1)
#       
#       
#       xo <- x
#       alpha_daily <- daily_constant*beta[j]*beta
#       
#       xo[runif(n) < alpha_daily] <- 1
#       xos <- which(xo == 1)
#       
#       infected_by_j <- unique(c(xhs, xns, xos))
#       
#       e_if_s[infected_by_j] <- (TRUE)
#     }
#     daily_infectees <- c()
#     x[x == 0 & e_if_s] <- 1
#     
#     S[i] <- sum(x == 0)
#     E[i] <- sum(x == 1)
#     I[i] <- sum(x == 2)
#     R[i] <- sum(x == 3)
#   }
#   
#   list(S=S, E=E, I=I, R=R, t=t)
# }

#-----------------------------
nseir <- function(beta, h, alink, alpha = c(.1, .01, .01), delta = .2, gamma = .4, 
                  nc = 15, nt = 100, pinf = .005) {
  
  n = length(beta)
  t = 1:nt
  
  initial_state <- c(0,2) # To start members of pop are only either suseptible or infected
  initial_prob <- c(1-pinf, pinf)
  x <- sample(initial_state, replace = TRUE, size = n, prob = initial_prob) # Randomly assigns initial states based on pinf
  
  beta_bar <- mean(beta)
  daily_constant <- (alpha[3] * nc)/((beta_bar^2)*(n-1))
  
  S <- E <- I <- R <- rep(0, nt)  # Defining vectors that will contain the pop in the states on each day
  S[1] <- sum(x == 0) # Initial suseptible population is all non-infected people
  I[1] <- sum(x == 2) # Initial infected population
  
  for (i in 2:nt) {
    u <- runif(n)
    v <- runif(n)
    prev_infectious <- which(x == 2)
    
    x[x == 2 & u < delta] <- 3
    x[x == 1 & u < gamma] <- 2
    
    if (length(prev_infectious) > 0) {
      
      # Household
      infected_hh_member_count <- tabulate(h[prev_infectious], nbins = max(h))
      alpha_by_hh <- 1 - (1 - alpha[1])^(infected_hh_member_count)
      alpha_hh <- alpha_by_hh[h]
      x[x == 0 & u < alpha_hh] <- 1
      
      # Network
      infected_network_ids <- unlist(lapply(alink[prev_infectious], as.numeric))
      # infected_network_ids <- infected_network_ids[!is.na(infected_network_ids)]
      if (length(infected_network_ids) > 0) {
        infected_by_network_count <- tabulate(infected_network_ids, nbins = n)
        alpha_net <- 1 - (1 - alpha[2])^(infected_by_network_count)
        x[x == 0 & u < alpha_net] <- 1
      }
      
      # Random
      alpha_daily <- daily_constant * outer(beta[prev_infectious], beta, FUN = "*")
      alpha_rm <- 1 - apply(1 - alpha_daily, 2, prod)
      x[x == 0 & v < alpha_rm] <- 1
    }
    
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
  }
  
  list(S = S, E = E, I = I, R = R, t = t)
}


epi = nseir(beta = runif(10000),h,alink = get.net(runif(10000),h,nc = 15),alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005)
epi


################## 
# 4: Write a function to nicely plot the dynamics of the simulated population states as returned by nseir (a better version of the plots in 6.2 in the notes)
################## 
par(mfcol=c(2,2),mar=c(4,4,1,1)) # set plot window up for multiple plots

plot.output <- function(beta,h,alink,title) {
  
  epi <- nseir(beta,h,alink,alpha=c(.1,.01,.01),delta=.2,gamma=.4,nc=15, nt = 100,pinf = .005)
  title = "title"
  plot(epi$S,ylim=c(0,max(epi$S)),main=title,xlab="day",ylab="N",las=1) ## S black
  points(epi$E,col=4);points(epi$I,col=2);points(epi$R,col=3) ## E (blue) and I (red)
  legend(x="right",legend = c("Suceptible", "Exposed", "Infected", "Recovered"),
         bty='n',
         pch=1, cex = .75,
         col = c("black","blue","red", "green"),
         text.col = c("black","blue","red", "green"))
  
}

# full model, random mixing, constant beta, random mixing and constant beta
################## 
# 5: Setting beta to a vector of U(0,1) random variables, use the model to compare 4 scenarios and plot them next to each other (suitably labelled).
# First: full model with default parameters
# Second: what happens when you remove the household and regular network structure, while keeping the average initial number of infectious contacts per day the same for each person
# by setting ah=ac=0 and ar=0.04
# Third: consider the full model but with the beta vector set to simply contain the average of the previous beta vector for every element
# Fourth: combine the previous two scenarios (constant beta, random mixing)
# Comment on the apparent effect of the household and network structure relative to random mixing
################## 

plot.output(beta = runif(10000),h,alink = get.net(runif(10000),h,nc = 15),"Full Model") #run 1: default params
plot.output(beta = runif(10000),h,alink = get.net(runif(10000),h,nc = 15),"Random Mixing") #run 2: remove household and regular network structure, while keeping the average initial number of infectious contacts per day the same for each person by setting ah=ac=0 and ar=0.04
plot.output(beta = runif(10000),h,alink = get.net(runif(10000),h,nc = 15),"Constant Beta") #run 3: consider the full model but with the beta vector set to simply contain the average of the previous beta vector for every element
plot.output(beta = runif(10000),h,alink = get.net(runif(10000),h,nc = 15),"Combined") #run 4: combine the previous two scenarios (constant beta, random mixing)


########### code given in section 6.2 that is referenced in the assignment instructions
# seir <- function(n=10000,ni=10,nt=100,gamma=1/3,delta=1/5,bmu=5e-5,bsc=1e-5) {
#   ## SEIR stochastic simulation model.
#   ## n = population size; ni = initially infective; nt = number of days
#   ## gamma = daily prob E -> I; delta = daily prob I -> R;
#   ## bmu = mean beta; bsc = var(beta) = bmu * bsc
#   x <- rep(0,n) ## initialize to susceptible state
#   beta <- rgamma(n,shape=bmu/bsc,scale=bsc) ## individual infection rates
#   x[1:ni] <- 2 ## create some infectives
#   S <- E <- I <- R <- rep(0,nt) ## set up storage for pop in each state
#   S[1] <- n-ni;I[1] <- ni ## initialize
#   for (i in 2:nt) { ## loop over days
#     u <- runif(n) ## uniform random deviates
#     x[x==2&u<delta] <- 3 ## I -> R with prob delta
#     x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
#     x[x==0&u<beta*I[i-1]] <- 1 ## S -> E with prob beta*I[i-1]
#     S[i] <- sum(x==0); E[i] <- sum(x==1)
#     I[i] <- sum(x==2); R[i] <- sum(x==3)
#   }
#   list(S=S,E=E,I=I,R=R,beta=beta)
# } ## seir
# 
# par(mfcol=c(2,3),mar=c(4,4,1,1)) ## set plot window up for multiple plots
# epi <- seir(bmu=7e-5,bsc=1e-7) ## run simulation
# hist(epi$beta,xlab="beta",main="") ## beta distribution
# plot(epi$S,ylim=c(0,max(epi$S)),xlab="day",ylab="N") ## S black
# points(epi$E,col=4);points(epi$I,col=2) ## E (blue) and I (red)