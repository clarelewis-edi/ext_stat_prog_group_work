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
n <- 1000

h <- sample(rep(1:n, sample(1:hmax, n, replace = TRUE))[1:n])

################## 
# 2: Write function get.net(beta,nc=15) where beta is the n vector of Bi value for each person. The function should return a list,
# the ith element of which is a vector of indices of the regular (non-household) contacts of a person i. Must be careful to implement
# the model properly so that you do not CREATE any links twice (links need to be RECORDED twice, however)
################## 


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


################## 
# 4: Write a function to nicely plot the dynamics of the simulated population states as returned by nseir (a better version of the plots in 6.2 in the notes)
################## 


################## 
# 5: Setting beta to a vector of U(0,1) random variables, use the model to compare 4 scenarios and plot them next to each other (suitably labelled).
# First: full model with default parameters
# Second: what happens when you remove the household and regular network structure, while keeping the average initial number of infectious contacts per day the same for each person
# by setting ah=ac=0 and ar=0.04
# Third: consider the full model but with the beta vector set to simply contain the average of the previous beta vector for every element
# Fourth: combine the previous two scenarios (constant beta, random mixing)
# Comment on the apparent effect of the household and network structure relative to random mixing
################## 


########### code given in section 6.2 that is referenced in the assignment instructions
seir <- function(n=10000,ni=10,nt=100,gamma=1/3,delta=1/5,bmu=5e-5,bsc=1e-5) {
  ## SEIR stochastic simulation model.
  ## n = population size; ni = initially infective; nt = number of days
  ## gamma = daily prob E -> I; delta = daily prob I -> R;
  ## bmu = mean beta; bsc = var(beta) = bmu * bsc
  x <- rep(0,n) ## initialize to susceptible state
  beta <- rgamma(n,shape=bmu/bsc,scale=bsc) ## individual infection rates
  x[1:ni] <- 2 ## create some infectives
  S <- E <- I <- R <- rep(0,nt) ## set up storage for pop in each state
  S[1] <- n-ni;I[1] <- ni ## initialize
  for (i in 2:nt) { ## loop over days
    u <- runif(n) ## uniform random deviates
    x[x==2&u<delta] <- 3 ## I -> R with prob delta
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma
    x[x==0&u<beta*I[i-1]] <- 1 ## S -> E with prob beta*I[i-1]
    S[i] <- sum(x==0); E[i] <- sum(x==1)
    I[i] <- sum(x==2); R[i] <- sum(x==3)
  }
  list(S=S,E=E,I=I,R=R,beta=beta)
} ## seir

par(mfcol=c(2,3),mar=c(4,4,1,1)) ## set plot window up for multiple plots
epi <- seir(bmu=7e-5,bsc=1e-7) ## run simulation
hist(epi$beta,xlab="beta",main="") ## beta distribution
plot(epi$S,ylim=c(0,max(epi$S)),xlab="day",ylab="N") ## S black
points(epi$E,col=4);points(epi$I,col=2) ## E (blue) and I (red)
