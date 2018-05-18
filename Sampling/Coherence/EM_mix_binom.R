####################
#### TITLE:     Estimate parameters of mixture binomial using EM
#### Contents: 	
#### 
#### Source Files: /FreddieFreeloader/Sampling/Coherence
#### First Modified: 15/05/2018
#### Notes: 
#################



##
###############
### Notes
###############
##

# Check calculations of EM algorithm to estimate mixture of two binomials


##
###############
### Preparation
###############
##

# Libraries
library(dplyr)

# Number of steps in algorithm for testing purpose
nstep <- 100

# Number of voxels for testing purpose
nvox <- 100000

# Number of independent maps (R) for testing purpose
R <- 10

######################
## Custom functions ##
######################

# Define function weights (using current step in algorithm)
weight <- function(lambda, pi_A1, pi_I1, R, Gi){
  # First calculate numinator: lambda times density of binomial
  num <- lambda * dbinom(x = Gi, size = R, prob = pi_A1)
  # Now calculate denominator
  denom <- num + ((1 - lambda) * dbinom(x = Gi, size = R, prob = pi_I1))
  # Return value
  return(num/denom)
}

# Maximalization function that returns updated 
#   value for pi_A1, given current values (_c)
update_pi_A <- function(lambda_c, pi_A1_c, pi_I1_c, R, Gi){
  # Weight of current step
  weight_c <- weight(lambda_c, pi_A1_c, pi_I1_c, R, Gi = Gi)
  # First calculate numinator
  num <- sum(weight_c * Gi)
  # Now calculate denominator
  denom <- sum(weight_c * R)
  # Return value
  return(num/denom)
}

# Same for pi_I1
update_pi_I <- function(lambda_c, pi_A1_c, pi_I1_c, R, Gi){
  # Weight of current step
  weight_c <- weight(lambda_c, pi_A1_c, pi_I1_c, R, Gi = Gi)
  # First calculate numinator
  num <- sum(((1 - weight_c) * Gi))
  # Now calculate denominator
  denom <- sum(((1 - weight_c) * R))
  # Return value
  return(num/denom)
}

# Same for lambda
update_lambda <- function(lambda_c, pi_A1_c, pi_I1_c, R, Gi){
  mean(weight(lambda_c, pi_A1_c, pi_I1_c, R, Gi))
}


##
###############
### Algorithm
###############
##

# Create test data
set.seed(pi*3.14)

## Data as mixture of binomials
true_PIA1 <- 0.75
true_PII1 <- 0.2
true_lambda <- 0.25

# Small simulation 
nsim <- 100

sim_PIA <- sim_PII <- sim_lambda <- c()

for(i in 1:nsim){
  # Create data
  # First generate the mixing parameter (Z)
  TDat_mix <- data.frame(Z = rbinom(n = nvox, size = 1, prob = true_lambda)) %>%
    # Mutate column with data, depending on mixing parameter
    mutate(G = ifelse(Z == 1,
              rbinom(n = 1, size = R, prob = true_PIA1),
              rbinom(n = 1, size = R, prob = true_PII1)))
  
  # Start for loop of algorithm steps
  for(v in 1:nstep){
    ################################
    ## 1: Initial starting values ##
    ################################
    if(v == 1){
      lambda <- true_lambda
      pi_A1 <- true_PIA1
      pi_I1 <- true_PII1
      #cat( "Iteration: 0", "| PI A1 =", pi_A1, "| PI I0 =", pi_I0, "| Lambda =", lambda,  "\n" )
    }
  
    ################################
    ## 2: Conditional expectation ##
    ################################
    
    # See custom functions and derivations in PDF
  
    #####################################
    ## 3: Maximise (update parameters) ##
    #####################################
    
    # Parameter: pi_A1
    pi_A_up <- update_pi_A(lambda, pi_A1, pi_I1, R, TDat_mix$G)
    # Parameter: pi_I1
    pi_I_up <- update_pi_I(lambda, pi_A1, pi_I1, R, TDat_mix$G)
    # Paramater: lambda
    lambda_up <- update_lambda(lambda, pi_A1, pi_I1, R, TDat_mix$G)
    
    # Now we replace the parameters
    lambda <- lambda_up
    pi_A1 <- pi_A_up
    pi_I1 <- pi_I_up
    
    #cat( "Iteration:", v, "| PI A1 =", pi_A1, "| PI I0 =", pi_I0, "| Lambda =", lambda,  "\n" )
  }

  # Save parameters
  sim_PIA <- c(sim_PIA, pi_A1)
  sim_PII <- c(sim_PII, pi_I1)
  sim_lambda <- c(sim_lambda, lambda)
}

true_PIA1
mean(sim_PIA)

true_PII1
mean(sim_PII)

true_lambda
mean(sim_lambda)


# Test function for R package
EM_binom <- function(Y, N, inLam, inPI1, inPI2, max.iter = 250){
  
  ##################
  # Local functions
  ##################
  
  # Define function weights (using current step in algorithm)
  weight <- function(lambda, pi_1, pi_2, N, Yi){
    # First calculate numinator: lambda times density of binomial
    num <- lambda * dbinom(x = Yi, size = N, prob = pi_1)
    # Now calculate denominator
    denom <- num + ((1 - lambda) * dbinom(x = Yi, size = N, prob = pi_2))
    # Return value
    return(num/denom)
  }
  
  # Maximalization function that returns updated
  #   value for pi_1, given current values (_c)
  update_pi_1 <- function(lambda_c, pi_1_c, pi_2_c, N, Yi){
    # Weight of current step
    weight_c <- weight(lambda_c, pi_1_c, pi_2_c, N, Yi = Yi)
    # First calculate numinator
    num <- sum(weight_c * Yi)
    # Now calculate denominator
    denom <- sum(weight_c * N)
    # Return value
    return(num/denom)
  }
  
  # Same for pi_2
  update_pi_2 <- function(lambda_c, pi_1_c, pi_2_c, N, Yi){
    # Weight of current step
    weight_c <- weight(lambda_c, pi_1_c, pi_2_c, R, Yi = Yi)
    # First calculate numinator
    num <- sum(((1 - weight_c) * Yi))
    # Now calculate denominator
    denom <- sum(((1 - weight_c) * N))
    # Return value
    return(num/denom)
  }
  
  # Same for lambda
  update_lambda <- function(lambda_c, pi_1_c, pi_2_c, N, Yi){
    mean(weight(lambda_c, pi_1_c, pi_2_c, N, Yi))
  }
  
  ##################
  # EM ALGORITHM  #
  ##################
  
  # Start for loop of algorithm steps
  for(v in 1:max.iter){
    ################################
    ## 1: Initial starting values ##
    ################################
    if(v == 1){
      lambda <- inLam
      pi_1 <- inPI1
      pi_2 <- inPI2
    }
    
    ################################
    ## 2: Conditional expectation ##
    ################################
    
    # See custom functions
    
    #####################################
    ## 3: Maximise (update parameters) ##
    #####################################
    
    # Parameter: pi_A1
    pi_1_up <- update_pi_1(lambda, pi_1, pi_2, N, Y)
    # Parameter: pi_I1
    pi_2_up <- update_pi_2(lambda, pi_1, pi_2, N, Y)
    # Paramater: lambda
    lambda_up <- update_lambda(lambda, pi_1, pi_2, N, Y)
    
    # Now we replace the parameters
    lambda <- lambda_up
    pi_1 <- pi_1_up
    pi_2 <- pi_2_up
  }
  
  # data frame with results
  res <- data.frame('lambda' = lambda, 'PI1' = pi_1, 'PI2' = pi_2)
  return(res)
}


EM_binom(Y = TDat_mix$G, N = R, inLam = true_lambda, inPI1 = true_PIA1, inPI2 = true_PII1, max.iter = 100)








