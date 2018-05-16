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

# Number of steps in algorithm for testing purpose
nstep <- 100

# Number of voxels for testing purpose
nvox <- 100

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
  # First calculate numinator
  num <- sum((weight(lambda_c, pi_A1_c, pi_I1_c, R, Gi = Gi) * Gi))
  # Now calculate denominator
  denom <- sum((weight(lambda_c, pi_A1_c, pi_I1_c, R, Gi = Gi) * R))
  # Return value
  return(num/denom)
}

# Same for pi_I1
update_pi_I <- function(lambda_c, pi_A1_c, pi_I1_c, R, Gi){
  # First calculate numinator
  num <- sum(((1 - weight(lambda_c, pi_A1_c, pi_I1_c, R, Gi = Gi)) * Gi))
  # Now calculate denominator
  denom <- sum(((1 - weight(lambda_c, pi_A1_c, pi_I1_c, R, Gi = Gi)) * R))
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

# Test weight function
weight(0.8,0.2,0.7,8,1:8)

# Create test data
set.seed(pi*3.14)

# No coherence if sampled at random (?)
TDat_R <- sample(x = 1:R, size = nvox, replace = TRUE)

# many null voxels, some random and some coherence (partial coherence)
TDat_PC <- c(rep(0, 60), sample(x = 1:R, size = 20, replace = TRUE), rep(9,20))

# data with only 0 and R (full coherence)
TDat_FC <- sample(c(0, R), size = nvox, replace = TRUE)

# Data with only R
TDat_FR <- rep(R, nvox)

# Data with only 0
TDat_F0 <- rep(0, nvox)

# Start for loop of algorithm steps
for(v in 1:nstep){
  ################################
  ## 1: Initial starting values ##
  ################################
  if(v == 1){
    lambda <- 0.5
    pi_A1 <- 0.1
    pi_I1 <- 0.1
  }

  ################################
  ## 2: Conditional expectation ##
  ################################
  
  # See custom functions and derivations in PDF

  #####################################
  ## 3: Maximise (update parameters) ##
  #####################################
  
  # Parameter: pi_A1
  pi_A_up <- update_pi_A(lambda, pi_A1, pi_I1, R, TDat_FC)
  # Parameter: pi_I1
  pi_I_up <- update_pi_I(lambda, pi_A1, pi_I1, R, TDat_FC)
  # Paramater: lambda
  lambda_up <- update_lambda(lambda, pi_A1, pi_I1, R, TDat_FC)
  
  # Now we replace the parameters
  lambda <- lambda_up
  pi_A1 <- pi_A_up
  pi_I1 <- pi_I_up
  
  cat( "Iteration:", v, "| PI A1 =", pi_A1, "| PI I1 =", pi_I1, "| Lambda =", lambda,  "\n" )
}


























