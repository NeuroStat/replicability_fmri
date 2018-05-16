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
nstep <- 1000

# Number of voxels for testing purpose
nvox <- 100

# Number of independent maps (R) for testing purpose
R <- 10

######################
## Custom functions ##
######################

# Define function weights (using current step in algorithm)
weight <- function(lambda, pi_A0, pi_I0, R, Gi){
  # First calculate numinator
  num <- lambda * (pi_A0^(R - Gi)*(1 - pi_A0)^Gi)
  # Now calculate denominator
  denom <- num + (1 - lambda)*(pi_I0^(R - Gi)*(1 - pi_I0)^Gi)
  # Return value
  return(num/denom)
}

# Maximalization function that returns updated 
#   value for pi_A0, given current values (_c)
update_pi_A <- function(lambda_c, pi_A0_c, pi_I0_c, R, Gi){
  # First calculate numinator
  num <- sum((weight(lambda_c, pi_A0_c, pi_I0_c, R, Gi = Gi) * R) -
                (weight(lambda_c, pi_A0_c, pi_I0_c, R, Gi = Gi) * Gi))
  # Now calculate denominator
  denom <- sum((weight(lambda_c, pi_A0_c, pi_I0_c, R, Gi = Gi) * R))
  # Return value
  return(num/denom)
}

# Same for pi_I0
update_pi_I <- function(lambda_c, pi_A0_c, pi_I0_c, R, Gi){
  # First calculate numinator
  num <- sum(((1 - weight(lambda_c, pi_A0_c, pi_I0_c, R, Gi = Gi)) * R) -
               ((1 - weight(lambda_c, pi_A0_c, pi_I0_c, R, Gi = Gi)) * Gi))
  # Now calculate denominator
  denom <- sum(((1 - weight(lambda_c, pi_A0_c, pi_I0_c, R, Gi = Gi)) * R))
  # Return value
  return(num/denom)
}

# Same for lambda
update_lambda <- function(lambda_c, pi_A0_c, pi_I0_c, R, Gi){
  # Note: length(Gi) equals amount of voxels (N)!
  sum(weight(lambda_c, pi_A0_c, pi_I0_c, R, Gi)) / length(Gi)
  #mean(weight(lambda_c, pi_A0_c, pi_I0_c, R, Gi))
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
    lambda <- 0.05
    pi_A0 <- 0.9
    pi_I0 <- 0.9
  }

  ################################
  ## 2: Conditional expectation ##
  ################################
  
  # See custom functions and derivations in PDF

  #####################################
  ## 3: Maximise (update parameters) ##
  #####################################
  
  # Parameter: pi_A0
  pi_A_up <- update_pi_A(lambda, pi_A0, pi_I0, R, TDat_FC)
  # Parameter: pi_I0
  pi_I_up <- update_pi_I(lambda, pi_A0, pi_I0, R, TDat_FC)
  # Paramater: lambda
  lambda_up <- update_lambda(lambda, pi_A0, pi_I0, R, TDat_FC)
  
  
  # Now we replace the parameters
  lambda <- lambda_up
  pi_A0 <- pi_A_up
  pi_I0 <- pi_I_up
  
  cat( "Iteration:", v, "| PI A0 =", pi_A0, "| PI I0 =", pi_I0, "| Lambda =", lambda,  "\n" )
}



# Using density function of binomial distribution
lambda <- 0.05
pi_A0 <- 0.1
pi_I0 <- 0.1

for( i in 1:nstep ) {
  
  T_1 = lambda * dbinom(TDat_FC, size = R, prob = pi_A0)
  T_2 = (1 - lambda) * dbinom(TDat_FC, size = R, prob = pi_I0)
  
  Weight_1 = T_1 / (T_1 + T_2)
  Weight_2 = T_2 / (T_1 + T_2)  #The same as 1-Weight_1
  
  p_1 = mean(Weight_1)
  p_2 = mean(Weight_2)
  
  pi_A0 <- sum(Weight_1 * TDat_FC ) / sum(Weight_1)
  pi_I0 <- sum(Weight_2 * TDat_FC ) / sum(Weight_2)
  
  sd_1 = sqrt(sum( Weight_1 * (x-mu_1)^2) / sum(Weight_1))
  sd_2 = sqrt(sum( Weight_2 * (x-mu_2)^2) / sum(Weight_2))
  
  
  ## print the current estimates
  
  cat( "Iteration: ", i, "| Mean", mu_1, mu_2, "| St. Dev", sd_1, sd_2,"| Prob", p_1, p_2,  "\n" )
  
}

























