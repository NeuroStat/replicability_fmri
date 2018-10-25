####################
#### TITLE:  Illustration of effect of threshold on overlap measure: data generation
#### Contents: 	
#### 
#### Source Files: /replicability_fmri/Analyses/_Figures/
#### First Modified: 18/09/2018
#### Notes: 
#################


##
###############
### Notes
###############
##

# Generating images with different proportion of truly activated areas,
# different sample sizes (one-sample t-tests), different thresholds for 
# activation and calculate the observed overlap.

##
###############
### Preparation
###############
##

# Take argument from master file
input <- commandArgs(TRUE)
# K'th simulation
hpcID <- try(as.numeric(as.character(input)[1]),silent=TRUE)
# Which machine
MACHINE <- try(as.character(input)[2],silent=TRUE)
# Source the paths
source(blind_illustration_overlap_HPC.R)
# If no machine is specified, then it has to be this machine!
if(is.na(MACHINE)){
  MACHINE <- 'MAC'
  hpcID <- 1
}

# Implement for loop over r iterations here: hpcID goes from 1 to 100 in master file
rIter <- 10
startIndex <- try(1 + rIter * (hpcID - 1), silent = TRUE)
endIndex <- try(startIndex + (rIter - 1), silent = TRUE)

# Set WD: this is location where results are written
if(MACHINE=='HPC'){
  wd <- wdHPC
}
if(MACHINE=='MAC'){
  wd <- wdMAC
}

##
###############
### Libraries
###############
##

# Load in libraries
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(magrittr)
library(RColorBrewer)

##
###############
### Custom functions
###############
##

# Custom function to calculate overlap
overl_calc <- function(bin_map1, bin_map2){
  # Check if maps are only binarized
  if(!all(bin_map1 %in% c(0,1))){
    stop('Input should be binary maps only!')
  }
  if(!all(bin_map2 %in% c(0,1))){
    stop('Input should be binary maps only!')
  }
  # Sum both maps
  SumMap <- bin_map1 + bin_map2
  
  # Now calculate the number of intersecting voxels when summing
  Vab <- sum(SumMap == 2)
  Va <- sum(bin_map1 == 1)
  Vb <- sum(bin_map2 == 1)
  
  # Now calculate overlap
  overlap <- Vab/(Va + Vb - Vab)
  return(overlap)
}

# Custom function for overlap based on one threshold
overlap_fu <- function(map1, map2, threshold){
  # Check dimensions
  if(!(any(dim(map1) == dim(map2)))){
    stop('Both maps should have the same dimension!')
  }
  
  # Check threshold (should be p-value)
  if((threshold > 1) || (threshold < 0)){
    stop('Threshold should be a p-value between 0 and 1')
  }
  
  # Local function to binarize maps, based on P-value!
  bin_fu <- function(map, threshold){
    # Binary map
    bin_map <- map
    bin_map[map <= threshold] <- 1
    bin_map[map > threshold] <- 0
    
    return(bin_map)
  }
  
  # First threshold & binarize the maps
  bin_map1 <- bin_fu(map1, threshold)
  bin_map2 <- bin_fu(map2, threshold)
  
  # Proportions declared significant
  propD_map1 <- prop.table(table(bin_map1))[length(prop.table(table(bin_map1)))]
  propD_map2 <- prop.table(table(bin_map2))[length(prop.table(table(bin_map2)))]
  
  # Now calculate the overlap using overl_calc
  overlap <- overl_calc(bin_map1, bin_map2)
  
  # The end
  return(data.frame('overlap' = overlap, 'P_threshold' = threshold,
                    'propSignA' = propD_map1, 'propSignB' = propD_map2))
}

# Custom function to generate P-values
gen_data <- function(TrueGrid, TrueES, trueSigma, N){
  # Custom function to return p-values when fitting GLM
  pval_FU <- function(data){
    summary(lm(data ~ 1))$coefficients[1,4]
  }
  # Dimension of data grid
  DIM2D <- dim(TrueGrid)
  # Generate NULL DATA for N subjects for all voxels (voxel x N)
  dat <- matrix(rnorm(n = prod(DIM2D) * N, mean = 0, sd = trueSigma),
                ncol = N, nrow = prod(DIM2D))
  
  # At the true location: add true effect size ==> shift values as data is normal
  #   Note: true value = true location (0/1) * true effect size
  dat_T <- apply(dat, 2, function(x,y){x+y}, 
                 y = matrix(TrueGrid*TrueES, ncol = 1))
  # Now fit the GLM and return the p-value
  PVals <- apply(dat_T, 1, pval_FU)
  return(PVals)
}

# Custom function generating the true grid
defineTrueLoc <- function(DIM2D, baseProp){
  # Check for proportion, only 3 options
  if(!baseProp %in% c(0,0.05,0.1,0.2,0.8,1)){
    stop('Baseline proportion of activated voxels can only be 0, 5, 10, 20, 80 or 100%!')
  }
  # True location
  if(baseProp == 0){
    TrueGrid <- array(0, dim = DIM2D)
  }
  if(baseProp == 0.05){
    TrueLoc <- list(c(16:25),c(16:25))
    TrueGrid <- array(0, dim = DIM2D)
    TrueGrid[TrueLoc[[1]],TrueLoc[[2]]] <- 1
    TrueGrid[16:25,26:27] <- 1
    TrueGrid[16:20,28] <- 1
  }
  if(baseProp == 0.10){
    TrueLoc <- list(c(11:25),c(11:26))
    TrueGrid <- array(0, dim = DIM2D)
    TrueGrid[TrueLoc[[1]],TrueLoc[[2]]] <- 1
    TrueGrid[11:20,27] <- 1
  }
  if(baseProp == 0.20){
    TrueLoc <- list(c(11:32),c(11:32))
    TrueGrid <- array(0, dim = DIM2D)
    TrueGrid[TrueLoc[[1]],TrueLoc[[2]]] <- 1
    TrueGrid[11:26,33] <- 1
  }
  if(baseProp == 0.8){
    TrueLoc <- list(c(3:46),c(3:46))
    TrueGrid <- array(0, dim = DIM2D)
    TrueGrid[TrueLoc[[1]],TrueLoc[[2]]] <- 1
    TrueGrid[3:46,47] <- 1
    TrueGrid[3:22,48] <- 1
  }
  if(baseProp == 1){
    TrueGrid <- array(1, dim = DIM2D)
  }
  
  # Check if proportion is correct
  if(!baseProp %in% c(0,1)){
    if(table(TrueGrid)[2]/(table(TrueGrid)[1] + table(TrueGrid)[2]) != baseProp){
      stop('Returned proportion of truly activated voxels does not match
           user requested baseline proportion!')
    }
  }
  if(baseProp == 0){
    if(any(TrueGrid != 0)){
      stop('Returned proportion of truly activated voxels does not match
             user requested baseline proportion!')
    }
  }
  if(baseProp == 1){
    if(any(TrueGrid != 1)){
      stop('Returned proportion of truly activated voxels does not match
             user requested baseline proportion!')
    }
  }
  return(TrueGrid)
}

# Adaptive thresholding: one image only, signLevel = starting value
# baseProp is the true baseline proportion of activated voxels!
# Target is the percentage activated voxels. 
AdapThresholdingI <- function(PVals, DIM2D, signLevel, baseProp, target = NULL){
  # DIM2D should be in 2D
  if(length(DIM2D) != 2){
    stop('DIM2D should be 2 dimensional!')
  }
  
  # If no target is set, take the baseline proportion.
  if(is.null(target)){
    # Baseline proportion activated voxels is target percentage.
    #   But when proportion == 0, then set target to 0.05
    if(baseProp == 0){
      target <- 0.05
    }else{
      target <- baseProp
    }
  }
  
  # Threshold at signLevel p-value: significant P-values get 1!
  idP <- PVals <= signLevel
  threshPval <- PVals
  threshPval[!idP] <- 0
  threshPval[idP] <- 1
  
  # Calculate starting (base) percentage of activated masked voxels
  Vt <- sum(idP,na.rm=TRUE)
  percentage <- Vt/prod(DIM2D)
  
  # Now if percentage is below (target - 0.01), increase threshold (by .001), 
  #     if above (target + 0.01); decrease by .01
  while(percentage < (target - 0.01) | percentage > (target + 0.01)){
    if(percentage < (target - 0.01)){
      if(signLevel >= 0.995){
        signLevel <- signLevel + 0.00001
      }else{
        signLevel <- signLevel + runif(1,min=.0001,max=0.005)
      }
    }
    if(percentage > (target + 0.01)){
      if(signLevel <= 0.005){
        signLevel <- signLevel - 0.00001
      }else{
        signLevel <- signLevel - runif(1,min=.0001,max=0.005)
      }
    }
    # Recalculate thresholded map and percentage
    idP <- PVals <= signLevel
    threshPval <- PVals
    threshPval[!idP] <- 0
    threshPval[idP] <- 1
    Vt <- sum(idP,na.rm=TRUE)
    percentage <- Vt/prod(DIM2D)
  }
  
  # Return percentage, signLevel and threshPval
  result <- list('signLevel' = signLevel, 'percentage' = percentage, 'SPM' = threshPval)
  return(result)
}


##
###############
### Preparation
###############
##

# Dimension of the 2D slice
DIM2D <- c(50,50)

# Define the true effect size
TrueES <- 1

# Sigma of data generated
trueSigma <- 2

# Proportion of image truly activated
baseProp <- c(0,0.05,0.20,0.80,1)

# Thresholds considered
PVal_thr <- seq(0.1, 0.0001, by = -0.01)

# Target percentages considered in the adaptive thresholding procedure.
#   For baseProp of 0, 0.05 and 1, we take the baseProp as the target.
#   For baseProp of .20 and .80, we also have two extra targets:
#   10% higher and 10% below as baseline
targets <- c(0, -.10, .10)

# Number of subjects in the analysis
nsub <- seq(10, 100, by = 10)

##
###############
### Data generation
###############
##

# Empty data frame
sim_data <- data.frame() %>% as_tibble()

# For loop over the simulations
for(ID in startIndex:endIndex){
  # Set seed
  set.seed(ID*pi)
  
  # For loop over the proportion of voxels with true activation
  for(l in 1:length(baseProp)){
    # Generate the true grid for this scenario
    TrueGrid_S <- defineTrueLoc(DIM2D, baseProp[l])
    # For loop over the amount of subjects
    for(s in 1:length(nsub)){
      # Generate p value data for both maps and bind to each other
      PVals <- data.frame(PVal1 = gen_data(TrueGrid_S, TrueES, trueSigma, N = nsub[s]),
                          PVal2 = gen_data(TrueGrid_S, TrueES, trueSigma, N = nsub[s]))
      
      ##########################################################################
      # First approach: for loop over the thresholds
      ##########################################################################
      for(j in 1:length(PVal_thr)){
        # Now send them through the overlap function with one threshold
        sim_data <- overlap_fu(map1 = array(PVals$PVal1, dim = DIM2D), 
                               map2 = array(PVals$PVal2, dim= DIM2D), 
                               threshold = PVal_thr[j]) %>%
          # Add info (no target percentage activated voxels in this case)
          mutate(N = nsub[s],
                 BaseProp = baseProp[l],
                 sigma = trueSigma,
                 target = NA,
                 sim = ID,
                 approach = 'normal') %>%
          # Bind to data frame
          bind_rows(sim_data, .)
      }

      ##########################################################################
      # Second approach: adaptive thresholding.
      ##########################################################################
      # For loop over the targets 
      for(j in 1:length(targets)){
        # Note: when baseProp %in% c(0, 0.05 or 1), then only have one target!
        if(baseProp[l] %in% c(0,0.05,1) && j > 1) next
        
        # Define target of this iteration: note, when baseline = 0 or 1, then
        # have either 0.05 or 0.95 (otherwise no voxels are ever selected),
        # or all
        if(baseProp[l] == 0){
          target_j <- 0.05
        }
        if(baseProp[l] == 1){
          target_j <- 0.95
        }
        if(!baseProp[l] %in% c(0,1)){
          target_j <- baseProp[l] + targets[j]
        }
          
        # Using the data, calculate the appropriate threshold.
        PVal_thr1 <- AdapThresholdingI(PVals = PVals$PVal1, DIM2D = DIM2D, 
                                       signLevel = 0.05, baseProp = baseProp[l],
                                       target = target_j)
        PVal_thr2 <- AdapThresholdingI(PVals = PVals$PVal2, DIM2D = DIM2D, 
                                       signLevel = 0.05, baseProp = baseProp[l],
                                       target = target_j)
        
        # Now calculate the overlap: threshold is the average used threshold in both maps
        sim_data <- data.frame('overlap' = overl_calc(PVal_thr1$SPM, PVal_thr2$SPM),
                               'P_threshold' = mean(PVal_thr1$signLevel,
                                                    PVal_thr2$signLevel),
                               'propSignA' = PVal_thr1$percentage,
                               'propSignB' = PVal_thr2$percentage) %>%
          # Add info
          mutate(N = nsub[s],
                 BaseProp = baseProp[l],
                 sigma = trueSigma,
                 target = target_j,
                 sim = ID,
                 approach = 'adaptive') %>%
          # Bind to data frame
          bind_rows(sim_data, .)
      }
    }
  }
}




##
###############
### Save results
###############
##



# Save R object
saveRDS(sim_data, file = paste0(wd, '/Results/illus_overl_', hpcID, '.rds'))




















