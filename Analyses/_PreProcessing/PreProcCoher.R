####################
#### TITLE:     Calculate coherence, after fitting a mixture of binomials using EM.
#### Contents: 	
#### 
#### Source Files: /FreddieFreeloader/Analyses/_PreProcessing/_Coherence
#### First Modified: 29/05/2018
#### Notes: 
#################



##
###############
### Notes
###############
##

# Data analysis:
# Make independent batches of subjects.
# Each run, we estimate the mixture of two binomials on the summed thresholded maps within a sample size.

# Here we read in the summed maps and then estimate the mixture using an EM algorithm.

# All contrasts can be selected here

##
###############
### Preparation
###############
##

# Source paths
source('blind_PreProcessing.R')

# Possible thresholding scenario
scenario_pos <- c('uncorrected', 'fdr')
scenario <- scenario_pos[2]

# Possible contrasts
contrast <- c('ML', 'Faces', 'Incentive', 'StopGo')
contr <- contrast[4]

# Location of raw data: not included in Github folder (too large)
# Also take the correct contrast
if(contr == 'ML'){
  if(scenario == 'uncorrected'){
    RawDat <- RawDatCohUnc
  }
  if(scenario == 'fdr'){
    RawDat <- RawDatCohFDR
  }
}
if(contr == 'Faces'){
  RawDat <- RawDatCohFDRF
}
if(contr == 'Incentive'){
  RawDat <- RawDatCohFDRI
}
if(contr == 'StopGo'){
  RawDat <- RawDatCohFDRS
}

# Load in libraries
library(tidyverse)
library(magrittr)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(NeuRRoStat)

#####################
#### Custom functions
#####################

# Function to create new values near given values, with restriction of the 
# updated values not going below 0 or above 1 (used in two-step EM).
nearVal <- function(values){
  add <- round(values, 4) + 0.025
  if(any(add <= 0 | add >= 1)){
    return(values)
  }else{
    return(add)
  }
}


# Dimension of the brains
DIM <- c(53,63,46)

# Number of runs
NRUNS <- 50

# Starting amount of subjects in group analysis
startingTotal <- 10

# Database total
DATTOTAL <- 1400

# Total amount of possible subjects if we have a maximum of 3 disjoint subsets
NTOT <- 460

# Steps in the sequence (sample sizes)
steps <- seq(startingTotal, NTOT, by = 10)

# Create the sequence
sequence <- data.frame()
for(k in 1:NRUNS){
  for(i in 1:length(steps)){
    # Number of disjoint sets in this sample size
    numDS <- floor(DATTOTAL/steps[i])
    for(s in 1:numDS){
      sequence <- rbind(sequence,
                        data.frame('run' = k, 'step' = i, 'set' = s))
    }
  }
}

# Dimension of this data frame
dimSeq <- dim(sequence)[1]

# Initial guess for starting values in the EM algorithm
# We want to estimate 
#     lambda: the proportion of truly active voxels, 
#     PI_A1: the probability of a truly active voxel being declared significant
#     PI_I1: the probability of a truly inactive voxel being incorrectly declared significant
ini_lambda <- 0.1
ini_PI_A1 <- 0.6
ini_PI_I1 <- 0.05

##
###############
### Estimation
###############
##

# We want to estimate 
#     lambda (the proportion of truly active voxels), 
#     PI_A1 (the probability of a voxel being truly active)
# and PI_I1 (the probability of a voxel being incorrectly activated)
# These parameters will then be used to calculate Cohen's kappa. 

# We do this for each run and for each sample size (step)

# Empty data frame
EM_params <- data.frame() %>% as_tibble()

# Start for loop over all runs 
for(k in 1:NRUNS){
  print(paste0('At run ', k))
  for(i in 1:length(steps)){
    print(i)
    # Number of disjoint sets in this sample size
    numDS <- floor(DATTOTAL/steps[i])
    
    # Read in the summed map and its mask
    sumMap <- readNIfTI(paste(RawDat, '/Run_', k, '/Step_', i, '/Gmap.nii.gz', sep = ''))[,,]
    maskStep <- readNIfTI(paste(RawDat, '/Run_', k, '/Step_', i, '/Mask_Gmap.nii.gz', sep = ''))[,,]
    
    # Put masked values to NA
    sumMap[maskStep == 0] <- NA
    
    # Go to vector and remove NA
    sumVec <- c(array(sumMap, dim = prod(DIM)))
    sumVec_val <- sumVec[!is.na(sumVec)]
    
    # Initiate EM algorithm (see NeuRRoStat package)
    # We first run the algorithm with sensible starting values for fMRI data.
    # Then we run it with values close to those obtained as a starting value.
    # CAREFUL: N IS NOT THE SAMPLE SIZE, BUT THE NUMBER OF DISJOINT SETS!!!
    iniGuess <- NeuRRoStat::EMbinom(Y = sumVec_val, N = numDS, iniL = ini_lambda,
                        iniPI1 = ini_PI_A1, iniPI2 = ini_PI_I1)
    nearValues <- nearVal(iniGuess[,-4])
    if(nearValues$PI1 == nearValues$PI2){
      EM_param <- iniGuess
    }else{
      EM_param <- NeuRRoStat::EMbinom(Y = sumVec_val, N = numDS, 
                          iniL = nearValues$lambda,
                          iniPI1 = nearValues$PI1, iniPI2 = nearValues$PI2)
    }
    
    # Save the parameters (of initial and final estimation) with info about step
    EM_params <- bind_rows(EM_params,
            data.frame(iniGuess, 'step' = i, 'run' = k, 'final' = FALSE),
            data.frame(EM_param, 'step' = i, 'run' = k, 'final' = TRUE))
    
    # Reset parameters
    rm(iniGuess, EM_param, sumVec_val, sumVec, sumMap, maskStep)
  }
}


##
###############
### Save processed data
###############
##

# ID for EM_params
if(scenario == 'uncorrected'){
  saveRDS(EM_params, 
          file = paste(SaveLoc, '/', contr, '/EM_params_unc.rda', sep = ''))  
}
if(scenario == 'fdr'){
  saveRDS(EM_params, 
          file = paste(SaveLoc, '/', contr, '/EM_params_fdr.rda', sep = ''))  
}





















