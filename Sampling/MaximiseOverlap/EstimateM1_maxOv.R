####################
#### TITLE:  Estimate number of true active voxels through maximizing the overlap
#### Contents: 	
#### 
#### Source Files: /replicability_fmri/Sampling/
#### First Modified: 25/09/2018
#### Notes: 
#################


##
###############
### Notes
###############
##

# Generating images with different proportion of truly activated areas,
# for N = 50. Then maximize the overlap of activation by changing the 
# thresholds for activation.
# Note: the proportion of truly active voxels = m1.


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
  wd <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/MaximOverl'
}
if(MACHINE=='MAC'){
  wd <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/MaximOverlap'
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
library(NeuRRoStat)


##
###############
### Custom functions
###############
##


# Custom function generating a true grid with a baseline proportion (m1) of truly activated voxels
defineTrueLocProp <- function(DIM2D, baseProp){
  # Total number of voxels
  total <- prod(DIM2D)
  
  # True location
  if(baseProp == 0){
    TrueGrid <- array(0, dim = DIM2D)
  }
  if(baseProp == 1){
    TrueGrid <- array(1, dim = DIM2D)
  }
  if(!baseProp %in% c(0,1)){
    TrueGrid <- array(0, dim = c(prod(DIM2D), 1))
    TrueGrid[1:floor(baseProp * total)] <- 1
    TrueGrid <- array(TrueGrid, dim = DIM2D)
  }
  return(TrueGrid)
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

# Percentile function, replacing the less precise adaptive thresholding function
PercThresholding <- function(pvals, idMask, target){
  # Under construction
}

# Adaptive Thresholding to be used in maximizing the overlap
AdapThresholdingM <- function(pvals, signLevel, idMask, target = NULL, 
                              seed_fu = 1112, MaxPValThr = 0.1){
  # Algorithm works with a random number generator (to avoid infinite loops).
  # However, I want to preserve seeds generated in global environment.
  # Hence, we get the global seed, save it, and set it back when exiting this function.
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
  
  # Now set the seed for this function
  set.seed(seed_fu)

  # If no target is set, take 20% as target (original value of first algorithm)
  if(is.null(target)){
    target <- 0.20
  }
  
  # Threshold at signLevel p-value: significant P-values get 1!
  idP <- pvals <= signLevel
  threshPval <- pvals
  threshPval[!idP] <- 0
  threshPval[idP] <- 1
  
  # Calculate starting (base) percentage of activated masked voxels
  # NOTE: THIS IS DIFFERENT THAN ORIGINAL AdapThresholdingM FUNCTION!
  maskedVox <- sum(idMask)		# idMask: TRUE for voxel inside of mask. 
  Vt <- sum(idP,na.rm=TRUE)
  percentage <- Vt/maskedVox
  
  # Now if percentage is below (target - 0.01), increase threshold (by .001), 
  #     if above (target + 0.01); decrease by .01
  while(percentage < (target - 0.01) | percentage > (target + 0.01)){
    if(percentage < (target - 0.01)){
      # Set a boundary on the significance threshold at MaxPValThr
      if(signLevel > MaxPValThr){
        break
      }else{
        # Using random values to increase the significance level
        # Random values are needed as algorithm gets stuck in loop sometimes.
        signLevel <- signLevel + runif(1,min=.0001,max=0.0005)
      }
    }
    if(percentage > (target + 0.01)){
      if(signLevel <= 0.005){
        signLevel <- signLevel - 0.00001
      }else{
        # Using random values to decrease the significance level
        # Random values are needed as algorithm gets stuck in loop sometimes.
        signLevel <- signLevel - runif(1,min=.0001,max=0.0005)
      }
    }
    # Recalculate thresholded map and percentage
    idP <- pvals <= signLevel
    threshPval <- pvals
    threshPval[!idP] <- 0
    threshPval[idP] <- 1
    Vt <- sum(idP,na.rm=TRUE)
    percentage <- Vt/maskedVox
  }
  
  # Return percentage, signLevel and threshPval
  result <- list('signLevel' = signLevel,'percentage' = percentage,'SPM' = threshPval)
  return(result)
}

# Function to maximise the overlap within a boundary of P-value threshold (i.e. MaxPValThr)
maximizeOverlap <- function(PVal_map1, PVal_map2, idMask_map1, idMask_map2,
                            signLevel, startingTarget, MaxPValThr){
  # Algorithm works with a random number generator (to avoid infinite loops).
  # However, I want to preserve seeds generated in global environment.
  # Hence, we get the global seed, save it, and set it back when exiting this function.
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
  
  # Check dimensions of P-values and masks 
  # First check if all are either in one vector, or in more dimensions
  if(any(is.null(dim(PVal_map1)) | is.null(dim(PVal_map2)) | 
       is.null(dim(idMask_map1)) | is.null(dim(idMask_map2)))){
    if(!(is.null(dim(PVal_map1)) & is.null(dim(PVal_map2)) & 
         is.null(dim(idMask_map1)) & is.null(dim(idMask_map2)))){
      stop('Please provide the P-values image and mask in same the dimension!')
    }
  }
  # If it is null, then it is a vector
  if(is.null(dim(PVal_map1)) | is.null(dim(PVal_map2))){
    if(length(PVal_map1) != length(PVal_map2)){
      stop('Dimensions of P-value images need to match!')
    }
    if(length(PVal_map1) != length(idMask_map1)){
      stop('Please provide the P-values image and mask in same the dimension!')
    }
    if(length(PVal_map2) != length(idMask_map2)){
      stop('Please provide the P-values image and mask in same the dimension!')
    }
  }
  # If not null, then it is >= 2D
  if(!(is.null(dim(PVal_map1)) | is.null(dim(PVal_map2)))){
    # All dimensions need to match
    if(!all(dim(PVal_map1) == dim(PVal_map2))){
      stop('Dimensions of P-value images need to match!')
    }
    if(!all(dim(PVal_map1) == dim(idMask_map1))){
      stop('Dimensions of PVal_map1 and idMask_map1 need to match!')
    }
    if(!all(dim(PVal_map2) == dim(idMask_map2))){
      stop('Dimensions of PVal_map2 and idMask_map2 need to match!')
    }
  }
  
  # Local function to return overlap based on two (adaptively) thresholded maps
  checkOverlap <- function(PVal_map1, PVal_map2, idMask_map1, idMask_map2,
                           signLevel, target_perc){
    
    # First do thresholding based on starting target
    adapG1 <- AdapThresholdingM(pvals = PVal_map1, signLevel = signLevel,
                                idMask = idMask_map1, target = target_perc,
                                seed_fu = 1112, MaxPValThr = MaxPValThr)
    adapG2 <- AdapThresholdingM(pvals = PVal_map2, signLevel = signLevel,
                                idMask = idMask_map2, target = target_perc,
                                seed_fu = 1112, MaxPValThr = MaxPValThr)
    
    # Now calculate the current overlap
    cur_overlap <- NeuRRoStat::overlapAct(bin_map1 = adapG1$SPM, 
                                          bin_map2 = adapG2$SPM)
    return(cur_overlap)
  }
  
  # Check the overlap for a range of target percentages
  targPerc <- seq(0,1,by=0.01)
  obsOverl <- c()
  for(i in 1:length(targPerc)){
    obsOverl <- c(obsOverl,
        checkOverlap(PVal_map1, PVal_map2, idMask_map1, idMask_map2,
                 signLevel, target_perc = targPerc[i]))
  }
  
  # Make it a data frame
  obsOverl <- data.frame('target_perc' = targPerc, 'obs_overlap' = obsOverl)
  return(obsOverl)
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
baseProp <- seq(0, 1, by = 0.1)

# Number of subjects in the analysis
nsub <- 50

# Baseline, starting significance level
signLevel <- 0.05

# Boundary on the P-value threshold
MaxPValThr <- 0.25

# Starting percentage in the algorithm
startingTarget <- 0.5



##
###############
### Data generation
###############
##

# Empty data frame with estimated target percentages
target_data <- data.frame() %>% as_tibble()

# Empty data frame with full data
full_data <- data.frame() %>% as_tibble()

# Common mask 
comMask <- array(1, dim = DIM2D)

# For loop over the simulations
for(ID in startIndex:endIndex){
  # Set seed
  set.seed(ID*pi)
  print(ID)
  
  # For loop over the proportion of voxels with true activation
  for(l in 1:length(baseProp)){
    # Generate the true grid for this scenario
    TrueGrid_S <- defineTrueLocProp(DIM2D, baseProp[l])
    
    # Generate p value data for both maps and bind to each other
    PVals <- data.frame(PVal1 = gen_data(TrueGrid_S, TrueES, trueSigma, N = nsub),
                        PVal2 = gen_data(TrueGrid_S, TrueES, trueSigma, N = nsub))
    
    ##########################################################################
    # Maximising overlap through adaptive thresholding
    ##########################################################################
    ObsOverl <- maximizeOverlap(PVal_map1 = array(PVals$PVal1, dim = DIM2D), 
                               PVal_map2 = array(PVals$PVal2, dim = DIM2D),
                               idMask_map1 = comMask, idMask_map2 = comMask,
                               signLevel = signLevel, startingTarget = startingTarget,
                               MaxPValThr = MaxPValThr) %>% as_tibble() %>%
      # Add info to raw data
      mutate(TruePerc = baseProp[l],
             sim = ID)
    # Bind to the data frame
    full_data <- bind_rows(full_data, ObsOverl)
    
    ##########################################################################
    # Return the estimated proportion of truyly active voxels
    ##########################################################################
    
    #ObsOverl %>% mutate(prev_overlap = lag(obs_overlap),
     #                   increase = obs_overlap >= prev_overlap) %>% View()
    
  }
}



# 
# 
# 
# data.frame('target' = targPerc, 'overlap' = obsOverl) %>%
#   ggplot(., aes(x = target, y = overlap)) + 
#   geom_line() + 
#   geom_vline(xintercept = 0.3)
# 
# 
# 
# 
# 
# ObsOverl %>% mutate(prev_overlap = lag(obs_overlap),
#                     increase = obs_overlap >= prev_overlap)








##
###############
### Save results
###############
##



# Save R object
saveRDS(full_data, file = paste0(wd, '/Results/full_data_', hpcID, '.rds'))









# Test code 
adapG1 <- AdapThresholdingM(pvals = PVal_map1, signLevel = signLevel,
                            idMask = idMask_map1, target = target_perc,
                            seed_fu = 1112, MaxPValThr = MaxPValThr)

adapG1$percentage
adapG1$signLevel

data.frame(pval = array(PVal_map1, dim = prod(dim(PVal_map1)))) %>% 
  as_tibble() %>% 
  arrange(pval) %>% 
  mutate(perc = percent_rank(pval)) %>% 
  filter(perc <= 0.25) %>% tail()
  
?dplyr::percent_rank
































