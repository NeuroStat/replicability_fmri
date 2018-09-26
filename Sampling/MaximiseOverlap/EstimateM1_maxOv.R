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
# Which algorithm do you want to try?
algorithm <- try(as.character(input)[3],silent=TRUE)

# If no machine is specified, then it has to be this machine!
if(is.na(MACHINE)){
  MACHINE <- 'MAC'
  hpcID <- 1
  algorithm <- 'percentile'
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
PercThresholding <- function(pvals, target, MaxPValThr = 1){
  # NOTE: I will use dplyr functions, but avoid the pipe in this function
  # In order not to break things in the future...
  
  # Set boundary on the target
  if(target >= MaxPValThr){
    target <- MaxPValThr
  }
  
  # Switch to array of 1 dimension
  pvals_1D <- data.frame(pvals = array(pvals, dim = prod(dim(pvals))))
  
  # Give every voxel an ID
  pvals_ID <- pvals_1D
  pvals_ID$ID <- 1:prod(dim(pvals))
  
  # Sort the data frame on the p-values
  sorted <- dplyr::arrange(pvals_ID, pvals)
  
  # Use percentile function
  perc <- dplyr::mutate(sorted, perc_pval = percent_rank(pvals))
  
  # Voxels under target get 1, 0 otherwise
  bin_pval <- dplyr::mutate(perc, sel_pval = 
                              ifelse(perc_pval <= target, 1, 0))
  
  # Now sort on the original ID value, to have voxels back in correct location
  orig_loc <- dplyr::arrange(bin_pval,ID)
  
  # Binarized image --> switch back to original dimension of p-values
  orDim_bin <- array(orig_loc$sel_pval, dim = dim(pvals))
  
  # The p-value that corresponds to the threshold
  selected <- dplyr::filter(bin_pval, perc_pval <= target) 
  thrPval <- selected[dim(selected)[1],'pvals']
  
  # Final percentile selected voxels
  percSel <- selected[dim(selected)[1],'perc_pval']
  
  # Return percentage, signLevel and thresholded image
  result <- list('signLevel' = thrPval,'percentage' = percSel,'ThrImage' = orDim_bin)
  return(result)
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
  result <- list('signLevel' = signLevel, 'percentage' = percentage, 
                 'ThrImage' = threshPval)
  return(result)
}

# Function to maximise the overlap
# Two approaches are possible:
# * adaptive thresholding: adapts the P-value threshold within a boundary 
#   (i.e. MaxPValThr) until target percentage is achieved
# * percentile thresholding: takes percentile lowest P-values corresponding to
#   target percentage
maximizeOverlap <- function(PVal_map1, PVal_map2, idMask_map1, idMask_map2,
                            signLevel, startingTarget, MaxPValThr, 
                            algorithm = c('percentile', 'adaptive')){
  # Adaptive algorithm works with a random number generator (to avoid infinite loops).
  # However, I want to preserve seeds generated in global environment.
  # Hence, we get the global seed, save it, and set it back when exiting this function.
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
  
  # Check algorithm options
  algorithm <- match.arg(algorithm)
  
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
  checkOverlap <- function(PVal_map1, PVal_map2, target_perc){
    
    # First do thresholding based on starting target
    # Two approaches: using adaptive thresholding or using a percentile approach
    if(algorithm == 'adaptive'){
      adapG1 <- AdapThresholdingM(pvals = PVal_map1, signLevel = signLevel,
                                 idMask = idMask_map1, target = target_perc,
                                seed_fu = 1112, MaxPValThr = MaxPValThr)
      adapG2 <- AdapThresholdingM(pvals = PVal_map2, signLevel = signLevel,
                                 idMask = idMask_map2, target = target_perc,
                                seed_fu = 1112, MaxPValThr = MaxPValThr)
    }
    if(algorithm == 'percentile'){
      adapG1 <- PercThresholding(pvals = PVal_map1, target = target_perc)
      adapG2 <- PercThresholding(pvals = PVal_map2, target = target_perc)
    }
    
    # Now calculate the current overlap
    cur_overlap <- NeuRRoStat::overlapAct(bin_map1 = adapG1$ThrImage, 
                                          bin_map2 = adapG2$ThrImage)
    return(cur_overlap)
  }
  
  # Check the overlap for a range of target percentages
  targPerc <- seq(0,1,by=0.001)
  obsOverl <- c()
  for(i in 1:length(targPerc)){
    obsOverl <- c(obsOverl,
        checkOverlap(PVal_map1, PVal_map2, target_perc = targPerc[i]))
  }
  
  # Make it a data frame
  obsOverl <- data.frame('target_perc' = targPerc, 'obs_overlap' = obsOverl)
  return(obsOverl)
}

# gridSearch is a function to obtain a local maximum for a function
# that is by default increasing to 1
gridSearch <- function(obsOverl, gapWidth, 
                       algorithm = c('percentile', 'adaptive')){
  # Check choice of algorithm
  algorithm <- match.arg(algorithm)
  
  # If percentile approach, need to do reverse direction search
  if(algorithm == 'percentile'){
    # Starting from target percentage = 1
    for(i in dim(obsOverl)[1]:1){
      # Select the other indices with some gap in between
      if(i != gapWidth){
        indices <- 1:(i - gapWidth)
        # If no value has been found outside this gap, then the estimate is
        # the maximum target percentage
      }else{
        estimateM1 <- max(obsOverl[,'target_perc'])
        break
      }
      # Is there a value lower than the one we observe at the moment?
      idTargets <- as.numeric(obsOverl[i,"obs_overlap"]) <= obsOverl[indices,"obs_overlap"]
      if(any(idTargets, na.rm = TRUE)){
        XcutOff <- obsOverl[i,"target_perc"]
        YcutOff <- obsOverl[i,"obs_overlap"]
        # Save the obtained target percentage that is in the remaining part of the graph
        remainder <- obsOverl[indices,]
        estimateM1 <- as.numeric(unlist(remainder[idTargets,'target_perc']))
        # If multiple values are found, then take the average
        if(length(estimateM1) != 1){
          estimateM1 <- mean(estimateM1, na.rm = TRUE) 
        }
        break
      }
    }
  }
  # Adaptive algorithm has for some reason a global maximum which is our point estimate
  if(algorithm == 'adaptive'){
    # Take the maximum overlap as the point estimate, there is no need to do
    # reverse search as the function is different.
    maxVal <- dplyr::filter(obsOverl, obs_overlap == max(obs_overlap, na.rm = TRUE))
    estimateM1 <- as.numeric(unlist(maxVal$target_perc))
    # If multiple values are found, then take the average
    if(length(estimateM1) != 1){
      estimateM1 <- mean(estimateM1, na.rm = TRUE) 
    }
  }
  # Return the estimate
  return(data.frame('m1' = estimateM1))
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

# Starting percentage in the adaptive algorithm
startingTarget <- 0.5


##
###############
### Data generation
###############
##

# Empty data frame with full data
full_data <- data.frame() %>% as_tibble()

# Common mask 
comMask <- array(1, dim = DIM2D)

# For loop over the simulations
for(ID in startIndex:endIndex){
  # Set seed
  set.seed(ID*pi)

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
                               signLevel = signLevel, 
                               startingTarget = startingTarget,
                               MaxPValThr = MaxPValThr, 
                               algorithm = algorithm) %>% as_tibble() %>%
      # Add info to raw data
      mutate(TruePerc = baseProp[l],
             sim = ID)
    
    ##########################################################################
    # Return the estimated proportion of truyly active voxels
    ##########################################################################
    # Using the function gridSearch and two widths.
    # This is only needed for the percentile approach:
      # We need to use the small width for proportions that are close to 1.
      # But we need a large gip if the proportion = 1.
    estimateM1 <- min(gridSearch(ObsOverl, gapWidth = 10, algorithm = algorithm), 
                    gridSearch(ObsOverl, gapWidth = 100, algorithm = algorithm))
    
    # Let us add this to the raw data frame
    ObsOverl <- mutate(ObsOverl, estM1 = estimateM1)
    
    ##########################################################################
    # Bind to data frame
    ##########################################################################
    full_data <- bind_rows(full_data, ObsOverl)
  }
}




##
###############
### Save results
###############
##

# Keywords for the algorithm choice
if(algorithm == 'adaptive'){
  keyAlg <- 'adap'
}
if(algorithm == 'percentile'){
  keyAlg <- 'perc'
}

# Save R object
saveRDS(full_data, file = paste0(wd, '/Results/full_data_', keyAlg, '_', hpcID, '.rds'))

















































