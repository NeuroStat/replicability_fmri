####################
#### TITLE:     Adaptive thresholding with maximizing overlap.
####            Which target percentages gives the highest overlap?
####            
#### Contents:
####
#### Source Files: FreddieFreeloader/Script/Analyses/AdaptiveThresholding
#### First Modified: 24/09/2018
#### Notes:
#################



##
###############
### Notes
###############
##


# Start at maximum sample size.
# We calculate the overlap using adaptive thresholding where we fix
# the percentage of declared significant voxels.
# Then we calculate the overlap between the independent images.
# Next we maximise the overlap by changing the target percentage.
# We asusme the function between target percentage and overlap has one maximum!

# Over all runs, we take the mean overlap as outcome measure!



##
###############
### Preparation
###############
##

# Location of the (raw) data
# This data is not stored at Github (too large)
dat <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/SplitSubjects/ScenarioA_50'

# Save intermediate results
SaveLoc <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Script.git/FreddieFreeloader/Analyses/_IntData'

# Load in libraries
library(tidyverse)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(NeuRRoStat)

# Dimension of the brains
DIM <- c(53,63,46)
# Number of runs
NRUNS <- 50
# Number of steps
NSTEP <- 70

# Seed
set.seed(11121990)

# Baseline, starting significance level
signLevel <- 0.05

# Starting target percentage
startTarget <- 0.20

# Boundary on the P-value threshold
MaxPValThr <- 0.1


##
###############
### Functions
###############
##

# Adaptive Thresholding
AdapThresholdingM <- function(pvals, signLevel, idMask, target = NULL, 
                              seed_fu = 1112, MaxPValThr = 0.1){
  # Algorithm works with a random number generator (to avoid infinite loops).
  # However, I want to preserve seeds generated in global environment.
  # Hence, we get the global seed, save it, and set it back when exiting this function.
  old <- .Random.seed
  on.exit( { .Random.seed <<- old } )
  
  # Now set the seed for this function
  set.seed(seed_fu)

  # No adjustment for multiple testing (not reallly needed, not of interest)
  adjPvals <- pvals
  
  # If no target is set, take 20% as target (original value of first graph)
  if(is.null(target)){
    target <- 0.20
  }
  
  # Threshold at signLevel p-value: significant P-values get 1!
  idP <- adjPvals <= signLevel
  threshPval <- adjPvals
  threshPval[!idP] <- 0
  threshPval[idP] <- 1
  
  # Calculate starting (base) percentage of activated masked voxels
  maskedVox <- sum(!idMask)		# idMask: TRUE for voxel outside of mask. Hence sum of reverse.
  Vt <- sum(idP,na.rm=TRUE)
  percentage <- Vt/maskedVox
  
  # Now if percentage is below (target - 0.01), increase threshold (by .001), 
  #     if above (target + 0.01); decrease by .01
  while(percentage < (target - 0.01) | percentage > (target + 0.01)){
    if(percentage < (target - 0.01)){
      # Set a boundary on the significance threshold at MaxPValThr
      if(signLevel > MaxPValThr){
        print(paste('Reached maximum threshold at P <= ', MaxPValThr, sep = ''))
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
    idP <- adjPvals <= signLevel
    threshPval <- adjPvals
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
  if(!(is.null(dim(a)) | is.null(dim(b)))){
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
    adapG2 <- AdapThresholdingM(pvals = idMask_map2, signLevel = signLevel,
                                idMask = idMaskG2, target = target_perc,
                                seed_fu = 1112, MaxPValThr = MaxPValThr)
    
    # Now calculate the current overlap
    cur_overlap <- NeuRRoStat::overlapAct(bin_map1 = adapG1$SPM, 
                                          bin_map2 = adapG2$SPM)
    return(cur_overlap)
  }
  
  # First get overlap based on starting target percentage
  curOverlap <- checkOverlap(PVal_map1, PVal_map2, idMask_map1, idMask_map2,
                             signLevel, target_perc = startingTarget)
  
  # Now have two directions of the 
  #   target to test whether overlap increases/decreases.
  upTarget <- startingTarget + 0.1
  upOverlap <- checkOverlap(PVal_map1, PVal_map2, idMask_map1, idMask_map2,
                            signLevel, target_perc = upTarget)
  
  # Now run the while loop
  while(curOverlap < (upOverlap - 0.005) | curOverlap > (upOverlap + 0.005)){
    if(curOverlap < (upOverlap - 0.005)){
      # Update the target: increase it
      upTarget <- upTarget + runif(n = 1, min = 0.020, max = 0.08)
    }
    if(curOverlap > (upOverlap + 0.005)){
      # Update the target: decrease it
      upTarget <- upTarget - runif(n = 1, min = 0.020, max = 0.08)
    }
    # Recalculate current and updated overlap
    print(paste0('Current overlap = ', curOverlap, ' while updated = ', upOverlap))
    curOverlap <- upOverlap
    upOverlap <- checkOverlap(PVal_map1, PVal_map2, idMask_map1, idMask_map2,
                              signLevel, target_perc = upTarget)    
  }
  
  # Now we have the maximised target percentage, within a given P-value
  # boundary.
  # Return the target percentage that maximises the overlap within this boundary
  return(upTarget)
}



##
###############
### Read data and calculate adaptive overlap
###############
##

# Empty data frame
MaxTargetPerc <- data.frame() %>% as_tibble()

# For loop over all runs
for(i in 1:NRUNS){
  print(paste('--------------- Run ', i,'--------------',sep=''))
  # Note: we can restrict to last step, assuming this gives the largest value of overlap
  
  # We have two images (G1 and G2) from group 1 and group 2: 0 is non significant and 1 is significant
  # We also need to check whether both of the thresholded maps are present (if not, skip step)
  
  ############################ READ IN IMAGES ##################################
  ####### IMAGE 1 #######
  imageG1 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',NSTEP,'/Group1','/stats/zstat1.nii',sep=''))[,,],silent=TRUE)
  if(!class(imageG1)=='array') next
  maskG1 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',NSTEP,'/Group1','/masked.nii.gz',sep=''))[,,],silent=TRUE)
  if(!class(maskG1)=='array') next
  # Create mask ID
  idMaskG1 <- maskG1==0
  # These voxels are excluded
  imageG1[idMaskG1] <- NA
  
  ####### IMAGE 2 #######
  imageG2 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',NSTEP,'/Group2','/stats/zstat1.nii',sep=''))[,,],silent=TRUE)
  if(!class(imageG2)=='array') next
  maskG2 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',NSTEP,'/Group2','/masked.nii.gz',sep=''))[,,],silent=TRUE)
  if(!class(maskG2)=='array') next
  # Create mask ID
  idMaskG2 <- maskG2==0
  # These voxels are excluded
  imageG2[idMaskG2] <- NA
  
  ############################## P-VALUES ######################################
  pvalsG1 <- 1-pnorm(imageG1)
  pvalsG2 <- 1-pnorm(imageG2)
  
  ########################## MAXIMISE OVERLAP ##################################
  MaxPerc <- maximizeOverlap(PVal_map1 = pvalsG1, PVal_map2 = pvalsG2,
                  idMask_map1 = idMaskG1, idMask_map2 = idMaskG2,
                  signLevel = signLevel, startingTarget = startingTarget,
                  MaxPValThr = MaxPValThr)
  
  # Save results
  MaxTargetPerc <- bind_rows(MaxTargetPerc,
                             data.frame('MaxPerc' = MaxPerc,
                                        'run' = i))
}

# Now retrieve summary on the percentage
summary(MaxTargetPerc$MaxPerc)
  #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  #  0.2516  0.3000  0.3764  0.3541  0.4045  0.4701 


ggplot(MaxTargetPerc, aes(x = run, y = MaxPerc)) +
  geom_line()






















