####################
#### TITLE:     HPC script for adaptive thresholding, FACES dataset
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Script/Analyses/Faces/Scenario9
#### First Modified: 13/11/2015
#### Notes: 
#################



##
###############
### Notes
###############
##

## Fixate percentage of activated voxels.
## We do this by starting with a starting value of significance level and then updating p-value by steps of .01.
## If percentage < lowerbound, then increase p-value. If > upperbound, decrease p-value.


## Data comes from scenario 6, FACES dataset

## HPC version to get parallel version on the runs for testing low percentages.


##
###############
### Preparation
###############
##

# Take arguments
args <- commandArgs(TRUE)

# Reset working directory
rm(list=ls())
gc(verbose = FALSE)


# Set WD
wd <- "/user/scratch/gent/gvo000/gvo00022/vsc40728/Freddie/FACES9"
setwd(wd)


# Load in libraries
library(oro.nifti)

# Location of the data
dat <- '/user/scratch/gent/gvo000/gvo00022/vsc40728/Freddie/FACES6/Results'

# Dimension of the brains
DIM <- c(53,63,46)

# Number of runs (so far)
NRUNS <- 42
# Number of steps
NSTEP <- 70

# Seed
set.seed(11121990)


# Baseline, starting significance level
signLevel <- 0.05


# Index
i <- as.numeric(as.character(args)[1])


##
###############
### Functions
###############
##

# Adaptive Thresholding
AdapThresholding <- function(pvals,DIM,signLevel,idMask,percentageToAch){
	# No adjustment for multiple testing with P-values with Benjamin & Hochberg procedure (built in R function)
	adjPvals <- pvals

	# Upper and lower bound based on percentage
	lower <- percentageToAch-.01
	upper <- percentageToAch+.01

	# Threshold at signLevel BH p-value: significant P-values get 1!
	idP <- adjPvals <= signLevel
	threshPval <- adjPvals
	threshPval[!idP] <- 0
	threshPval[idP] <- 1

	# Calculate starting (base) percentage of activated masked voxels
	maskedVox <- sum(!idMask)		# idMask: TRUE for voxel outside of mask. Hence sum of reverse. 
	Vt <- sum(idP,na.rm=TRUE)
	percentage <- Vt/maskedVox

	# Now if percentage is below lowerlimit, increase threshold (by .00001), if above upperlimit; decrease by .00001
	while(percentage < lower | percentage > upper){
		if(percentage < lower){
			if(signLevel >= 0.995){
				signLevel <- signLevel + 0.00001
			}else{
				signLevel <- signLevel + runif(1,min=.0001,max=0.005)
			}
		}
		if(percentage > upper){
			if(signLevel <= 0.005){
				signLevel <- signLevel - 0.00001	
			}else{
				signLevel <- signLevel - runif(1,min=.0001,max=0.005)
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


##
###############
### For loop: start with group A, sample size S within a RUN, calculate percentage, adapt if necessary, 
### move to group B, calculate overlap between the two images, then move to sample size S+1
###############
##


# Matrix were overlap results will get into
MatrixOverlapAdap <- array(NA,dim=c(NSTEP,1))

# Vector where we will store the final percentage of activated masked voxels
	# As we have two images per step and run, we will take the average of the percentage. 
PercActAdap <- array(NA,dim=c(NSTEP,1))

# Vector to store average of percentage activated masked voxels in each step and run
SignLevels <- array(NA,dim=c(NSTEP,1))

# Vector for activation map of first run only
SPMRun1G1 <- array(NA,dim=c(prod(DIM),NSTEP))
SPMRun1G2 <- array(NA,dim=c(prod(DIM),NSTEP))

# What percentage do we want to achieve?
PercToAchieve <- 0.10


# For loop over all runs
print(paste('--------------- Run ', i,'--------------',sep=''))
# For loop over all steps
for(j in 1:NSTEP){
	print(paste('Step ',j,sep=''))
	###################################################################
	# We have two images (G1 and G2) from group 1 and group 2: 0 is non significant and 1 is significant
		# We also need to check whether both of the thresholded maps are present (if not, skip step)
	imageG1 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group1','/stats/zstat1.nii',sep=''))[,,],silent=TRUE)
		if(!class(imageG1)=='array') next
		maskG1 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group1','/masked.nii.gz',sep=''))[,,],silent=TRUE)
			if(!class(maskG1)=='array') next
		idMaskG1 <- maskG1==0
		imageG1[idMaskG1] <- NA
		# Check second image (there is no point in calculating percentages and stuff if not both images are available)
	imageG2 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group2','/stats/zstat1.nii',sep=''))[,,],silent=TRUE)
		if(!class(imageG2)=='array') next
		maskG2 <- try(readNIfTI(paste(dat,'/Run_',i,'/Step_',j,'/Group2','/masked.nii.gz',sep=''))[,,],silent=TRUE)
			if(!class(maskG2)=='array') next
		idMaskG2 <- maskG2==0
		imageG2[idMaskG2] <- NA

	###################################################################
	# Now that we have z-statistic images with its masks, calculate BF-adjusted p-values using baseline significance level
	pvalsG1 <- 1-pnorm(imageG1)
	pvalsG2 <- 1-pnorm(imageG2)

	# Function AdapThresholding: calculates adjusted p-values, then finds correct significance level according to percentage.
	adapG1 <- AdapThresholding(pvalsG1,DIM,signLevel,idMaskG1,PercToAchieve)
		signLevel <- adapG1$signLevel
		perc <- adapG1$percentage

	# Image G2
	adapG2 <- AdapThresholding(pvalsG2,DIM,signLevel,idMaskG2,PercToAchieve)
		signLevel <- round(mean(c(signLevel,adapG2$signLevel)),6)
		perc <- round(mean(c(perc,adapG2$percentage)),4)


	###################################################################
	# Calculate overlap: summing image K and K-1 to know the voxels in both maps (1+1 = 2)
	sumMap <- adapG1$SPM+adapG2$SPM
	# Minus map: know the voxels different in both images
	minusMap <- adapG1$SPM-adapG2$SPM
		# Mask this image
	idMask <- maskG1*maskG2
		minusMap[idMask] <- NA

	# Union of activated voxels, voxels in image G1 and voxels in image G2
	Vjt <- length(sumMap[which(sumMap==2)])
	Vt <- length(adapG1$SPM[which(adapG1$SPM==1)])
	Vj <- length(adapG2$SPM[which(adapG2$SPM==1)])
	# Voxels different from both images
	VjtS <- length(minusMap[which(minusMap==0)])


	# Put overlap, percentage and signLevel in matrix
	MatrixOverlapAdap[j,i] <- round((Vjt)/(Vj + Vt - Vjt),6)
	PercActAdap[j,i] <- perc
	SignLevels[j,i] <- signLevel


	# Remove objects
	rm(Vjt,Vt,Vj,imageG1,imageG2,sumMap,maskG1,idMaskG1,maskG2,idMaskG2,adapG1,adapG2)
	gc()
}



# Transform NaN values to 0
MatrixOverlapAdap[is.nan(MatrixOverlapAdap)] <- 0



## Save objects
save(MatrixOverlapAdap,file=paste(wd,'/MatrixOverlapAdapFace',i,sep=''))
save(PercActAdap,file=paste(wd,'/PercActSplitAdapFace',i,sep=''))
save(SignLevels,file=paste(wd,'/SignLevelsAdapFace',i,sep=''))







