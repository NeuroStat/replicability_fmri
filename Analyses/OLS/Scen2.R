####################
#### TITLE:     Different thresholds for calculating within samplesize stability for increasing sample sizes in group study. Scenario 2. OLS Pooling.
#### Contents: 	
#### 
#### Source Files: //FreddieFreeloader/Script.git/Sampling/Combinations/
#### First Modified: 05/10/2015
#### Notes: 
#################


##
###############
### Notes
###############
##

# The idea is to stepwise increase the sample size for a group study. 
# We repeat (resample) this process several times. Then we can compare all the group studies with a given sample size [k(k-1)/2] combinations.

# This R file is with all steps combined (going from sample size 10 to 740, total of 74 steps)

# We want to see the effect when using a different thresholding level

# Based on OLS pooling of subjects

# Note, some objects are saved for scenario 3 as well

##
###############
### Preparation
###############
##

rm(list=ls())

# Set WD: location of data. Same data as scenario 1!
wd <- "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/OLS/Scenario1"
setwd(wd)



# Libraries
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(fslr)
source('~/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')


# Dimension of the brains
DIM <- c(53,63,46)

# Number of runs
NRUNS <- 15
# Number of steps
NSTEP <- 74





##
###############
### Re-threshold the Z-statistics using voxelwise FDR.
###############
##

# Create list (each list element is a run, the fourth dimension is a step)
AllData <- NULL

Data <- array(NA,dim=c(DIM,NSTEP))
for(k in 1:NRUNS){
	AllData <- c(AllData, k = list(Data))
}

# Create list for all the masks
AllMasks <- NULL

Mask <- array(NA,dim=c(DIM,NSTEP))
for(k in 1:NRUNS){
	AllMasks <- c(AllMasks, k = list(Mask))
}

# Significance level
signLevel <- 0.001

# Start for loop over all runs and then steps
for(i in 1:NRUNS){
	print(paste('At run ',i,sep=''))
	for(j in 1:NSTEP){
		## Preparation
		# Temporary wd
		wd.tmp <- paste(wd,'/Run_',i,'/Step_',j,sep='')
		# Get df
		df <- as.numeric(read.table(paste(wd.tmp,'/stats/dof',sep='')))
		# tmap
		tmap <- readNIfTI(paste(wd.tmp,"/stats/tstat1.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
		# Mask
		mask <- readNIfTI(paste(wd.tmp,"/masked.nii.gz",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,]
			# Put the mask in the list
			AllMasks[[i]][,,,j] <- mask

		## FDR correction
		# Dimension of the tmap
		DIM <- dim(tmap)
		# Mask the t-values map
		idMask <- mask==0
		MTmap <- tmap
		MTmap[idMask] <- NA
		# Calculate P-values
		pvals <- 1-pt(MTmap,df)	

		# Adjust P-values with Benjamin & Hochberg procedure
		adjPvals <- array(p.adjust(pvals,method="BH"),dim=DIM)

		## Threshold at signLevel BH p-value: significant P-values get 1!
		idP <- adjPvals <= signLevel
		threshPval <- adjPvals
		threshPval[!idP] <- 0
		threshPval[idP] <- 1

		# Put it in the list
		AllData[[i]][,,,j] <- threshPval

		# Remove tmp objects
		rm(wd.tmp,df,tmap,mask,DIM,idMask,MTmap,pvals,adjPvals,idP,threshPval)
	}
}


# Save R Object
save(AllData,file=paste(wd,'/AllDataScen2',sep=''))
save(AllMasks, file=paste(wd,'/AllMasksScen2',sep=''))



##
###############
### Calculate within sample size overlap
###############
##


# All combinations of the runs
combinations <- t(combn(NRUNS,2))
NCOMB <- dim(combinations)[1]

# Global array of overlap
GlobalOverlap <- array(0,dim=c(NSTEP,1))
GlobalVjt <- array(0,dim=c(NSTEP,1))
GlobalVt <- array(0,dim=c(NSTEP,1))
GlobalVj <- array(0,dim=c(NSTEP,1))
# Of percentage based overlap
PercOverlap <- array(0,dim=c(NSTEP,1))

# Matrix of overlap and percentage based overlap
MatrixOverlap <- array(0,dim=c(NSTEP,NCOMB))
MatrixPercOverlap <- array(0,dim=c(NSTEP,NCOMB))

# Print statement
PriST <- (c(1:NSTEP)/NSTEP)[seq(1,NSTEP,length.out=10)][-10]

# For loop over all steps
for(i in 1:NSTEP){
	IDprint <- c(i/NSTEP)==PriST
	if(any(IDprint)){
		print(paste(round(PriST,2)[IDprint]*100, "% done"))
	}
	# For loop over all combinations
	for(j in 1:NCOMB){
		# Read in the data from image K and K-1: 0 is non significant and 1 is significant
		imageK <- AllData[[combinations[j,1]]][,,,i]
		imageK_1 <- AllData[[combinations[j,2]]][,,,i]
		
		# Get the appropriate mask: union of both masks
		mask <- AllMasks[[combinations[j,1]]][,,,i]+AllMasks[[combinations[j,2]]][,,,i]
			idMask <- mask==2
			NumMasked <- sum(idMask)
		# Summing image K and K-1 to know the voxels in both maps (1+1 = 2)
		sumMap <- imageK+imageK_1
		# Minus map: know the voxels different in both images
		minusMap <- imageK-imageK_1
			# Mask this image
			minusMap[!idMask] <- NA
	
		# Union of activated voxels, voxels in image K and voxels in image K-1
		Vjt <- length(sumMap[which(sumMap==2)])
		Vt <- length(imageK[which(imageK==1)])
		Vj <- length(imageK_1[which(imageK_1==1)])
		# Voxels different from both images
		VjtS <- length(minusMap[which(minusMap==0)])
		# Now calculate overlap
		GlobalOverlap[i] <- sum(GlobalOverlap[i],(round((Vjt)/(Vj + Vt - Vjt),6)),na.rm=TRUE)

		# Percentage overlap
		PercOverlap[i] <- sum(PercOverlap[i],(round((VjtS/NumMasked),6)),na.rm=TRUE)
		# We want to know if Vjt, Vt or Vj is maybe responsible for the observed effect
		GlobalVjt[i] <- GlobalVjt[i]+Vjt
		GlobalVt[i] <- GlobalVt[i]+Vt
		GlobalVj[i] <- GlobalVj[i]+Vj

		# We want to have the SD of the overlap, so calculate a matrix as well
		MatrixOverlap[i,j] <- round((Vjt)/(Vj + Vt - Vjt),6)
		MatrixPercOverlap[i,j] <- round((VjtS/NumMasked),6)

		# Remove objects
		rm(Vjt,Vt,Vj,imageK,imageK_1,sumMap)
	}
if(i==NSTEP) print("100% done")
}

# Now divide by NCOMB to get the average overlap
AverageOverlap <- GlobalOverlap/NCOMB
	# Same for Vjt, Vt and Vj
	AverageVjt <- GlobalVjt/NCOMB
	AverageVt <- GlobalVt/NCOMB
	AverageVj <- GlobalVj/NCOMB
	# Same for percentage based overlap
	AveragePercOverlap <- PercOverlap/NCOMB

# Transpose matrix of overlap and percentage based overlap
MatrixOverlap <- t(MatrixOverlap)
MatrixPercOverlap <- t(MatrixPercOverlap)


# Save those objects
save(AverageOverlap,file=paste(wd,'/AverageOverlapT001Scen2',sep=''))
save(AveragePercOverlap,file=paste(wd,'/AveragePercOverlapT001Scen2',sep=''))
save(MatrixOverlap,file=paste(wd,'/MatrixOverlapT001Scen2',sep=''))
save(MatrixPercOverlap,file=paste(wd,'/MatrixPercOverlapT001Scen2',sep=''))



















