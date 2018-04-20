####################
#### TITLE:     Calculate within sample size overlap (Dice only) with increasing sample sizes
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/SplitSubjects
#### First Modified: 20/08/2015
#### Notes: 
#################



##
###############
### Notes
###############
##

# Scenario:
# Completely split up subjects.
# Each run we have an analysis between group 1 and 2 in a specific sample size.
# The subjects in group 1 are completely different from those in group 2 (sole restriction)

# Then we calculate the Maitra overlap measure (adapted Dice overlap).


##
###############
### Preparation
###############
##

# Location of raw data: not included in Github folder (too large)
RawDat <- "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/SplitSubjects/ScenarioA_50"

# Save intermediate results
SaveLoc <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Script.git/FreddieFreeloader/Analyses/_IntData'

# Load in libraries
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
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


##
###############
### For loop: start with sample size S, calculate overlap between all images, then move to sample size S+1
###############
##

# Matrix were data will get into
MatrixOverlap <- array(NA,dim=c(NSTEP,NRUNS))

# Vector where we will store the percentage of activated masked voxels
# As we have two images per step and run, we will take the average of the percentage. 
PercAct <- array(NA,dim=c(NSTEP,NRUNS))

# Print statement
PriST <- (c(1:NRUNS)/NRUNS)[seq(1,NRUNS,length.out=10)][-10]

# For loop over all runs
for(i in 1:NRUNS){
  IDprint <- c(i/NRUNS)==PriST
  if(any(IDprint)){
    print(paste(round(PriST,2)[IDprint]*100, "% done"))
  }
  # For loop over all steps
  for(j in 1:NSTEP){
    # We have two images (G1 and G2) from group 1 and group 2: 0 is non significant and 1 is significant
    # We also need to check whether both of the thresholded maps are present (if not, skip step)
    imageG1 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group1','/thresh_zstat1.nii',sep=''))[,,],silent=TRUE)
    if(!class(imageG1)=='array') next
    maskG1 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group1','/masked.nii.gz',sep=''))[,,],silent=TRUE)
    if(!class(maskG1)=='array') next
    idMaskG1 <- maskG1==0
    imageG1[idMaskG1] <- NA
    imageG2 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group2','/thresh_zstat1.nii',sep=''))[,,],silent=TRUE)
    if(!class(imageG2)=='array') next
    maskG2 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group2','/masked.nii.gz',sep=''))[,,],silent=TRUE)
    if(!class(maskG2)=='array') next
    idMaskG2 <- maskG2==0
    imageG2[idMaskG2] <- NA
    # Summing image K and K-1 to know the voxels in both maps (1+1 = 2)
    sumMap <- imageG1+imageG2
    # Minus map: know the voxels different in both images
    minusMap <- imageG1-imageG2
    # Mask this image
    idMask <- maskG1*maskG2
    minusMap[idMask] <- NA
    
    # Union of activated voxels, voxels in image G1 and voxels in image G2
    Vjt <- length(sumMap[which(sumMap==2)])
    Vt <- length(imageG1[which(imageG1==1)])
    Vj <- length(imageG2[which(imageG2==1)])
    # Voxels different from both images
    VjtS <- length(minusMap[which(minusMap==0)])
    
    
    # Put overlap in matrix
    MatrixOverlap[j,i] <- round((Vjt)/(Vj + Vt - Vjt),6)
    
    # Now calculate percentage of activated masked voxels in imageG1 and imageG2
    baseG1 <- sum(maskG1)
    percG1 <- Vt/baseG1
    
    baseG2 <- sum(maskG2)
    percG2 <- Vt/baseG2
    PercAct[j,i] <- round(mean(c(percG1,percG2)),4)
    
    # Remove objects
    rm(Vjt,Vt,Vj,imageG1,imageG2,sumMap)
  }
  if(i==NRUNS) print("100% done")
}

# Transform NaN values to 0
MatrixOverlap[is.nan(MatrixOverlap)] <- 0

# Check for missing values (TRUE = missing somewhere)
complete.cases(t(MatrixOverlap))
MissingValues <- function(...){
  any(is.na(...))
}
apply(MatrixOverlap,c(2),MissingValues)


##
###############
### Save intermediate results
###############
##

## Save R objects
saveRDS(MatrixOverlap, file=paste(SaveLoc,'/MaitraOverlap.rda',sep=''))
saveRDS(PercAct,file=paste(SaveLoc,'/PercActMaitraOverlap.rda',sep=''))
