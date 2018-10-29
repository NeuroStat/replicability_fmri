####################
#### TITLE:     Calculate within sample size overlap (Dice only) with increasing sample sizes
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/I****ATA/Data/FreddieFreeloader/SplitSubjects
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

# We also calculate the average significance threshold corresponding to 
# FDR 0.05.


##
###############
### Preparation
###############
##

# Source paths
source(blind_PreProcessing.R)

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


########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

##
###############
### Measure P-value threshold
###############
##

# Empty data frame
PvalThr <- array(NA,dim=c(NSTEP,NRUNS))

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
    # We have two images Z-statistic images
    imageZ1 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group1/stats/zstat1.nii',sep=''))[,,],silent=TRUE)
    if(!class(imageZ1)=='array') next
    maskG1 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group1','/masked.nii.gz',sep=''))[,,],silent=TRUE)
    if(!class(maskG1)=='array') next
    idMaskG1 <- maskG1==0
    imageZ1[idMaskG1] <- NA
    imageZ2 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group2/stats/zstat1.nii',sep=''))[,,],silent=TRUE)
    if(!class(imageZ2)=='array') next
    maskG2 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group2','/masked.nii.gz',sep=''))[,,],silent=TRUE)
    if(!class(maskG2)=='array') next
    idMaskG2 <- maskG2==0
    imageZ2[idMaskG2] <- NA
    
    # Now convert both images to P-values
    imageP1 <- pnorm(q = imageZ1, 0, 1, lower.tail = FALSE)
    imageP2 <- pnorm(q = imageZ2, 0, 1, lower.tail = FALSE)
    
    # Adjust the P-values using B&H approach
    imageP1_adj <- p.adjust(imageP1, method = 'BH')
    imageP2_adj <- p.adjust(imageP2, method = 'BH')
    
    # Now match the original P-values with the adjusted ones
    matching1 <- data.frame('orig' = array(imageP1, dim = prod(DIM)),
               'adj' = array(imageP1_adj, dim = prod(DIM)))
    matching2 <- data.frame('orig' = array(imageP2, dim = prod(DIM)),
                            'adj' = array(imageP2_adj, dim = prod(DIM)))
    # Sort the original P-values
    matching1 <- matching1[order(matching1$orig, na.last = TRUE, decreasing = FALSE),]
    matching2 <- matching2[order(matching2$orig, na.last = TRUE, decreasing = FALSE),]
    
    # Check the P-value at corrected 0.05
    Porig1 <- matching1[matching1$adj <= 0.05,'orig'][1]
    Porig2 <- matching2[matching2$adj <= 0.05,'orig'][1]
    
    # Take the average
    Porig <- mean(c(Porig1, Porig2), na.rm = TRUE)
    
    # Add to data frame
    PvalThr[j,i] <- Porig
  }
  if(i==NRUNS) print("100% done")
}

# Average over all runs
AvgRun <- apply(PvalThr, 1, mean, na.rm = TRUE)
# Now over all steps
mean(AvgRun, na.rm = TRUE)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################


##
###############
### Repeat the code from above (I KNOW I COPIED IT, SORRY GUYS) with other levels for FDR 
###############
##

# Thresholds considered
FDR_Pth <- 0.2
FDR_Pth <- 0.1
FDR_Pth <- 0.01
FDR_Pth <- 0.001

# Matrix were data will get into
MatrixOverlapFDR <- array(NA,dim=c(NSTEP,NRUNS))

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
    # We read in non-thresholded z-images
    imageG1Z <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group1','/stats/zstat1.nii',sep=''))[,,], silent=TRUE)
    if(!class(imageG1Z)=='array') next
    # Read in the mask
    maskG1 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group1','/masked.nii.gz',sep=''))[,,],silent=TRUE)
    idMaskG1 <- maskG1==0
    imageG1Z[idMaskG1] <- NA
    # Convert to P-values
    imageG1P <- pnorm(q = imageG1Z, lower.tail = FALSE)
    # Adjust using Benjaminig and Hochberg approach
    imageG1aP <- array(p.adjust(p = imageG1P, method = 'BH'), dim = DIM)
    # Threshold at P_fdr < FDR_Pth
    imageG1t <- imageG1aP
    imageG1t[imageG1aP <= FDR_Pth] <- 1
    imageG1t[imageG1aP > FDR_Pth] <- 0
    
    # Repeat with the replication
    imageG2Z <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group2','/stats/zstat1.nii',sep=''))[,,], silent=TRUE)
    if(!class(imageG2Z)=='array') next
    # Read in the mask
    maskG2 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group2','/masked.nii.gz',sep=''))[,,],silent=TRUE)
    idMaskG2 <- maskG2==0
    imageG2Z[idMaskG2] <- NA
    # Convert to P-values
    imageG2P <- pnorm(q = imageG2Z, lower.tail = FALSE)
    # Adjust using Benjaminig and Hochberg approach
    imageG2aP <- array(p.adjust(p = imageG2P, method = 'BH'), dim = DIM)
    # Threshold at P_fdr < 0.1
    imageG2t <- imageG2aP
    imageG2t[imageG2aP <= FDR_Pth] <- 1
    imageG2t[imageG2aP > FDR_Pth] <- 0
    
    # Summing image K and K-1 to know the voxels in both maps (1+1 = 2)
    sumMap <- imageG1t+imageG2t
    # Minus map: know the voxels different in both images
    minusMap <- imageG1t-imageG2t
    # Mask this image
    idMask <- maskG1*maskG2
    minusMap[idMask] <- NA
    
    # Union of activated voxels, voxels in image G1 and voxels in image G2
    Vjt <- length(sumMap[which(sumMap==2)])
    Vt <- length(imageG1t[which(imageG1t==1)])
    Vj <- length(imageG2t[which(imageG2t==1)])
    # Voxels different from both images
    VjtS <- length(minusMap[which(minusMap==0)])

    # Put overlap in matrix
    MatrixOverlapFDR[j,i] <- round((Vjt)/(Vj + Vt - Vjt),6)
    
    # Now calculate percentage of activated masked voxels in imageG1t and imageG2t
    baseG1 <- sum(maskG1)
    percG1 <- Vt/baseG1
    
    baseG2 <- sum(maskG2)
    percG2 <- Vt/baseG2
    #PercAct[j,i] <- round(mean(c(percG1,percG2)),4)
    
    # Remove objects
    rm(Vjt,Vt,Vj,imageG1t,imageG2t,sumMap)
  }
  if(i==NRUNS) print("100% done")
}

# Transform NaN values to 0
MatrixOverlapFDR[is.nan(MatrixOverlapFDR)] <- 0

# Check for missing values (TRUE = missing somewhere)
complete.cases(t(MatrixOverlapFDR))
MissingValues <- function(...){
  any(is.na(...))
}
apply(MatrixOverlapFDR,c(2),MissingValues)

# Identifier for the files
IDTh <- sub(pattern = '.', replacement = '_', x = FDR_Pth, fixed = TRUE)

## Save R objects
saveRDS(MatrixOverlapFDR, file=paste(SaveLoc,'/MaitraOverlapFDR',IDTh, '.rda',sep=''))














