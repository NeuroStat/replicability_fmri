####################
#### TITLE:     Calculate within samplesize stability with increasing sample sizes in group study: sample size from 10-100 controlling for the center. Scenario 4: OLS POOLING.
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Center
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

# Sampling proccess is controlled within centers: scenario 4.

# Based on OLS pooling of subjects.


##
###############
### Preparation
###############
##

rm(list=ls())

# Set WD
wd <- "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/OLS/Scenario4"
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
	# Print which libraries are needed
	print("Need packages: oro.nifti, fslr, lattice, ggplot2, reshape2, RColorBrewer, gridExtra and the functions.R file")


# Dimension of the brains
DIM <- c(53,63,46)

# Number of runs
NRUNS <- 21
# Number of steps
NSTEP <- 10



# Load Center information
load('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/IMAGEN/IMAGEN/IMAGEN/OurIDs')
head(OurIDs)


##
###############
### For loop: start with sample size S, calculate overlap between all images, then move to sample size S+1
###############
##

# All combinations of the runs
combinations <- t(combn(NRUNS,2))
NCOMB <- dim(combinations)[1]

# Centers: when both same number, then they come from same center.
centers <- cbind(combinations[,1]%%7,combinations[,2]%%7)
	# Function to use comparing the two numbers
	equal <- function(...){
		input <- c(...)
		if(length(input)!=2) stop('Function equal only accepts two numbers to compare')
		ifelse(input[1]==input[2],ind <- 1,ind <- 0) 
		return(ind)
	}
IndCenter <- apply(centers,c(1),equal)



# Matrix of overlap and percentage based overlap
MatrixOverlap <- array(0,dim=c(NSTEP,NCOMB))


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
		imageK <- readNIfTI(paste(wd,'/Run_',combinations[j,1],'/Step_',i,'/thresh_zstat1.nii',sep=''))[,,]
			# Mask image
			maskK1 <- readNIfTI(paste(wd,'/Run_',combinations[j,1],'/Step_',i,'/masked.nii.gz',sep=''))[,,]
			idK1 <- maskK1==0
			imageK[idK1] <- NA

		imageK_1 <- readNIfTI(paste(wd,'/Run_',combinations[j,2],'/Step_',i,'/thresh_zstat1.nii',sep=''))[,,]
			# Mask image
			maskK_1 <- readNIfTI(paste(wd,'/Run_',combinations[j,2],'/Step_',i,'/masked.nii.gz',sep=''))[,,]
			idK_1 <- maskK_1==0
			imageK_1[idK_1] <- NA		

		# Summing image K and K-1 to know the voxels in both maps (1+1 = 2)
		sumMap <- imageK+imageK_1
		# Minus map: know the voxels different in both images
		minusMap <- imageK-imageK_1

	
		# Union of activated voxels, voxels in image K and voxels in image K-1
		Vjt <- length(sumMap[which(sumMap==2)])
		Vt <- length(imageK[which(imageK==1)])
		Vj <- length(imageK_1[which(imageK_1==1)])
		# Voxels different from both images
		VjtS <- length(minusMap[which(minusMap==0)])
		
		# We want to have the SD of the overlap, so calculate a matrix as well
		MatrixOverlap[i,j] <- round((Vjt)/(Vj + Vt - Vjt),6)

		# Remove objects
		rm(Vjt,Vt,Vj,imageK,imageK_1,sumMap,maskK1,maskK_1,idK1,idK_1)
	}
if(i==NSTEP) print("100% done")
}

# Transpose matrix of overlap
MatrixOverlap <- t(MatrixOverlap)


# Save object
save(MatrixOverlap,file=paste(wd,'/MatOverlapCenterScen4',sep=''))

# Load object
load(paste(wd,'/MatOverlapCenterScen4',sep=''))



##
###############
### Create plots
###############
##

quartz.options(width=15,height=12)


### Histogram of centers
summary(OurIDs$ScanID)
hist(as.numeric(OurIDs$ScanID))
	sumCent <- as.numeric(OurIDs$ScanID)
	table(sumCent)

ggplot(OurIDs,aes(ScanID)) + geom_bar(width=.4,fill='#2b8cbe') + 
	scale_x_discrete(name='Scanner',labels=c(1:7)) +
	scale_y_continuous(name='Count') +
	geom_hline(yintercept=min(table(sumCent))) +
	stat_bin(geom="text", aes(label=..count.., vjust=-1)) +
	ggtitle('Number of subjects in each center of the IMAGEN database.')




#### Start with Jaccard overlap index
# Prepare the data: all points of MatrixOverlap. Add the center ID to it!
subjBreak <- seq(10,100,by=10)
BoxOverlap <- matrix(MatrixOverlap,ncol=1)
	boxSubj <- rep(seq(10,100,by=10),each=NCOMB)
boxOverlap <- data.frame(BoxOverlap,boxSubj,rep(IndCenter,NSTEP))
	names(boxOverlap) <- c("Overlap", 'SubjectSize','EqualCenter')
	boxOverlap$EqualCenter <- factor(boxOverlap$EqualCenter)


# Actual plot of all individual points
ggplot(boxOverlap, aes(x=factor(SubjectSize),y=Overlap)) + 
	geom_point(aes(color=EqualCenter),position='jitter') +
	scale_colour_manual(values = c('#542788', '#b35806'),name="Between (0) or within center (1)") +
	ggtitle(paste("Within sample size overlap. Controlled for center. Resampling runs = ",NRUNS, ". Combinations = ",NCOMB,sep="")) +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) +
	  			scale_x_discrete(breaks=subjBreak, name="Sample size")


# Or in boxplot format
ggplot(boxOverlap, aes(x=factor(SubjectSize),y=Overlap))+ 
	geom_boxplot(aes(color=EqualCenter)) + 
	scale_color_manual(values = c('#542788', '#b35806'),name="Between (0) or within center (1)") +
	ggtitle(paste("Within sample size overlap. Controlled for center. Resampling runs = ",NRUNS, ". Combinations = ",NCOMB,sep="")) +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) +
	  			scale_x_discrete(breaks=subjBreak, name="Sample size")


IDequal <- boxOverlap[,'EqualCenter']==1
	equalcenters <- boxOverlap[IDequal,'Overlap']
IDnonequal <- boxOverlap[,'EqualCenter']==0
	nonEqualcenters <- boxOverlap[IDnonequal,'Overlap']

fitMain <- lm(Overlap~SubjectSize+EqualCenter,data=boxOverlap)
	summary(fitMain)
fitInt <- lm(Overlap~SubjectSize*EqualCenter,data=boxOverlap)
	summary(fitInt)


# Average within versus between center
AvgOverlapEq <- mean(equalcenters,na.rm=TRUE)
AvgOverlapNonEq <- mean(nonEqualcenters,na.rm=TRUE)
AvgOverlap <- aggregate(boxOverlap$Overlap, by = list(size = boxOverlap$SubjectSize, center = boxOverlap$EqualCenter), FUN = mean,na.rm=TRUE)

	ggplot(AvgOverlap, aes(x=factor(size),y=x)) + 
	geom_point(aes(color=center),position='identity') +
	scale_colour_manual(values = c('#542788', '#b35806'),name="Between (0) or within center (1)")	+
	ggtitle(paste("Average within sample size overlap. Controlled for center. Resampling runs = ",NRUNS, ". Combinations = ",NCOMB,sep="")) +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) +
	  			scale_x_discrete(breaks=subjBreak, name="Sample size")





















