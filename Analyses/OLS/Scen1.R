####################
#### TITLE:     Calculate within samplesize stability for increasing sample sizes in group study: scenario 1 - OLS pooling.
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

# This is all based on OLS pooling of subjects, scenario 1.

# Note, some objects are saved for scenario 3 as well


##
###############
### Preparation
###############
##

rm(list=ls())

# Set WD: location of data
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
	# Print which libraries are needed
	print("Need packages: oro.nifti, fslr, lattice, ggplot2, reshape2, RColorBrewer, gridExtra and the functions.R file")


# Dimension of the brains
DIM <- c(53,63,46)

# Number of runs
NRUNS <- 15
# Number of steps
NSTEP <- 74




##
###############
### For loop: start with sample size S, calculate overlap between all images, then move to sample size S+1
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
		# Number of masked voxels
		BothMasks <- maskK1*maskK_1
			NumMasked <- sum(BothMasks)
	
		# Union of activated voxels, voxels in image K and voxels in image K-1
		Vjt <- length(sumMap[which(sumMap==2)])
		Vt <- length(imageK[which(imageK==1)])
		Vj <- length(imageK_1[which(imageK_1==1)])
		# Voxels different from both images
		VjtS <- length(minusMap[which(minusMap==0)])
		# Now calculate overlap
		GlobalOverlap[i] <- sum(GlobalOverlap[i],(round((Vjt)/(Vj + Vt - Vjt),6)),na.rm=TRUE)
		# Percentage overlap: over number of masked voxels
		PercOverlap[i] <- sum(PercOverlap[i],(round((VjtS/NumMasked),6)),na.rm=TRUE)
		# We want to know if Vjt, Vt or Vj is maybe responsible for the observed effect
		GlobalVjt[i] <- GlobalVjt[i]+Vjt
		GlobalVt[i] <- GlobalVt[i]+Vt
		GlobalVj[i] <- GlobalVj[i]+Vj

		# We want to have the SD of the overlap, so calculate a matrix as well
		MatrixOverlap[i,j] <- round((Vjt)/(Vj + Vt - Vjt),6)
		MatrixPercOverlap[i,j] <- round((VjtS/NumMasked),6)

		# Remove objects
		rm(Vjt,Vt,Vj,imageK,imageK_1,sumMap,idK1,maskK1,maskK_1,idK_1,NumMasked)
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

# Save objects
save(AverageOverlap,file=paste(wd,'/AvgOverlapScen1',sep=''))
save(AveragePercOverlap,file=paste(wd,'/PercOverlapScen1',sep=''))
save(MatrixOverlap,file=paste(wd,'/MatOverlapScen1',sep=''))
save(MatrixPercOverlap,file=paste(wd,'/PercMatOverlapScen1',sep=''))







 










##
###############
### For loop to calculate average percentage of voxels active in each sample size
###############
##

# Array of percentage of significant voxels
PercSignVox <- array(NA,dim=c(NSTEP,NRUNS))

for(i in 1:NRUNS){
	print(paste("At run ",i,sep=""))
	for(j in 1:NSTEP){
		print(paste("Step ",j,sep=""))
		# Read in the data from image K: 0 is non significant and 1 is significant
		imageK <- readNIfTI(paste(wd,'/Run_',i,'/Step_',j,'/thresh_zstat1.nii',sep=''))[,,]
		# Read in the mask
		mask <- readNIfTI(paste(wd,'/Run_',i,'/Step_',j,'/mask.nii',sep=''))[,,,1]

		# Voxels active in image K and mask
		Vt <- length(imageK[which(imageK==1)])
		SumMask <- sum(mask)

		# Percentage of masked voxels active
		PercSignVox[j,i] <- Vt/SumMask

		# Remove objects
		rm(imageK,mask,Vt,SumMask)
	}
}




## Calculate a rough number for percentage of significant voxels
PercSignVox <- ((AverageVj+AverageVt)/2)/sum(mask)


##
###############
### Save the objects created so far
###############
##

save(MatrixOverlap, file=paste(wd,'/MatrixOverlap',sep=''))
save(PercSignVox, file=paste(wd,'/PercSignVox',sep=''))

load(paste(wd,'/MatrixOverlap',sep=''))
load(paste(wd,'/PercSignVox',sep=''))


##
###############
### Create plots
###############
##


#### Start with Jaccard overlap index
# Make qqplot
subjects <- seq(10,740,by=10)
overlap <- data.frame(subjects,AverageOverlap,PercSignVox)
	names(overlap) <- c("SampleSize", "Overlap","PercentageSign")


subjBreak <- seq(10,740,by=50)
qplot(SampleSize, Overlap, data=overlap, geom = c("point","smooth"),
	xlab = "Sample size", ylab = "Overlap",
  main = "Within sample size overlap for increasing sample size.") + scale_x_continuous(breaks=subjBreak)


# Or with ggplot
scatPlot <- ggplot(overlap, aes(x = SampleSize, y = Overlap, color = Overlap > .8)) + geom_point() + geom_smooth(aes(group = 1),colour='#f7f7f7',method='loess') +
    scale_color_manual(values = c('#542788', '#b35806'), breaks = c(TRUE, FALSE), labels = c('> .8', '< .8'),name="Average overlap")+
    	ggtitle("Within sample size average overlap for increasing sample size.") +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) +
	  			scale_x_continuous(breaks=subjBreak, name="Sample size") +
	  			scale_y_continuous(name="Average overlap") +
	  			annotation_custom(grob=textGrob(paste("Based on ",NCOMB," within \n sample size images",sep="")),
	  				xmin=832,xmax=842,ymin=.25,ymax=.4)
# Code to override clipping
gt <- ggplot_gtable(ggplot_build(scatPlot))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)



# Make boxplots of MatrixOverlap
BoxOverlap <- matrix(MatrixOverlap,ncol=1)
	boxSubj <- rep(seq(10,740,by=10),each=NCOMB)
boxOverlap <- data.frame(BoxOverlap,boxSubj,rep(AverageOverlap,each=NCOMB))
	names(boxOverlap) <- c("Overlap", 'SubjectSize','AvgOverlap')


Plot <- ggplot(boxOverlap, aes(x=factor(SubjectSize),y=Overlap,fill=AvgOverlap>.8)) + 
	geom_boxplot() + geom_smooth(aes(group = 1),colour='#800026',method='loess') +
	scale_fill_manual(values = c('#2b8cbe', '#016c59'), breaks = c(TRUE, FALSE), labels = c('> .8', '< .8'),name="Average overlap")+
	ggtitle(paste("Within sample size overlap. Number of resampling runs = ",NRUNS,sep="")) +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) +
	  			scale_x_discrete(breaks=subjBreak, name="Sample size") +
	  			annotation_custom(grob=textGrob(paste("Based on ",NCOMB," within \n sample size images",sep="")),
	  				xmin=75,xmax=85,ymin=.25,ymax=.4)

# Code to override clipping
gt <- ggplot_gtable(ggplot_build(Plot))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)



# With percentage of voxels on top of it
ggplot() + geom_point(aes(x=SampleSize, y=Overlap, colour=Overlap>.8), overlap) + 
	geom_smooth(aes(x=SampleSize, y=Overlap,group = 1),data=overlap,colour='#f7f7f7',method='loess') +
  		geom_line(aes(x=SampleSize, y=PercentageSign), overlap) +
  		scale_color_manual(values = c('#542788', '#b35806'), breaks = c(TRUE, FALSE), labels = c('> .8', '< .8'),name="Average overlap") 






#### Percentage based overlap
# Make qqplot
subjects <- seq(10,740,by=10)
PercOverlap <- data.frame(subjects,AveragePercOverlap)
	names(PercOverlap) <- c("SampleSize", "PercentageOverlap")


subjBreak <- seq(10,740,by=50)
qplot(SampleSize, PercentageOverlap, data=PercOverlap, geom = c("point","smooth"),
	xlab = "Sample size", ylab = "Percentage Overlap based on Masked Images",
  main = "Within sample size overlap for increasing sample size.") + scale_x_continuous(breaks=subjBreak)


# Or with ggplot
scatPlot <- ggplot(PercOverlap, aes(x = SampleSize, y = PercentageOverlap, color = PercentageOverlap > .96)) + geom_point() + geom_smooth(aes(group = 1),colour='#f7f7f7',method='loess') +
    scale_color_manual(values = c('#542788', '#b35806'), breaks = c(TRUE, FALSE), labels = c('> .96', '< .96'),name="Average overlap")+
    	ggtitle("Within sample size average overlap for increasing sample size.") +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) +
	  			scale_x_continuous(breaks=subjBreak, name="Sample size") +
	  			scale_y_continuous(name="Average overlap") +
	  			annotation_custom(grob=textGrob(paste("Based on ",NCOMB," within \n sample size images",sep="")),
	  				xmin=832,xmax=842,ymin=.25,ymax=.4)
# Code to override clipping
gt <- ggplot_gtable(ggplot_build(scatPlot))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


# Make boxplots of MatrixOverlap
BoxPercOverlap <- matrix(MatrixPercOverlap,ncol=1)
	boxSubj <- rep(seq(10,740,by=10),each=NCOMB)
boxPercOverlap <- data.frame(BoxPercOverlap,boxSubj,rep(AveragePercOverlap,each=NCOMB))
	names(boxPercOverlap) <- c("PercOverlap", 'SubjectSize','AvgOverlap')


Plot <- ggplot(boxPercOverlap, aes(x=factor(SubjectSize),y=PercOverlap,fill=AvgOverlap>.96)) + 
	geom_boxplot() + geom_smooth(aes(group = 1),colour='#800026',method='loess') +
	scale_fill_manual(values = c('#2b8cbe', '#016c59'), breaks = c(TRUE, FALSE), labels = c('> .96', '< .96'),name="Average overlap")+
	ggtitle(paste("Within sample size overlap. Number of resampling runs = ",NRUNS,sep="")) +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) +
	  			scale_x_discrete(breaks=subjBreak, name="Sample size") +
	  			annotation_custom(grob=textGrob(paste("Based on ",NCOMB," within \n sample size images",sep="")),
	  				xmin=75,xmax=85,ymin=.25,ymax=.4)

# Code to override clipping
gt <- ggplot_gtable(ggplot_build(Plot))
gt$layout$clip[gt$layout$name=="panel"] <- "off"
grid.draw(gt)


















