####################
#### TITLE:     Calculate within sample size stability with increasing sample sizes: scenario B. OLS - POOLING
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/OLS
#### First Modified: 25/09/2015
#### Notes: 
#################



##
###############
### Notes
###############
##

# Pooling based on OLS!!

# We need to load in the unthresholded map, because we initially created thresholded maps at P < 0.001 with FDR.
# Now we would like to do this with P < 0.05

# Scenario B:
	# Completely split up subjects.
	# Each run we have an analysis between group 1 and 2 in a specific sample size.
	# The subjects in group 1 are completely different from those in group 2
	# They have to be from a different center





##
###############
### Preparation
###############
##

# Reset working directory
rm(list=ls())

# Set WD
wd <- "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/OLS/Scenario7"
setwd(wd)


# Load in libraries
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

# Number of runs (so far)
NRUNS <- 42
# Number of steps
NSTEP <- 14



# Load Center information
load('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/IMAGEN/IMAGEN/IMAGEN/OurIDs')
head(OurIDs)



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

# New significance level
signLevel <- 0.05

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
		#-------------------------------------------------------------------------------------------------------------#
		imageG1 <- try(readNIfTI(paste(wd,'/Run_',i,'/Step_',j,'/Group1','/stats/tstat1.nii',sep=''))[,,],silent=TRUE)
			if(!class(imageG1)=='array') next
			maskG1 <- try(readNIfTI(paste(wd,'/Run_',i,'/Step_',j,'/Group1','/masked.nii.gz',sep=''))[,,],silent=TRUE)
				if(!class(maskG1)=='array') next
			# Mask the image
			idMaskG1 <- maskG1==0
			imageG1[idMaskG1] <- NA
			# Now re-threshold: 1 = significant, 0 otherwise
				# Get df
				df <- as.numeric(read.table(paste(wd,'/Run_',i,'/Step_',j,'/Group1','/stats/dof',sep='')))
			MTmap1 <- imageG1
			# Calculate P-values
			pvals1 <- 1-pt(MTmap1,df)	
			# Adjust P-values with Benjamin & Hochberg procedure
			adjPvals1 <- array(p.adjust(pvals1,method="BH"),dim=DIM)

			## Threshold at signLevel BH p-value: significant P-values get 1!
			idP1 <- adjPvals1 <= signLevel
			threshPval1 <- adjPvals1
			threshPval1[!idP1] <- 0
			threshPval1[idP1] <- 1
		# Re-name threshPval as imageG1
		imageG1 <- threshPval1

		#-------------------------------------------------------------------------------------------------------------#
		# Group 2	
		imageG2 <- try(readNIfTI(paste(wd,'/Run_',i,'/Step_',j,'/Group2','/stats/tstat1.nii',sep=''))[,,],silent=TRUE)
			if(!class(imageG2)=='array') next
			maskG2 <- try(readNIfTI(paste(wd,'/Run_',i,'/Step_',j,'/Group2','/masked.nii.gz',sep=''))[,,],silent=TRUE)
				if(!class(maskG2)=='array') next
			idMaskG2 <- maskG2==0
			imageG2[idMaskG2] <- NA
			# Now re-threshold: 1 = significant, 0 otherwise
			MTmap2 <- imageG2
			# Calculate P-values
			pvals2 <- 1-pt(MTmap2,df)	
			# Adjust P-values with Benjamin & Hochberg procedure
			adjPvals2 <- array(p.adjust(pvals2,method="BH"),dim=DIM)
			## Threshold at signLevel BH p-value: significant P-values get 1!
			idP2 <- adjPvals2 <= signLevel
			threshPval2 <- adjPvals2
			threshPval2[!idP2] <- 0
			threshPval2[idP2] <- 1
		# Re-name threshPval as imageG1
		imageG2 <- threshPval2

		#-------------------------------------------------------------------------------------------------------------#
		# Summing image K and K-1 to know the voxels in both maps (1+1 = 2)
		sumMap <- imageG1+imageG2
		# Minus map: know the voxels different in both images
		minusMap <- imageG1-imageG2
			# Mask this image
		idMask <- maskG1*maskG2
			minusMap[idMask] <- NA

		#-------------------------------------------------------------------------------------------------------------#	
		# Union of activated voxels, voxels in image G1 and voxels in image G2
		Vjt <- length(sumMap[which(sumMap==2)])
		Vt <- length(imageG1[which(imageG1==1)])
		Vj <- length(imageG2[which(imageG2==1)])
		# Voxels different from both images
		VjtS <- length(minusMap[which(minusMap==0)])

		#-------------------------------------------------------------------------------------------------------------#
		# Put overlap in matrix
		MatrixOverlap[j,i] <- round((Vjt)/(Vj + Vt - Vjt),6)

		#-------------------------------------------------------------------------------------------------------------#
		# Now calculate percentage of activated masked voxels in imageG1 and imageG2
		baseG1 <- sum(maskG1)
		percG1 <- Vt/baseG1

		baseG2 <- sum(maskG2)
		percG2 <- Vt/baseG2
		PercAct[j,i] <- round(mean(c(percG1,percG2)),4)

		#-------------------------------------------------------------------------------------------------------------#
		# Remove objects
		rm(Vjt,Vt,Vj,imageG1,imageG2,sumMap,threshPval1,threshPval2,idP1,idP2,adjPvals1,adjPvals2,pvals1,pvals2,MTmap1,MTmap2,df)
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


## Save objects
save(MatrixOverlap,file=paste(wd,'/MatrixOverlapOLSScen7',sep=''))
save(PercAct,file=paste(wd,'/PercActOLSScen7',sep=''))

## Or load objects
load(paste(wd,'/MatrixOverlapOLSScen7',sep=''))
load(paste(wd,'/PercActOLSScen7',sep=''))

##
###############
### Create plots
###############
##

quartz.options(width=18,height=12)


## Prepare matrix
Overlap.tmp <- matrix(MatrixOverlap,ncol=1)
	sampleSize <- rep(seq(10,140,by=10),NRUNS)
Overlap <- data.frame('overlap' = Overlap.tmp, 'size' = sampleSize)

# Variables for plotting
subjBreak <- seq(0,140,by=10)

#################
## ggplot: points
#################
ggplot(Overlap, aes(x=factor(sampleSize),y=overlap)) + 
	geom_point() +
	scale_x_discrete(breaks=subjBreak, name="Sample size")


#################
## ggplot: boxplot
#################
ggplot(Overlap, aes(x=factor(size),y=overlap))+ 
	geom_boxplot(colour = '#1f18b4',fill='#1f78b4',width=.5) + 
	ggtitle('Boxplot of scenario B') +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) +
	  			scale_x_discrete(breaks=subjBreak, name="Sample size")



#################
## ggplot: plot the amount of observations/sample size
#################
# how many observations are there?
observations <- apply(MatrixOverlap,c(1),function(x){length(na.omit(x))})
	# put in data frame
	obs <- data.frame('count' = observations, 'size' = seq(10,700,by=10))

ggplot(obs, aes(x=factor(size), y=count))+
	geom_bar(stat='identity', colour='#1f78b4',fill='#a6cee3') +
	scale_x_discrete(breaks=seq(0,700,by=50), name="sample size") +
	ggtitle('Incomplete data: amount of comparissons in each sample size.') +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) + 
	  		scale_y_continuous(name='Amount of data points')+
	  		theme_bw()





#################
## plot the percentage of masked voxels being activated
#################


# Data frame with overlap and percentage
	perc <- matrix(PercAct,ncol=1)
OverPerc <- data.frame('Value' = rbind(Overlap.tmp,perc),
		'Size' = rep(sampleSize,2), 
		'Type' = c(rep('Overlap',588),rep('Percentage',588)))
	OverPerc$Type <- as.factor(OverPerc$Type)
	OverPerc$Value <- as.numeric(OverPerc$Value)

# Correlation
corr <- round(cor(na.omit(perc),na.omit(Overlap.tmp)),2)

ggplot(OverPerc, aes(x=factor(Size),y=Value)) + 
	geom_point(aes(colour=Type),position='identity',size=1.5) +
	scale_x_discrete(breaks=subjBreak, name="Sample size") +
	scale_y_continuous(name='Overlap') +
	theme(plot.title = element_text(lineheight=.2, face="bold")) +
	ggtitle('Overlap and percentage of masked voxels active.') +
	annotate("text", y = .55, x = 2, label = paste('Correlation = ', corr,sep=''),size=5)



# Calculate averages
tmp <- aggregate(OverPerc,by=list(OverPerc$Size,OverPerc$Type),mean,na.rm=TRUE)
AvgOverPerc <- tmp[,c(2:4)]
	names(AvgOverPerc) <- c('Type', 'Value', 'Size')

ggplot(AvgOverPerc, aes(x=factor(Size),y=Value,group=Type)) + 
	geom_line(aes(colour=Type)) +
	geom_point(aes(colour=Type),position='identity',size=1.5) +
	scale_x_discrete(breaks=subjBreak, name="Sample size") +
	scale_y_continuous(name='Overlap') +
	theme(plot.title = element_text(lineheight=.2, face="bold")) +
	ggtitle('Average overlap and percentage of masked voxels active.')+
	annotate("text", y = .55, x = 2, label = paste('Correlation = ', corr,sep=''),size=5)

































































