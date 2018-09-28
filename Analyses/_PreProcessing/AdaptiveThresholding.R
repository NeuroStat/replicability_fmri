####################
#### TITLE:     Adaptive thresholding: analysis in which significance 
####            thresholding level is adapted in order to keep percentage 
####            of voxels at constant level.
#### Contents:
####
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Script/Analyses/AdaptiveThresholding
#### First Modified: 31/08/2015
#### Notes:
#################



##
###############
### Notes
###############
##

## Fixate percentage of activated voxels at around target of 20% (range is between 19 and 21).
## We do this by starting with a starting value of significance level and then updating p-value by steps of .01.
## If percentage < (target - 0.01) percentage, then increase p-value. If > (target + 0.01), decrease p-value.

## The data comes from the SplitSubjects script.

## We will directly save the intermediate results in this script!


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

# Seed
set.seed(11121990)

# Baseline, starting significance level
signLevel <- 0.05

##
###############
### Functions
###############
##

# Adaptive Thresholding
AdapThresholding <- function(pvals,DIM,signLevel,idMask, target = NULL){
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
			if(signLevel >= 0.995){
				signLevel <- signLevel + 0.00001
			}else{
			  # Using random values to increase the significance level
				signLevel <- signLevel + runif(1,min=.0001,max=0.005)
			}
		}
		if(percentage > (target + 0.01)){
			if(signLevel <= 0.005){
				signLevel <- signLevel - 0.00001
			}else{
			  # Using random values to decrease the significance level
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
MatrixOverlapAdap <- array(NA,dim=c(NSTEP,NRUNS))

# Vector where we will store the final percentage of activated masked voxels
	# As we have two images per step and run, we will take the average of the percentage.
PercActAdap <- array(NA,dim=c(NSTEP,NRUNS))

# Vector to store average of percentage activated masked voxels in each step and run
SignLevels <- array(NA,dim=c(NSTEP,NRUNS))

# Vector for activation map of first run only
SPMRun1G1 <- array(NA,dim=c(prod(DIM),NSTEP))
SPMRun1G2 <- array(NA,dim=c(prod(DIM),NSTEP))

# For loop over all runs
for(i in 1:NRUNS){
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
		# Now that we have z-statistic images with its masks, 
			# calculate p-values using baseline significance level
		pvalsG1 <- 1-pnorm(imageG1)
		pvalsG2 <- 1-pnorm(imageG2)

		# Function AdapThresholding: calculates adjusted p-values, then finds 
		# correct significance level according to target percentage.
		# The latter might be obtained in the MaxAdapThresh.R script, which
		# maximises the overlap given the target percentage.
		adapG1 <- AdapThresholding(pvalsG1,DIM,signLevel,idMaskG1, target = 0.3541)
			signLevel <- adapG1$signLevel
			perc <- adapG1$percentage

		# Image G2
		adapG2 <- AdapThresholding(pvalsG2,DIM,signLevel,idMaskG2, target = 0.3541)
			signLevel <- round(mean(c(signLevel,adapG2$signLevel)),6)
			perc <- round(mean(c(perc,adapG2$percentage)),4)


		###################################################################
		# Calculate overlap: summing images to know the voxels in both maps (1+1 = 2)
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
if(i==NRUNS) print("100% done")
}


##
###############
### Save intermediate results
###############
##

# Save matrix with overlap using adaptive thresholding
#saveRDS(MatrixOverlapAdap, paste(SaveLoc,'/MatrixOverlapAdap.rda', sep = ''))

# Save the percentage of activated voxels
#saveRDS(PercActAdap, paste(SaveLoc,'/PercActAdap.rda', sep = ''))

# Save the significance thresholding levels
#saveRDS(SignLevels, paste(SaveLoc,'/SignLevels.rda', sep = ''))



##
###############
### Create plots
###############
##

# Some plots, just for fun
quartz.options(width=18,height=12)

## Prepare matrix
Overlap.tmp <- matrix(MatrixOverlapAdap,ncol=1)
	sampleSize <- rep(seq(10,700,by=10),NRUNS)
Overlap <- data.frame('overlap' = Overlap.tmp, 'size' = sampleSize)

# Variables for plotting
subjBreak <- c(seq(0,100,by=20),seq(100,700,by=50))

#################
## plot the percentage of masked voxels being activated
#################

# Data frame with overlap and percentage
	perc <- matrix(PercActAdap,ncol=1)
OverPerc <- data.frame('Value' = rbind(Overlap.tmp,perc),
		'Size' = rep(sampleSize,2),
		'Type' = c(rep('Overlap',3500),rep('Percentage',3500)))
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
	annotate("text", y = .1, x = 60, label = paste('Correlation = ', corr,sep=''),size=5)

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
	annotate("text", y = .1, x = 60, label = paste('Correlation = ', corr,sep=''),size=5)


#################
## plot percentage, overlap and significance level
#################
	sign <- matrix(SignLevels,ncol=1)
	SignLevels <- data.frame(	'Value' = sign,
								'Size' = sampleSize,
								'Type' = rep('SignLevel', 3500))
		SignLevels$Type <- as.factor(SignLevels$Type)
		SignLevels$Value <- as.numeric(SignLevels$Value)
tmp <- aggregate(SignLevels,by=list(SignLevels$Size,SignLevels$Type),mean,na.rm=TRUE)
TMP <- tmp[,c(2:4)]
	names(TMP) <- c('Type', 'Value', 'Size')
AvgOverPercSign <- rbind(AvgOverPerc,TMP)


colours <- c('#8c510a', '#4393c3','#35978f')

ggplot(AvgOverPercSign, aes(x=factor(Size),y=Value,group=Type)) +
	geom_line(aes(colour=Type),size=1.3) +
	scale_colour_manual(values=colours,labels=c('Overlap','Percentage','Significance level')) +
	scale_x_discrete(breaks=subjBreak, name="Sample size") +
	scale_y_continuous(name='Value') +
	theme(plot.title = element_text(lineheight=.2, face="bold")) +
	ggtitle('Average overlap, percentage and significance level (uncorrected p-values).')


#################
## Heatmap of all steps in run 1: group 1.
#################

# Calculate amount of voxels activated over all sample sizes
heat.tmp <- apply(SPMRun1G1,c(1),sum,na.rm=TRUE)
	HeatMap <- array(heat.tmp,dim=DIM)

maskHeat <- readNIfTI(paste(dat,'/Run_1/Step_70/Group1/masked.nii.gz',sep=''))[,,]

# Histogram of all values
hist(HeatMap,breaks=70,freq=FALSE)
	# Histogram with log y-axis
	hist.data = hist(HeatMap, plot=F)
	hist.data$counts = log10(hist.data$counts)
	plot(hist.data, ylab='log10(Frequency)')

	# Histogram without 0
	hist(HeatMap[HeatMap>0],breaks=68,xlim=c(1,69),xlab='Amount of voxels over 70 sample sizes',
		main='Distribution of amount of times a voxel is activated over all sample sizes')

# Levelplot
levelplot(HeatMap)

# Try with levelplot
HeatMap[maskHeat==0] <- -1

colours <- colorRampPalette(c('#fff7fb','#ece2f0', '#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636'))(n = 70)

levelplot(HeatMap,main="Heatmap going from sample size 10 to 700: how many times is a voxel activated?", pretty=FALSE,strip=FALSE,
	col.regions=c("#f0f0f0",'#d9d9d9',colours),xlab='X',ylab='Y',
	xlim=c(0,DIM[1]),ylim=c(0,DIM[2]),at=c(-2:70),tick.number=71)





##
###############
### JSM 2016: plot
###############
##

# Might want to reload data!
sign <- matrix(SignLevels,ncol=1)
SignLevels <- data.frame(	'Value' = sign,
							'Size' = sampleSize,
							'Type' = rep('SignLevel', 3500))
	SignLevels$Type <- as.factor(SignLevels$Type)
	SignLevels$Value <- as.numeric(SignLevels$Value)
tmp <- aggregate(SignLevels,by=list(SignLevels$Size,SignLevels$Type),mean,na.rm=TRUE)
TMP <- tmp[,c(2:4)]
names(TMP) <- c('Type', 'Value', 'Size')
AvgOverPercSign <- rbind(AvgOverPerc,TMP)

colours <- c('#8c510a', '#4393c3','#35978f')

ggplot(AvgOverPercSign, aes(x=factor(Size),y=Value,group=Type)) +
geom_line(aes(colour=Type),size=1.8) +
scale_colour_manual(values=colours,labels=c('Modified Dice index','Percentage signficant voxels','Significance level')) +
scale_x_discrete(breaks=subjBreak, name="Sample size") +
scale_y_continuous(name='Value') +
theme(plot.title = element_text(lineheight=.2, face="bold"),
			axis.text=element_text(size=20, face="bold"),
			axis.title=element_text(size=24),
			legend.text=element_text(size=18),
			legend.key.size=unit(1.5, "cm"),
			legend.title=element_text(size=18))






















