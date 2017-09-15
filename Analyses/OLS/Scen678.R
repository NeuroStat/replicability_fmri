####################
#### TITLE:     Calculate within sample size stability with increasing sample sizes: scenario 6, 7 and 8. OLS - POOLING
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/SplitSubjects
#### First Modified: 25/09/2015
#### Notes: 
#################



##
###############
### Notes
###############
##


# OLS pooling of subjects


##
###############
### Preparation
###############
##

# Reset working directory
rm(list=ls())

# Set WD
wd <- "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/OLS"
setwd(wd)


# Load in libraries
library(lattice)
library(gridExtra)
library(oro.nifti)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
source('~/Dropbox/PhD/PhDWork/Meta\ Analysis/R\ Code/Studie_FixRan/FixRanStudyGit.git/Development/functions.R')
	# Print which libraries are needed
	print("Need packages: oro.nifti, lattice, ggplot2, reshape2, RColorBrewer, gridExtra and the functions.R file")


# Dimension of the brains
DIM <- c(53,63,46)

# Information about the number of runs
NRUNSA <- 42
NRUNSB <- 42
NRUNSC <- 42

# Load Center information
load('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/IMAGEN/IMAGEN/IMAGEN/OurIDs')
head(OurIDs)



##
###############
### Load in all overlap data and create data frame
###############
##


load(paste(wd,'/Scenario6/MatrixOverlapOLSA',sep=''))
	OverlapA <- matrix(MatrixOverlap,ncol=1)
	scenA <- data.frame('Scenario' = rep('A',length(OverlapA)))
	subjA <- data.frame('Size' = rep(seq(10,700,by=10),NRUNSA))
		rm(MatrixOverlap)
load(paste(wd,'/Scenario7/MatrixOverlapOLSB',sep=''))
	OverlapB <- matrix(MatrixOverlap,ncol=1)
	scenB <- data.frame('Scenario' = rep('B',length(OverlapB)))
	subjB <- data.frame('Size' = rep(seq(10,140,by=10),NRUNSB))
		rm(MatrixOverlap)
load(paste(wd,'/Scenario8/MatrixOverlapOLSC',sep=''))
	OverlapC <- matrix(MatrixOverlap,ncol=1)
	scenC <- data.frame('Scenario' = rep('C',length(OverlapC)))
	subjC <- data.frame('Size' = rep(seq(10,70,by=10),NRUNSC))
		rm(MatrixOverlap)

AllScen <- data.frame('Overlap' = rbind(OverlapA,OverlapB,OverlapC), 
	'Scenario' = rbind(scenA,scenB,scenC),
	'Size' = rbind(subjA,subjB,subjC))



##
###############
### Load in data for percentages and create data frame
###############
##


load(paste(wd,'/Scenario6/PercActOLSA',sep=''))
	PercA <- matrix(PercAct,ncol=1)
	scenA <- data.frame('Scenario' = rep('A',length(PercA)))
	subjA <- data.frame('Size' = rep(seq(10,700,by=10),NRUNSA))
		rm(PercAct)
load(paste(wd,'/Scenario7/PercActOLSB',sep=''))
	PercB <- matrix(PercAct,ncol=1)
	scenB <- data.frame('Scenario' = rep('B',length(PercB)))
	subjB <- data.frame('Size' = rep(seq(10,140,by=10),NRUNSB))
		rm(PercAct)
load(paste(wd,'/Scenario8/PercActOLSC',sep=''))
	PercC <- matrix(PercAct,ncol=1)
	scenC <- data.frame('Scenario' = rep('C',length(PercC)))
	subjC <- data.frame('Size' = rep(seq(10,70,by=10),NRUNSC))
		rm(PercAct)

PercAllScen <- data.frame('Percentage' = rbind(PercA,PercB,PercC), 
	'Scenario' = rbind(scenA,scenB,scenC),
	'Size' = rbind(subjA,subjB,subjC))






##
###############
### Create plots
###############
##

quartz.options(width=18,height=12)


# Variables for plotting
subjBreak <- c(seq(0,100,by=20),seq(100,700,by=50))

# Labels
labls <- c('A: no restrictions','B: different centra','C: same centra')

# Colours
colour <- c('#377eb8','#4daf4a','#984ea3')




#################
## ggplot: points
#################
ggplot(AllScen, aes(x=factor(Size),y=Overlap)) + 
	geom_point(aes(colour=Scenario, alpha=Scenario),position='jitter') +
	scale_alpha_manual(values = c(.8,.6,.8),labels=labls) +
	scale_colour_manual(values = colour, labels=labls) +
	scale_x_discrete(breaks=subjBreak, name="Sample size") +
	ggtitle('Overlap for scenario A, B and C') +
	theme(plot.title = element_text(lineheight=.2, face="bold"))


# Calculate averages
tmp <- aggregate(AllScen,by=list(AllScen$Size,AllScen$Scenario),mean,na.rm=TRUE)
AvgAllScen <- tmp[,c(1:3)]
	names(AvgAllScen) <- c('Size', 'Scenario', 'Overlap')

ggplot(AvgAllScen, aes(x=factor(Size),y=Overlap, group=Scenario)) + 
	geom_point(aes(colour=Scenario)) +
	geom_line(aes(colour=Scenario),size=.2) +
	scale_colour_manual(values = colour, labels=labls) +
	scale_x_discrete(breaks=subjBreak, name="Sample size") +
	ggtitle('Overlap for scenario A, B and C') +
	theme(plot.title = element_text(lineheight=.2, face="bold"))




#################
## ggplot: boxplot
#################
ggplot(AllScen, aes(x=factor(Size),y=Overlap))+ 
	geom_boxplot(aes(color=Scenario),outlier.size = 1,lwd=.9) + 
	scale_color_manual(values = colour,name="Scenario", labels=labls) +
	ggtitle('Overlap for scenario A, B and C') +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) +
	  			scale_x_discrete(breaks=subjBreak, name="Sample size")


	# Take only sample sizes < 200
	ID <- AllScen$Size < 200
	AllScen200 <- AllScen[ID,]
	ggplot(AllScen200, aes(x=factor(Size),y=Overlap))+ 
		geom_boxplot(aes(color=Scenario),outlier.size = 1,lwd=1) + 
		scale_color_manual(values = colour,name="Scenario", labels=labls) +
		ggtitle('Overlap for scenario A, B and C') +
		  		theme(plot.title = element_text(lineheight=.2, face="bold")) +
		  			scale_x_discrete(breaks=seq(0,200,20), name="Sample size")





#################
## ggplot: plot the amount of observations/sample size
#################
# Counting the observations
AllScenObs <- AllScen
AllScenObs[,1] <- !is.na(AllScenObs[,1])
AllScenObs <- aggregate(AllScenObs$Overlap, by=list(AllScenObs$Size,AllScenObs$Scenario), sum)
	names(AllScenObs) <- c('Size', 'Scenario', 'Observations')

# Plot of each sample size
ggplot(AllScenObs, aes(x=factor(Size), y=Observations))+
	geom_bar(stat='identity', colour='#1f78b4',fill='#a6cee3') +
	facet_grid(~ Scenario) +
	scale_x_discrete(breaks=seq(0,700,by=50), name="Sample Size") +
	scale_y_continuous(name='Frequency') +
	ggtitle('Amount of observations in each scenario') +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) + 
	  				theme_bw()

# At least one observation:
ggplot(AllScen200, aes(x=factor(Size)))+
	geom_histogram(colour='#1f78b4',fill='#a6cee3') +
	facet_grid(~ Scenario) +
	scale_x_discrete(breaks=seq(0,200,by=20), name="Sample Size") +
	scale_y_continuous(name='Frequency') +
	ggtitle('Amount of observations in each scenario') +
	  		theme(plot.title = element_text(lineheight=.2, face="bold")) + 
	  				theme_bw()



#################
## plot the percentage of masked voxels being activated
#################

# Calculate averages
tmp <- aggregate(PercAllScen,by=list(PercAllScen$Size,PercAllScen$Scenario),mean,na.rm=TRUE)
AvgPercAllScen <- tmp[,c(1:3)]
	names(AvgPercAllScen) <- c('Size', 'Scenario', 'Percentage')

ggplot(AvgPercAllScen, aes(x=factor(Size),y=Percentage, group=Scenario)) + 
	geom_point(aes(colour=Scenario)) +
	geom_line(aes(colour=Scenario),size=.2) +
	scale_colour_manual(values = colour, labels=labls) +
	scale_x_discrete(breaks=subjBreak, name="Sample size") +
	ggtitle('Percentage activated masked voxels for scenario A, B and C') +
	theme(plot.title = element_text(lineheight=.2, face="bold"))





#################
## plot the percentage and overlap
#################

# First take only subject sizes < 200
ID <- AvgPercAllScen$Size < 200
AllAvgScenPerc200 <- AvgPercAllScen[ID,]
	names(AllAvgScenPerc200) <- c('Size', 'Scenario', 'Value')
AvgAllScen200 <- AvgAllScen[ID,]
	names(AvgAllScen200) <- c('Size', 'Scenario', 'Value')

	# Now bind them together
	BothScen200 <- data.frame(rbind(AvgAllScen200,AllAvgScenPerc200),
		'Type' = rep(c(1,2),each=length(AvgAllScen200$Value)))
	BothScen200$Type <- factor(BothScen200$Type, labels=c('Overlap', 'Percentage'))


# Plot
ggplot(BothScen200, aes(x=factor(Size),y=Value, group=Scenario)) + 
	geom_point(aes(colour=Scenario)) +
	geom_line(aes(colour=Scenario),size=.5) +
	facet_grid(~ Type) +
	scale_colour_manual(values = colour, labels=labls) +
	scale_x_discrete(breaks=seq(0,200,20), name="Sample size") +
	ggtitle('Overlap and percentage activated masked voxels for scenario A, B and C') +
	theme(plot.title = element_text(lineheight=.2, face="bold"))





















































































