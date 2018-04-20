####################
#### TITLE:     Plot within sample size overlap (Dice + Pearson) with increasing sample sizes
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/SplitSubjects
#### First Modified: 19/04/2018
#### Notes: 
#################



##
###############
### Notes
###############
##


# Data analysis:
# Completely split up subjects.
# Each run we have an analysis between group 1 and 2 in a specific sample size.
# The subjects in group 1 are completely different from those in group 2

# Here we read in the preprocessed data and then plot:
# Dice overlap (adapation from Maitra)
# Pearson correlation
# Overlap while fixing the percentage of significant voxels 

##
###############
### Preparation
###############
##

# Location of intermediate results
LocIntRes <- '../_IntData/'

# Load in libraries
library(tidyverse)
library(magrittr)
library(cowplot)
library(lattice)
library(gridExtra)
library(oro.nifti)
library(reshape2)
library(RColorBrewer)
library(Hmisc)
library(NeuRRoStat)

# Dimension of the brains
DIM <- c(53,63,46)

# Number of runs (so far)
NRUNS <- 42
# Number of steps
NSTEP <- 70

##
###############
### Read in data
###############
##

# Data with overlap values and percentage of activated voxels
MatrixOverlap <- readRDS(paste(LocIntRes,'/MaitraOverlap.rda',sep=''))
PercAct <- readRDS(paste(LocIntRes,'/PercActMaitraOverlap.rda',sep=''))


##
###############
### Create plots
###############
##

quartz.options(width=18,height=12)


## Prepare matrix
Overlap.tmp <- matrix(MatrixOverlap,ncol=1)
sampleSize <- rep(seq(10,700,by=10),NRUNS)
Overlap <- data.frame('overlap' = Overlap.tmp, 'size' = sampleSize)

# Variables for plotting
subjBreak <- c(seq(0,100,by=20),seq(100,700,by=50))

#################
## ggplot: points
#################
ggplot(Overlap, aes(x=factor(sampleSize),y=overlap)) + 
  geom_point() +
  scale_x_discrete(breaks=subjBreak, name="Sample size")


# Check with earlier results: load in, rename and reload this version
load("/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Script.git/Analyses/MatrixOverlap")
MatrixOverlapComb <- MatrixOverlap
# Reload results from subjectsplit
load(paste(wd,'/MatrixOverlapSplitSubjA',sep=''))


# Prepare the data, previous data had larger sample size, so we get rid of these
# Then we add the sample sizes
# we make a data frame and bind the data from this version to it (with splitted subjects)
# We add a variable with 1 for this version and 2 for previous version
MatrixOverlapCombn <- MatrixOverlapComb[,-c(71:74)]
MatrixOverlapCombn <- t(MatrixOverlapCombn)
ArrayComb <- matrix(MatrixOverlapCombn,ncol=1)
combSampleSize <- rep(seq(10,700,by=10),dim(MatrixOverlapCombn)[2])
combOverlap <- data.frame('overlap' = ArrayComb, 'size' = combSampleSize)

BothOverlap <- rbind(Overlap,combOverlap)
originData <- c(rep(2,length(sampleSize)),rep(1,length(combSampleSize)))
BothOverlap <- data.frame(BothOverlap, 'origin' = originData)
BothOverlap$origin <- factor(BothOverlap$origin)



## ggplot2
ggplot(BothOverlap, aes(x=factor(size),y=overlap)) + 
  geom_point(aes(colour=origin, alpha=origin),position='identity') +
  scale_x_discrete(breaks=subjBreak, name="Sample size") +
  scale_y_continuous(name='Overlap') +
  scale_alpha_discrete(name='Restrictions', labels = c('No restriction (prev)','No overlapping subjects (curr)')) +
  scale_colour_discrete(name='Restrictions', labels = c('No restriction (prev)','No overlapping subjects (curr)')) +
  theme(plot.title = element_text(lineheight=.2, face="bold")) +
  ggtitle('Comparing previous with current sampling frame')



#################
## ggplot: boxplot
#################
ggplot(BothOverlap, aes(x=factor(size),y=overlap))+ 
  geom_boxplot(aes(color=origin),outlier.size = 1,lwd=.9) + 
  scale_color_manual(values = c('#1f78b4', '#b2df8a'),name="Restrictions", labels=c('None', 'No overlapping subjects')) +
  ggtitle('Comparing previous with current sampling frame') +
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
                       'Type' = c(rep('Overlap',2940),rep('Percentage',2940)))
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

