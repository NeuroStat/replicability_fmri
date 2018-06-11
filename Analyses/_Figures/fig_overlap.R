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
# Overlap while fixing the percentage of significant voxels (called adaptive thresholding)

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

# Number of runs
NRUNS <- 50

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

# Data with overlap using adaptive thresholding, its percentage of activated voxels
# and the significance thresholding levels
MatrixOverlapAdap <- readRDS(paste(LocIntRes,'/MatrixOverlapAdap.rda', sep = ''))
PercActAdap <- readRDS(paste(LocIntRes,'/PercActAdap.rda', sep = ''))
SignLevels <- readRDS(paste(LocIntRes,'/SignLevels.rda', sep = ''))

# Data with the Pearson correlation 
MatrixCorrelation <- readRDS(paste(LocIntRes,'/MatrixCorrelation.rda', sep = ''))

# Vector with all sample sizes
sampleSize <- rep(seq(10,700,by=10), NRUNS)

# Variables for plotting
subjBreak <- c(seq(10,110,by=20), seq(150,700, by=50))

##
###############
### Create plots
###############
##

# Set window 
quartz.options(width=18,height=12)

#################
## Points for overlap with smoothed regression line
#################

# Data frame and plot
Overlap <- data.frame('overlap' = matrix(MatrixOverlap,ncol=1),
                      'sampleSize' = sampleSize)
# Put NaN to zero
Overlap$overlap[is.na(Overlap$overlap)] <- 0

# Create plot
overlapPlot <- ggplot(Overlap, aes(x = sampleSize, y = overlap)) + 
  geom_point(size = 0.6) +
  geom_smooth(aes(x = sampleSize, y = overlap),
          method = 'loess', 
          formula = y ~ x,
          colour = '#fc8d59') +
  scale_x_continuous(breaks=subjBreak, name="Sample size") +
  scale_y_continuous(name='Overlap') +
  theme_bw()
overlapPlot

# Or using a boxplot
overlapBoxPlot <- ggplot(Overlap, aes(x=factor(sampleSize), y = overlap)) + 
  geom_boxplot(outlier.size = .7) +
  scale_x_discrete(breaks = subjBreak, name="Sample size") +
  scale_y_continuous(name=expression(Overlap~~ (omega))) +
  theme_bw()
overlapBoxPlot

# Version for OHBM 2018
subjBreak <- c(seq(10,110,by=30), seq(150,700, by=50))
overlapBoxPlot <- ggplot(Overlap, aes(x=factor(sampleSize), y = overlap)) + 
  geom_boxplot(outlier.size = .7, outlier.color = 'orange') +
  scale_x_discrete(breaks = subjBreak, name="Sample size") +
  scale_y_continuous(name=expression(Overlap~~(omega))) +
  labs(subtitle = 'FDR = 0.05') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
overlapBoxPlot


#################
## Points for overlap and ADAPTIVE THRESHOLDING 
#################

# Prepare the data frame for plotting
OverlapAdapt <- data.frame('Value' = c(matrix(MatrixOverlapAdap, ncol = 1),
                         matrix(PercActAdap, ncol = 1),
                         matrix(SignLevels, ncol = 1)),
             'Type' = c(rep('Overlap', NRUNS * NSTEP),
                        rep('Percentage', NRUNS * NSTEP),
                        rep('SignLevels', NRUNS * NSTEP)),
             'size' = rep(sampleSize, times = 3)) %>% as.tibble()

# Plot overlap, % of masked voxels being activated and significance level
AdaptOverlap <- OverlapAdapt %>%
  mutate(SampleSize = size) %>% 
ggplot(., aes(x = SampleSize, y = Value, group = Type)) + 
  geom_point(aes(colour = Type), size = 0.6, alpha = 0.75) +
  stat_summary(aes(group = Type), fun.y = mean, 
               geom="line", colour = 'black', size = .6) +
  scale_color_manual('', values = c('#1b9e77','#d95f02','#7570b3'),
                     labels = c('Overlap', '% of activated voxels', 'Significance level')) +
  scale_y_continuous('') + 
  scale_x_continuous('Sample size', breaks = subjBreak) +
  theme_bw() +
  # Increase size in legend
  guides(colour = guide_legend(override.aes = list(size=3)),
         alpha = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = 'bottom')
AdaptOverlap  


#################
## Points for Pearson correlation
#################

## Prepare matrix
Correlation <- data.frame('PearsonCorr' = matrix(MatrixCorrelation,ncol=1), 
                          'sampleSize' = sampleSize)

corrPlot <- ggplot(Correlation, aes(x = sampleSize, y = PearsonCorr)) +
  geom_point(colour='black',size = 0.6) +
  geom_smooth(aes(x = sampleSize, y = PearsonCorr),
              method = 'loess', 
              formula = y ~ x,
              colour = '#fc8d59') +
  scale_x_continuous(breaks = subjBreak, name="Sample size") +
  scale_y_continuous(name='Pearson product-moment correlation coefficient') +
  theme_bw()
corrPlot

# Or using a boxplot
corrBoxPlot <- ggplot(Correlation, aes(x=factor(sampleSize), y = PearsonCorr)) + 
  geom_boxplot(outlier.size = .7) +
  scale_x_discrete(breaks = subjBreak, name="Sample size") +
  scale_y_continuous(name='Pearson product-moment correlation coefficient') +
  theme_bw()
corrBoxPlot

# Version for OHBM 2018
subjBreak <- c(seq(10,110,by=30), seq(150,700, by=50))
corrBoxPlot <- ggplot(Correlation, aes(x=factor(sampleSize), y = PearsonCorr)) + 
  geom_boxplot(outlier.size = .7, outlier.colour = 'orange') +
  scale_x_discrete(breaks = subjBreak, name="Sample size") +
  scale_y_continuous(name=expression(Pearson~ product~-~moment~ correlation~ coefficient~~(rho))) +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
corrBoxPlot

##
###############
### Save plots
###############
##

# Overlap plot: points
ggsave(filename = paste0(getwd(), '/overlapPlot.png'),
       plot = overlapPlot,
       width = 20, height = 14, units = 'cm', scale = 0.9)

# Overlap plot: boxplots
ggsave(filename = paste0(getwd(), '/overlapBoxPlot.png'),
       plot = overlapBoxPlot,
       width = 20, height = 14, units = 'cm', scale = 0.9)

# Adaptive overlap plot
ggsave(filename = paste0(getwd(), '/AdaptOverlap.png'),
       plot = AdaptOverlap,
       width = 20, height = 14, units = 'cm', scale = 0.9)

# Pearson product moment correlation coefficient
ggsave(filename = paste0(getwd(), '/corrPlot.png'),
       plot = corrPlot,
       width = 20, height = 14, units = 'cm', scale = 0.9)

# Correlation: boxplot
ggsave(filename = paste0(getwd(), '/corrBoxPlot.png'),
       plot = corrBoxPlot,
       width = 20, height = 14, units = 'cm', scale = 0.9)



##
###############
### Some extra plots
###############
##

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
                       'Type' = c(rep('Overlap',3500),rep('Percentage',3500)))
OverPerc$Type <- as.factor(OverPerc$Type)
OverPerc$Value <- as.numeric(OverPerc$Value)

# Correlation
corr <- round(cor(perc, Overlap.tmp, use = "complete.obs"),2)
ggplot(OverPerc, aes(x=factor(Size),y=Value)) + 
  geom_point(aes(colour=Type),position='identity',size=1.5) +
  scale_x_discrete(breaks=subjBreak, name="Sample size") +
  scale_y_continuous(name='Overlap') +
  theme(plot.title = element_text(lineheight=.2, face="bold")) +
  ggtitle('Overlap and percentage of masked voxels active.') +
  annotate("text", y = .1, x = 50, label = paste('Correlation = ', corr,sep=''),size=5)

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
  annotate("text", y = .1, x = 50, label = paste('Correlation = ', corr,sep=''),size=5)


