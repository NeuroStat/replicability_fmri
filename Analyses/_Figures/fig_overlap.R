####################
#### TITLE:     Plot within sample size overlap (Dice + Pearson) with increasing sample sizes
#### Contents: 	
#### 
#### Source Files: /replicability_fmri/SplitSubjects
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
subjBreak <- c(seq(10,110,by=30), seq(150,700, by=50))

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

# Version for OHBM 2018 (and paper!)
subjBreak <- c(seq(10,110,by=30), seq(150,700, by=50))
overlapBoxPlot <- ggplot(Overlap, aes(x=factor(sampleSize), y = overlap)) + 
  geom_boxplot(outlier.size = .7, outlier.color = 'orange') +
  scale_x_discrete(breaks = subjBreak, name="Sample size") +
  scale_y_continuous(name=expression(Overlap~~(omega))) +
  labs(title = 'Conditional test-retest reliability',
       subtitle = 'FDR = 0.05') +
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

# Median overlap at N = 200
Overlap %>%
  filter(sampleSize == 200) %>%
  summarise(med = median(overlap))

# Maximum overlap
Overlap %>% as_tibble() %>%
  group_by(sampleSize) %>%
  summarise(AvOver = mean(overlap)) %>%
  filter(AvOver == max(AvOver))

# Median overlap at N = 30
Overlap %>%
  filter(sampleSize == 30) %>%
  summarise(med = median(overlap))



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

# What sample size corresponds to P < 0.001?
SSize001 <- OverlapAdapt %>%
  filter(Type == 'SignLevels') %>%
  group_by(size) %>%
  summarise(AvSignL = mean(Value)) %>%
  filter(AvSignL <= 0.001) %>%
  dplyr::slice(1) %>%
  dplyr::select(size) %>% as.numeric()
SSize001

# What is the median overlap associated with this sample size?
MedOv001 <- OverlapAdapt %>%
  filter(size == SSize001) %>%
  filter(Type == 'Overlap') %>%
  summarise(MedOv = median(Value)) %>%
  dplyr::select(MedOv) %>% as.numeric() %>% round(.,2)
MedOv001

# Median overlap at N = 20 and median P value
OverlapAdapt %>% 
  group_by(size, Type) %>%
  summarise(MedOv = median(Value)) %>%
  filter(size == 20)

# Maximum overlap and significant level
OverlapAdapt %>% 
  filter(Type == 'Overlap') %>%
  group_by(size, Type) %>%
  summarise(MedOv = median(Value)) %>%
  ungroup() %>%
  filter(MedOv == max(MedOv))
OverlapAdapt %>% 
  group_by(size, Type) %>%
  summarise(MedOv = median(Value)) %>%
  filter(size == 700)

# Plot overlap, % of masked voxels being activated and significance level
AdaptOverlap <- OverlapAdapt %>%
  mutate(SampleSize = size) %>% 
ggplot(., aes(x = SampleSize, y = Value, group = Type)) + 
  geom_point(aes(shape = Type, colour = Type), size = 0.8, alpha = 0.75) +
  geom_vline(xintercept = SSize001, linetype = 'dashed') +
  geom_hline(yintercept = MedOv001, linetype = 'dashed') +
  stat_summary(aes(group = Type), fun.y = mean, 
               geom="line", colour = 'black', size = .6) +
#  scale_color_manual('', values = c('#000000','#525252','#969696'),
#                    labels = c('Overlap', 'x100 % of activated voxels', 'Significance level')) +
  scale_color_brewer('', type = 'qual', palette = 2,
                    labels = c('Overlap', 'x100 % of activated voxels', 'Significance level')) +
  scale_shape_manual('', values = c(17,18,1),
                     labels = c('Overlap', 'x100 % of activated voxels', 'Significance level')) +
  scale_y_continuous('', breaks = sort(c(seq(0,0.8,by = 0.2), MedOv001))) + 
  scale_x_continuous('Sample size', breaks = subjBreak) +
  theme_bw() +
  # Increase size in legend
  guides(colour = guide_legend(override.aes = list(size=3)),
         alpha = guide_legend(override.aes = list(alpha = 1))) +
  labs(title = 'Conditional test-retest reliability',
       subtitle = 'Adaptive thresholding') +
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

# Version for OHBM 2018 (and paper!)
subjBreak <- c(seq(10,110,by=30), seq(150,700, by=50))
corrBoxPlot <- ggplot(Correlation, aes(x=factor(sampleSize), y = PearsonCorr)) + 
  geom_boxplot(outlier.size = .7, outlier.colour = 'orange') +
  scale_x_discrete(breaks = subjBreak, name="Sample size") +
  scale_y_continuous(name=expression(Pearson~ product~-~moment~ correlation~ coefficient~~(rho))) +
  theme_classic() +
  labs(title = 'Unconditional test-retest reliability') +
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

# Maximum value
Correlation %>%
  filter(PearsonCorr == max(PearsonCorr))

# Maximum median value
Correlation %>% as_tibble() %>%
  group_by(sampleSize) %>%
  summarise(MedPearson = round(median(PearsonCorr),3)) %>%
  filter(MedPearson == max(MedPearson))

# At N = 700?
Correlation %>% as_tibble() %>%
  group_by(sampleSize) %>%
  summarise(MedPearson = median(PearsonCorr)) %>%
  filter(sampleSize == 700)

# When do we have a median rho of 0.80?
Correlation %>% as_tibble() %>%
  group_by(sampleSize) %>%
  summarise(MedPearson = median(PearsonCorr)) %>%
  filter(MedPearson >= 0.80)

# Median rho at N = 30
Correlation %>% as_tibble() %>%
  group_by(sampleSize) %>%
  summarise(MedPearson = median(PearsonCorr)) %>%
  filter(sampleSize == 30)

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
### APPENDIX
###############
##

# Different FDR levels: average overlap

# These are the extra FDR levels considered:
FDRlevels <- c(0.001, 0.01, 0.1, 0.2)

# Read the data in and add to a new data frame
OverlapFDR <- data.frame('overlap' = matrix(MatrixOverlap,ncol=1),
                         'sampleSize' = sampleSize,
                         'FDR' = 0.05) %>% as_tibble()

# Loop over the levels
for(i in 1:length(FDRlevels)){
  # Identifier
  IDFDR <- sub(pattern = '.', replacement = '_', x = FDRlevels[i], fixed = TRUE)
  
  # Read in the data
  MatrixOverlapFDRtmp <- readRDS(paste(LocIntRes,'/MaitraOverlapFDR',IDFDR,'.rda',sep=''))
  
  # Convert to data frame and bind to final data frame
  OverlapFDR <- bind_rows(
    OverlapFDR, 
    data.frame('overlap' = matrix(MatrixOverlapFDRtmp, ncol = 1),
               'sampleSize' = sampleSize) %>% as_tibble() %>%
      mutate('FDR' = FDRlevels[i])
  )
}
# Change FDR to factor
OverlapFDR$FDR <- factor(OverlapFDR$FDR)

# Average overlap
OverlapFDR %>% 
  group_by(sampleSize, FDR) %>%
  summarise(AvOverl = mean(overlap)) %>%
  ungroup() %>%
  ggplot(., aes(x = sampleSize, y = AvOverl)) +
  geom_line(aes(colour = FDR), size = 1) + 
  scale_x_continuous(breaks = subjBreak, name="Sample size") +
  scale_y_continuous(name=expression(Overlap~~(omega))) +
  labs(title = 'Conditional test-retest reliability') +
  theme_classic()

# Plot using smoother
ggplot(OverlapFDR, aes(x=factor(sampleSize), y = overlap, group = factor(FDR))) + 
geom_smooth(aes(colour = factor(FDR))) +
scale_x_discrete(breaks = subjBreak, name="Sample size") +
scale_y_continuous(name=expression(Overlap~~(omega))) +
labs(title = 'Conditional test-retest reliability') +
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

# Boxplots
overlapBoxFDRs <- ggplot(OverlapFDR, aes(x=factor(sampleSize), y = overlap)) + 
  geom_boxplot(aes(fill = factor(FDR)), outlier.size = .5, 
               outlier.color = 'orange', width = .90) +
  scale_x_discrete(breaks = subjBreak, name="Sample size") +
  scale_y_continuous(name=expression(Overlap~~(omega))) +
  labs(title = 'Conditional test-retest reliability') +
  scale_fill_brewer('FDR control at ', type = 'qual', palette = 6) +
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
overlapBoxFDRs

# Overlap plot: boxplots
ggsave(filename = paste0(getwd(), '/overlapBoxFDRs.png'),
       plot = overlapBoxFDRs,
       width = 28, height = 18, units = 'cm', scale = 1.5)


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
OverPerc <- data.frame('Value' = c(Overlap$overlap,perc),
                       'Size' = rep(sampleSize,2), 
                       'Type' = c(rep('Overlap',3500),rep('Percentage',3500)))
OverPerc$Type <- as.factor(OverPerc$Type)
OverPerc$Value <- as.numeric(OverPerc$Value)

# Correlation
corr <- round(cor(perc, Overlap$overlap, use = "complete.obs"),2)
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

# Mimic figure from paper
subjBreak <- c(seq(10,110,by=30), seq(150,700, by=50))
overlap_percBoxPlot <- ggplot(OverPerc, aes(x=factor(Size), y = Value)) + 
  geom_boxplot(aes(fill = Type), outlier.size = .7, outlier.color = 'orange') +
  scale_x_discrete(breaks = subjBreak, name="Sample size") +
  scale_y_continuous('Value') +
  scale_fill_brewer('Type of value', labels = c('overlap', 'proportion'),
                    type = 'qual', palette = 7) +
  labs(title = 'Overlap and proportion of significant voxels',
       subtitle = 'FDR = 0.05') +
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
overlap_percBoxPlot







