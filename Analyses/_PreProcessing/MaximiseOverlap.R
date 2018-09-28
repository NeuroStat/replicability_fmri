####################
#### TITLE:  Run the estimation algorithm to obtain m1 and plot some figures
#### Contents: 	
#### 
#### Source Files: /replicability_fmri/Analyses/_PreProcessing/
#### First Modified: 26/09/2018
#### Notes: 
#################


##
###############
### Notes
###############
##

# In the file: EstimateM1_maxOv.R, we created an algorithm to estimate
# m1 --> the proportion of truly active voxels.
# We tested the algorithm on simulated data (small grid of 50x50 voxels).
# Here we will read in the results and plot the functions.
# Then we need to obtain the maxima of the functions 
# and check whether results are unbiased.



##
###############
### Preparation
###############
##

# Location of raw data (not available on Github)
LocRawData <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/MaximOverlap'

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

# Dimension of the 2D slice
DIM2D <- c(50,50)

# Number of batches (each batch contains an amount of simulations)
nbatch <- 100


##
###############
### Read in data
###############
##

# Empty data frame
full_data <- data.frame() %>% as_tibble()

# Loop over the batches
for(i in 1:nbatch){
  full_data <- readRDS(file = paste(LocRawData, '/BackUPData/full_data_', i, '.rds', sep = '')) %>%
    bind_rows(full_data, .)
}


##
###############
### Plot the functions of the target percentages
###############
##

# Quick plot
full_data %>%
  filter(TruePerc != 0) %>%
  group_by(target_perc, TruePerc) %>%
  summarise(AvgObs_overlap = mean(obs_overlap)) %>%
  ungroup() %>%
  ggplot(., aes(x = target_perc, y = AvgObs_overlap)) +
    geom_line(size = .75) +
    geom_vline(aes(xintercept = TruePerc), size = .75) +
    facet_wrap( ~ TruePerc, ncol = 5) +
  scale_x_continuous('target % activated voxels in the algorithm') +
  scale_y_continuous('average observed overlap of activation') +
  labs(caption = 'vertical line = true % truly activated voxels in the image') +
    theme_classic() +
    theme(panel.grid.major = element_line(size = 0.8),
          panel.grid.minor = element_line(size = 0.8),
          axis.title.x = element_text(face = 'plain'),
          axis.title.y = element_text(face = 'plain'),
          axis.text = element_text(size = 11, face = 'plain'),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks = element_line(size = 1.3),
          axis.ticks.length=unit(.20, "cm"),
          axis.line = element_line(size = .75),
          title = element_text(face = 'plain'),
          plot.title = element_text(hjust = 0.5),
          strip.background = element_rect(colour="white", fill="white"),
          strip.text = element_text(size = 11),
          plot.subtitle = element_text(hjust = 0, vjust = -2),
          legend.key.size = unit(.85,'cm'),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 13))


##
###############
### Estimate m1
###############
##

# Bias of m1
full_data %>%
  filter(!TruePerc %in% c(0,1)) %>%
  # Remove outliers
  filter(estM1 > 0.05) %>%
  group_by(TruePerc) %>%
  summarise(AvgEstM1 = mean(estM1))

# In plot
full_data %>%
  filter(!TruePerc %in% c(0,1)) %>%
  # Remove outliers
  #filter(estM1 > 0.05) %>%
  # Get rid of the intermediate steps (target percentages)
  group_by(TruePerc, sim) %>%
  summarise(estM1 = mean(estM1)) %>%
  # Now plot
  ggplot(., aes(x = factor(TruePerc), y = estM1)) + 
  geom_violin(aes(fill = factor(TruePerc))) +
  scale_x_discrete('Percentage truly activated voxels') + 
  scale_y_continuous('Estimated proportion of truly activated voxels',
                     breaks = seq(0.1,0.9,by=0.1)) +
  scale_fill_brewer('', type = 'qual', palette = 3) +
  guides(fill=FALSE) + 
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 11),
        plot.subtitle = element_text(hjust = 0, vjust = -2),
        legend.key.size = unit(.85,'cm'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))
  
# Using the adaptive algorithm
full_data %>%
  filter(TruePerc != 0) %>% 
  group_by(TruePerc, sim) %>%
  filter(obs_overlap == max(obs_overlap, na.rm = TRUE)) %>%
  mutate(estM1 = target_perc) %>% 
  dplyr::select(-target_perc, -obs_overlap) %>%
  ungroup() %>%
  # Now plot
  ggplot(., aes(x = factor(TruePerc), y = estM1)) + 
  geom_boxplot(aes(fill = factor(TruePerc))) +
  scale_x_discrete('Percentage truly activated voxels') + 
  scale_y_continuous('Estimated proportion of truly activated voxels',
                     breaks = seq(0.1,0.9,by=0.1)) +
  scale_fill_brewer('', type = 'qual', palette = 3) +
  guides(fill=FALSE) + 
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length=unit(.20, "cm"),
        axis.line = element_line(size = .75),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 11),
        plot.subtitle = element_text(hjust = 0, vjust = -2),
        legend.key.size = unit(.85,'cm'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13))
# Bias of m1
full_data %>%
  filter(TruePerc != 0) %>% 
  group_by(TruePerc, sim) %>%
  filter(obs_overlap == max(obs_overlap, na.rm = TRUE)) %>%
  mutate(estM1 = target_perc) %>% 
  dplyr::select(-target_perc, -obs_overlap) %>%
  ungroup() %>%
  group_by(TruePerc) %>%
  summarise(AvgEstM1 = mean(estM1))
  

















































