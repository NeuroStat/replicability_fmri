####################
#### TITLE:  Illustration of effect of threshold on overlap measure: plot figure
#### Contents: 	
#### 
#### Source Files: /replicability_fmri/Analyses/_Figures/
#### First Modified: 18/09/2018
#### Notes: 
#################


##
###############
### Notes
###############
##

# First we generated data in illustration_overlap_HPC.R.
# Here we gather the results and plot into figure. 


##
###############
### Preparation
###############
##

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

# Which sigma do you want?
sigma_opt <- c(1,2)
sigma <- sigma_opt[2]

# Location of results (relative to this file)
rel_loc_res <- paste0('Results/Sigma_', sigma)

# Number of independent batches on HPC
if(sigma == 1){
  nID <- 100
}
if(sigma == 2){
  nID <- 100
}

##
###############
### Gather results
###############
##

# Empty data frame
sim_data <- data.frame() %>% as_tibble()

# Get the data, which are simulated using the HPC in nID batches
for(i in 1:nID){
  sim_data <- bind_rows(sim_data,
    readRDS(paste(rel_loc_res,'/illus_overl_', i, '.rds', sep = ''))
    )
}

# Create labels for the baseline proportion and make target a factor
sim_data <- sim_data %>% 
  mutate(BasePropL = paste('Prop. true active voxels = ', BaseProp * 100, '%', sep='')) %>%
  # Order the facets, based on proportion
  mutate(BasePropLOr = reorder(BasePropL, BaseProp)) %>%
  mutate(targetF = factor(target)) %>%
  # Drop BasePropL (not needed in the figures) %>%
  dplyr::select(-BasePropL)
  

##
###############
### Plot figure
###############
##


# First plotting when fixing P-threshold.
# Variable that for facetting: baseProp
fixP <- sim_data %>% 
  filter(approach == 'normal') %>%
  # Drop unused variables
  dplyr::select(-target) %>%
  # Now group and summarise
  group_by(BasePropLOr, P_threshold, N) %>%
  summarise(avgOv = mean(overlap, na.rm = TRUE)) %>%
  ggplot(., aes(x = N, y = avgOv)) + 
  geom_line(aes(colour = P_threshold), size = .75) +
  facet_wrap(~BasePropLOr, as.table = TRUE, dir = 'h') + 
  scale_color_continuous('Threshold (P)') + 
  scale_y_continuous(expression(Overlap~(omega))) +
  labs(title = 'Average overlap of activation - numerical simulation') + 
       # subtitle = expression(X %~% N(mu,sigma^2))) +
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
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 11),
        plot.subtitle = element_text(hjust = 0, vjust = -2),
        legend.key.size = unit(.85,'cm'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.position = c(0.86, 0.22))
fixP


# Now plotting for adaptive thresholding.
# Here we have an extra variable: target.
# Use as extra lines.
adapP <- sim_data %>% 
  filter(approach == 'adaptive') %>% 
  # Remove unused variables
  dplyr::select(-P_threshold, -sigma, -target) %>%
  group_by(BaseProp, BasePropLOr, targetF, N) %>%
  summarise(avgOv = mean(overlap, na.rm = TRUE)) %>%
  ggplot(., aes(x = N, y = avgOv)) + 
  geom_line(aes(colour = targetF), size = 1) +
  facet_wrap(~BasePropLOr) + 
  scale_color_brewer('Algorithm target proportion of \n activated voxels', type = 'qual', palette = 7) +
  #scale_linetype_discrete('Target % \n activated voxels') + 
  scale_y_continuous(expression(Overlap~(omega))) +
  labs(title = 'Average overlap of activation - numerical simulation',
       subtitle = 'Adaptive thresholding') +
       #subtitle = expression(X %~% N(mu,sigma^2))) +
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
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 11),
        plot.subtitle = element_text(hjust = 0, vjust = -2),
        legend.key.size = unit(.75,'cm'),
        legend.key.width = unit(1,'cm'),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.position = c(0.86, 0.22))
adapP




##
###############
### Save figures
###############
##

# Fixed thresholds
ggsave(filename = paste(getwd(), '/Illus_overl_fixP.png', sep = ''),
       plot = fixP,
       width = 32, height = 22, units = 'cm', scale = .75)


# Adaptive thresholding
ggsave(filename = paste(getwd(), '/Illus_overl_adapP.png', sep = ''),
       plot = adapP,
       width = 32, height = 22, units = 'cm', scale = 1)

























