####################
#### TITLE:     Plot within sample size stability with increasing sample sizes.
#### Contents: 	
#### 
#### Source Files: /FreddieFreeloader/Analyses/_PreProcessing/_Stability
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
# Number of clusters
# Average cluster size
# Variance of cluster sizes
# Out of the X runs, how many clusters are at least one voxel overlapping?
# Average fraction of the clusters that have overlapping voxels.


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
# Number of groups in each step
NGROUP <- 2


##
###############
### Read in data
###############
##

# Data frame with cluster sizes
ClustSize <- readRDS(paste(LocIntRes, 'ClustSize.rda', sep = ''))

# Data frame with number of unique and total clusters
numUniqClust <- readRDS(paste(LocIntRes, 'numUniqClust.rda', sep = ''))

# Data frame with proportion of overlapping voxels in the clusters
propOverVox <- readRDS(paste(LocIntRes, 'propOverVox.rda', sep = ''))


##
###############
### Process and plot various measurements
###############
##

#########################
##### CLUSTER SIZES #####
#########################

# Average cluster size within one study
# For each sample size, I calculate the average cluster size within one study
# We plot points for each simulation and group within the sample size.
# Then add the best fitting smoothed regression line
AvgClustS <- ClustSize %>%
  group_by(step, run, group) %>%
  # Calculate average within cluster size
  summarise(AWCS = mean(size)) %>%
  mutate(SampleSize = step * 10) %>%
  ungroup() %>%
  ggplot(., aes(x = SampleSize, y = AWCS)) + 
  geom_point(size = 0.75) +
  geom_smooth(method = 'loess') + 
  scale_x_continuous('Sample size') + 
  scale_y_continuous('Average cluster size within one fMRI study') +
  theme_bw()
AvgClustS


# Average cluster size over both replications of largest cluster: percentage of masked voxels
# Here we focus on the largest cluster over both groups.
# We measure its size and then calculate the proportion it 
# overlaps with in the masked brain.
PropLargGroupC <- ClustSize %>%
  # First select largest cluster per group (and step, run)
  group_by(step, run, group, NumMask) %>%
  top_n(n = 1, size) %>%
  ungroup() %>%
  # Now group by step and run
  group_by(step, run, NumMask) %>%
  # Calculate average over the two independent replications of the cluster size
  # Hence we have half of the dataset remaining
  summarise(AWCS = mean(size)) %>%
  # Add sample size
  mutate(SampleSize = step * 10) %>%
  # Percentage of masked voxels
  mutate(PercClust = AWCS/NumMask) %>%
  ungroup() %>%
  ggplot(., aes(x = SampleSize, y = PercClust)) + 
  geom_point(size = 0.75) +
  geom_smooth(method = 'loess') + 
  scale_x_continuous('Sample size') + 
  scale_y_continuous('Proportion of largest cluster over all masked voxels ') +
  ggtitle('Average over both replications') +
  theme_bw()
PropLargGroupC


# Average cluster size for each replications of largest cluster: percentage of masked voxels
# Here we focus on the largest cluster in each group.
# We measure its size and then calculate the proportion it 
# overlaps with in the masked brain.
PropLargC <- ClustSize %>%
  # First select largest cluster per group (and step, run)
  group_by(step, run, group, NumMask) %>%
  top_n(n = 1, size) %>%
  ungroup() %>%
  # Add sample size
  mutate(SampleSize = step * 10) %>%
  # Percentage of masked voxels
  mutate(PercClust = size/NumMask) %>%
  ungroup() %>%
  ggplot(., aes(x = SampleSize, y = PercClust, group = factor(group))) + 
  geom_point(aes(colour = factor(group)), size = 0.75) +
  geom_smooth(method = 'loess') + 
  scale_x_continuous('Sample size') + 
  scale_y_continuous('Proportion of all masked voxels in largest cluster') +
  scale_colour_manual('Replication', labels = c('1','2'), values = c('#045a8d',
                        '#74a9cf')) +
  theme_bw() + 
  # Increase size in legend
  guides(colour = guide_legend(override.aes = list(size=6))) +
  theme(legend.position = 'bottom')
PropLargC


# Variance (SD) in terms of amount of voxels of largest cluster
# Calculated over runs and groups
SDClustSize <- ClustSize %>%
  group_by(step, run, group) %>%
  # Filter only largest cluster
  top_n(n=1, size) %>%
  ungroup() %>%
  # now calculate variance
  group_by(step) %>%
  summarise(amount = max(index),
            avgSize = mean(size),
            varSize = var(size),
            sdSize = sd(size)) %>%
  mutate(SampleSize = step * 10) %>%
  ggplot(., aes(x = SampleSize, y = sdSize)) + 
  geom_line(size = 0.7) + 
  scale_x_continuous('Sample size') + 
  scale_y_continuous('Standard deviation of size of largest clusters') +
  theme_bw()
SDClustSize



##########################################
##### UNIQUE VS OVERLAPPING CLUSTERS #####
##########################################

# Here we calculate the number of overlapping clusters
# We look at each test and retest and compare those two images.
# If one voxel of the cluster is found in the other cluster, then they are
# overlapping!
# We average over all runs for each sample size to obtain the total amount 
# of clusters in both images (tests) and the number of overlapping clusters.

OverlClust1Vox <- numUniqClust %>% 
  # Take overlapping clusters instead of unique clusters
  mutate(OverlClust = TotClus - UniClus) %>%
  # Gather clusters in one column
  gather(key = 'cluster', value = 'count', 1,2,5) %>%
  mutate(SampleSize = step * 10) %>%
  filter(cluster != 'UniClus') %>%
  group_by(cluster) %>% 
  # Plot
  ggplot(., aes(x = SampleSize, y = count, colour = cluster)) +
  geom_smooth(aes(colour = cluster), size = 1.1) +
  scale_x_continuous('Sample Size') +
  scale_y_continuous('Count (clusters)') +
  scale_color_manual('', values = c('#1b9e77','#d95f02'),
                     labels = c('Overlapping clusters',
                                'Total amount of clusters')) + 
  theme_bw() +
  theme(legend.position = 'bottom') 
OverlClust1Vox


###############################################
##### PROPORTION OVERLAP BETWEEN CLUSTERS #####
###############################################


# Here we calculate the amount of voxels within each cluster that are also
# found in a cluster in the retest-image.
# So this is acutally the union of clustered voxels. 
# We then fit a smoothed regression over all runs and sample sizes
UnionCluster <- propOverVox %>%
  # Mutate proportion
  mutate(propOver = OverVox/TotVox,
         SampleSize = step * 10) %>%
  group_by(step, run) %>%
  ungroup() %>%
  ggplot(., aes(x = SampleSize, y = propOver)) + 
  geom_point(size = 0.75) +
  geom_smooth(method = 'loess') + 
  scale_x_continuous('Sample size') + 
  scale_y_continuous('Union of overlapping clustered voxels') +
  theme_bw()







