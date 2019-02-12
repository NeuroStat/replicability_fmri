####################
#### TITLE:     Plot within sample size stability with increasing sample sizes.
#### Contents: 	
#### 
#### Source Files: /replicability_fmri/Analyses/_Figures/
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

# Possible contrasts: default = MATH > LANGUAGE
contrast <- c('ML', 'Faces', 'Incentive', 'StopGo')
contr <- contrast[4]

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
# Number of groups in each step
NGROUP <- 2

# Variables for plotting
subjBreak <- c(10, seq(100,700, by=100))


##
###############
### Read in data
###############
##

# Data frame with cluster sizes
ClustSize <- readRDS(paste(LocIntRes, contr, '/', 'ClustSize.rda', sep = ''))

# Data frame with number of unique and total clusters
numUniqClust <- readRDS(paste(LocIntRes, contr, '/', 'numUniqClust.rda', sep = ''))

# Data frame with number of unique and total clusters when using a cut-off percentage
numPercClust <- readRDS(paste(LocIntRes, contr, '/', 'numPercClust.rda', sep = ''))

# Data frame with proportion of overlapping voxels in the clusters
propOverVox <- readRDS(paste(LocIntRes, contr, '/', 'propOverVox.rda', sep = ''))


##
###############
### Process and plot various measurements
###############
##

####################################
##### CLUSTER SIZES AND COUNTS #####
####################################

# Average cluster size for N = 700
ClustSize %>%
  mutate(SampleSize = step * 10) %>%
  dplyr::select(index, size, run, group, SampleSize) %>%
  group_by(SampleSize) %>%
  filter(SampleSize == 700) %>% 
  # Calculate average cluster size
  summarise(avSize = mean(size))

# Proportion of clusters with cluster size > 10e+04 at N = 100
ClustSize %>%
  mutate(SampleSize = step * 10) %>%
  dplyr::select(index, size, run, group, SampleSize) %>%
  # Calculate average within-study cluster size
  group_by(SampleSize, run, group) %>%
  summarise(avCS = mean(size)) %>%
  # Filter at N = 100
  filter(SampleSize == 100) %>% 
  ungroup() %>%
  mutate(PercRank = cume_dist(desc(avCS))) %>%
  # Sort
  arrange(avCS) %>%
  # Cut-off
  filter(avCS >= 10^4)

# Same calculation:
ClustSize %>%
  mutate(SampleSize = step * 10) %>%
  dplyr::select(index, size, run, group, SampleSize) %>%
  # Calculate average within-study cluster size
  group_by(SampleSize, run, group) %>%
  summarise(avCS = mean(size)) %>%
  # Filter at N = 100
  filter(SampleSize == 100) %>% 
  ungroup() %>%
  # Sort
  arrange(avCS) %>%
  mutate(TotalOb = n()) %>%
  # Cut-off
  filter(avCS >= 10^4) %>%
  # Remaining number of clusters
  mutate(RemOb = n()) %>%
  # Proportion
  mutate(ProbLarge = RemOb/TotalOb)


# Proportion of masked voxels in a cluster
# --> note index 2 is also the large cluster, but in two analyses there is 
#     a smaller cluster, then with index = 1.
ClustSize %>%
  # Create sample size
  mutate(SampleSize = step * 10) %>%
  # Select largest sample size
  filter(SampleSize == 700) %>%
  dplyr::select(index, size, group, NumMask) %>%
  # Average per index and group
  group_by(index, group) %>%
  summarise(AvInSizeG = mean(size),
            AvMask = mean(NumMask)) %>%
  # Average over both groups
  group_by(index) %>%
  summarise(AvInSize = mean(AvInSizeG),
            AvMask = mean(AvMask)) %>%
  # Sum both indices
  ungroup() %>%
  summarise(SumInd = sum(AvInSize),
            AvMask = mean(AvMask)) %>%
  # Calculate proportion
  mutate(TotProp = SumInd/AvMask)

# Proportion of masked voxels in a cluster
# --> but now averaging over the index
ClustSize %>%
  # Create sample size
  mutate(SampleSize = step * 10) %>%
  # Select largest sample size
  filter(SampleSize == 700) %>%
  dplyr::select(index, size, group, NumMask) %>%
  # Average over both groups and index
  summarise(AvSize = mean(size),
            AvMask = mean(NumMask)) %>%
  # Calculate proportion
  mutate(TotProp = AvSize/AvMask)

# Proportion of masked voxels in a cluster, separated by the two groups at N = 700
ClustSize %>%
  # Create sample size
  mutate(SampleSize = step * 10) %>%
  # Select largest sample size
  filter(SampleSize == 700) %>%
  # Only take the largest cluster (there are two analyses with an extra small cluster)
  group_by(step, run, group, SampleSize) %>% 
  top_n(1, size) %>% 
  ungroup() %>%
  dplyr::select(size, group, run, NumMask) %>%
  # Average per group
  group_by(group) %>%
  summarise(AvInSize = mean(size),
            AvMask = mean(NumMask)) %>%
  # Calculate proportion
  mutate(TotProp = AvInSize/AvMask)

# Now let's check the convergence into the amount of clusters.
ClustSize %>%
  group_by(step, run, group) %>%
  # Select number of clusters per analyses (i.e. highest index)
  top_n(n=1, index) %>%
  # Calculate average over runs
  ungroup() %>%
  group_by(step) %>%
  summarise(AvgCount = mean(index)) %>%
  mutate(SampleSize = step * 10) %>%
  filter(SampleSize == 700)

# Average cluster size within one study
# For each sample size, I calculate the average cluster size within one study
# We plot points for each simulation and group within the sample size.
# Then add the best fitting smoothed regression line
AvgClustS <- ClustSize %>%
  group_by(step, run, group) %>%
  # Calculate average within-study cluster size
  summarise(AWCS = mean(size)) %>%
  mutate(SampleSize = step * 10) %>%
  ungroup() %>%
  ggplot(., aes(x = SampleSize, y = AWCS)) + 
  geom_point(size = 0.75) +
  geom_smooth(method = 'loess', colour = '#8856a7') + 
  scale_x_continuous(breaks = subjBreak, name="Sample size") + 
  scale_y_continuous('Number of voxels', labels = scales::scientific) +
  labs(title = 'Average cluster size in an fMRI group analysis',
       subtitle = 'Z = 2.3 and FWER = 0.05') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length = unit(.20, "cm"),
        axis.line = element_line(size = .75),
        panel.grid = element_line(linetype = 'dotted'),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
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
PropLargC <- 
  ClustSize %>%
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
  geom_smooth(method = 'loess', colour = '#8856a7') + 
  scale_x_continuous(breaks = subjBreak, name="Sample size") + 
  scale_y_continuous('Proportion') +
  scale_colour_manual('Replication', labels = c('1','2'), values = c('#045a8d',
                        '#74a9cf')) +
  # Title
  labs(title = 'Proportion of all masked voxels within largest cluster',
       subtitle = 'Z = 2.3 and FWER = 0.05') +
  # Increase size in legend
  guides(colour = guide_legend(override.aes = list(size=4))) +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain', size = 9),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length = unit(.20, "cm"),
        axis.line = element_line(size = .75),
        panel.grid = element_line(linetype = 'dotted'),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
PropLargC

# Variability of number of selected clusters
# We do not care about the two groups in this case, that is we calculate SD
# over both groups. Hence we have NRUNS * 2 amount of analyses per sample size.
ClustSize %>%
  group_by(step, run, group) %>%
  # Select number of clusters per analyses (i.e. highest index)
  top_n(n=1, index) %>%
  ungroup() %>%
  mutate(SampleSize = step * 10) %>%
  ggplot(., aes(x = factor(SampleSize), y = index)) +
  geom_boxplot(outlier.size = .7) +
  scale_x_discrete(breaks = subjBreak, name="Sample size") +
  scale_y_continuous(name='Count') +
  facet_grid( ~ group) +
  theme_bw()

# Draw line for the variability
SDClustCount <- 
  ClustSize %>%
  group_by(step, run, group) %>%
  # Select number of clusters per analyses (i.e. highest index)
  top_n(n=1, index) %>%
  ungroup() %>%
  # Calculate variability over runs and group
  group_by(step) %>%
  summarise(SDcount = sd(index)) %>%
  mutate(SampleSize = step * 10) %>%
  # Draw line
  ggplot(., aes(x = SampleSize, y = SDcount)) +
  geom_line(size = 0.9) + 
  scale_x_continuous(breaks = subjBreak, 'Sample size') + 
  scale_y_continuous('Standard deviation on number of clusters') +
  labs(title = 'Variability of number of clusters',
       subtitle = 'Z = 2.3 and FWER = 0.05') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length = unit(.20, "cm"),
        axis.line = element_line(size = .75),
        panel.grid = element_line(linetype = 'dotted'),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
SDClustCount

# Variance (SD) in terms of amount of voxels of largest cluster
# Calculated over runs and groups
SDClustSize <- 
  ClustSize %>%
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
  geom_line(size = 0.9) + 
  scale_x_continuous(breaks = subjBreak, 'Sample size') +
  scale_y_continuous('Standard deviation on number of voxels', 
                     labels = scales::scientific) +
  labs(title = 'Variability of size of largest cluster',
       subtitle = 'Z = 2.3 and FWER = 0.05') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length = unit(.20, "cm"),
        axis.line = element_line(size = .75),
        panel.grid = element_line(linetype = 'dotted'),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
SDClustSize


##########################################
##### UNIQUE VS OVERLAPPING CLUSTERS #####
##########################################


# Here we calculate the number of overlapping clusters
# We look at each test and retest and compare those two images.
# If one voxel of the cluster is found in the other cluster, then they are
# overlapping!
# We average over all runs for each sample size and then average over the test and replication
# to obtain the average number of overlapping vs non-overlapping (and total number of) clusters.
OverlClust1Vox <- numUniqClust %>% 
  # Gather clusters in one column
  gather(key = 'cluster', value = 'count', 1,2,3,4,5,6) %>%
  mutate(SampleSize = step * 10) %>%
  # Filter out the average unique, overlapping and total number of clusters
  filter(cluster %in% c('AvUniClus', 'AvOverlCluster', 'AvTotClus')) %>%
  group_by(cluster) %>% 
  # Plot
  ggplot(., aes(x = SampleSize, y = count, colour = 
                # Recode cluster
                factor(cluster, levels = 
                    c('AvTotClus', 'AvOverlCluster', 'AvUniClus'), 
                labels = 
                    c('AvTotClus', 'AvOverlCluster', 'AvUniClus')))) +
  geom_smooth(aes(colour = 
                # Recode cluster
                factor(cluster, levels = 
                         c('AvTotClus', 'AvOverlCluster', 'AvUniClus'), 
                       labels = 
                         c('AvTotClus', 'AvOverlCluster', 'AvUniClus'))),
              size = 1.7,
              method = 'gam',
              formula = y ~ s(x, bs = "cs")) +
  scale_x_continuous(breaks = subjBreak, 'Sample size') +
  scale_y_continuous('Average number of clusters') +
  scale_color_manual('', values = c('#1b9e77','#984ea3', '#d95f02'),
                     labels = c('Total count',
                                'Overlapping',
                                'Non-overlapping')) + 
  labs(title = 'Total, overlapping & non-overlapping clusters',
       subtitle = 'Z = 2.3 and FWER = 0.05') +
  guides(colour = guide_legend(order = 3)) +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length = unit(.20, "cm"),
        axis.line = element_line(size = .75),
        panel.grid = element_line(linetype = 'dotted'),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = 'bottom')
OverlClust1Vox


# When is average non-overlapping more or less zero?
numUniqClust %>% 
  # Gather clusters in one column
  gather(key = 'cluster', value = 'count', 1,2,3,4,5,6) %>%
  mutate(SampleSize = step * 10) %>%
  # Filter out the average unique, overlapping and total number of clusters
  filter(cluster == 'AvUniClus') %>%
  # Summarise over sample sizes
  group_by(SampleSize) %>%
  summarise(AvUniClus = mean(count)) %>% 
  # Sample size when average < 0.05
  filter(AvUniClus <= 0.05)



#######################################################
##### UNIQUE VS OVERLAPPING CLUSTERS WITH CUT-OFF #####
#######################################################

# Here we calculate the number of overlapping clusters
# We look at each test and retest and compare those two images.
# If 50% of the voxels of a cluster are found in the other cluster, then they are
# considered replicated!
# We average over all runs for each sample size and then average over the test and replication
# to obtain the average number of overlapping vs non-overlapping (and total number of) clusters.
OverlClustCutOff <- numPercClust %>% 
  # Gather clusters in one column
  gather(key = 'cluster', value = 'count', 1,2,3,4) %>%
  mutate(SampleSize = step * 10) %>%
  # Filter out the average unique, overlapping and total number of clusters
  filter(cluster %in% c('UniPercClust', 'OverlPercClust', 'AvTotClus')) %>%
  group_by(cluster)  %>% 
  # Plot
  ggplot(., aes(x = SampleSize, y = count, colour = 
                  # Recode cluster
                  factor(cluster, levels = 
                           c('AvTotClus', 'OverlPercClust', 'UniPercClust'), 
                         labels = 
                           c('AvTotClus', 'OverlPercClust', 'UniPercClust')))) +
  geom_smooth(aes(colour = 
                    # Recode cluster
                    factor(cluster, levels = 
                             c('AvTotClus', 'OverlPercClust', 'UniPercClust'), 
                           labels = 
                             c('AvTotClus', 'OverlPercClust', 'UniPercClust'))),
              size = 1.7,
              method = 'gam',
              formula = y ~ s(x, bs = "cs")) +
  scale_x_continuous(breaks = subjBreak, 'Sample size') +
  scale_y_continuous('Average number of clusters') +
  scale_color_manual('', values = c('#1b9e77','#984ea3', '#d95f02'),
                     labels = c('Total count',
                                'Overlapping',
                                'Non-overlapping')) + 
  labs(title = 'Total, overlapping & non-overlapping clusters',
       subtitle = 'Z = 2.3 and FWER = 0.05') +
  guides(colour = guide_legend(order = 3)) +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length = unit(.20, "cm"),
        axis.line = element_line(size = .75),
        panel.grid = element_line(linetype = 'dotted'),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = 'bottom')
OverlClustCutOff

# When is average non-overlapping more or less zero?
numPercClust %>% 
  # Gather clusters in one column
  gather(key = 'cluster', value = 'count', 1,2,3,4) %>%
  mutate(SampleSize = step * 10) %>%
  # Filter out the average unique, overlapping and total number of clusters
  filter(cluster == 'UniPercClust') %>%
  # Summarise over sample sizes
  group_by(SampleSize) %>%
  summarise(AvUniClus = mean(count)) %>% 
  # Sample size when average < 0.05
  filter(AvUniClus <= 0.05)

###################################################################
##### UNIQUE VS OVERLAPPING CLUSTERS WITH AND WITHOUT CUT-OFF #####
###################################################################


OverlClust1VoxCutOff <- numUniqClust %>% 
  # Gather clusters in one column
  gather(key = 'cluster', value = 'count', 1,2,3,4,5,6) %>%
  mutate(SampleSize = step * 10) %>%
  # Filter out the average unique, overlapping and total number of clusters
  filter(cluster %in% c('AvUniClus', 'AvOverlCluster', 'AvTotClus')) %>%
  mutate(CutOff = 'A: at least one voxel overlapping') %>% 
  bind_rows(.,
      numPercClust %>% 
        # Gather clusters in one column
        gather(key = 'cluster', value = 'count', 1,2,3,4) %>%
        mutate(SampleSize = step * 10) %>%
        # Filter out the average unique, overlapping and total number of clusters
        filter(cluster %in% c('UniPercClust', 'OverlPercClust', 'AvTotClus')) %>%
        # recode cluster
        mutate(cluster = recode(cluster,  
                'AvTotClus' = 'AvTotClus',
               'OverlPercClust' = 'AvOverlCluster', 
               'UniPercClust' = 'AvUniClus')) %>%
        mutate(CutOff = 'B: 50% of the cluster overlapping')) %>%
  group_by(cluster) %>% 
  # Plot
  ggplot(., aes(x = SampleSize, y = count, colour = 
                  # Recode cluster
                  factor(cluster, levels = 
                           c('AvTotClus', 'AvOverlCluster', 'AvUniClus'), 
                         labels = 
                           c('AvTotClus', 'AvOverlCluster', 'AvUniClus')))) +
  geom_smooth(aes(colour = 
                    # Recode cluster
                    factor(cluster, levels = 
                             c('AvTotClus', 'AvOverlCluster', 'AvUniClus'), 
                           labels = 
                             c('AvTotClus', 'AvOverlCluster', 'AvUniClus'))),
              size = 1.7,
              method = 'gam',
              formula = y ~ s(x, bs = "cs")) +
  facet_wrap(~CutOff) +
  scale_x_continuous(breaks = subjBreak, 'Sample size') +
  scale_y_continuous('Average number of clusters') +
  scale_color_manual('', values = c('#1b9e77','#984ea3', '#d95f02'),
                     labels = c('Total count',
                                'Overlapping',
                                'Non-overlapping')) + 
  labs(title = 'Total, overlapping & non-overlapping clusters',
       subtitle = 'Z = 2.3 and FWER = 0.05') +
  guides(colour = guide_legend(order = 3)) +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length = unit(.20, "cm"),
        axis.line = element_line(size = .75),
        panel.grid = element_line(linetype = 'dotted'),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = 'bottom')
OverlClust1VoxCutOff


  
###############################################
##### PROPORTION OVERLAP BETWEEN CLUSTERS #####
###############################################


# Here we calculate the amount of voxels within each cluster that are also
# found in a cluster in the retest-image.
# So this is acutally the intersection of clustered voxels. 
# We then fit a smoothed regression over all runs and sample sizes
IntsCluster <- propOverVox %>%
  # Mutate proportion
  mutate(propOver = OverVox/TotVox,
         SampleSize = step * 10) %>%
  group_by(step, run) %>%
  ungroup() %>%
  ggplot(., aes(x = SampleSize, y = propOver)) + 
  geom_point(size = 0.75) +
  geom_smooth(method = 'loess', colour = '#8856a7') + 
  scale_x_continuous(breaks = subjBreak, 'Sample size') + 
  scale_y_continuous('Proportion') +
  labs(title = 'Proportion of overlapping voxels in a cluster',
     subtitle = 'Z = 2.3 and FWER = 0.05') +
  theme_classic() +
  theme(panel.grid.major = element_line(size = 0.8),
        panel.grid.minor = element_line(size = 0.8),
        axis.title.x = element_text(face = 'plain'),
        axis.title.y = element_text(face = 'plain'),
        axis.text = element_text(size = 11, face = 'plain'),
        axis.ticks = element_line(size = 1.3),
        axis.ticks.length = unit(.20, "cm"),
        axis.line = element_line(size = .75),
        panel.grid = element_line(linetype = 'dotted'),
        title = element_text(face = 'plain'),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'bottom')
IntsCluster

# Median proportion at N = 200
propOverVox %>%
  # Mutate proportion
  mutate(propOver = OverVox/TotVox,
         SampleSize = step * 10) %>%
  # Summarise
  group_by(SampleSize) %>%
  summarise(MedProp = median(propOver, na.rm = TRUE)) %>%
  # Sample size = 200
  filter(SampleSize >= 200)

# Maximum median proportion
propOverVox %>%
  # Mutate proportion
  mutate(propOver = OverVox/TotVox,
         SampleSize = step * 10) %>%
  # Summarise
  group_by(SampleSize) %>%
  summarise(MedProp = median(propOver, na.rm = TRUE)) %>%
  # Sample size
  filter(SampleSize == max(SampleSize))



##
###############
### Measurements at N = 30
###############
##

# Let us limit to overlapping vs non-overlapping (and total clusters) and the intersection
numUniqClust %>% 
  # Gather clusters in one column
  gather(key = 'cluster', value = 'count', 1,2,3,4,5,6) %>%
  mutate(SampleSize = step * 10) %>%
  # Filter out the average unique, overlapping and total number of clusters
  filter(cluster %in% c('AvUniClus', 'AvOverlCluster', 'AvTotClus')) %>%
  # Summarise over sample sizes
  group_by(SampleSize, cluster) %>%
  summarise(MedClust = median(count)) %>% 
  # N = 30
  filter(SampleSize == 30)

# When using conservative definition
numPercClust %>% 
  # Gather clusters in one column
  gather(key = 'cluster', value = 'count', 1,2,3,4) %>%
  mutate(SampleSize = step * 10) %>%
  # Filter out the average unique, overlapping and total number of clusters
  filter(cluster %in% c('UniPercClust', 'OverlPercClust', 'AvTotClus')) %>%
  group_by(SampleSize, cluster) %>%
  summarise(MedClust = median(count)) %>% 
  # N = 30
  filter(SampleSize == 30)

# Proportion intersecting
propOverVox %>%
  # Mutate proportion
  mutate(propOver = OverVox/TotVox,
         SampleSize = step * 10) %>%
  # Summarise
  group_by(SampleSize) %>%
  summarise(MedProp = median(propOver, na.rm = TRUE)) %>%
  # N = 30
  filter(SampleSize == 30)


##
###############
### Combine plots
###############
##


# Combine AvgClustS and PropLargC
plot_grid(AvgClustS, PropLargC, labels = c("A", "B"), nrow = 1, align = "h",
          axis = 'b')
ggsave(filename = paste0(getwd(), '/clusterSizes.png'),
       plot = ggplot2::last_plot(),
       width = 20, height = 14, units = 'cm', scale = 0.9)


# Combine AvgClustS, PropLargC, OverlClust1Vox and IntsCluster
plot_grid(AvgClustS, PropLargC, OverlClust1Vox, IntsCluster,
          labels = c("A", "B", "C", "D"), nrow = 2, align = "hv",
          axis = 'tblr')
ggsave(filename = paste0(getwd(), '/clusterStab.png'),
       plot = ggplot2::last_plot(),
       width = 20, height = 24, units = 'cm', scale = 1)

# Save the SD of amount of voxels in largest cluster
ggsave(filename = paste0(getwd(), '/clusterSD.png'),
       plot = SDClustSize, scale = 0.3)
      # width = 12, height = 10, units = 'cm', scale = 0.7)

# Stability: SD of cluster count and cluster size in largest cluster
plot_grid(SDClustCount, SDClustSize, nrow = 1, align = 'hv', axis = 'tblr')
ggsave(filename = paste0(getwd(), '/stabilitySD.png'),
       plot = ggplot2::last_plot(), 
       width = 12, height = 6, units = 'cm')

# Combine AvgClustS, PropLargC, OverlClust1Vox, IntsCluster, SDClustCount and SDClustSize
# --> function unstable, sometimes need to re-run if crashes!
plot_grid(AvgClustS, PropLargC, 
          SDClustCount, SDClustSize, 
          OverlClust1Vox, IntsCluster, 
          labels = c("A", "B", "C", "D", "E", "F"), nrow = 3, align = "hv",
          axis = 'tblr')
ggsave(filename = paste0(getwd(), '/FullStability.png'),
       plot = ggplot2::last_plot(),
       width = 26, height = 32, units = 'cm', scale = 1)


# Combine AvgClustS, PropLargC, OverlClust1VoxCutOff, IntsCluster, SDClustCount and SDClustSize
# --> function unstable, sometimes need to re-run if crashes!
plot_grid(AvgClustS, PropLargC, 
          SDClustCount, SDClustSize, 
          OverlClust1VoxCutOff, IntsCluster, 
          labels = c("A", "B", "C", "D", "E", "F"), nrow = 3, align = "hv",
          axis = 'tblr')
ggsave(filename = paste0(getwd(), '/FullStabilityCutOff_',contr,'.png'),
       plot = ggplot2::last_plot(),
       width = 26, height = 32, units = 'cm', scale = 1)



















