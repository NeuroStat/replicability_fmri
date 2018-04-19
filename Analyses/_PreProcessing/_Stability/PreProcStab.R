####################
#### TITLE:     Calculate within sample size stability with increasing sample sizes.
#### Contents: 	
#### 
#### Source Files: /FreddieFreeloader/Analyses/_PreProcessing/_Stability
#### First Modified: 18/04/2018
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

# Here we read in the raw data and then calculate:
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

# Location of raw data: not included in Github folder (too large)
RawDat <- "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Data/Stability"

# Save intermediate results
SaveLoc <- '/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Script.git/FreddieFreeloader/Analyses/_IntData'

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


########################
### CUSTOM FUNCTIONS ###
########################

# Function to read in clusters:
  # ARGUMENTS: location of data file, ID of run, step and group
ReadClusters <- function(location, run, step, group){
  # Read in the cluster_zstat1_std.txt file
  clust <- read.delim(file = paste(location,'/Run_',run,'/Step_',step,
                      '/Group', group, '/cluster_zstat1_std.txt', sep = ''),
                       header = TRUE)
  # If no clusters are found, then add row with elements 0
  if(dim(clust)[1] == 0){
    clust[1,] <- rep(0, length(clust))
  }
  names(clust) <- NULL
  # Add to data frame
  ClustSize <- data.frame('index' = clust[1],
                'size' = clust[2],
                'step' = as.integer(step),
                'run' = as.integer(run),
                'group' = as.integer(group))
  return(ClustSize)
}


##
###############
### Read in data
###############
##

# NOTE:
# For loop: start with sample size S, read in cluster files, 
#   then move to sample size S+1

# Data frame with the cluster size
ClustSize <- data.frame('index' = as.integer(),
          'size' = as.numeric(),
          'step' = as.integer(),
          'run' = as.integer(),
          'group' = as.integer()) %>% as.tibble()

# Data frame with the number of unique clusters
numUniqClust <- data.frame('UniClus' = as.integer(),
          'TotClus' = as.integer(),
          'step' = as.integer(),
          'run' = as.integer()) %>% as.tibble()

# Data frame with proportion of overlapping voxels in each image
propOverVox <- data.frame('OverVox' = as.integer(),
          'TotVox' = as.integer(),
          'step' = as.integer(),
          'run' = as.integer()) %>% as.tibble()

# Print statement
PriST <- (c(1:NRUNS)/NRUNS)[seq(1,NRUNS,length.out=10)][-10]

# For loop over all runs
for(i in 1:NRUNS){
  IDprint <- c(i/NRUNS)==PriST
  if(any(IDprint)){
    print(paste(round(PriST,2)[IDprint]*100, "% done"))
  }
  # For loop over all steps
  for(j in 1:NSTEP){
    #########################
    ##### CLUSTER SIZES #####
    #########################
    # We have data files from group 1 and group 2: read in cluster index + size
    # then bind to ClustSize data frame
    ClustSize %<>% bind_rows(
      ReadClusters(location = RawDat, run = i, step = j, group = 1),
      ReadClusters(location = RawDat, run = i, step = j, group = 2))

    #########################
    ##### NIFTI IMAGES ######
    #########################
    # Next we read in the nifti file with the clusters
    imageG1 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group1','/clusteroutput.nii.gz',
                                   sep=''))[,,], silent=TRUE)
    imageG2 <- try(readNIfTI(paste(RawDat,'/Run_',i,'/Step_',j,'/Group2','/clusteroutput.nii.gz',
                                   sep=''))[,,], silent=TRUE)

    ###########################
    ##### UNIQUE CLUSTERS #####
    ###########################    
    # First we select the image with the most clusters
    IDmax <- which.max(c(max(imageG1), max(imageG2)))
    # The other image
    IDmin <- which.min(c(max(imageG1), max(imageG2)))
    # If both are equal, set manually
    if(IDmax == IDmin){
      IDmax <- 1
      IDmin <- 2
    }
    # Switch to 2D matrix
    clusters <- matrix(c(array(imageG1, dim = prod(DIM)),
                       array(imageG2, dim = prod(DIM))),
                       ncol = 2)

    # Vector with common clusters in other image
    comClustMin <- comClustMax <- c(0)
    
    # Loop over the clusters
    for(r in 1:max(clusters[,IDmax])){
      # Select the voxels of r^th cluster
      IDr <- clusters[,IDmax] == r
      # Take matching clusters in other image
      matchC <- clusters[IDr,IDmin]
      # Cluster IDs that are found in other image
      comClustMin <- c(comClustMin, unique(matchC))
      # Overlapping clusters (if they exist) from this image that are found in other image
      if(any(matchC > 0)) comClustMax <- c(comClustMax, r)
    }
    
    # Remove 0 from both vectors
    comClustMin <- comClustMin[comClustMin != 0]
    comClustMax <- comClustMax[comClustMax != 0]
    # Remaining clusters: the ones that are not overlapping
      # sum over both images
    uniCluster <- sum(!1:max(clusters[,IDmin]) %in% comClustMin) + 
                    sum(!1:max(clusters[,IDmax]) %in% comClustMax)
    # Add to data frame
    numUniqClust <- bind_rows(numUniqClust,
          data.frame('UniClus' = uniCluster,
                     'TotClus' = as.integer(c(max(imageG1) + max(imageG2))),
                     'step' = as.integer(j),
                     'run' = as.integer(i)))
     
    ######################################################
    ###### PROPORTION WITHIN CLUSTERS OVERLAPPING ########
    ######################################################
    
    # First create binary images
    BinImageG1 <- imageG1
    BinImageG2 <- imageG2
    # Number of voxels in clusters
    BinImageG1[BinImageG1 > 0] <- 1
    BinImageG2[BinImageG2 > 0] <- 1
    # Sum both images and take numerator and denumerator of union
    numUnion <- sum((BinImageG1 + BinImageG2) == 2)
    denumUnion <- sum((BinImageG1 + BinImageG2) == 1)
    
    # Add to data frame
    propOverVox <- bind_rows(propOverVox,
      data.frame('OverVox' = numUnion,
                 'TotVox' = numUnion + denumUnion,
                 'step' = as.integer(j),
                 'run' = as.integer(i)))
  }
}


# I want to know how many masked voxels there are.
# Therefore, I read in one masked nifti image and sum the amount of 1's
NumMask <- sum(readNIfTI(paste(RawDat, '/Run_1/Step_1/Group1/masked.nii.gz', sep = ''))[,,])

# Add this to ClustSize dataframe
ClustSize %<>% mutate(NumMask = NumMask)


##
###############
### Save intermediate results
###############
##


saveRDS(ClustSize, paste(SaveLoc,'/ClustSize.rda', sep = ''))
saveRDS(numUniqClust, paste(SaveLoc,'/numUniqClust.rda', sep = ''))
saveRDS(propOverVox, paste(SaveLoc,'/propOverVox.rda', sep = ''))


# Analysis
ClustSize %>%
  group_by(step) %>%
  summarise(amount = max(index),
            avgSize = mean(size),
            varSize = var(size),
            sdSize = sd(size)) %>%
  mutate(SampleSize = step * 10)

# Variance (SD) of largest cluster
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
  geom_line() + 
  scale_x_continuous('Sample size') + 
  scale_y_continuous('Standard deviation of largest clusters') +
  theme_bw()


# Variance of largest clusters, but first with group mean centering
ClustSize %>%
  group_by(step, run, group) %>%
  # Filter only largest cluster
  top_n(n=1, size) %>%
  ungroup() %>%
  group_by(step) %>%
  # Average cluster size per sample size (over runs and groups)
  summarise(avgSize = mean(size)) %>% 
  # bind to data again
  right_join(., ClustSize, by = 'step') %>% 
  # again, filter on largest clusters
  group_by(step, run, group) %>%
  top_n(n=1, size) %>%
  # Now do group mean centering
  mutate(GrCent = size - avgSize) %>% 
  ungroup() %>%
  # group by step
  group_by(step, NumMask) %>% 
  # calculate variance of cluster size
  summarise(varSize = var(GrCent),
            sdSize = sd(GrCent)) %>%
  # Add sample size info
  mutate(SampleSize = step * 10) %>%
  ungroup() %>%
  ggplot(., aes(x = SampleSize, y = sdSize)) + 
  geom_line()



# Average cluster size within one study
ClustSize %>%
  group_by(step, run) %>%
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

# Average cluster size within one study: percentage of masked voxels
ClustSize %>%
  group_by(step, run) %>%
  # Calculate average within cluster size
  summarise(AWCS = mean(size)) %>%
  mutate(SampleSize = step * 10) %>%
  mutate(PercClust = AWCS/NumMask) %>%
  ungroup() %>%
  ggplot(., aes(x = SampleSize, y = PercClust)) + 
  geom_point(size = 0.75) +
  geom_smooth(method = 'loess') + 
  scale_x_continuous('Sample size') + 
  scale_y_continuous('Average percentage of masked voxels in cluster within one fMRI study') +
  theme_bw()

# Average cluster size of largest cluster: percentage of masked voxels
ClustSize %>%
  # First select largest cluster per group (and step, run)
  group_by(step, run, group) %>%
  summarise(MSize = max(size)) %>% 
  ungroup() %>%
  group_by(step, run) %>%
  # Calculate average within cluster size
  summarise(AWCS = mean(MSize)) %>%
  mutate(SampleSize = step * 10) %>%
  mutate(PercClust = AWCS/NumMask) %>%
  ungroup() %>%
  ggplot(., aes(x = SampleSize, y = PercClust)) + 
  geom_point(size = 0.75) +
  geom_smooth(method = 'loess') + 
  scale_x_continuous('Sample size') + 
  scale_y_continuous('Average percentage of masked voxels in largest cluster within one fMRI study') +
  theme_bw()

# Average number of unique clusters
numUniqClust %>% 
  # Gather nunique and total clusters in one column
  gather(key = 'cluster', value = 'count', 1:2) %>%
  mutate(SampleSize = step * 10) %>%
  group_by(cluster) %>% 
  # Plot
  ggplot(., aes(x = SampleSize, y = count, colour = cluster)) +
  geom_smooth(aes(colour = cluster))


# Average number of overlapping clusters
numUniqClust %>% 
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



# Average proportion of overlapping voxels within the clusters
propOverVox %>%
  # Mutate proportion
  mutate(propOver = OverVox/TotVox,
         SampleSize = step * 10) %>%
  group_by(step, run) %>%
  ungroup() %>%
  ggplot(., aes(x = SampleSize, y = propOver)) + 
  geom_point(size = 0.75) +
  geom_smooth(method = 'loess') + 
  scale_x_continuous('Sample size') + 
  scale_y_continuous('Average proportion of overlapping voxels') +
  theme_bw()


