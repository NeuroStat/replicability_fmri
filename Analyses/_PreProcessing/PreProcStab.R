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
  # Out of all analyses, how many clusters are at least one voxel overlapping with the replication analysis?
  # --> Also take the average over test and replication of the latter.
  # Average fraction of the clusters that have overlapping voxels.



##
###############
### Preparation
###############
##

# Source paths
source(blind_PreProcessing.R)
RawDat <- RawDatStab

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
          'OverlCluster' = as.integer(),
          'TotClus' = as.integer(),
          'AvUniClus'= as.integer(),
          'AvOverlCluster'= as.integer(),
          'AvTotClus'= as.integer(),
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
    
    # The total amount of clusters
      # --> NOTE: this is the total over BOTH images (so not the total per analysis)
    TotClus <- as.integer(c(max(imageG1) + max(imageG2)))
    
    # The number of overlapping clusters --> count the cluster in both images
      # that overlaps (so you have 2 clusters overlapping, one in each analysis, 
      # per overlap event). 
    OverlCluster <- TotClus - uniCluster
    
    # NOTE: we will also add the average number of unique, overlapping and total clusters.
    #   This makes more sense as we then can speak about the individual analysis (on average).
    #   While the other numbers are about both analyses together. 
    AvUniClus <- uniCluster/2
    AvOverlCluster <- OverlCluster/2
    AvTotClus <- TotClus/2

    # Add to data frame
    numUniqClust <- bind_rows(numUniqClust,
          data.frame('UniClus' = uniCluster,
                     'OverlCluster' = OverlCluster,
                     'TotClus' = TotClus,
                     'AvUniClus' = AvUniClus,
                     'AvOverlCluster' = AvOverlCluster,
                     'AvTotClus' = AvTotClus,
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



