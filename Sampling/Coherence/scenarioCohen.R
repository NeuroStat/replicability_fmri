####################
#### TITLE:     Sample in steps of 10 subjects from the total pool of all 
####            subjects, without controlling for external factors and create
####            independent disjoint sets.
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Script.git/Sampling/Coherence
#### First Modified: 23/05/2018
#### Notes: 
#################

##
###############
### Notes
###############
##

# SCENARIO: 
# We need to divide the total pool of subjects PER sample size (in steps of 10)
# in independent subsets (disjoints) to calculate the measure of coherence.
# The sampling is complete at random (sole restriction is that the subjects are different).
# The same is done for sample size 20, 30,..., 460.

# However, this might take a very long time to do sequentially. Hence, we can use the HPC to work on in parallel. 
# We already sample the subjects.
# Then each analysis can be executed independent of the previous one. 
# Later on we can calculate the coherence.


##
###############
### Preparation
###############
##

# Reset wm
rm(list=ls())

# Provide working directory
wd <- "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/SampleSizes/Coherence/"


# Set the starting seed
StartingSeed <- 12
set.seed(StartingSeed)

# Load the file with information about dataset
  # This file contains subject number without and with leading zero (1,2), 
  # then center ID (3) and the number I have of the subject in my folder.
load('/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/IMAGEN/IMAGEN/IMAGEN/OurIDs')
# Some columns are obsolote (refer to research questions earlier)
OurIDs <- OurIDs[,-c(5:11)]


# Libraries
library(RColorBrewer)
library(ggplot2)
library(grid)


# Global variables
# As information, the next subjects are deleted from the sampling process:
# We do not remove sujbects here anymore, this is already done in the OurIDs file!
# I just let them here as a reference!!!!
remove <- c(1284,680,814,1460,1470,788,																															# Corrupted files, visual inspection
40, 54, 61, 69, 155, 349, 384, 591, 869, 1219, 																											# Less amount of total masked voxels
43,77,89,169,180,235,375,394,471,489,544,637,744,787,903,957,1043,1065,1129,1168,1232,1324,1335,1351,																					# Less amount of total masked voxels
10, 166, 245, 373, 428, 474, 482, 491, 536, 565, 588, 600, 613, 618, 662, 668, 685, 699, 722, 747, 799, 829, 890, 987, 997, 1013, 1078, 1097,	# Less amount of total masked voxels
1140, 1212, 1217, 1228, 1302, 1310, 1368, 1403, 1406, 1462, 1478,																		# No activity in upper brain regions
68, 250, 328, 475, 551, 757, 1085, 1420)	


# Starting amount of subjects in group analysis
startingTotal <- 10

# Number of resampling runs
NRUNS <- 50

# Total amount of possible subjects if we have a maximum of 3 disjoint subsets
NTOT <- 460


##
###############
### Sampling
###############
##


# We add sample sizes starting from 10 to 460.
steps <- seq(startingTotal, NTOT, by = 10)

# Data frame with all subjects
AllSub <- c()

# Vector with group number
Group <- c()

# Vector with sample size
SSize <- c()

# Vector with run
run <- c()

# Iterating the sampling procedure
for(k in 1:NRUNS){
	print(paste('In run ', k, sep=''))
	seed <- StartingSeed*k
	set.seed(seed)
	# Start sequence of increasing sample sizes (steps)
	for(i in 1:length(steps)){
	  # Number of disjoint sets in this sample size
	  numDS <- floor(dim(OurIDs)[1]/steps[i])

	  # Starting set of subjects in this sample size (step)
	  ToSample <- OurIDs[,'OwnID']
	  # For loop over number of disjoint sets
	  for(s in 1:numDS){
      # Sample subjects
	    SampleSet <- sample(ToSample, steps[i], replace=FALSE)
	    # Print those subjects to the file
	    cat(SampleSet, sep='\n', file=paste(wd, 'Run_', k, '/SubjectsStep', i, '_', s, '.txt',sep=''))
	    
	    # Remove those subjects from total pool
	    IDSet <- ToSample %in% SampleSet
	    ToSample <- ToSample[!IDSet]
	    
      # Add the selected subjects to a data frame
	    AllSub <- c(AllSub,SampleSet)
  		# Save set number
  		Group <- c(Group,rep(s,steps[i]))
  		# Save sample size vector
  		SSize <- c(SSize, rep(steps[i],length(SampleSet)))
  		# Save run
  		run.tmp <- rep(k,length(SampleSet))
  		run <- c(run,run.tmp)
	  }
	}
}

# Now create data frame
AllData <- data.frame('Subject' = AllSub, 'Run' = run,'Group' = Group, 'Size' = SSize)
















