####################
#### TITLE:     Sample in steps of 10 subjects from the total pool of all 
####            subjects, without controlling for external factors and create
####            independent disjoint sets.
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/I****DATA/Data/FreddieFreeloader/Script.git/Sampling/Faces/Coherence
#### First Modified: 23/05/2018
#### Notes: 
#################

##
###############
### Notes
###############
##

# FACES dataset

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

# Source information not provided on Github (paths and IDs)
source('blind_scenarioCohen.R')

# Set the starting seed
StartingSeed <- 12
set.seed(StartingSeed)

# Libraries
library(RColorBrewer)
library(ggplot2)
library(grid)

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
















