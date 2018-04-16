####################
#### TITLE:     Sample stepwise 10subjects from the total pool of all 
####            subjects, without controlling for external factors
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Script.git/Sampling/SplitSubjects
#### First Modified: 17/08/2015
#### Notes: 
#################

##
###############
### Notes
###############
##

# SCENARIO A: 
# The idea is to take 10 subjects to start with and sample 10 OTHER subjects to compare with.
# The sampling is complete at random (sole restriction is that the subjects are different).
# This is repeated for 50 iterations.
# The same is done for sample size 20, 30,..., 700.

# However, this might take a very long time to do sequentially. Hence, we can use the HPC to work on in parallel. 
# We already sample the subjects.
# Then each analysis can be executed independent of the previous one. 
# Later on we can calculate the overlap.


##
###############
### Preparation
###############
##

# Reset wm
rm(list=ls())

# Provide working directory
wd <- "/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/SampleSizes/SplitSubjects/ScenarioA_redo"


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

# Total amount of subjects
NTOT <- 700


##
###############
### Sampling
###############
##


# We add sample sizes starting from 10 to 700.
steps <- seq(startingTotal,NTOT,by=10)

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
	print(paste('In run ',k,sep=''))
	seed <- StartingSeed*k
	set.seed(seed)
	# Start sequence of increasing sample sizes (steps)
	for(i in 1:length(steps)){
		# Possible subjects for group A
		ToSample <- OurIDs[,'OwnID']
		SampleA <- sample(ToSample,steps[i],replace=FALSE)
			IDA <- OurIDs[,'OwnID'] %in% SampleA
		# Now have subjects for group B (other subjects)
		ToSample <- OurIDs[!IDA,'OwnID']
		SampleB <- sample(ToSample,steps[i],replace=FALSE)

		# Print them to files: A and B are called _1 and _2!
		cat(SampleA,sep='\n',file=paste(wd,'/Run_',k,'/SubjectsStep',i,'_1.txt',sep=''))
		cat(SampleB,sep='\n',file=paste(wd,'/Run_',k,'/SubjectsStep',i,'_2.txt',sep=''))

		# Now have them in a dataframe
		AllSub <- c(AllSub,SampleA,SampleB)
		# Save group letter (1:A, 2:B)
		Group <- c(Group,rep(1,steps[i]),rep(2,steps[i]))

		# Save sample size vector
		SSize <- c(SSize, rep(steps[i],length(c(SampleA,SampleB))))

		# Save run
		run.tmp <- rep(k,length(c(SampleA,SampleB)))
		run <- c(run,run.tmp)
	}
}

# Now create data frame
AllData <- data.frame('Subject' = AllSub, 'Run' = run,'Group' = Group, 'Size' = SSize)
















