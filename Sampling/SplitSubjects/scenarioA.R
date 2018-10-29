####################
#### TITLE:     Sample stepwise 10subjects from the total pool of all 
####            subjects
#### Contents: 	
#### 
#### Source Files: /Volumes/2_TB_WD_Elements_10B8_Han/PhD/I****DATA/Data/FreddieFreeloader/Script.git/Sampling/SplitSubjects
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

# Source information not provided on Github (paths and IDs)
source(blind_scenarioA.R)

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
















