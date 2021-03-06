####################
#### TITLE:     Give back the correspondig run and step from a sequence of numbers
#### Contents: 	
#### 
#### Source Files: \\FreddieFreeloader/Script.git/Sampling/SplitSubjects
#### First Modified: 26/06/2015
#### Notes: 
#################

##
###############
### Notes
###############
##

# INCENTIVE dataset

##
###############
### Preparation
###############
##

# Take arguments from master file
args <- commandArgs(TRUE)

# Which index are we in?
index <- as.numeric(as.character(args[1]))

# Which scenario is this?
scenario <- as.character(args[2])

##
###############
### RUN and STEP
###############
##

# There are in scenario A:
	# A) 50 runs 
	# B) 70 steps in each run
	# C) 2 groups in each step within a run
	if(scenario=='A'){NRUNS <- 50;NSTEP <- 70}

# Provide the sequenced numbers
RUNsequence <- rep(c(1:NRUNS),each=(NSTEP*2))
	run <- RUNsequence[index]
STEPsequence <- rep(rep(c(1:NSTEP),each=2),NRUNS)
	step <- STEPsequence[index]
GROUPsequence <- rep(c(1,2),(NSTEP*NRUNS))
	group <- GROUPsequence[index]
cat(c(run,step,group))






