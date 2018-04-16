####################
#### TITLE:     Give back the correspondig run and step from a sequence of numbers for experiment controlling for centers
#### Contents: 	
#### 
#### Source Files: \\FreddieFreeloader/Script.git/Sampling/CenterEffect
#### First Modified: 26/06/2015
#### Notes: 
#################

##
###############
### Notes
###############
##


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
	# A) 42 runs 
	# B) 70 steps in each run
	# C) 2 groups in each step within a run
	if(scenario=='A'){NRUNS <- 42;NSTEP <- 70}


# There are in scenario B:
	# A) 42 runs 
	# B) 14 steps in each run
	# C) 2 groups in each step within a run
	if(scenario=='B'){NRUNS <- 42;NSTEP <- 13}


# There are in scenario C:
	# A) 42 runs 
	# B) 7 steps in each run
	# C) 2 groups in each step within a run
	if(scenario=='C'){NRUNS <- 42;NSTEP <- 6}


RUNsequence <- rep(c(1:NRUNS),each=(NSTEP*2))
	run <- RUNsequence[index]
STEPsequence <- rep(rep(c(1:NSTEP),each=2),NRUNS)
	step <- STEPsequence[index]
GROUPsequence <- rep(c(1,2),(NSTEP*NRUNS))
	group <- GROUPsequence[index]
cat(c(run,step,group))






