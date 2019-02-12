####################
#### TITLE:     Give back the correspondig run and step from a sequence of numbers
#### Contents: 	
#### 
#### Source Files: \\FreddieFreeloader/Script.git/Sampling/Coherence
#### First Modified: 24/05/2018
#### Notes: 
#################

##
###############
### Notes
###############
##

# STOPGO dataset

# There are in this scenario:
# A) 50 runs 
# B) 46 steps in each run
# C) varying amount of sets in each step within a run 


##
###############
### Preparation
###############
##

# Take arguments from master file
args <- commandArgs(TRUE)

# Which index are we in?
index <- as.numeric(as.character(args[1]))


##
###############
### RUN and STEP
###############
##

# Number of runs
NRUNS <- 50
# Starting amount of subjects in group analysis
startingTotal <- 10
# Database total
DATTOTAL <- 1400
# Total amount of possible subjects if we have a maximum of 3 disjoint subsets
NTOT <- 460
steps <- seq(startingTotal, NTOT, by = 10)

# Create the sequence
sequence <- data.frame()
for(k in 1:NRUNS){
  for(i in 1:length(steps)){
    # Number of disjoint sets in this sample size
    numDS <- floor(DATTOTAL/steps[i])
    for(s in 1:numDS){
      sequence <- rbind(sequence,
          data.frame('run' = k, 'step' = i, 'set' = s))
    }
  }
}

# Dimension of this data frame
dimSeq <- dim(sequence)[1]

# Now we print the run, step and set to the console, which will be read in by 
# master file
cat(unlist(c(sequence[index,'run'], sequence[index,'step'],sequence[index,'set'])))

