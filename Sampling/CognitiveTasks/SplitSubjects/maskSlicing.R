####################
#### TITLE:     Process the masks (only keep slice 1)
#### Contents: 	
#### 
#### Source Files: ///Data/FreddieFreeloader/Script.git/Sampling/SplitSubjects
#### First Modified: 26/20/2015
#### Notes: 
#################

##
###############
### Notes
###############
##

# We only want to keep slice 1 of the masks. 
# For large group analyses, it is not conventient to have lot's of the same slices. 


##
###############
### Preparation
###############
##

# Library
suppressMessages(library(oro.nifti))
	# Print which libraries are needed
	print("Need package oro.nifti")

# Take arguments from master file
args <- commandArgs(TRUE)

# Set working directory
wd <- as.character(args)[1]
setwd(wd)



##
###############
### Keeping slice 1: going from 4D to 3D
###############
##


mask <- readNIfTI(paste(wd,"/mask.nii",sep=""), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,,1]
	# Dimension of the mask
	DIM <- dim(mask)

# Writing the new image
img <- nifti(img=mask,dim=DIM)
writeNIfTI(img,filename=paste(wd,'/masked',sep=''),gzipped=TRUE)
	file.remove(paste(wd,'/mask.nii',sep=''))




