####################
#### TITLE:     Create the design.mat, design.con and design.grp files
#### Contents: 	
#### 
#### Source Files: //Meta\ Analyis/R\ Code/Studie_FixRan/FixRanStudy.git/Imagen
#### First Modified: 14/11/2014
#### Notes: 
#################

##
###############
### Notes
###############
##

# We also need a mask, so we read in a contrast file and put all non zero elements to 1

##
###############
### Preparation
###############
##

# Library
library(oro.nifti)
	# Print which libraries are needed
	print("Need package: AnalyzeFMRI")

# Take arguments from master file
args <- commandArgs(TRUE)

# Set working directory
wd <- as.character(args)[1]
setwd(wd)
# Number of subjects in the study
NumSub <- as.numeric(as.character(args)[2])


##
###############
### Write design.mat
###############
##

# File connection
fileCon <- paste(wd,"design.mat",sep="")
# Text to be written to the file
cat('/NumWaves\t1
/NumPoints\t',paste(NumSub,sep=''),'
/PPheights\t\t1.000000e+00

/Matrix
',rep("1.000000e+00\n",NumSub),file=fileCon)


##
###############
### Write design.con
###############
##

fileCon <- file("design.con")
	writeLines('/ContrastName1	Group Average 
/NumWaves	1
/NumContrasts	1
/PPheights		1.000000e+00
/RequiredEffect		5.034

/Matrix
1.000000e+00 
',fileCon)
close(fileCon)


##
###############
### Write design.grp
###############
##

# File connection
fileCon <- paste(wd,"design.grp",sep="")
# Text to be written to the file
cat('/NumWaves\t1
/NumPoints\t',paste(NumSub,sep=''),'

/Matrix
',rep("1\n",NumSub),file=fileCon)


# ##
# ###############
# ### Create masks
# ###############
# ##

# # Read in the cope image: all copes from each subject are read in.
# 	# Then we take the union of the voxels being activated in all subjects. 
# 		# First process the cope images such that:
# 			# Voxels with any value get 1.
# 			# Voxels with NaN get 0.
# 	cope <- readNIfTI(paste(wd,'cope.nii',sep=''), verbose=FALSE, warn=-1, reorient=TRUE, call=NULL)[,,,]
# 	cope[which(cope!=0)] <- 1
# 	cope[which(is.nan(cope)==TRUE)] <- 0

# # Now create a mask file
# mask <- array(0,dim=dim(cope[,,,1]))							# A mask of same dimension as one cope. All with zeroes.
# mask[which(apply(cope,c(1,2,3),sum)==dim(cope)[4])] <- 1 		# Only the voxels that sum to the total amount of subjects (dim(cope)[4]), get 1. Apply is used to sum over all subjects (4the dimension)

# # Now write the amount of masks according to the amount of subjects.
# for(j in 1:NumSub){
# 	niftiimage <- nifti(img=mask,dim=dim(cope[,,,1]))
# 	writeNIfTI(niftiimage,paste(wd,'mask_',j,sep=''),gzipped=FALSE)
# }








