####################
#### TITLE:     Script where we stepwise increase the sample size of the 
####            group analysis (the GT) to see when stability is achieved, 
####            splitting up subjects.
####
#### Contents:
####
#### Source Files: \\FreddieFreeloader/Script.git/Sampling/
#### First Modified: 18/08/2015
#### Notes:
#################

# Make sure you have all the sample sizes for each step (obtained by scenarioCohen.R)!
# So you need to manually create Results folder and then create SampleSizes in this folder.

#### Structure of folders (NOTE: structure is slightly different on local machine):
# Scripts (location of the main program and all scripts it uses)
#--# Results ====>>> this is the working directory
#-----# SampleSizes (this contains already all the sample sizes so we can work parallel)
#---------# e.g. SubjectsStep1.txt (file with the subjects needed for that specific step)
#-----# Steps
#---------# e.g. Step_1 (results of Step 1)
#-----------# Group
#---------------# e.g. Group1 (analysis of set 1)

# INCENTIVE dataset!

# Actual data can be in any other location!

# NOTE: Group and Set are synonyms in this script...

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pre Step: Global variables
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	# In which index are we?
	index=$1
	# Which is your vsc number
	vsc=$2
	# Which thresholding scenario (uncorrected or fdr)
	scenario=$3
	# Location of scripts (only to source paths)
	loc_scrpt=$4
	# Which computer: HPC or MAC
	COMPUTER=HPC
	# Source the paths
	source "$loc_scrpt"/blind_StepWiseGTDisJoint.sh "$vsc"
	# Location of the data
	if [ $COMPUTER = MAC ] ; then
		data=$dataMAC
	fi
	if [ $COMPUTER = HPC ] ; then
		data=$dataHPC
	fi
	# Significance level in group analysis: uncorrected = 0.001, BH FDR = 0.05
	if [ $scenario = uncorrected ] ; then
	  signLevel=0.001
	fi
		if [ $scenario = fdr ] ; then
	  signLevel=0.05
	fi
	echo .... INDEX "$index" ....


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step one: Preparation
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Location of scripts folder
if [ $COMPUTER = MAC ] ; then
	SCRPT=$SCRPTmac
fi
if [ $COMPUTER = HPC ] ; then
	SCRPT=$SCRPThpc
fi
cd "${SCRPT}"

# In which RUN are we? Get this from the RUNSTEP.R file, based on the index given.
RUN=($(Rscript RUNSTEP.R ${index} | awk '{ print $1 }'))
# In which STEP are we? Get this from the RUNSTEP.R file, based on the index given.
step=($(Rscript RUNSTEP.R ${index} | awk '{ print $2 }'))
# Which group/set do we analyze?
group=($(Rscript RUNSTEP.R ${index} | awk '{ print $3 }'))


# Define workspace: based on run
WD=$SCRPT/Results/Run_$RUN
cd "${WD}"

# Create folder Step_$step for each step. This will be the step working directory
mkdir Step_$step
STEPWD=$WD/Step_$step
cd "${STEPWD}"

# Create folder Group$group for each set. This will be the group/set working directory
mkdir Group$group
GROUPWD=$STEPWD/Group$group


# Define working directory of sample sizes (on local machine, this is other location because it would otherwise appear in GIT folder)
if [ $COMPUTER = MAC ] ; then
	SAMPLEWD=$SAMPLEWDmac/Run_$RUN
fi
if [ $COMPUTER = HPC ] ; then
	SAMPLEWD=$SAMPLEWDhpc/Run_$RUN
fi

# Need to have a temporary folder, one for each run, step and group so there is no overlap in calculations (parallel processing)!
	# Create temp folder
	mkdir $GROUPWD/temp
	temp=$GROUPWD/temp


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step two: Copy and merge all the subjects for this step
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo .... Copy the data ....
# Subglobal variables
subjects=0
cd "${SAMPLEWD}"

# Read the subjects in the .txt file
	IFS=$'\n' read -d '' -r -a subjects < SubjectsStep"$step"_"$group".txt


# Go to the temp foler located at the data folder. Copy the data from the subjects to it.
	cd "${temp}"
	NumSub=${#subjects[@]}
	NumSub=$(($NumSub - 1))

	for j in $(eval echo "{0..$NumSub}"); do
		SubID=${subjects[$j]}
		# Copy contrast map
		scp -r $data/$SubID/EPI_short_MID/con_0032_feedback_hit_largewin_-_smallwin.nii.gz .
		# Rename with subject ID (otherwise we have duplicates in the temp folder)
		mv con_0032_feedback_hit_largewin_-_smallwin.nii.gz con_0032_feedback_hit_largewin_-_smallwin_$SubID.nii.gz
		# Extract
		gunzip con_0032_feedback_hit_largewin_-_smallwin_$SubID.nii.gz
		# Same with variance
		scp -r $data/$SubID/EPI_short_MID/varcon_0032_feedback_hit_largewin_-_smallwin.nii.gz .
		mv varcon_0032_feedback_hit_largewin_-_smallwin.nii.gz varcon_0032_feedback_hit_largewin_-_smallwin_$SubID.nii.gz
		gunzip varcon_0032_feedback_hit_largewin_-_smallwin_$SubID.nii.gz
		# Copy the universal mask to the temp folder, rename it with a subject index and then send it to the GROUPWD
		scp -r $data/AllSubMaskInc.nii .
		mv AllSubMaskInc.nii AllSubMaskInc_$SubID.nii
		mv AllSubMaskInc_$SubID.nii "$GROUPWD/"
		# Array with ConFile and VarConFile
		ConFile=("${ConFile[@]}" con_0032_feedback_hit_largewin_-_smallwin_$SubID.nii)
		VarConFile=("${VarConFile[@]}" varcon_0032_feedback_hit_largewin_-_smallwin_$SubID.nii)
	done

	# Use fslmerge with output again to group analysis folder and then extract the files.
	$FSLDIR/bin/fslmerge -t "$GROUPWD/cope" ${ConFile[@]}
	$FSLDIR/bin/fslmerge -t "$GROUPWD/varcope" ${VarConFile[@]}
	cd "${GROUPWD}"
		gunzip cope.nii.gz
		gunzip varcope.nii.gz

	# Now create the needed files for flameo (design.mat, design.con and design.grp) through R file designCrossVal.R
	Rscript "${SCRPT}"/designCrossVal.R "$GROUPWD/" "$(($NumSub + 1))" &> ROutput_designCrossVal.txt

	# Merge the masks: 4D mask
	$FSLDIR/bin/fslmerge -t "${GROUPWD}"/mask.nii AllSubMaskInc*.nii

	# Remove individual masks
	rm -r AllSubMaskInc*.nii

	# Now perform flameo to pool subjects: cope is a 4D data structure (4th dimension are the subjects)
	echo .... Pooling subjects ....
	$FSLDIR/bin/flameo --cope=cope --vc=varcope --mask=mask --ld=stats --dm=design.mat --cs=design.grp --tc=design.con --runmode=flame1 &> pooling.txt


# Clean up temp folder
cd "${temp}"
rm -r ../temp
# Re-initialize variables
unset -v 'ConFile'
unset -v 'VarConFile'
unset -v 'subjects'

cd "${WD}"



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step three: inference (voxelwise Uncorrected/FDR)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo .... Voxel wise inference "for" group analysis, "(un)corrected" "for" multiple testing"!!!" ....
cd "${GROUPWD}"

# Some preparation
/bin/rm -f stats/zem* stats/zols* stats/mask* # zem and zols are unknown actually... [taken from Nichols blog, so I leave it here]
gunzip mask.nii.gz

# Inference: depending on scenario either uncorrected P < 0.001 or BH FDR < 0.05
Rscript "${SCRPT}"/voxelInference.R "$GROUPWD" "$GROUPWD/stats" "$signLevel" "$scenario" &> ROutput_voxelInference.txt



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step four: copy header information
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo .... Copy Header Information ....
# Copy header information from subject 33 (random subject): same header, assume safe to copy.
scp -r $data/33/EPI_short_MID/con_0032_feedback_hit_largewin_-_smallwin.nii.gz .
gunzip con_0032_feedback_hit_largewin_-_smallwin.nii.gz 
$FSLDIR/bin/fslcpgeom con_0032_feedback_hit_largewin_-_smallwin.nii thresh_zstat1.nii
# Also for tstat1.nii, because we use it later on in transformations
		cd stats
		scp -r $data/33/EPI_short_MID/con_0032_feedback_hit_largewin_-_smallwin.nii.gz .
		gunzip con_0032_feedback_hit_largewin_-_smallwin.nii.gz
		$FSLDIR/bin/fslcpgeom con_0032_feedback_hit_largewin_-_smallwin.nii tstat1.nii.gz
		$FSLDIR/bin/fslcpgeom con_0032_feedback_hit_largewin_-_smallwin.nii zstat1.nii.gz
		rm -r con_0032_feedback_hit_largewin_-_smallwin.nii
		cd ..
rm -r con_0032_feedback_hit_largewin_-_smallwin.nii



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step five: remove cope and varcope (too much data)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo "Removing COPE, VARCOPE and weights1"
cd "${GROUPWD}"

rm -r cope.nii
rm -r varcope.nii
rm -r stats/weights1.nii.gz
rm -r stats/cope1.nii.gz
rm -r stats/mean_random_effects_var1.nii.gz
rm -r stats/pe1.nii.gz
rm -r stats/res4d.nii.gz
rm -r stats/tdof_t1.nii.gz
rm -r stats/varcope1.nii.gz
rm -r stats/zflame1lowertstat1.nii.gz
rm -r stats/zflame1uppertstat1.nii.gz



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step six: process the mask (only keep first slice otherwise too much data!)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Rscript "${SCRPT}"/maskSlicing.R "$GROUPWD" &> SlicingMask.txt
