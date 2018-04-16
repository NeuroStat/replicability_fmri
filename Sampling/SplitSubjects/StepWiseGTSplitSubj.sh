#!/bin/sh

####################
#### TITLE:     SCENARIO A, B and C: Program to stepwise increase the sample size of the group analysis (the GT) to see when stability is achieved, controlling for centers and splitting up subjects.
####			
#### Contents: 	
#### 
#### Source Files: \\FreddieFreeloader/Script.git/Sampling/SplitSubjects
#### First Modified: 18/08/2015
#### Notes: 
#################

# Make sure you have all the sample sizes for each step (obtained by scenarioA.R)!
# So you need to manually create Results folder and then create SampleSizes in this folder.

#### Structure of folders (NOTE: structure is slightly different on local machine):
# Scripts (location of the main program and all scripts it uses)
#--# Results ====>>> this is the working directory
#-----# SampleSizes (this contains already all the sample sizes so we can work parallel)
#---------# e.g. SubjectsStep1.txt (file with the subjects needed for that specific step)
#-----# Steps
#---------# e.g. Step_1 (results of Step 1)
#-----------# Group
#---------------# e.g. Group1 (comparing group 1 vs group 2)



# Actual IMAGEN data can be in any other location!


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Pre Step: Global variables
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	# In which index are we?
	index=$1
	# Which is your vsc number
	vsc=$2
	# Which scenario: A, B or C
	scenario=$3
	# Which computer: HPC or MAC
	COMPUTER=HPC
	# Location of the data
	if [ $COMPUTER = MAC ] ; then
		data=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/IMAGEN/IMAGEN/IMAGEN/Renamed/processed/data
	fi
	if [ $COMPUTER = HPC ] ; then
		data=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/Imagen/data
	fi
	# Significance level in group analysis
	signLevel=0.05
	echo .... INDEX "$index" ....
	
	
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step one: Preparation
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Location of scripts folder
if [ $COMPUTER = MAC ] ; then
	SCRPT=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/Script.git/Sampling/SplitSubjects
fi
if [ $COMPUTER = HPC ] ; then
	SCRPT=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/SplitSubjects
	#SCRPT=/user/scratch/gent/vsc407/vsc"$vsc"/Freddie/SplitSubjects"$scenario"
fi
cd "${SCRPT}"

# In which RUN are we? Get this from the RUNSTEP.R file, based on the index given.
RUN=($(Rscript RUNSTEP.R ${index} "$scenario" | awk '{ print $1 }'))
# In which STEP are we? Get this from the RUNSTEP.R file, based on the index given.
step=($(Rscript RUNSTEP.R ${index} "$scenario" | awk '{ print $2 }'))
# Which group do we analyze? Group A or B?
group=($(Rscript RUNSTEP.R ${index} "$scenario" | awk '{ print $3 }'))


# Define workspace: based on run 
WD=$SCRPT/Results/Run_$RUN
cd "${WD}"

# Create folder Step_$step for each step. This will be the step working directory
mkdir Step_$step
STEPWD=$WD/Step_$step
cd "${STEPWD}"

# Create folder Group$group for each of the two groups. This will be the group working directory
mkdir Group$group
GROUPWD=$STEPWD/Group$group


# Define working directory of sample sizes (on local machine, this is other location because it would otherwise appear in GIT folder)
if [ $COMPUTER = MAC ] ; then
	SAMPLEWD=/Volumes/2_TB_WD_Elements_10B8_Han/PhD/IMAGENDATA/Data/FreddieFreeloader/SampleSizes/SplitSubjects/Scenario"$scenario"/Run_$RUN
fi
if [ $COMPUTER = HPC ] ; then
	SAMPLEWD=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/SplitSubjects/Results/SampleSizes/Scenario"$scenario"/Run_$RUN
	#SAMPLEWD=/user/scratch/gent/vsc407/vsc"$vsc"/Freddie/SplitSubjects"$scenario"/Results/SampleSizes/Scenario"$scenario"/Run_$RUN
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
		scp -r $data/$SubID/SessionB/EPI_global/swea/con_0023_$SubID.nii .
		scp -r $data/$SubID/SessionB/EPI_global/swea/varcon_0023_$SubID.nii .
		# Copy the universal mask to the temp folder, rename it with a subject index and then send it to the GROUPWD
		scp -r $data/AllSubMask.nii .
		mv AllSubMask.nii AllSubMask_$SubID.nii
		mv AllSubMask_$SubID.nii "$GROUPWD/"
		# Array with ConFile and VarConFile
		ConFile=("${ConFile[@]}" con_0023_$SubID.nii)
		VarConFile=("${VarConFile[@]}" varcon_0023_$SubID.nii)
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
	$FSLDIR/bin/fslmerge -t "${GROUPWD}"/mask.nii AllSubMask*.nii
	
	# Remove individual masks
	rm -r AllSubMask*.nii

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
# Step three: inference (voxelwise FDR)
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

echo .... Voxel wise inference "for" group analysis, corrected "for" multiple testing"!!!" ....
cd "${GROUPWD}"

# Some preparation
/bin/rm -f stats/zem* stats/zols* stats/mask* # zem and zols are unknown actually... [taken from Nichols blog, so I leave it here]
gunzip mask.nii.gz

Rscript "${SCRPT}"/voxelInference.R "$GROUPWD" "$GROUPWD/stats" "$signLevel" &> ROutput_voxelInference.txt



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step four: copy header information
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
echo .... Copy Header Information ....
# Copy header information from subject 33 (random subject): same header, assume safe to copy.
scp -r $data/33/SessionB/EPI_global/swea/con_0023_33.nii .
$FSLDIR/bin/fslcpgeom con_0023_33.nii thresh_zstat1.nii
# Also for tstat1.nii, because we use it later on in transformations
		cd stats
		scp -r $data/33/SessionB/EPI_global/swea/con_0023_33.nii .
		$FSLDIR/bin/fslcpgeom con_0023_33.nii tstat1.nii.gz
		$FSLDIR/bin/fslcpgeom con_0023_33.nii zstat1.nii.gz
		rm -r con_0023_33.nii		
		cd ..
rm -r con_0023_33.nii



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



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Step seven: check if there are results
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ! [ -f "GroupAnalysis.pdf" ]
	then
	echo "0 ---ERROR--- NO THRESHOLDED RESULT ON NUMBER $index"
else 
	echo " 1 ----- FINISHED -----"
fi



















