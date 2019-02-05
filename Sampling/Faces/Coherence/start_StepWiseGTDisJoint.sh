#!/bin/sh
#
#
#PBS -N FreddieCohFac
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=03:00:00
#


#----------------------------------------------------#
# SWAP TO CLUSTER DELCATTY
module swap cluster/delcatty
#----------------------------------------------------#


module load R/3.1.0-ictce-5.5.0
module load FSL/5.0.6-ictce-5.5.0
module load fmri/1.5-0-ictce-5.5.0-R-3.1.0
. $FSLDIR/etc/fslconf/fsl.sh

#----------------------------------------------------#
# CHANGE YOUR VSC NUMBER HERE AND GOD WILL DO THE REST
vsc=40728
#----------------------------------------------------#


#----------------------------------------------------#
# WHICH SCENARIO DO YOU WANT TO RUN (uncorrected or fdr)?
scenario=fdr
#----------------------------------------------------#

# Location of scripts 
srcdir=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/Freddie/Faces/Coherence
cd $srcdir

# There are 30050 jobs to be submitted!
# We need to submit in batches so that we can consolidate jobs.
# To do so, we take the array ID and loop 10 batches based on this ID
hpcID=${PBS_ARRAYID}
startIndex=$(( 1+10*($hpcID-1) ))
endIndex=$(( $startIndex+(10-1) ))

# This loop takes i which goes from startIndex to endIndex based on PBS_ARRAYID
for i in $(eval echo "{$startIndex..$endIndex}")
do
	./StepWiseGTDisJoint.sh $i $vsc $scenario $srcdir
done
echo "job finished"
