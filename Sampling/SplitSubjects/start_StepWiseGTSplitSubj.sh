#!/bin/sh
#
#
#PBS -N FreddieFSCNA
#PBS -o output/output.file
#PBS -e error/error.file
#PBS -m a
#PBS -l walltime=11:25:00
#PBS -l vmem=30GB
#

module load R/3.1.0-ictce-5.5.0
module load FSL/5.0.6-ictce-5.5.0
module load fmri/1.5-0-ictce-5.5.0-R-3.1.0
. $FSLDIR/etc/fslconf/fsl.sh

#----------------------------------------------------#
# CHANGE YOUR VSC NUMBER HERE AND GOD WILL DO THE REST
vsc=40728
#----------------------------------------------------#


#----------------------------------------------------#
# WHICH SCENARIO DO YOU WANT TO RUN?
scenario=A
#----------------------------------------------------#


srcdir=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/SplitSubjects
#srcdir=/user/scratch/gent/vsc407/vsc"$vsc"/Freddie/SplitSubjects"$scenario"
cd $srcdir

./StepWiseGTSplitSubj.sh ${PBS_ARRAYID} $vsc $scenario

echo "job finished"