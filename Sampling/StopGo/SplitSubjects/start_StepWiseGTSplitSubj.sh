#!/bin/sh
#
#
#PBS -N FreOveSG
#PBS -o output/
#PBS -e error/
#PBS -m a
#PBS -l walltime=04:25:00
#

#----------------------------------------------------#
# SWAP TO CLUSTER DELCATTY
# module swap cluster/delcatty

# SWAP TO CLUSTER GOLETT
module swap cluster/golett
#----------------------------------------------------#

# New modules
module load FSL/5.0.9-intel-2016a-Mesa-11.2.1
module load R/3.2.3-intel-2016a
. $FSLDIR/etc/fslconf/fsl.sh

#module load R/3.1.0-ictce-5.5.0
#module load FSL/5.0.6-ictce-5.5.0
#module load fmri/1.5-0-ictce-5.5.0-R-3.1.0


#----------------------------------------------------#
# CHANGE YOUR VSC NUMBER HERE AND GOD WILL DO THE REST
vsc=40728
#----------------------------------------------------#


#----------------------------------------------------#
# WHICH SCENARIO DO YOU WANT TO RUN?
scenario=A
#----------------------------------------------------#


srcdir=/user/scratch/gent/gvo000/gvo00022/vsc"$vsc"/Freddie/StopGo/SplitSubjects
cd $srcdir

# There are 7000 jobs to be submitted!
./StepWiseGTSplitSubj.sh ${PBS_ARRAYID} $vsc $scenario $srcdir

echo "job finished"