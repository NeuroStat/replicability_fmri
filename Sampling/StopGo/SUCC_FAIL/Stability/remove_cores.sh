#!/bin/sh


# Number of runs
NRUN=50
# Number of steps
NSTEP=70
# Number of groups
NGROUP=2

# Run over the runs
for i in $(eval echo "{1..$NRUN}"); do
	cd Run_"$i"/
	# Run over the steps
	for j in $(eval echo "{1..$NSTEP}"); do
		cd Step_"$j"/
		# Run over the groups
		for l in $(eval echo "{1..$NGROUP}"); do
			# Go to the folder
			cd Group"$l"
			# Check if core file exists
			if test -n "$(find . -maxdepth 1 -name 'core.*' -print -quit)"
			then
				# If so remove
				rm core.*
			fi
			cd ..
		done
		cd ..
	done
	cd ..
done



echo "A job well done!"		

