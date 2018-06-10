#! /bin/bash

# Bash script for sumbitting matlab jobs to an SGE cluster queue.
# This is the script you run from the command line to submit the jobs.
# By David Black-Schaffer, June 2007.
# Permission to use and modify this script is granted.
# I am not responsible for any errors in this script, so be forewarned!


# Modify this to set the number of jobs you want to run
# qsub -t 1,2:1 job.sh


# Array is not good because I can't avoid they cross two nodes
# qsub -q fast.q -t 1,$1:1 -v TotalNodes=$1 job_array.sh # Job array to three nodes

# Instead, by generating independent jobs, I can use ' -l h="xxx" ' to constrain them.
# Call this by  ./submit_job node1 node2 node3...
i=0
for node # Equivalent to "for node in "$@"" which loops over input variables
do
	((i++))
	qsub -l h="clc00$node" -v ThisNode=$i,TotalNodes=$# job_ind.sh # Independent jobs
	echo "Submitted ($i/$#)th part to clc00$node"
done
qstat
