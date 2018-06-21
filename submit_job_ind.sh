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

#i=0
#for node # Equivalent to "for node in "$@"" which loops over input variables
#do
#	((i++))
#	qsub -l h="clc00$node" -v ThisNode=$i,TotalNodes=$# job_ind.sh # Independent jobs
#	echo "Submitted ($i/$#)th part to clc00$node"
#done

# 20180621 Update: allow manually assignment of different CPUs for different nodes
# Syntax: ./submit_job_ind.sh 1 10 2 10 3 10 4 5 5 5  --> 10 CPUs for 1-3 nodes, and 5 CPUs for 4-5 nodes

config=("$@")   # Convert input to array, length should be an even number
nNode=$[${#config[@]}/2] # Get node number
totalNode=0

for ii in `seq 1 $nNode`
do
	thisNode=${config[(ii-1)*2]}   # Get this node name
	thisCPU=${config[ii*2-1]}      # Get CPU requested for this node
    
	echo " "  
	echo $ii". Submitting to Node clc00"$thisNode", requesting "$thisCPU" CPUs"
	# Pass all parameters including config to job_ind.sh
        qsub -l h="clc00$thisNode" -pe openmpi $thisCPU -v n=$ii,config="${config[*]}" job_ind.sh

	totalNode=$[totalNode+thisCPU]
done
echo "Total number of nodes: "$totalNode
echo ""
qstat
