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

# 20180621 Update: allow manually assignment of different CPUs for different nodes
# Syntax: ./submit_job_ind.sh 1 10 2 10 3 10 4 5 5 5  --> 10 CPUs for 1-3 nodes, and 5 CPUs for 4-5 nodes

# 20180627 Update: Acquire number of CPUs available automatically and allocate tasks
# if no  CPU available for all nodes now, check status, run submit_job_ind.sh and manually input #nodes and #CPUs as above  

totalNode=0 # to count available nodes
totalCPU=0 # to count available CPUs

if [ $# -gt 0 ]  # Manually control
then

	config=("$@")   # Convert input to array, length should be an even number
	nNode=$[${#config[@]}/2] # Get node number
	totalCPU=0

	for ii in `seq 1 $nNode`
	do
		thisNode=${config[(ii-1)*2]}   # Get this node name
		thisCPU=${config[ii*2-1]}      # Get CPU requested for this node
	    
		echo " "  
		echo $ii". Submitting to Node clc00"$thisNode", requesting "$thisCPU" CPUs"
		# Pass all parameters including config to job_ind.sh
		qsub -l h="clc00$thisNode" -pe openmpi $thisCPU -v n=$ii,config="${config[*]}" job_ind.sh

		totalCPU=$[totalCPU+thisCPU]
	done

	totalNode=$nNode

else

	for ii in `seq 1 5`
	do
	  thisNode=$ii   # Get this node name
          
	  # get number of CPU running of this node 
	  CPURunThisNode=`qstat -q "*@clc00$thisNode" -u "*" -f | grep "  r  " \
		| awk '{print $8}' | awk '{total = total + $1}END{print total}'` 
	  
	  # calculate the number of CPU available of this node  loadThisNode=$[${loadThisNode%.*}+1]
	  CPUThisNode=$[24-$CPURunThisNode-4] 
	   
	  if [ $CPUThisNode -gt 0 ] # 
	  then
	    totalNode=`expr $totalNode + 1`
	    thisCPU=$CPUThisNode  
	    config[$totalNode*2-1]=$thisNode     # name of this available node
	    config[$totalNode*2]=$thisCPU  # number of CPUs requested for this node
	  fi
	  
	  totalCPU=$[totalCPU+thisCPU]
	  
	done

	for ii in `seq 1 $totalNode`
	do
	    thisNode=${config[$ii*2-1]}
	    thisCPU=${config[$ii*2]}
	    echo " "  
	    echo ". Submitting to Node clc00"$thisNode", requesting "$thisCPU" CPUs"
	    # Pass all parameters including config to job_ind.sh
	    qsub -l h="clc00$thisNode" -pe openmpi $thisCPU -v n=$ii,config="${config[*]}" job_ind.sh
	done

	if [ $totalNode -eq 0 ]  #  No CPU available for all nodes now, print status
	then
	  echo " "
	  echo "No CPU available now"
	  echo " "
	  echo "check status..."
	  echo " "
	  qstat -u "*" -f
	  echo " "
	  echo "Input nodes and CPUS manually:"
	  echo " "
	  read -a config
	  totalNode=$[${#config[@]}/2] # Get node number
	  totalCPU=0 # to count available CPUs
	  
	  for ii in `seq 1 $totalNode`
	  do
	    thisNode=${config[(ii-1)*2]}   # Get this node name
	    thisCPU=${config[ii*2-1]}      # Get CPU requested for this node

	    echo " "
	    echo $ii". Submitting to Node clc00"$thisNode", requesting "$thisCPU" CPUs"

	    # Pass all parameters including config to job_ind.sh
	      qsub -l h="clc00$thisNode" -pe openmpi $thisCPU -v n=$ii,config="${config[*]}" job_ind.sh
	    
	    totalCPU=$[totalCPU+thisCPU]
	    
	  done
	    
	fi

fi

echo " "  
echo "Total number of nodes: "$totalNode  
echo "Total number of CPUs: "$totalCPU
echo " "  
qstat
