#!/bin/bash
exp_name=$(basename $(dirname $(dirname $(pwd))))_$(basename $(pwd))

tmpq=($(bjobs -noheader -J $exp_name -o 'exec_cwd jobid' 2>/dev/null))

jobsfound=true

if [ "${#tmpq[*]}" -eq "0" ]; then
	jobsfound=false
fi

status="RUNNING"
if [ "${tmpq[0]}" == "-" ]; then
	status="PENDING"
fi

if [ $jobsfound == "true" ]; then
  echo
  echo  .....................................................................
  echo "       FOUND $status JOBS FOR JOBNAME $exp_name .."
  echo  .....................................................................
  echo  .......NOT RECOMENDED TO PROCEED WITHOUT KILLING THE JOB.............
  echo  .....................................................................
  echo
  echo "JOB DIRECTORY IS: " ${tmpq[0]}
  echo "JOB ID IS: " ${tmpq[1]}
	echo "ENTER 'yes' FOR KILLING THE JOB (else press ctrl+c):"  
  read tmp
	if [ $tmp == "yes" ]; then
  	bkill ${tmpq[1]}
  fi
  exit 1
fi


exp_name=${exp_name}.post

tmpq=($(bjobs -noheader -J $exp_name -o 'exec_cwd jobid' 2>/dev/null))

jobsfound=true

if [ "${#tmpq[*]}" -eq "0" ]; then
	jobsfound=false
fi

status="RUNNING"
if [ "${tmpq[0]}" == "-" ]; then
	status="PENDING"
fi

if [ $jobsfound == "true" ]; then
  echo
  echo  .....................................................................
  echo "       FOUND $status POSTPROCESSING JOBS FOR JOBNAME $exp_name .."
  echo  .....................................................................
  echo  .......NOT RECOMENDED TO PROCEED WITHOUT JOB END....... .............
  echo  .....................................................................
  echo
  echo "JOB DIRECTORY IS: " ${tmpq[0]}
  echo "JOB ID IS: " ${tmpq[1]}
  exit 1
fi

echo
echo  .....................................................................
echo "...........Could not Find any job for current directory.............."
echo  .....................................................................
echo  .........................SAFE TO PROCEED.............................
echo  .....................................................................
echo
exit 0
