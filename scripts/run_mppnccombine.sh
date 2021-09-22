#!/bin/bash
#USER DEFINED PARAMETERS

# queue - Queue name
queue="cccr"

# WCLOCK - Wall clock limit (in hours)
WCLOCK="240"

# END OF USER DEFINED PARAMETERS

#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------- 



# nproc_combine - number of processors for postproc
nproc=36
while getopts 'p:n:' flag; do
    case "${flag}" in
    p) jobid=$OPTARG ;;
    n) nproc=$OPTARG ;;
    esac
done


#!/bin/bash
while IFS= read -r line; do
	calendar_type=
done < "INPUT/odtm.res"


JOBNAME=_EXPNAME_
RUNNCCP2R=_ROOTDIR_/exec/run_mppnccp2r/run_mppnccp2r
WCLOCK=$((WCLOCK*2))
ppn=36
nnodes=$((nproc_combine/ppn))
rem=$((nproc_combine%ppn))
if [ "$rem" -gt "0" ]; then
	nnodes=$((nnodes+1))
fi

tfile=$(mktemp)

child_run=0
if [ ! -z "$jobid"]; then
	child_run=1
	COND="PBS -W depend=after:$jobid"
fi

cat <<EOF > 'mppncc.nml'
&opts_nml 
 removein=1, 
 child_run=$child_run 
/
EOF 


cat <<EOF > $tfile
#!/bin/bash
#PBS -q $queue
#PBS -l nodes=$nnodes:ppn=$ppn
#PBS -l walltime=${WCLOCK}:00:00
#PBS -N ${JOBNAME}_post
#PBS -j oe

ulimit -c unlimited
set -xu
cd \$PBS_O_WORKDIR
source _ROOTDIR_/bin/env.pratyush_intel

aprun -n $nproc_combine $RUNNCCP2R &> ${JOBNAME}.mppncc.out

EOF

qsub < $tfile
