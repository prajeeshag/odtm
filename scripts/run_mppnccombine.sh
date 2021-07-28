#!/bin/bash
#USER DEFINED PARAMETERS

startdate=2004,01,01,00,00,00
calendar_type=4
enddate=2013,12,31,00,00,00

# queue - Queue name
queue="cccr"

# WCLOCK - Wall clock limit (in hours)
WCLOCK="240"

# nproc_combine - number of processors for postproc
nproc_combine=36





# END OF USER DEFINED PARAMETERS

#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#--------------------------------------------------------------------------------   
#-------------------------------------------------------------------------------- 

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

aprun -n $nproc_combine $RUNNCCP2R \
<<< "&opts_nml removein=1, atmpes=$npes, child_run=$combine_child_run, startdate=$startdate, calendar_type=$calendar_type, enddate=$enddate /" &> ${JOBNAME}.post.out
EOF

qsub < $tfile
