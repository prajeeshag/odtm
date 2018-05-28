#!/bin/bash
#########################################
# USER DEFINED PARAMETER SECTION STARTS #
set -ea

NPROCS=16
curr_date=19950101000000
clock_lim=24

exp_name=$(basename $(dirname $(dirname $(pwd))))_$(basename $(pwd))

export NPROCS

function nkjobid { output=$($@); echo $output | head -n1 | cut -d'<' -f2 | cut -d'>' -f1; }

export f77_dump_flag=TRUE
export MALLOC_CHECK_=2

jobid=$(nkjobid bsub -q "cccr-res" -n $NPROCS -J $exp_name -W ${clock_lim}:00 -o stdlog.out -e stdlog.err < odtm_submit.sh)

clc_lim=$((clock_lim*2))

bsub -w "started($jobid)" -q "cccr-res" -W ${clc_lim}:00 -o mppnccombine.out -e mppnccombine.out -J ${exp_name}.post  "./domppnccombine.sh -s $curr_date -p $jobid"

