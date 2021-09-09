#!/bin/bash
set -e
odtmroot=/scratch/cccr/prajeesh/odtm/

$odtmroot/make_odtm.sh

cd $odtmroot/work/test32/
rm -f RESTART/*
rm test32.*
qsub odtm_submit.pbs

cd $odtmroot/work/test64/
rm -f RESTART/*
rm test64.*
qsub odtm_submit.pbs

#cdo -sub test32/RESTART/20120101.003000.odtm_restart.nc test64/RESTART/20120101.003000.odtm_restart.nc out1.nc
#cdo -abs out1.nc out1abs.nc
#cdo -fldsum out1abs.nc out1_1.nc
#cdo -vertsum out1_1.nc out1_2.nc
#
#cdo -sub test32/RESTART/20120101.010000.odtm_restart.nc test64/RESTART/20120101.010000.odtm_restart.nc out2.nc
#cdo -abs out2.nc out2abs.nc
#cdo -fldsum out2abs.nc out2_1.nc
#cdo -vertsum out2_1.nc out2_2.nc

cdo -sub test32/RESTART/odtm_restart.nc test64/RESTART/odtm_restart.nc outf.nc
cdo -abs outf.nc outfabs.nc
cdo -fldsum outfabs.nc outf_1.nc
cdo -vertsum outf_1.nc outf_2.nc
