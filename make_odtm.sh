#!/bin/bash
set -e

export srcdir=$(pwd)

thisdir=$(pwd)

EXE="odtm.exe"

execdir="$thisdir/exec"
mkmf="$thisdir/bin/mkmf"

mkmftemplate="$thisdir/bin/mkmf.template"

FMS_UTILS=$thisdir/src/fms_shared

FMS_UTILITIES="$FMS_UTILS/include $FMS_UTILS/platform $FMS_UTILS/constants $FMS_UTILS/fms $FMS_UTILS/time_manager $FMS_UTILS/mpp $FMS_UTILS/diag_manager  $FMS_UTILS/memutils $FMS_UTILS/constants $FMS_UTILS/mpp/include"

paths="$thisdir/src/odtm"

mkdir -p $execdir

cd $execdir

$mkmf -f -p $EXE -t $mkmftemplate $paths $FMS_UTILITIES

make $@

