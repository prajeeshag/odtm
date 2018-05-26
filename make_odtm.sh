
#!/bin/bash

set -e

export srcdir=$(pwd)

thisdir=$(pwd)

while getopts 'd:' flag; do
    case "${flag}" in
    d) debug='.debug' ;;
    e) enddate="$OPTARG" ;;
    n) npes="$OPTARG" ;;
    p) pjobid="$OPTARG" ;;
    *) error "Unexpected option ${flag}" ;;
    esac
done

EXE="odtm.exe"

execdir="$thisdir/exec"
mkmf="$thisdir/bin/mkmf"

mkmftemplate="$thisdir/bin/mkmf.template"

FMS_UTILS=$thisdir/src/fms_shared

FMS_UTILITIES="$FMS_UTILS/include \
			   $FMS_UTILS/platform \
               $FMS_UTILS/constants \
			   $FMS_UTILS/fms \
			   $FMS_UTILS/time_manager \
				$FMS_UTILS/mpp \
				$FMS_UTILS/diag_manager  \
				$FMS_UTILS/memutils \
				$FMS_UTILS/constants \
				$FMS_UTILS/mpp/include \
				$FMS_UTILS/data_override \
				$FMS_UTILS/horiz_interp \
				$FMS_UTILS/time_interp \
				$FMS_UTILS/axis_utils \
				$FMS_UTILS/mosaic"

paths="$thisdir/src/odtm"

mkdir -p $execdir/lib_fms


cd $execdir/lib_fms

$mkmf -f -p lib_fms.a -t $mkmftemplate $FMS_UTILITIES

make -j 16

mkdir -p $execdir/odtm

cd $execdir/odtm

$mkmf -f -p $EXE -t $mkmftemplate -o "-I$execdir/lib_fms" -l "$execdir/lib_fms/lib_fms.a" $paths

make $@

