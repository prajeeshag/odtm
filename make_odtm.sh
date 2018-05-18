#!/bin/bash
set -e

export srcdir=$(pwd)

thisdir=$(pwd)

EXE="odtm.exe"

execdir="$thisdir/exec"
mkmf="$thisdir/bin/mkmf"

mkmftemplate="$thisdir/bin/mkmf.template"

paths="$thisdir/src"

mkdir -p $execdir

cd $execdir

$mkmf -f -p $EXE -t $mkmftemplate $paths

make $@

