#!/bin/bash
set -e

mach=generic_intel

while getopts 'rm:' flag; do
    case "${flag}" in
    m) mach="$OPTARG" ;;
    *)
        echo "error"
        exit 1
        ;;
    esac
done

if [ -f .env ]; then
	echo ".env file already exist!"
	exit
fi

rootdir=$(pwd)

avmach=''
for f in bin/env.*; do
	mc=$(echo $f | sed 's/bin\/env.//g')
	avmach="$avmach $mc"
done

if [[ "$mach" != "none" ]] && [[ "$avmach" == *"$mach"* ]]; then
	echo "export MACH=$mach" >> .env
	echo "Setting Machine as : $mach"
else
	echo "Specify a machine using the option -m"
	echo "For eg. $0 -m machine_name"
	echo 
	echo "Available machines are: " $avmach
	exit
fi

	
echo "export rootdir=$rootdir" >> .env


