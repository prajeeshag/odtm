#!/bin/bash
#set -x

datefunc="./datefunc"
diag_table="diag_table"
mppnccombinecmd="./mppnccombine -n4 -c -r"
outdir="OUTPUT"

if [[ ! -e $datefunc ]]; then
    echo Error: $datefunc do not exist.
    exit 1
fi	

export TZ=UTC


npes=0
pjobid=-1

# Function for submitting parallel combine jobs
maxnpjob=8
parjobs=""
function submitpar {
    job=$@

    while [ "true" == "true" ]; do
	njobs=0
	for id in $parjobs; do
	    if kill -0 $id > /dev/null 2>&1; then
		njobs=$((njobs+1))
	    else
		parjobs=${parjobs/$id}
		parjobs=$(echo $parjobs | tr -s " ")
	    fi
	done

	if [ "$njobs" -lt "$maxnpjob" ]; then
	    eval $job &
	    parjobs="$parjobs $!"
	    echo submited job $job
	    break
	else
	    echo waiting for submitting $job
	    sleep 1
	    continue
	fi
    done
}
# End of Function for submitting parrallel combine jobs

while getopts 's:e:n:p:' flag; do
    case "${flag}" in
	s) startdate="$OPTARG" ;;
	e) enddate="$OPTARG" ;;
	n) npes="$OPTARG" ;;
	p) pjobid="$OPTARG" ;;
	*) error "Unexpected option ${flag}" ;;
    esac
done

if [[ -z $startdate ]] || \
       [[ -z $npes ]] || [[ -z $pjobid ]]; then
    echo Usage:
    echo $0 -s startdate -n ocnstrtnpes -p parrentjobid
    exit 1
fi

mkdir -p $outdir

nfile=0

count=0
while read LINE
do
    if [[ -z $LINE ]]; then
	continue
    fi
    
    firstltr=$(echo $LINE | head -c 1)
    if [ $firstltr == "#" ]; then
	continue
    fi

    let count++
    if [ $count -lt 3 ]; then
	continue
    fi

    col2=$(echo $LINE | awk -F"," '{print $2}')
    if ! [ "$col2" -eq "$col2" ] 2>/dev/null; then
  	break
    fi
    
    nfile=$((nfile+1))
    
    tmpfmt=$(echo "$LINE" | awk -F"," '{print $1}')
    tmpfmt=${tmpfmt//\"/}
    
    pos=$(($(expr index $tmpfmt "%") - 1 ))
    bname=${tmpfmt:0:$pos}
    tmpfmt=${tmpfmt:$pos}
    
    tmpfmt=${tmpfmt//\%/_%0}
    tmpfmt=${tmpfmt/yr/d}
    tmpfmt=${tmpfmt/mo/d}
    tmpfmt=${tmpfmt/dy/d}
    tmpfmt=${tmpfmt/hr/d}
    tmpfmt=${tmpfmt/mi/d}
    tmpfmt=${tmpfmt/sc/d}
    tmpfmt=${tmpfmt}.nc
    
    tmpfmt=$bname$tmpfmt
    
    tmpordr=$(echo "$LINE" | awk -F"," '{print $1}')
    tmpordr=${tmpordr//\"/}
    tmpordr=${tmpordr//\%/ }
    tmpordr=$(echo $tmpordr | sed 's/[0-9]/$/g')
    tmpordr=$(echo $tmpordr | awk '{$1=""; print $0}') 
    
    fileprintfcmdl[$nfile]=$(echo printf \"$tmpfmt\" $tmpordr)

    tmp_intr=$(echo "$LINE" | awk -F"," '{print $7}')
    tmp_intrunit=$(echo "$LINE" | awk -F"," '{print $8}')
    
    if [[ -z $tmp_intr || -z $tmp_intrunit ]]; then
	echo Error file reading diag_table
	exit 1
    fi

    
    tmp_intr=${tmp_intr//\"/}
    tmp_intrunit=${tmp_intrunit//\"/}
    tmp_intrunit=${tmp_intrunit,,}
    
    intrl[$nfile]=$tmp_intr
    intrunitl[$nfile]=$tmp_intrunit

    fname=$(echo "$LINE" | awk -F"," '{print $1}')
    fname=${fname//\"/}

    mix=$(awk -F"=" '{ if (tolower($1)~"mix_snapshot_average_field") print $2}' input.nml)
    x=$(awk -v fname="$fname" -F "," '{ if (substr($1,1,1) != "#" && $4 ~ fname) print $6 }' diag_table )

    if [[ $mix == *".true."* ]]
    then
	x=".false."
    fi
    
    if [[ $x == *".false."* ]]
    then
	inst[$nfile]=true
    else

	inst[$nfile]=false
    fi
    
done < $diag_table

i=1
filelist=""
incunit=("years" "months" "days" "hours" "minutes" "seconds")

for i in $(seq 1 $nfile); do 
    dateoffile[$i]=$startdate
    nfiledone[$i]=0
done

chcknxt=2

alldone=false

donefiles=""
firsttime=true

while [ $alldone == false ]; do
    alldone=true
    for i in $(seq 1 $nfile); do
	#check if all files have been done for a particular file type
  	if [ ${nfiledone[$i]} == 0 ]; then
	    alldone=false
	else
	    continue
	fi
	
	fileprintfcmd=${fileprintfcmdl[$i]}
  	inc=(0 0 0 0 0 0)
  	incchk=(0 0 0 0 0 0)
  	intr=${intrl[$i]}
  	intrunit=${intrunitl[$i]}
  	intrunit=${intrunit,,}
	icdt=${dateoffile[$i]}

	# Setting up increment
  	n=0
  	while [ $n -le "5" ]; do
	    if [ $intrunit == ${incunit[$n]} ]; then
		inc[$n]=$intr
  		incchk[$n]=$((intr*chcknxt))
		break
	    fi
  	    n=$((n+1))
	done

	n=0; flag=false
	while [ $n -le "5" ]; do
	    if [ ${inc[$n]} -ne "0" ]; then
		flag=true
  		break
	    fi
	    n=$((n+1))
	done
	if [ $flag != true ]; then
	    echo increment 0 for $fileprintfcmd
  	    exit 1
	fi

  	inc=${inc[@]}
  	inc=${inc// /,}
  	incchk=${incchk[@]}
  	incchk=${incchk// /,}

	if [[ ${inst[$i]} == "false" ]]; then
	    cdt=$($datefunc incdatemid "julian $icdt ${inc}")
	    cdt=${cdt// /}
	else
	    cdt=$($datefunc incdate "julian $icdt ${inc}")
	    cdt=${cdt// /}
	fi
	

	# Done Setting up increment

	# Setting up current filename 
  	yr=${cdt:0:4}; if [ "$yr" -ne "0" ]; then yr=${yr#0}; fi
  	mo=${cdt:4:2}; if [ "$mo" -ne "0" ]; then mo=${mo#0}; fi
  	dy=${cdt:6:2}; if [ "$dy" -ne "0" ]; then dy=${dy#0}; fi
  	hr=${cdt:8:2}; if [ "$hr" -ne "0" ]; then hr=${hr#0}; fi 
  	mi=${cdt:10:2}; if [ "$mi" -ne "0" ]; then mi=${mi#0}; fi
  	sc=${cdt:12:2}; if [ "$sc" -ne "0" ]; then sc=${sc#0}; fi
  	ffdfile=$(eval ${fileprintfcmd})
	# Done Setting up current filename 
	
	# Setting up check filename 
	cdtchk=$($datefunc incdate "julian $cdt ${incchk}")
	cdtchk=${cdtchk// /}

  	yr=${cdtchk:0:4}; if [ "$yr" -ne "0" ]; then yr=${yr#0}; fi
  	mo=${cdtchk:4:2}; if [ "$mo" -ne "0" ]; then mo=${mo#0}; fi
  	dy=${cdtchk:6:2}; if [ "$dy" -ne "0" ]; then dy=${dy#0}; fi
  	hr=${cdtchk:8:2}; if [ "$hr" -ne "0" ]; then hr=${hr#0}; fi 
  	mi=${cdtchk:10:2}; if [ "$mi" -ne "0" ]; then mi=${mi#0}; fi
  	sc=${cdtchk:12:2}; if [ "$sc" -ne "0" ]; then sc=${sc#0}; fi
  	fcffile=$(eval ${fileprintfcmd})
	# Done Setting up check filename 
	
	if [ $pjobid -gt "0" ]; then 
    	    stat=$(bjobs -o -noheader stat $pjobid)
	else
	    stat="DONE"
	fi

	if [ $stat == "EXIT" ] || \
	       [ $stat == "DONE" ]; then
	    fcffile=$ffdfile
	fi

	fdfile=$ffdfile
	cffile=$fcffile

	if [ $npes -ne 0 ] && [ ${fdfile:0:5} == "ocean" ]; then
      	    cnpes=$(printf "%04d" $npes)
            fdfile=${fdfile}.$cnpes
            cffile=${cffile}.$cnpes
            combinecmd="$mppnccombinecmd -n $npes"
	else
	    fdfile=${fdfile}.0000
            cffile=${cffile}.0000
            combinecmd=$mppnccombinecmd
	fi

	if [[ -e $cffile ]] && [[ ! -e $fdfile ]]; then
	    echo skipping $ffdfile
    	    icdt=$($datefunc incdate "julian $icdt ${inc}")
    	    dateoffile[$i]=${icdt// /}
	    continue
	fi
	
	if [[ -e $fdfile ]] && [[ -e $cffile ]]; then
	    echo $ffdfile submitting for combine 
    	    submitpar "$combinecmd $ffdfile && mv $ffdfile $outdir"
	    donefiles="$donefiles $ffdfile"
    	    icdt=$($datefunc incdate "julian $icdt ${inc}")
    	    dateoffile[$i]=${icdt// /}
	    continue
	fi

	if [[ ! -e $cffile ]] && [[ ! -e $fdfile ]]; then
	    #echo files $fdfile and $cffile not found!
	    if [ $stat == "EXIT" ] || \
		   [ $stat == "DONE" ]; then
		nfiledone[$i]=-1
		continue
	    fi
	fi

	if [[ ! -e $cffile ]] && [[ -e $fdfile ]]; then
	    #echo $fdfile found but $cffile not found!!!
	    continue
	fi
	
    done
    
    if [ $stat == "EXIT" ] || [ $stat == "DONE" ]; then
	continue
    else
	sleep 0.5
    fi
done

wait


