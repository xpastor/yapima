#!/bin/bash

while getopts 'c:vCFA' OPTION
do
    case $OPTION in
	c) cflag=1
	   CONFIG_FILE="$OPTARG"
	   ;;
	v) vflag=1 # verbose
	   ;;
	\?) printf "Usage: %s -c CONFIG_FILE [-v]\n" $(basename $0) >&2
		exit 2
		;;
	esac
done

shift $(($OPTIND - 1))
PIDS=$*

if [[ -z $CONFIG_FILE || ! -f $CONFIG_FILE ]] 
then
	printf "Error: Config file %s not found. Please specify the absolute path to the config file.\n" $CONFIG_FILE >&2
	exit 2
fi
source $CONFIG_FILE

if [[ "$vflag" == 1 ]]
then
	grep -v '###' $CONFIG_FILE   ### FIX ME
	printf "%s\n" "`RSCRIPT_BIN | head -n1`"
fi

if [ ! -d $PIPELINE_DIR ]
then
	printf "Analysis pipeline directory not found. Exiting...\n" >&2
	exit 2
fi

if [[ ! -d $CLUSTER_EO ]]; then mkdir $CLUSTER_EO; fi

if [[ ! -d $OUTDIR ]]
then
	mkdir -p $OUTDIR
fi
if [[ ! -d $OUTDIR ]]
then
	echo "Not possible to create directory $OUTDIR."
	exit 2
fi

rscript=`qsub -o $CLUSTER_EO -j oe -M $EMAIL -N yapima -l $PBS_RESOURCES -v CONFIG_FILE=$CONFIG_FILE $PIPELINE_DIR/process450k.sh | cut -d '.' -f 1`
echo "yapima submitted, job ID $rscript"
