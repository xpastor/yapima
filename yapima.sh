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

if [[ -z $CONFIG_FILE || ! -f $CONFIG_FILE ]] 
then
	printf "Error: Config file %s not found. Please specify the absolute path to the config file.\n" $CONFIG_FILE >&2
	exit 2
fi

set -a
source $CONFIG_FILE
set +a

if [[ ! -z $CONDA_PREFIX && ! -z $CONDA_ENV ]]
then
  . activate $CONDA_ENV
fi

PIPELINE_DIR=$(dirname `readlink -f $0`)

# Redefine the 0/1 switch to the R boolean values
RUN_CNV=`echo $RUN_CNV | tr 01 FT`
RUN_DIFFERENTIAL_METHYLATION=`echo $RUN_DIFFERENTIAL_METHYLATION | tr 01 FT`
REMOVE_SNPS=`echo $REMOVE_SNPS | tr 01 FT`
USE_PREDICTED_SEX=`echo $USE_PREDICTED_SEX | tr 01 FT`

Rscript $PIPELINE_DIR/yapima.R $PIPELINE_DIR/config_yapima.R

if [[ $? == 0 ]]
then
	cp $CONFIG_FILE $OUTDIR
	cp $PIPELINE_DIR/*rda $OUTDIR
	$PIPELINE_DIR/script_analysis.sh > $OUTDIR/script.R
fi
