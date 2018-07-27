#!/bin/bash

source $CONFIG_FILE
#set -x

# Redefine the 0/1 switch to the R boolean values
RUN_CNV=`echo $RUN_CNV | tr 01 FT`
RUN_DIFFERENTIAL_METHYLATION=`echo $RUN_DIFFERENTIAL_METHYLATION | tr 01 FT`
REMOVE_SNPS=`echo $REMOVE_SNPS | tr 01 FT`
USE_PREDICTED_SEX=`echo $USE_PREDICTED_SEX | tr 01 FT`

# Export environment variables for R execution
export NCORES
export PIPELINE_DIR
export RUN_CNV
export RUN_DIFFERENTIAL_METHYLATION
export IDAT_DIR
export SAMPLE_ANNOTATION
export BLACKLIST
export OUTDIR
export REMOVE_SNPS
export POPULATION
export USE_PREDICTED_SEX
export BATCH_VARS
export SEED

Rscript $PIPELINE_DIR/yapima.R $PIPELINE_DIR/config_yapima.R

if [[ $? == 0 ]]
then
	cp $CONFIG_FILE $OUTDIR
	$PIPELINE_DIR/script_analysis.sh > $OUTDIR/script.R
fi
