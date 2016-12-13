#!/bin/bash

source $CONFIG_FILE
#set -x

# Redefine the 0/1 switch to the R boolean values
RUN_BATCH_CORRECTION=`echo $RUN_BATCH_CORRECTION | tr 01 FT`
RUN_CNV=`echo $RUN_CNV | tr 01 FT`
RUN_PROBE_SELECTION=`echo $RUN_PROBE_SELECTION | tr 01 FT`
RUN_DIFFERENTIAL_METHYLATION=`echo $RUN_DIFFERENTIAL_METHYLATION | tr 01 FT`
CORRECT_BACKGROUND=`echo $CORRECT_BACKGROUND | tr 01 FT`
NORMALIZE=`echo $NORMALIZE | tr 01 FT`
REMOVE_EUROPEAN_SNPS=`echo $REMOVE_EUROPEAN_SNPS | tr 01 FT`

# Export environment variables for R execution
export NCORES
export PIPELINE_DIR
export RSCRIPT_BIN
export RUN_BATCH_CORRECTION
export RUN_CNV
export RUN_PROBE_SELECTION
export RUN_DIFFERENTIAL_METHYLATION
export IDAT_DIR
export SAMPLE_ANNOTATION
export BLACKLIST
export OUTDIR
export REMOVE_EUROPEAN_SNPS
export BATCH_VARS
export SEED

$RSCRIPT_BIN $PIPELINE_DIR/run_process450k.R $PIPELINE_DIR/config_yapima.R

if [[ $? == 0 ]]
then
	cp $CONFIG_FILE $OUTDIR
	$PIPELINE_DIR/script_analysis.sh > $OUTDIR/script.R
	echo -e "\n#yapima `cd $PIPELINE_DIR/;git tag | tail -1`" >> $OUTDIR/methods.txt
	echo -e "\n#commit version: `tail -1 $PIPELINE_DIR/.git/logs/HEAD | cut -d' ' -f2`" >> $OUTDIR/methods.txt
	echo -e "\n#commit version: `tail -1 $PIPELINE_DIR/.git/logs/HEAD | cut -d' ' -f2 | cut -c-6`" >> $OUTDIR/methods.txt
fi
