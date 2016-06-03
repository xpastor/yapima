#!/bin/bash

isOn () {
	status=0
	case $1 in
		0 | 'F' | 'FALSE' ) status=1;;
		1 | 'T' | 'TRUE' ) status=0;;
	esac
	return $status
}

echo -e "#!/usr/bin/env $RSCRIPT_BIN\n"
echo -e "# Variables initialization"
echo -e "idat_dir <- '$IDAT_DIR'"
echo -e "sample.annotation <- '$SAMPLE_ANNOTATION'"
echo -e "wd <- '${OUTDIR}_rep'"
echo -e "pipeline_dir <- '$PIPELINE_DIR'"
echo -e "seed <- $SEED"
echo -e "non_specific_cg <- '$OUTDIR/$(basename $NON_SPECIFIC_CG)'"
echo -e "non_specific_ch <- '$OUTDIR/$(basename $NON_SPECIFIC_CH)'"
if [[ -n $BLACKLIST ]]
then
	echo -e "blacklist <- '$BLACKLIST'"
fi
echo -e "batch.vars <- '$BATCH_VARS'"

if isOn $RUN_PROBE_SELECTION
then
	echo -e "ncores <- $NCORES"
	gawk '/Clustering #/,/#o#/' $PIPELINE_DIR/functions.R
fi

#if isOn $RUN_DIFFERENTIAL_METHYLATION
#then
#	gawk '/variables #/,/#o#/' $PIPELINE_DIR/functions.R
#fi

gawk '/# Create output/,/#o#/' $PIPELINE_DIR/run_process450k.R
gawk '/# Define/,/#o#/' $PIPELINE_DIR/run_process450k.R
gawk '/# Extract/,/#o#/' $PIPELINE_DIR/run_process450k.R
if [[ -n $BATCH_VARS ]]
then
	gawk '/# Produce/,/#o#/' $PIPELINE_DIR/run_process450k.R
fi

gawk '/# Load libraries/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if isOn $CORRECT_BACKGROUND
then
	gawk '/^# Remove background/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R
else
	gawk '/# Produce raw/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R
fi

if isOn $NORMALIZE
then
	gawk '/# Normalization/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R
fi

gawk '/# Filter/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if [[ -n $BLACKLIST ]]
then
	gawk '/# Load blacklist/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R | sed 's/^\s*//g'
fi

if isOn $REMOVE_EUROPEAN_SNPS
then
	echo -e "polymorphic <- '$OUTDIR/$(basename $POLYMORPHIC)'"
	gawk '/# Process SNPs/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R
fi

gawk '/# Exclude probes/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if isOn $RUN_BATCH_CORRECTION
then
	echo -e "\n#### Batch correction ####"
	cat $PIPELINE_DIR/batch_correction.R
fi

if isOn $RUN_PROBE_SELECTION
then
	gawk '/# Probe selection/ || /# Output/,/#o#/' $PIPELINE_DIR/probe_selection.R
fi

if isOn $RUN_DIFFERENTIAL_METHYLATION
then
	gawk '/# Differential methylation analysis #/,/#o#/' $PIPELINE_DIR/differential_methylation.R	
	if isOn $SURROGATE_CORRECTION
	then
		gawk '/# Selection/,/#o#/' $PIPELINE_DIR/differential_methylation.R
	else
		gawk '/# Limma/,/#o#/' $PIPELINE_DIR/differential_methylation.R
	fi
	echo -e "\t}\n}"
fi