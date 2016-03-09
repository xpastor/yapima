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
echo -e "wd <- '$OUTDIR'"
echo -e "pipeline_dir <- '$PIPELINE_DIR'"
echo -e "seed <- $SEED"
if [[ -n $BLACKLIST ]]
then
	echo -e "blacklist <- '$BLACKLIST'"
fi
if isOn $RUN_BATCH_CORRECTION
then
	echo -e "batch.vars <- '$BATCH_VARS'"
fi

if isOn $RUN_DIFFERENTIAL_METHYLATION
then
	gawk '/variables #/,/#o#/' $PIPELINE_DIR/qc_functions.R
fi

gawk '/# Create output/,/#o#/' $PIPELINE_DIR/run_process450k.R
gawk '/# Define/,/#o#/' $PIPELINE_DIR/run_process450k.R

gawk '/# Load libraries/ || /# Reading in data #/ || /# Remove ambiguous/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if [[ -n $BLACKLIST ]]
then
	gawk '/# Load blacklist/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R | sed 's/^\s*//g'
fi

gawk '/# Exclude probes/ || /# Output raw tables #/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if isOn $CORRECT_BACKGROUND
then
	gawk '/# Remove background/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R
else
	gawk '/# Produce raw/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R
fi

if isOn $NORMALIZE
then
	gawk '/# Normalization/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R
fi

gawk '/# Extract/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if isOn $REMOVE_EUROPEAN_SNPS
then
	gawk '/# Process SNPs/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R
fi

gawk '/# Mask probes/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if isOn $RUN_BATCH_CORRECTION
then
	echo -e "\n#### Batch correction ####"
	cat $PIPELINE_DIR/batch_correction.R
fi

if isOn $RUN_PROBE_SELECTION
then
	gawk '/# Probe selection/,/#o#/' $PIPELINE_DIR/probe_selection.R
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
fi
