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
echo -e "ncores <- $NCORES"
echo -e "usePredictedSex <- $USE_PREDICTED_SEX"

if [[ -n $BLACKLIST ]]
then
	echo -e "blacklist <- '$BLACKLIST'"
fi
echo -e "batch.vars <- '$BATCH_VARS'"

if isOn $RUN_PROBE_SELECTION
then
	gawk '/Clustering #/,/#o#/' $PIPELINE_DIR/functions.R
fi

if isOn $RUN_DIFFERENTIAL_METHYLATION
then
	gawk '/BED-like/,/#o#/' $PIPELINE_DIR/functions.R
	cat $PIPELINE_DIR/extractCoords.R
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

if isOn $REMOVE_SNPS
then
	gawk '/# Process SNPs/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R
fi

gawk '/# Load crossreactive/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if [[ -n $BLACKLIST ]]
then
	gawk '/# Load blacklist/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R | sed 's/^\s*//g'
fi

gawk '/# Output raw/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if isOn $RUN_PROBE_SELECTION
then
	gawk '/# Probe selection/ || /# Output/,/#o#/' $PIPELINE_DIR/probe_selection.R
fi

if isOn $RUN_DIFFERENTIAL_METHYLATION
then
	gawk '/# Differential methylation/,/#o#/' $PIPELINE_DIR/differential_methylation.R
	gawk '/# Limma/,/#o#/' $PIPELINE_DIR/differential_methylation.R
	gawk '/# DMR/, /#o#/' $PIPELINE_DIR/differential_methylation.R
	echo -e "\n\t}\n}"
fi
