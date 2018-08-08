#!/bin/bash

isOn () {
	status=0
	case $1 in
		0 | 'F' | 'FALSE' ) status=1;;
		1 | 'T' | 'TRUE' ) status=0;;
	esac
	return $status
}

PIPELINE_DIR=$(dirname $0)

echo -e "#!/usr/bin/env Rscript\n"
echo -e "# Variables initialization"
echo -e "idat_dir <- '$IDAT_DIR'"
echo -e "sample.annotation <- '$SAMPLE_ANNOTATION'"
echo -e "wd <- '${OUTDIR}_rep'"
#echo -e "pipeline_dir <- '$PIPELINE_DIR'"
echo -e "seed <- $SEED"
echo -e "ncores <- $NCORES"
echo -e "usePredictedSex <- $USE_PREDICTED_SEX"

if [[ -n $BLACKLIST ]]
then
	echo -e "blacklist <- '$BLACKLIST'"
fi
echo -e "batch.vars <- '$BATCH_VARS'"

if isOn $RUN_DIFFERENTIAL_METHYLATION
then
	gawk '/BED-like/,/#o#/' $PIPELINE_DIR/functions.R
	cat $PIPELINE_DIR/extractCoords.R
fi

gawk '/# Create output/,/#o#/' $PIPELINE_DIR/yapima.R
gawk '/# Define/,/#o#/' $PIPELINE_DIR/yapima.R
gawk '/# Extract/,/#o#/' $PIPELINE_DIR/yapima.R
if [[ -n $BATCH_VARS ]]
then
	gawk '/# Produce/,/#o#/' $PIPELINE_DIR/yapima.R
fi

gawk '/# Load libraries/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if isOn $REMOVE_SNPS
then
	echo -e "load(file.path('$OUTDIR', 'polymorphic.rda'))"
	echo -e "population <- '$POPULATION'"
	gawk '/# Process SNPs/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R
fi

echo -e "load(file.path('$OUTDIR', 'crossreactive.rda'))"
gawk '/# Process crossreactive/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if [[ -n $BLACKLIST ]]
then
	gawk '/# Load blacklist/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R | sed 's/^\s*//g'
fi

gawk '/# Output raw/,/#o#/' $PIPELINE_DIR/methylation_preprocessing.R

if isOn $RUN_DIFFERENTIAL_METHYLATION
then
	gawk '/# Differential methylation/,/#o#/' $PIPELINE_DIR/differential_methylation.R
	gawk '/# Limma/,/#o#/' $PIPELINE_DIR/differential_methylation.R
	gawk '/# GO analysis/,/#o#/' $PIPELINE_DIR/differential_methylation.R
	echo -e "\n\t}\n}" 
fi
