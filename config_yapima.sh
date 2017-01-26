#!/usr/bin/bash

NCORES=1
PIPELINE_DIR=$HOME/pipelines/devel/yapima

### Cluster-related parameters ###
EMAIL=x.pastorhostench@dkfz-heidelberg.de
PBS_RESOURCES="walltime=3:00:00,nodes=1:ppn=$NCORES,mem=10g"
#PBS_RESOURCES="walltime=10:00:00,nodes=1:ppn=$NCORES,mem=121g -q highmem" # for jobs with more than 120 gb
CLUSTER_EO=/ibios/co02/xavier/eo/yapima

### Binaries ###
RSCRIPT_BIN=/ibios/tbi_cluster/13.1/x86_64/R/R-3.3.1/bin/Rscript

### Steps ###
# The variables must take an R boolean value (T, TRUE, F or FALSE) #
RUN_PROBE_SELECTION=F
RUN_CNV=F
RUN_DIFFERENTIAL_METHYLATION=F

### Input data ###
IDAT_DIR=/icgc/dkfzlsdf/analysis/hipo/hipo_054/infinium450k/idat
SAMPLE_ANNOTATION=/icgc/dkfzlsdf/analysis/hipo/hipo_054/user_folders/pastor/methylation/sample_sheet_H054.csv # Illumina's sample sheet with additional columns for sample annotation
BLACKLIST='' # file with additional probes to be discarded from any analysis; '' to not to remove additional probes

### Output directory ###
OUTDIR=/icgc/dkfzlsdf/analysis/hipo/hipo_054/user_folders/pastor/test_yapima_filters_noob_swan

### Params ###
REMOVE_EUROPEAN_SNPS=T # T or F; remove SNPs that may be present in at least 1 sample
USE_PREDICTED_SEX=T
BATCH_VARS='' # comma separated list of batch variables present in $SAMPLE_ANNOTATION
SEED=11
