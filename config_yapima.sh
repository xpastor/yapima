#!/usr/bin/bash

NCORES=1
PIPELINE_DIR=$HOME/pipelines/yapima

### Cluster-related parameters ###
EMAIL=x.pastorhostench@dkfz-heidelberg.de
PBS_RESOURCES="walltime=1:00:00,nodes=1:ppn=$NCORES,mem=6g"
CLUSTER_EO=/ibios/co02/xavier/eo/yapima

### Binaries ###
RSCRIPT_BIN=/ibios/tbi_cluster/13.1/x86_64/R/R-3.2.0/bin/Rscript

### Steps ###
# The variables must take an R boolean value (T, TRUE, F or FALSE) #
RUN_BATCH_CORRECTION=F
RUN_CNV=F
RUN_PROBE_SELECTION=F
RUN_DIFFERENTIAL_METHYLATION=F

### Input data ###
IDAT_DIR=/icgc/dkfzlsdf/analysis/hipo/hipo_054/infinium450k/idat
SAMPLE_ANNOTATION=/icgc/dkfzlsdf/analysis/hipo/hipo_054/user_folders/pastor/methylation/sample_sheet_H054.csv # Illumina's sample sheet with additional columns for sample annotation
BLACKLIST='' # file with additional probes to be discarded from any analysis; '' to not to remove additional probes

### Output directory ###
OUTDIR=/icgc/dkfzlsdf/analysis/hipo/hipo_054/user_folders/pastor/test_pipeline

### Params ###
CORRECT_BACKGROUND=T # T or F; Noob correction from the methylumi package
NORMALIZE=F # T or F; SWAN correction from the minfi package
REMOVE_EUROPEAN_SNPS=T # T or F; remove SNPs that may be present in at least 1 sample
SURROGATE_CORRECTION=F
BATCH_VARS='' # comma separated list of batch variables present in $SAMPLE_ANNOTATION
VARIANCE_PROPORTION=0.01 # proportion of variance that surrogate variables should explain
SEED=11
