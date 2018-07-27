#!/usr/bin/bash

NCORES=1
PIPELINE_DIR=$HOME/pipelines/devel/yapima

### Cluster-related parameters ###
EMAIL=x.pastorhostench@dkfz-heidelberg.de
PBS_RESOURCES="walltime=3:00:00,nodes=1:ppn=$NCORES,mem=10g"
#PBS_RESOURCES="walltime=10:00:00,nodes=1:ppn=$NCORES,mem=121g -q highmem" # for jobs with more than 120 gb
CLUSTER_EO=/ibios/co02/xavier/eo/yapima

### Modules ###
R_MOD=R/3.5.0

### Steps ###
# The variables must take an R boolean value (T, TRUE, F or FALSE) #
RUN_CNV=F
RUN_DIFFERENTIAL_METHYLATION=F

### Input data ###
IDAT_DIR=/icgc/dkfzlsdf/analysis/hipo/hipo_054/infinium450k/idat
SAMPLE_ANNOTATION=/icgc/dkfzlsdf/analysis/hipo/hipo_054/user_folders/pastor/methylation/sample_sheet_H054.csv # Illumina's sample sheet with additional columns for sample annotation
BLACKLIST='' # file with additional probes to be discarded from any analysis; '' to not to remove additional probes

### Output directory ###
OUTDIR=/icgc/dkfzlsdf/analysis/hipo/hipo_054/user_folders/pastor/test_yapima_filters_noob_swan

### Params ###
REMOVE_SNPS=T # T or F; remove SNPs that may be present in at least 1 sample
POPULATION='EUR' # '' for global populations, 'AFR', 'AMR', 'ASN' only in 450k, 'EAS' only in EPIC, 'EUR', 'SAS' only in EPIC
USE_PREDICTED_SEX=T
BATCH_VARS='' # comma separated list of batch variables present in $SAMPLE_ANNOTATION
SEED=11
