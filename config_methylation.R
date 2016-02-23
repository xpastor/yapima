#!/usr/bin/env Rscript-3.2.0

pipeline_dir <- '/home/pastor/ngs2/branches/pipelines/yapima'

### Input data ###
idat_dir <- file.path('/icgc/dkfzlsdf/analysis/hipo/hipo_054/infinium450k/idat')
sample.annotation <- '/icgc/dkfzlsdf/analysis/hipo/hipo_054/user_folders/pastor/methylation/sample_sheet_H054.csv' # Illumina's sample sheet with additional columns for sample annotation
blacklist <- '' # file with additional probes to be discarded from any analysis; '' to not to remove additional probes

### Output directory ###
wd <- '/icgc/dkfzlsdf/analysis/hipo/hipo_054/user_folders/pastor/test_RnBeads'

### Steps ###
batchCorrection <- F
surrogateCorrection <- T
runQC <- F
probeSelection <- F
diffMeth <- F

### Params ###
backgroundCorrection <- T # Noob correction from the methylumi package
normalization <- F # SWAN correction from the minfi package
removeEuropeanSNPs <- T
batch.vars <- '' # comma separated list of batch variables present in sample.annotation
varianceProportion <- 0.01 # proportion of variance that surrogate variables should explain

seed <- 11
ncores <- 4
