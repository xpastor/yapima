#!/usr/bin/env Rscript-3.2.0

pipeline_dir <- Sys.getenv("PIPELINE_DIR")

### Input data ###
idat_dir <- Sys.getenv("IDAT_DIR")
sample.annotation <- Sys.getenv("SAMPLE_ANNOTATION") # Illumina's sample sheet with additional columns for sample annotation
blacklist <- Sys.getenv("BLACKLIST") # file with additional probes to be discarded from any analysis; '' to not to remove additional probes

### Output directory ###
wd <- Sys.getenv("OUTDIR")

### Steps ###
batchCorrection <- as.logical(Sys.getenv("RUN_BATCH_CORRECTION"))
surrogateCorrection <- as.logical(Sys.getenv("RUN_SURROGATE_CORRECTION"))
runQC <- as.logical(Sys.getenv("RUN_QC"))
probeSelection <- as.logical(Sys.getenv("RUN_PROBE_SELECTION"))
diffMeth <- as.logical(Sys.getenv("RUN_DIFFERENTIAL_METHYLATION"))

### Params ###
backgroundCorrection <- as.logical(Sys.getenv("CORRECT_BACKGROUND")) # Noob correction from the methylumi package
normalization <- as.logical(Sys.getenv("NORMALIZE")) # SWAN correction from the minfi package
removeEuropeanSNPs <- as.logical(Sys.getenv("REMOVE_EUROPEAN_SNPS"))
batch.vars <- Sys.getenv("BATCH_VARS") # comma separated list of batch variables present in sample.annotation
varianceProportion <- Sys.getenv("VARIANCE_PROPORTION") # proportion of variance that surrogate variables should explain

seed <- Sys.getenv("SEED")
ncores <- Sys.getenv("NCORES")
