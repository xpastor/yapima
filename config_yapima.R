#!/usr/bin/env Rscript-3.2.0

### Input data ###
idat_dir <- Sys.getenv("IDAT_DIR")
sample.annotation <- Sys.getenv("SAMPLE_ANNOTATION") # Illumina's sample sheet with additional columns for sample annotation
blacklist <- Sys.getenv("BLACKLIST") # file with additional probes to be discarded from any analysis; '' to not to remove additional probes

### Output directory ###
wd <- Sys.getenv("OUTDIR")

### Steps ###
runCNV <- as.logical(Sys.getenv("RUN_CNV"))
diffMeth <- as.logical(Sys.getenv("RUN_DIFFERENTIAL_METHYLATION"))

### Params ###
removeSNPs <- as.logical(Sys.getenv("REMOVE_SNPS"))
population <- Sys.getenv("POPULATION") # '' for global populations, 'AFR', 'AMR', 'ASN' only in 450k, 'EAS' only in EPIC, 'EUR', 'SAS' only in EPIC
usePredictedSex <- as.logical(Sys.getenv("USE_PREDICTED_SEX"))
#surrogateCorrection <- as.logical(Sys.getenv("SURROGATE_CORRECTION"))
batch.vars <- Sys.getenv("BATCH_VARS") # comma separated list of batch variables present in sample.annotation
#varianceProportion <- Sys.getenv("VARIANCE_PROPORTION") # proportion of variance that surrogate variables should explain

seed <- as.integer(Sys.getenv("SEED"))
ncores <- as.integer(Sys.getenv("NCORES"))
