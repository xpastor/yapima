#!/usr/bin/env Rscript

params <- commandArgs(F)
trail_params <- grep('--args', params)

if (length(params) == trail_params) {
	stop("Error:\n\tYou need to specify a configuration file.")
} else if (length(params) > trail_params + 1) {
	message("Warning:\n\tYou specified too many parameters, only one is allowed. The first one will be taken as the configuration file.")
}

config <- params[trail_params + 1]
if (!file.exists(config)) {
	stop(paste0("\n\tThe config file cannot be accessed:\n\t", config))
}
pipeline_dir <- gsub('--file=', '', grep('--file', params, value=T))
pipeline_dir <- dirname(pipeline_dir)

source(config)

# Create output directory
#now <- strftime(Sys.time(), '%y%m%d%H%M%S')
#wd <- file.path(wd, paste0('yapima_', now))
if (! dir.exists(wd)) {
	success <- try(dir.create(wd, recursive=T, mode='0770'), T)
	if (! success) {
		stop(paste0("\n\tThe output directory could not be created.\n\t", wd))
	}
} else {
	stop(paste0('The directory \'', wd, '\' already exists. Specify a non existing directory.'))
}
#message(paste0("The results will be stored in:\n\t",wd))
#o#

if (! dir.exists(idat_dir)) {
	stop(paste0("\n\tThe directory with the IDAT files can not be accessed.\n\t", idat_dir))
}

if (! dir.exists(pipeline_dir)) {
	stop(paste0("\n\tThe directory with the pipeline scripts can not be accessed.\n\t", pipeline_dir))
}

pipeline_scripts <- list.files(pipeline_dir)
if (! 'methylation_preprocessing.R' %in% pipeline_scripts) {
	stop(paste0("\n\tThe script 'methylation_preprocessing.R' is not present in ", pipeline_dir))
}
if (! 'methylation_qc.R' %in% pipeline_scripts) {
	stop(paste0("\n\tThe script 'methylation_qc.R' is not present in ", pipeline_dir))
}
if ('functions.R' %in% pipeline_scripts) {
	source(file.path(pipeline_dir, 'functions.R'))
} else {
	stop(paste0("\n\tThe script 'functions.R' is not present in ", pipeline_dir))
}

if (! file.exists(sample.annotation) | file.access(sample.annotation) == -1) {
	stop(paste0("\n\tThe sample sheet could not be accessed.\n\t", sample.annotation))
}

if ( blacklist != '' & (! file.exists(blacklist) | file.access(blacklist) == -1)) {
	stop(paste0("\n\tThe file with blacklisted probes could not be accessed.\n\t", blacklist))
}

if (!is.logical(runCNV)) {
	stop("\n\t'runQC' must be a valid R boolean: T, F, TRUE or FALSE.")
}

if (!is.logical(diffMeth)) {
	stop("\n\t'diffMeth' must be a valid R boolean: T, F, TRUE or FALSE.")
}

if (!is.logical(removeSNPs)) {
	stop("\n\t'removeSNPs' must be a valid R boolean: T, F, TRUE or FALSE.")
}
if (removeSNPs) {
	if (is.null(population)) population <- ''
	populations <- c('', 'AFR', 'AMR', 'ASN', 'EAS', 'EUR', 'SAS')
	if (!(population %in% populations)) {
		stop(paste0(population, ' is not a valid population. Specify one of ', paste0(populations, collapse=', '), '.'))
	}
}
	
if (is.na(ncores)) {
	stop("'ncores' must be integer.")
}

if (is.na(seed)) {
	stop("'seed' must be integer.")
}

qcdir <- file.path(wd, 'qc')
dir.create(qcdir, recursive=T)

# Define the seed
set.seed(seed)
#o#

# Extract variables of interest from sample sheet
library(minfi)
targets <- read.metharray.sheet(dirname(sample.annotation), paste0('^', basename(sample.annotation), '$'), recursive=F)
header <- colnames(targets)
#header <- readLines(sample.annotation, n=1)
#header <- unlist(strsplit(header, ','))
illumina.vars <- c('Sample_Name', 'Sample_Well', 'Sample_Plate', 'Sentrix_ID', 'Sentrix_Position', 'Slide', 'Array', 'Basename')

interest.vars <- header[! header %in% illumina.vars]
if (length(interest.vars) != 0) interest.vars <- interest.vars[!apply(is.na(targets[,interest.vars,drop=F]), 2, all)]

# Produce vector with batch variables
batch.vars <- unlist(strsplit(batch.vars, ','))
batch.vars <- gsub('Sentrix_ID', 'Slide', batch.vars)
batch.vars <- gsub('Sentrix_Position', 'Array', batch.vars)
# Remove batch variables from list of variables for analysis
interest.vars <- interest.vars[! interest.vars %in% batch.vars]
#o#

methods <- file.path(wd, 'methods.txt')
citations.txt <- file.path(wd, 'citations.txt')

library(tools)
source(file.path(pipeline_dir, 'methylation_preprocessing.R'))
source(file.path(pipeline_dir, 'methylation_qc.R'))
if (runCNV) source(file.path(pipeline_dir, 'methylation_CNV.R'))
if (diffMeth) {
	if (isEmpty(interest.vars)) {
	    message("There's no factor eligible for a differential methylation analysis and it will be skipped.")
	} else {
		source(file.path(pipeline_dir, 'differential_methylation.R'))
	}
}
source(file.path(pipeline_dir, 'document_run.R'))
