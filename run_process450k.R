#!/usr/bin/env Rscript

params <- commandArgs(T)

if (length(params) == 0) {
	stop("Error:\n\tYou need to specify a configuration file.")
} else if (length(params) > 1) {
	message("Warning:\n\tYou specified too many parameters, only one is allowed. The first one will be taken as the configuration file.")
}

config <- params[1]
if (!file.exists(config)) {
	stop(paste0("\n\tThe config file cannot be accessed:\n\t", config))
}

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
if (! 'probe_selection.R' %in% pipeline_scripts) {
	stop(paste0("\n\tThe script 'probe_selection.R' is not present in ", pipeline_dir))
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

if (!is.logical(probeSelection)) {
	stop("\n\t'probeSelection' must be a valid R boolean: T, F, TRUE or FALSE.")
}

if (!is.logical(diffMeth)) {
	stop("\n\t'diffMeth' must be a valid R boolean: T, F, TRUE or FALSE.")
}

if (!is.logical(removeEuropeanSNPs)) {
	stop("\n\t'removeEuropeanSNPs' must be a valid R boolean: T, F, TRUE or FALSE.")
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
header <- readLines(sample.annotation, n=1)
header <- unlist(strsplit(header, ','))
illumina.vars <- c('Sample_Name', 'Sample_Well', 'Sample_Plate', 'Sample_Group', 'Pool_ID', 'Sentrix_ID', 'Sentrix_Position')
interest.vars <- header[! header %in% illumina.vars]

# Produce vector with batch variables
batch.vars <- unlist(strsplit(batch.vars, ','))
# Remove batch variables from list of variables for analysis
interest.vars <- interest.vars[! interest.vars %in% batch.vars]
#o#

methods <- file.path(wd, 'methods.txt')
citations.txt <- file.path(wd, 'citations.txt')

library(tools)
source(file.path(pipeline_dir, 'methylation_preprocessing.R'))
source(file.path(pipeline_dir, 'methylation_qc.R'))
if (probeSelection) source(file.path(pipeline_dir, 'probe_selection.R'))
if (runCNV) source(file.path(pipeline_dir, 'methylation_CNV.R'))
if (diffMeth) {
	if (isEmpty(interest.vars)) {
	    message("There's no factor eligible for a differential methylation analysis and it will be skipped.")
	} else {
		source(file.path(pipeline_dir, 'differential_methylation.R'))
	}
}
source(file.path(pipeline_dir, 'methods.R'))
source(file.path(pipeline_dir, 'citations.R'))

### sessionInfo() ###
source('http://bioconductor.org/biocLite.R')
session.tex <- file.path(wd, 'session.tex')
write(paste0('\\documentclass{report}\n\\title{\'yapima\' sessionInfo}\n\n\\usepackage{hyperref}\n\n\\begin{document}\n\\section*{\\centerline{\'yapima\' sessionInfo}}\n\\center{\\today}\n\n\\begin{itemize}\\raggedright\n  \\item Bioconductor version ', biocVersion(), '\n\\end{itemize}'),file=session.tex)
write(toLatex(sessionInfo()), file=session.tex, append=T)
write('\n\\end{document}', file=session.tex, append=T)
setwd(wd)
texi2pdf(session.tex, clean=T)
