text <- paste0('---
title: "yapima - Methods, Session & References"
author: Xavier Pastor

date: "`r Sys.Date()`"
output:
  word_document:
    output_file: ',
file.path(wd, 'methods_references.doc'),
'\n    toc: true
    number_sections: true
link-citations: true
bibliography: ',
file.path(pipeline_dir, 'REFERENCES.bib'),
'\n---

# Methods

## Preprocessing
')

flag.ref <- NULL
if (array.type == 'IlluminaHumanMethylation450k') {
	flag.ref <- '[@450k]'
} else if (array.type == 'IlluminaHumanMethylationEPIC') {
	flag.ref <- '[@epic]'
}

text <- paste0(text,
	'The array data were read into the R environment (R version ',
	R.version$major, '.', R.version$minor, ' \'', R.version$nickname, '\', ', R.version$year, '-', R.version$month, '-', R.version$day, ') using the *minfi* package [@minfi].')

text <- paste(text, 'The crossreactive probes', flag.ref, 'were flagged.')

 if (blacklist != '') {
 	text <- paste0(text,
		'A list of blacklisted probes was also flagged (describe origin of blacklisted probes).')
 }

if (removeSNPs) {
	if (is.null(population) | population == '') population <- 'all'
	text <- paste0(text,
		' Probes with the single base extension (SBE) position overlaping SNPs with allele frequency higher than 0.005 in ', population,' populations ', flag.ref, ' were also flagged (n=', num_snps, ').')
}

text <- paste(text,
	'Measures with a detection P-value higher than 0.01, as estimated by *minfi* and tipically of low quality, were masked. Probes with low quality in at least 50% of the samples were also flagged.')

text <- c(text, '\n## Processing',
	'The intensities were adjusted with the *ENmix* package [@ENmix] using combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type, and enabling *RELIC* dye bias correction [@RELIC]. Methylation values were normalized applying the *RCP* normalization mehtod also available in *ENmix* [@RCP].')

if (diffMeth) {
	if (! isEmpty(interest.vars)) {
		text <- c(text, '\n## Differential Methylation Analysis\n')
		limmaURL <- 'http://bioconductor.org/packages/3.1/bioc/vignettes/limma/inst/doc/usersguide.pdf'
		dmrcateURL <- 'https://www.bioconductor.org/packages/release/bioc/vignettes/DMRcate/inst/doc/DMRcate.pdf'
		text <- c(text,
			'The analysis of differentially methylated positions (DMP) was done on the M-values of the reliable probes (i.e. flag=0) using the biconductor *limma* package [@limma] and multiple testing correction was applied [@fdr] (protocol described in ', limmaURL, ').')
		text <- c(text,
			'The detection of differentially methylated regions (DMR) was done using the bioconductor *DMRcate* package [@dmr]. For two groups comparisons, the t statistics from the DMP analysis were used and the beta log fold change was computed running the standard *limma* workflow on the beta values. For comparisons with more than two groups the squared F statistics from the DMP analysis were provided and the beta log fold change was set to 0. All the other parameters were left as default.')
		text <- c(text, '\n',
			'When necessary, the beta-values were derived from the M-values as described by [@mval].')
	}
}

source('http://bioconductor.org/biocLite.R')
repo <- repository(pipeline_dir)
commit <- revparse_single(repo, "HEAD")
repo_status <- do.call(c, status(repo))
sha <- ifelse(any(grepl('\\.modified', names(repo_status))), '', commit@sha)
if (sha_ini != sha) sha <- ''

session <- c('```{r echo=F}',
"message(paste0('Bioconductor version ', biocVersion()))",
'sessionInfo()',
"message(paste('commit: ', sha))",
'```')
text <- c(text, '\n# Session Info', session)

text <- c(text, '\n# References')

rmd <- file.path(wd, 'methods_session_references.Rmd')
writeLines(text, rmd)
#rmarkdown::render(file.path(rmd), rmarkdown::word_document(), output_file='methods_session_references.doc', output_dir=wd)
rmarkdown::render(file.path(rmd), rmarkdown::word_document())
