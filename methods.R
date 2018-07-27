text <- "The Illumina Infinium methylation array analysis was done using the \'yapima\' pipeline developed by Xavier Pastor (xavier.pastor@compbio-dev) at DKFZ.\n\n"

.short_citation <- function(package, idx=1)
{
	cite <- citation(package)[idx]
	first <- cite$author[1]
	first <- gsub(' <.*', '', first)
	short_cite <- gsub('.* ', '', first)
	year <- cite$year
	if (! is.null(year)) short_cite <- paste(short_cite, year, sep=', ')
	return(short_cite)
}

.cite_package <- function(package, idx=1)
{
	short_cite <- .short_citation(package, idx)
	version <- paste0('v', packageVersion(package))
	cite_package <- paste0('(', short_cite, '; ', version, ')')
	return(cite_package)
}

### Preprocessing ###

flag.ref <- NULL
if (array.type == 'IlluminaHumanMethylation450k') {
	flag.ref <- '(Chen, 2013)'
} else if (array.type == 'IlluminaHumanMethylationEPIC') {
	flag.ref <- '(McCartney, 2016)'
}

text <- paste0(text,
	'The array data were read into the R environment (R version ',
	R.version$major, '.', R.version$minor, ' \'', R.version$nickname, '\', ', R.version$year, '-', R.version$month, '-', R.version$day, ') using the \'minfi\' package ', .cite_package('minfi'), '.')

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
	'The intensities were adjusted with the \'ENmix\' package (', .short_citation('ENmix', 1), ') using combined methylated and unmethylated intensities to estimate background distribution parameters separately for each color channel and each probe type, and enabling \'RELIC\' dye bias correction (', .short_citation('ENmix', 4), '). Methylation values were normalized applying the \'RCP\' normalization mehtod also available in \'ENmix\' (', .short_citation('ENmix', 2), ').')
#norm.meth <- preprocessENmix(rgset, bpParaEst='oob', dyeCorr=T, QCinfo=NULL, exQCsample=F, exQCcpg=F, exSample=NULL, exCpG=NULL, nCores=ncores)

text <- paste(text,
	'Measures with a detection P-value higher than 0.01, as estimated by \'minfi\' and tipically of low quality, were masked. Probes with low quality in at least 50% of the samples were also flagged.')

### Differential methylation analysis ###
if (diffMeth) {
	if (! isEmpty(interest.vars)) {
		limmaURL <- 'http://bioconductor.org/packages/3.1/bioc/vignettes/limma/inst/doc/usersguide.pdf'
		dmrcateURL <- 'https://www.bioconductor.org/packages/release/bioc/vignettes/DMRcate/inst/doc/DMRcate.pdf'
		text <- paste0(text, '\n',
			'The analysis of differentially methylated positions (DMP) was done on the M-values of the reliable probes (i.e. flag=0) using the biconductor \'limma\' package ', .cite_package('limma'), ' and multiple testing correction was applied (Benjamini, 1995) (protocol described in ', limmaURL, ').')
		text <- paste0(text, '\n',
			'The detection of differentially methylated regions (DMR) was done using the bioconductor \'DMRcate\' package ', .cite_package('DMRcate'), '. For two groups comparisons, the t statistics from the DMP analysis were used and the beta log fold change was computed running the standard \'limma\' workflow on the beta values. For comparisons with more than two groups the squared F statistics from the DMP analysis were provided and the beta log fold change was set to 0. All the other parameters were left as default.')
	}
	text <- paste0(text, '\n\n',
		'When necessary, the beta-values were derived from the M-values as described by Du et. al. (Du, 2010).')
}

write(text, file=file.path(wd, 'methods.txt'))
