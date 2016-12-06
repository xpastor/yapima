text <- "The Illumina Infinium methylation array analysis was done using the \'yapima\' pipeline developed by Xavier Pastor (x.pastorhostench@dkfz-heidelberg.de) at DKFZ.\n\n"

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


if (removeEuropeanSNPs) {
	text <- paste0(text,
		' Probes with the single base extension (SBE) position overlaping SNPs with allele frequency higher than ', round(1/nrow(targets), digits=2),' in european populations ', flag.ref, ' were also flagged (n=', length(european.snps), ').')
}

text <- paste(text,
	'The intensities were adjusted with the \'ENmix\' background correction available in the \'ENmix\' package, using the out-of-band probes to estimate the background and enabling dye bias correction and normalized applying the \'SWAN\' normalization mehtod (', .short_citation('minfi', 2), ').')
#norm.meth <- preprocessENmix(rgset, bpParaEst='oob', dyeCorr=T, QCinfo=NULL, exQCsample=F, exQCcpg=F, exSample=NULL, exCpG=NULL, nCores=ncores)

text <- paste(text,
	'Measures with a detection P-value higher than 0.01, as estimated by \'minfi\' and tipically of low quality, were masked.')

### ComBat batch effect correction ###
if (batchCorrection) {
	text <- paste0(text, '\n',
		'The M-values were corrected for possible batch effects for the following variables using the \'ComBat\' (Johnson, 2007) implementation in the \'sva\' package ', .cite_package('sva'), ' from bioconductor: ', paste(batch.vars, collapse=', '), '.')
}

### Probe selection ###
if (probeSelection) {
	text <- paste0(text, '\n',
 		'To find the most stable unsupervised clustering of the samples, 12 sets of probes were analysed. Initially, the flagged probes were removed and from the remaining, the ones with the highest variation in their beta-values were chosen for the 12 sets, from 1000 to 12000 probes in steps of 1000 probes. The bootstrap clustering approach implemented by the \'pvclust\' package ', .cite_package('pvclust'), ' from bioconductor was applied, with 10000 iterations for each using the euclidean distance measure. Finally, to find the most stable set of probes a score was defined as the sum of the pvalue of the  top ', round((nrow(targets)-1)*0.25), ' edges multiplied by the height of the edge, and the set of probes with the highest score was chosen.')
}

### Differential methylation analysis ###
if (diffMeth) {
	if (! isEmpty(interest.vars)) {
		limmaURL <- 'http://bioconductor.org/packages/3.1/bioc/vignettes/limma/inst/doc/usersguide.pdf'
		text <- paste0(text, '\n',
			'Differential methylation analysis was done on the M-values using the biconductor \'limma\' package ', .cite_package('limma'), ' and multiple testing correction was applied (Benjamini, 1995) (protocol described in ', limmaURL, ').')
	}
}

if (batchCorrection | diffMeth) {
	text <- paste0(text, '\n\n',
		'When necessary, the beta-values were derived from the M-values as described by Du et. al. (Du, 2010).')
}

write(text, file=file.path(wd, 'methods.txt'))
