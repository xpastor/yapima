# Load libraries
#library(minfi)
library(ENmix)
library(GenomicRanges)

# Adjust number of usable cores
if (Sys.getenv('PBS_NUM_PPN') != '') {
	ncores <- min(ncores, as.integer(Sys.getenv('PBS_NUM_PPN')))
} else {
	ncores <- min(ncores, detectCores()-1)
}
ncores <- max(1, ncores)

#### Reading in data ####
### Preparing targets data frame ###
targets$Basename <- paste(targets$Slide, targets$Array, sep='_')
row.names(targets) <- targets$Basename

not.unique <- colnames(targets)[apply(targets, 2, function(x) length(unique(x)) > 1)]
interest.vars <- interest.vars[interest.vars %in% not.unique]

### Reading methylation data ###
message("Reading in methylation files...")
rgset <- read.metharray.exp(idat_dir, targets, extended=T, recursive=T, force=TRUE)
sampleNames(rgset) <- pData(rgset)$Sample_Name
message("Data read.")

### Fetching array annotation ###
array.type <- annotation(rgset)['array'] # IlluminaHumanMethylation450k / IlluminaHumanMethylationEPIC
array.annot <- getAnnotation(rgset)
array.annot.gr <- GRanges(array.annot$chr, ranges=IRanges(array.annot$pos, array.annot$pos+1), strand=NULL, array.annot[,! colnames(array.annot) %in% c('chr', 'pos', 'strand')], seqinfo=Seqinfo(paste0('chr', c(1:22, 'X', 'Y'))))
annot.bed <- data.frame(chrom=array.annot$chr, chromStart=array.annot$pos-1, chromEnd=array.annot$pos, name=array.annot$Name, score=rep(0, nrow(array.annot)), strand=array.annot$strand, stringsAsFactors=F, row.names=array.annot$Name)
annot.bed <- annot.bed[with(annot.bed, order(chrom, chromStart)),]
colnames(annot.bed)[1] <- paste0('#', colnames(annot.bed)[1])

## Output annotation ##
write.table(array.annot, file.path(wd, 'annotation.txt'), sep="\t", quote=F, row.names=T)

### Filter data ###

## Flag low quality probes ##
lowQ <- detectionP(rgset) > 0.01
lowQ.probes <- rowSums(lowQ)/ncol(lowQ) >= 0.5
annot.bed$score[lowQ.probes] <- annot.bed$score[lowQ.probes] + 8

## Load crossreactive probes ##
crossreactive <- NULL
if (array.type == 'IlluminaHumanMethylation450k') {
#	source(file.path(pipeline_dir, 'xlsxToR.R'))
	source('https://gist.githubusercontent.com/schaunwheeler/5825002/raw/a7e1844d2abcb6c51f7d2479d5d9f64a473cb50f/xlsxToR.r')
	download.file('http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48639-non-specific-probes-Illumina450k.xlsx', destfile=file.path(wd, 'crossreactive.xlsx'))
	crossreactive <- xlsxToR(file.path(wd, 'crossreactive.xlsx'), keep_sheets=c('nonspecific cg probes', 'nonspecific ch probes'), header=T)
	crossreactive <- do.call(rbind, crossreactive)
	crossreactive <- crossreactive[,1]
} else if (array.type == 'IlluminaHumanMethylationEPIC') {
	crossreactive <- scan('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909830/bin/mmc2.txt', what='character')
	crossreactive <- c(crossreactive, scan('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909830/bin/mmc3.txt', what='character'))
}

crossreactive <- crossreactive[crossreactive %in% annot.bed$name]
exclude <- crossreactive
annot.bed[crossreactive, 'score'] <- annot.bed[crossreactive, 'score'] + 2
#o#

## Flag blacklisted probes ##
if (blacklist != '') {
# Load blacklist
	remove <- read.delim(blacklist, sep='\t', stringsAsFactors=F, colClasses=c('character', 'character'), header=F)
	colnames(remove) <- c('probeID', 'sample')
	remove <- remove[remove$probeID %in% annot.bed$name,]
	flag <- unique(remove[remove$sample=='',1])
	annot.bed[flag, 'score'] <- annot.bed[flag, 'score'] + 4
	remove <- remove[remove$sample!='',]
	if (nrow(remove) > 0) {
		remove <- remove[remove$sample %in% targets$Sample_Name,]
		remove$value <- TRUE
		library(reshape2)
		remove <- dcast(probeID ~ sample, data=remove, fill=FALSE, value.var='value')
		row.names(remove) <- remove$probeID
		remove <- remove[,-1]
		remove <- as.matrix(remove)
		colnames(remove) <- rownames(targets)[match(colnames(remove), targets$Sample_Name)]
		remove <- remove[row.names(remove) %in% row.names(lowQ),]
		lowQ[row.names(remove), colnames(remove)] <- lowQ[row.names(remove), colnames(remove)] | remove
	}
#o#
	rm('flag', 'remove')
}

### Flag polymorphic probes ###
if (removeEuropeanSNPs) {
# Process SNPs
	if (array.type == 'IlluminaHumanMethylation450k') {
		download.file('http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48640-polymorphic-CpGs-Illumina450k.xlsx', destfile=file.path(wd, 'polymorphic.xlsx'))
		polymorphic <- xlsxToR(file.path(wd, 'polymorphic.xlsx'), keep_sheets='Polymorphic CpGs & SNPs at SBE', header=T)
		polymorphic[,grep('AF', colnames(polymorphic))] <- apply(polymorphic[,grep('AF', colnames(polymorphic))], 2, as.numeric)
	} else if (array.type == 'IlluminaHumanMethylationEPIC') {
		download.file('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4909830/bin/mmc1.txt', destfile=file.path(wd, 'polymorphic.txt'))
		polymorphic <- read.delim(file.path(wd, 'polymorphic.txt'), header=T, stringsAsFactors=F)
	}
	af <- 1/ncol(rgset)
	european.snps <- polymorphic[polymorphic[,'EUR_AF'] > af,1]
	european.snps <- european.snps[!is.na(european.snps)]
	european.snps <- european.snps[european.snps %in% annot.bed$name]
	annot.bed[european.snps, 'score'] <- annot.bed[european.snps, 'score'] + 1
#o#
	exclude <- c(exclude, european.snps)
}

### Output raw tables ###
raw.betas <- getBeta(preprocessRaw(rgset))
#colnames(raw.betas) <- targets[colnames(raw.betas), 'Sample_Name']
raw.betas.bed <- data.frame(annot.bed[rownames(raw.betas),], raw.betas, stringsAsFactors=F, check.names=F)

save(rgset, file=file.path(wd, 'rgset.RData'))
gz <- gzfile(file.path(wd, 'raw_betas.bed.gz'), 'w', compression=9)
write.table(raw.betas.bed, gz, sep="\t", quote=F, row.names=F)
close(gz)

#### Process Data ####

### Remove background ###
norm.meth <- NULL
message("Correcting the data...")
norm.meth <- preprocessENmix(rgset, bgParaEst='est', dyeCorr='RELIC', QCinfo=NULL, exQCsample=F, exQCcpg=F, exSample=NULL, exCpG=NULL, nCores=ncores)
#norm.meth <- preprocessNoob(rgset, offset=15, dyeCorr=T, dyeMethod='single')
#norm.meth <- preprocessSWAN(rgset, norm.meth)
processed.betas <- rcp(norm.meth)
message("Data corrected.")

### Mask probes
#message("Masking unreliable values...")

## Remove measures with low detection ###
#message("Masking done.")

### Extract genotyping probes ###
genotype.betas <- getSnpBeta(rgset)
#colnames(genotype.betas) <- targets[colnames(genotype.betas), 'Sample_Name']

### Sex determination ###
ratio.meth <- mapToGenome(norm.meth, mergeManifest=T)
gender <- getSex(ratio.meth)
if (usePredictedSex) {
	pData(norm.meth)$predictedSex <- factor(gender$predictedSex)
	interest.vars <- c(interest.vars, 'predictedSex')
}

### Output preprocessed data ###
message("Writing output...")
#processed.betas <- getBeta(norm.meth)
processed.betas[lowQ] <- NA
#colnames(processed.betas) <- targets[colnames(processed.betas), 'Sample_Name']
betas.bed <- data.frame(annot.bed[rownames(processed.betas),], processed.betas, stringsAsFactors=F, check.names=F)
#processed.mval <- getM(norm.meth)
processed.mval <- logit2(processed.betas)
processed.mval[lowQ] <- NA
#colnames(processed.mval) <- targets[colnames(processed.mval), 'Sample_Name']
mval.bed <- data.frame(annot.bed[rownames(processed.mval),], processed.mval, stringsAsFactors=F, check.names=F)
save(norm.meth, file=file.path(wd, 'bgCorr_meth.RData'))
gz <- gzfile(file.path(wd, 'normalized_betas.bed.gz'), 'w', compression=9)
write.table(betas.bed, gz, sep="\t", quote=F, row.names=F)
close(gz)
gz <- gzfile(file.path(wd, 'normalized_M.bed.gz'), 'w', compression=9)
write.table(mval.bed, gz, sep="\t", quote=F, row.names=F)
close(gz)
message("Data ready in the output folder.")

pdata <- pData(norm.meth)
#rownames(pdata) <- pdata$Sample_Name
pdata.out <- pdata
if (!usePredictedSex) {
	pdata.out$predictedSex <- factor(gender$predictedSex)
}
write.csv(pdata.out, file.path(wd, 'extended_sample_sheet.csv'), quote=F, row.names=F)
#o#

rm('array.annot', 'betas.bed', 'crossreactive', 'lowQ', 'lowQ.probes', 'mval.bed', 'ratio.meth', 'raw.betas.bed')
gc()
