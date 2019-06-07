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
exclude <- NULL

lowQ <- detectionP(rgset) > 0.01
lowQ.probes <- rowSums(lowQ)/ncol(lowQ) >= 0.5
annot.bed$score[lowQ.probes] <- annot.bed$score[lowQ.probes] + 8

#o#

## Flag polymorphic probes ##
if (removeSNPs) {
	load(file.path(pipeline_dir, 'polymorphic.rda'))
	# Process SNPs
	polymorphic <- polymorphic[[array.type]]$CpG_SBE
	available <- grep('AF$', colnames(polymorphic), value=T)
	available <- gsub('_?AF$', '', available)
	if (! population %in% available) {
		available <- gsub('(.*)', '\'\\1\'', available)
		available <- paste(available, collapse=', ')
		stop(paste0('\'',population, '\' population is not available for \'', array.type,'\'. Available populations are: ', available, '.'))
	}
	af <- 0.005
	polymorphic <- polymorphic$PROBE[!is.na(polymorphic[polymorphic[,paste(c(population, 'AF'), collapse="_")] > af,1])]
	polymorphic <- polymorphic[polymorphic %in% annot.bed$name]
	num_snps <- length(polymorphic)
	annot.bed[polymorphic, 'score'] <- annot.bed[polymorphic, 'score'] + 1
#o#
	exclude <- c(exclude, polymorphic)
}

## Load crossreactive probes ##
load(file.path(pipeline_dir, 'crossreactive.rda'))
#load('https://github.com/xpastor/yapima/crossreactive.rda')
# Process crossreactive
crossreactive <- crossreactive[[array.type]]
crossreactive <- do.call(c, crossreactive)
crossreactive <- crossreactive[crossreactive %in% annot.bed$name]
exclude <- unique(c(exclude, crossreactive))
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
#norm.meth <- preprocessENmix(rgset, bgParaEst='est', dyeCorr='RELIC', QCinfo=NULL, exQCsample=F, exQCcpg=F, exSample=NULL, exCpG=NULL, nCores=ncores)
norm.meth <- preprocessNoob(rgset, offset=15, dyeCorr=T, dyeMethod='single')
norm.meth <- preprocessSWAN(rgset, norm.meth)
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
pdata$Slide <- as.character(pdata$Slide)
pdata$ArrayRow <- gsub('C..', '', pdata$Array)
pdata$ArrayColumn <- gsub('R..', '', pdata$Array)
#rownames(pdata) <- pdata$Sample_Name
pdata.out <- pdata
if (!usePredictedSex) {
	pdata.out$predictedSex <- factor(gender$predictedSex)
}
write.csv(pdata.out, file.path(wd, 'extended_sample_sheet.csv'), quote=F, row.names=F)
#o#

rm('array.annot', 'betas.bed', 'crossreactive', 'lowQ', 'lowQ.probes', 'mval.bed', 'ratio.meth', 'raw.betas.bed')
gc()
