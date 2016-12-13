# Load libraries
library(minfi)
library(ENmix)
library(GenomicRanges)

#### Reading in data ####
### Preparing targets data frame ###
#header <- readLines(sample.annotation, n=1)
#header <- unlist(strsplit(header, ','))
colClasses <- data.frame(t(rep('character', length(header))), stringsAsFactors=F)
colnames(colClasses) <- header
#interest.vars <- variablesOfInterest(colClasses, batch.vars)
colClasses[,batch.vars] <- 'factor'
colClasses[,interest.vars] <- 'factor'
targets <- read.csv(sample.annotation, colClasses=colClasses, stringsAsFactors=F)
colnames(targets)[colnames(targets)=='Sentrix_ID'] <- 'Slide'
colnames(targets)[colnames(targets)=='Sentrix_Position'] <- 'Array'
targets$Basename <- paste(targets$Slide, targets$Array, sep='_')
row.names(targets) <- targets$Basename

### Reading methylation data ###
message("Reading in methylation files...")
#raw.meth <- read.450k.exp(idat_dir, targets, extended=T, recursive=T)
rgset <- read.metharray.exp(idat_dir, targets, extended=T, recursive=T)
message("Data read.")

### Fetching array annotation ###
array.type <- annotation(rgset)['array'] # IlluminaHumanMethylation450k / IlluminaHumanMethylationEPIC
array.annot <- getAnnotation(rgset)
array.annot.gr <- GRanges(array.annot$chr, ranges=IRanges(array.annot$pos, array.annot$pos+1), strand=NULL, array.annot[,! colnames(array.annot) %in% c('chr', 'pos', 'strand')])
annot.bed <- data.frame(chrom=array.annot$chr, chromStart=array.annot$pos-1, chromEnd=array.annot$pos, name=array.annot$Name, score=rep(0, nrow(array.annot)), strand=array.annot$strand, stringsAsFactors=F, row.names=array.annot$Name)
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
	crossreactive <- scan('http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S221359601630071X/1-s2.0-S221359601630071X-mmc2.txt/283447/html/S221359601630071X/d122a9289e5bc82609233a2cdeac8fa4/mmc2.txt', what='character')
	ch_crossreactive <- scan('http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S221359601630071X/1-s2.0-S221359601630071X-mmc3.txt/283447/html/S221359601630071X/daff1031095143b45eeb59c65f1dd940/mmc3.txt', what='character')
	crossreactive <- c(crossreactive, ch_crossreactive)
}

crossreactive <- crossreactive[crossreactive %in% annot.bed$name]
exclude <- crossreactive
annot.bed[crossreactive, 'score'] <- annot.bed[crossreactive, 'score'] + 2
#o#

## Flag blacklisted probes ##
if (blacklist != '') {
# Load blacklist
	remove <- scan(blacklist, what='character')
	remove <- remove[remove %in% annot.bed$name]
	exclude <- c(exclude, remove)
	annot.bed[remove, 'score'] <- annot.bed[remove, 'score'] + 4
#o#
}

### Flag polymorphic probes ###
if (removeEuropeanSNPs) {
# Process SNPs
	if (array.type == 'IlluminaHumanMethylation450k') {
		download.file('http://www.sickkids.ca/MS-Office-Files/Research/Weksberg%20Lab/48640-polymorphic-CpGs-Illumina450k.xlsx', destfile=file.path(wd, 'polymorphic.xlsx'))
		polymorphic <- xlsxToR(file.path(wd, 'polymorphic.xlsx'), keep_sheets='Polymorphic CpGs & SNPs at SBE', header=T)
		polymorphic[,grep('AF', colnames(polymorphic))] <- apply(polymorphic[,grep('AF', colnames(polymorphic))], 2, as.numeric)
	} else if (array.type == 'IlluminaHumanMethylationEPIC') {
		polymorphic <- read.delim('http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S221359601630071X/1-s2.0-S221359601630071X-mmc1.txt/283447/html/S221359601630071X/69a1cbea394507150394a9faaef9d876/mmc1.txt', header=T, stringsAsFactors=F)
	}
	af <- 1/ncol(rgset)
	european.snps <- polymorphic[polymorphic[,'EUR_AF'] > af,1]
	european.snps <- european.snps[!is.na(european.snps)]
	european.snps <- european.snps[european.snps %in% annot.bed$name]
	annot.bed[european.snps, 'score'] <- annot.bed[european.snps, 'score'] + 1
#o#
	exclude <- c(exclude, european.snps)
}

# Exclude probes
#exclude <- unique(exclude[!is.na(exclude)])
#exclude <- exclude[exclude %in% array.annot$Name]

### Output raw tables ###
raw.betas <- getBeta(rgset)
colnames(raw.betas) <- targets[colnames(raw.betas), 'Sample_Name']
raw.betas.bed <- data.frame(annot.bed[rownames(raw.betas),], raw.betas, stringsAsFactors=F, check.names=F)

save(rgset, file=file.path(wd, 'rgset.RData'))
gz <- gzfile(file.path(wd, 'raw_betas.bed.gz'), 'w', compression=9)
write.table(raw.betas.bed, gz, sep="\t", quote=F, row.names=F)
close(gz)

#### Process Data ####

### Remove background ###
norm.meth <- NULL
raw.meth <- preprocessRaw(rgset) 
message("Correcting the data...")
norm.meth <- preprocessENmix(rgset, bgParaEst='oob', dyeCorr=T, QCinfo=NULL, exQCsample=F, exQCcpg=F, exSample=NULL, exCpG=NULL, nCores=ncores)
norm.meth <- preprocessSWAN(rgset, norm.meth)
message("Data corrected.")

### Mask probes
message("Masking unreliable values...")
filtered.norm.meth <- norm.meth
#assayDataElement(filtered.norm.meth, 'Meth')[exclude, ] <- NA
#assayDataElement(filtered.norm.meth, 'Unmeth')[exclude, ] <- NA

## Remove measures with low detection ###
assayDataElement(filtered.norm.meth, 'Meth')[lowQ] <- NA
assayDataElement(filtered.norm.meth, 'Unmeth')[lowQ] <- NA
message("Masking done.")

### Extract genotyping probes ###
genotype.betas <- getSnpBeta(rgset)
colnames(genotype.betas) <- targets[colnames(genotype.betas), 'Sample_Name']

## Output betas of genotyping probes ##
write.table(genotype.betas, file.path(wd, 'genotyping_betas.txt'), sep="\t", quote=F, row.names=T)

### Sex determination ###
ratio.meth <- mapToGenome(filtered.norm.meth, mergeManifest=T)
gender <- getSex(ratio.meth)
pData(filtered.norm.meth)$predictedSex <- factor(gender$predictedSex)

### Output preprocessed data ###
message("Writing output...")
processed.betas <- getBeta(filtered.norm.meth)
colnames(processed.betas) <- targets[colnames(processed.betas), 'Sample_Name']
betas.bed <- data.frame(annot.bed[rownames(processed.betas),], processed.betas, stringsAsFactors=F, check.names=F)
processed.mval <- getM(filtered.norm.meth)
colnames(processed.mval) <- targets[colnames(processed.mval), 'Sample_Name']
mval.bed <- data.frame(annot.bed[rownames(processed.mval),], processed.mval, stringsAsFactors=F, check.names=F)
save(filtered.norm.meth, file=file.path(wd, 'filtered_normalized_meth.RData'))
gz <- gzfile(file.path(wd, 'filtered_normalized_betas.bed.gz'), 'w', compression=9)
write.table(betas.bed, gz, sep="\t", quote=F, row.names=F)
close(gz)
gz <- gzfile(file.path(wd, 'filtered_normalized_M.bed.gz'), 'w', compression=9)
write.table(mval.bed, gz, sep="\t", quote=F, row.names=F)
close(gz)
message("Data ready in the output folder.")

pdata <- pData(filtered.norm.meth)
rownames(pdata) <- pdata$Sample_Name
write.csv(pdata, file.path(wd, 'extended_sample_sheet.csv'), quote=F, row.names=F)
#o#
