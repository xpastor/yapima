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
array.annot <- getAnnotation(rgset)
array.annot.gr <- GRanges(array.annot$chr, ranges=IRanges(array.annot$pos, array.annot$pos+1), mcols=array.annot[,! colnames(array.annot) %in% c('chr', 'pos')])
annot.bed <- data.frame(chrom=array.annot$chr, chromStart=array.annot$pos-1, chromEnd=array.annot$pos, name=array.annot$Name, score=rep(0, nrow(array.annot)), strand=array.annot$strand, stringsAsFactors=F, row.names=array.annot$Name)
colnames(annot.bed)[1] <- paste0('#', colnames(annot.bed)[1])

## Output annotation ##
write.table(array.annot, file.path(wd, 'annotation.txt'), sep="\t", quote=F, row.names=T)

### Filter data ###

## Load ambiguous probes ##
exclude <- read.csv(non_specific_cg, quote='', header=T, stringsAsFactors=F)
exclude <- exclude$TargetID
exclude2 <- read.csv(non_specific_ch, quote='', header=T, stringsAsFactors=F)
exclude <- c(exclude, exclude2$TargetID)
exclude <- exclude[exclude %in% annot.bed$name]
annot.bed[exclude, 'score'] <- annot.bed[exclude, 'score'] + 2
#o#

## Remove blacklisted probes ##
if (blacklist != '') {
# Load blacklist
	remove <- scan(blacklist, what='character')
	remove <- remove[remove %in% annot.bed$name]
	exclude <- c(exclude, remove)
	annot.bed[remove, 'score'] <- annot.bed[remove, 'score'] + 4
#o#
}

### Remove probes with interrogated CpGs mapping SNPs ###
if (removeEuropeanSNPs) {
# Process SNPs
	cpg.snps <- read.csv(polymorphic, stringsAsFactors=F)
	af <- 1/ncol(rgset)
	european.snps <- unique(cpg.snps$PROBE[cpg.snps[,'EUR_AF'] > af])
	european.snps <- european.snps[!is.na(european.snps)]
	european.snps <- european.snps[european.snps %in% annot.bed$name]
	annot.bed[european.snps, 'score'] <- annot.bed[european.snps, 'score'] + 1
	exclude <- c(exclude, european.snps)
#o#
}

# Exclude probes
exclude <- unique(exclude[!is.na(exclude)])
exclude <- exclude[exclude %in% array.annot$Name]

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
lowQ <- detectionP(rgset) > 0.01
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
