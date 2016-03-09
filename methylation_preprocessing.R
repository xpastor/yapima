# Load libraries
library(minfi)
#o#

#### Reading in data ####
### Preparing targets data frame ###
targets <- read.csv(sample.annotation, stringsAsFactors=F)
colnames(targets)[colnames(targets)=='Sentrix_ID'] <- 'Slide'
colnames(targets)[colnames(targets)=='Sentrix_Position'] <- 'Array'
targets$Basename <- paste(targets$Slide, targets$Array, sep='_')
row.names(targets) <- targets$Basename

### Reading methylation data ###
message("Reading in methylation files...")
raw.meth <- read.450k.exp(idat_dir, targets, extended=T)
message("Data read.")

#### Fetching array annotation ####
array.annot <- getAnnotation(raw.meth)
array.annot.gr <- GRanges(array.annot$chr, ranges=IRanges(array.annot$pos, array.annot$pos+1), mcols=array.annot[,! colnames(array.annot) %in% c('chr', 'pos')])

## Output annotation ##
write.table(array.annot, file.path(wd, '450k_annotation.txt'), sep="\t", quote=F, row.names=T)

#### Filter data ####

### Remove ambiguous probes ###
exclude <- scan(file.path(pipeline_dir, 'data', 'non-specific-probes-Illumina450k.csv'), what='character')
#o#

### Remove blacklisted probes ###
if (blacklist != '') {
# Load blacklist
	remove <- scan(blacklist, what='character')
	exclude <- c(exclude, remove)
#o#
}
# Exclude probes
exclude <- unique(exclude[!is.na(exclude)])

excluded <- c(array.annot[exclude, 'AddressA'], array.annot[exclude, 'AddressB'])
excluded <- excluded[excluded != '']
filtered.raw.meth <- raw.meth
sampleNames(filtered.raw.meth) <- targets[sampleNames(filtered.raw.meth), 'Sample_Name']
assayDataElement(filtered.raw.meth, 'Green')[excluded,] <- NA
assayDataElement(filtered.raw.meth, 'Red')[excluded,] <- NA
#o#

### Output raw tables ###
filtered.raw.betas <- getBeta(filtered.raw.meth)

#save(filtered.raw.meth, file=file.path(wd, 'filtered_raw_meth.RData'))
write.table(filtered.raw.betas, file.path(wd, 'filtered_raw_betas.txt'), sep="\t", quote=F, row.names=T)
#o#

#### Process Data ####

### Remove background ###
norm.meth <- NULL
#o#
if (! backgroundCorrection) {
# Produce raw objects
	norm.meth <- preprocessRaw(filtered.raw.meth) 
#o#
} else {
# Remove background
	message("Removing the background...")
	norm.meth <- preprocessNoob(filtered.raw.meth)
	message("Background corrected.")
#o#
}

### Normalize data ###
if (normalization) {
# Normalization
	message("Normalizing data...")
	norm.meth <- preprocessSWAN(filtered.raw.meth, mSet=norm.meth)
	message("Data normalized.")
#o#
}

### Extract genotyping probes ###
genotype.betas <- getSnpBeta(filtered.raw.meth)

## Output betas of genotyping probes ##
write.table(genotype.betas, file.path(wd, 'genotyping_betas.txt'), sep="\t", quote=F, row.names=T)

### Sex determination ###
ratio.meth <- mapToGenome(norm.meth, mergeManifest=T)
gender <- getSex(ratio.meth)
pData(norm.meth)$predictedSex <- gender$predictedSex

### Remove probes with interrogated CpGs mapping SNPs ###
remove.probes <- NULL
#o#
if (removeEuropeanSNPs) {
# Process SNPs
	cpg.snps <- read.csv(file.path(pipeline_dir, 'data', 'polymorphic-CpGs-SNPs-Illumina450k.csv'), stringsAsFactors=F)
	af <- 1/ncol(norm.meth)
	european.snps <- unique(cpg.snps$PROBE[cpg.snps$EUR_AF > af])
	european.snps <- european.snps[!is.na(european.snps)]
	remove.probes <- c(remove.probes, european.snps)
#o#
}

# Mask probes
filtered.norm.meth <- norm.meth
assayDataElement(filtered.norm.meth, 'Meth')[remove.probes, ] <- NA
assayDataElement(filtered.norm.meth, 'Unmeth')[remove.probes, ] <- NA

### Remove measures with low detection ###
lowQ <- detectionP(filtered.raw.meth) > 0.01
assayDataElement(filtered.norm.meth, 'Meth')[lowQ] <- NA
assayDataElement(filtered.norm.meth, 'Unmeth')[lowQ] <- NA

### Output preprocessed data ###
message("Writing output...")
processed.betas <- getBeta(filtered.norm.meth)
processed.mval <- getM(filtered.norm.meth)
save(filtered.norm.meth, file=file.path(wd, 'filtered_normalized_meth.RData'))
write.table(processed.betas, file.path(wd, 'filtered_normalized_betas.txt'), sep="\t", quote=F, row.names=T)
write.table(processed.mval, file.path(wd, 'filtered_normalized_M.txt'), sep="\t", quote=F, row.names=T)
message("Data ready in the output folder.")

pdata <- pData(filtered.norm.meth)
write.table(pdata, file.path(wd, 'extended_sample_sheet.csv'), quote=F, row.names=F)
#o#
