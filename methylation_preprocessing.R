# Load libraries
library(minfi)

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
raw.meth <- read.metharray.exp(idat_dir, targets, extended=T, recursive=T)
message("Data read.")

#### Fetching array annotation ####
array.annot <- getAnnotation(raw.meth)
array.annot.gr <- GRanges(array.annot$chr, ranges=IRanges(array.annot$pos, array.annot$pos+1), mcols=array.annot[,! colnames(array.annot) %in% c('chr', 'pos')])

## Output annotation ##
write.table(array.annot, file.path(wd, '450k_annotation.txt'), sep="\t", quote=F, row.names=T)

### Output raw tables ###
raw.betas <- getBeta(raw.meth)
colnames(raw.betas) <- targets[colnames(raw.betas), 'Sample_Name']

#save(raw.meth, file=file.path(wd, 'filtered_raw_meth.RData'))
gz <- gzfile(file.path(wd, 'raw_betas.gz'), 'w', compression=9)
write.table(raw.betas, gz, sep="\t", quote=F, row.names=T)
close(gz)
#o#

#### Process Data ####

### Remove background ###
norm.meth <- NULL
if (! backgroundCorrection) {
# Produce raw objects
	norm.meth <- preprocessRaw(raw.meth) 
#o#
} else {
# Remove background
	message("Removing the background...")
	norm.meth <- preprocessNoob(raw.meth)
	message("Background corrected.")
#o#
}

### Normalize data ###
if (normalization) {
# Normalization
	message("Normalizing data...")
	norm.meth <- preprocessSWAN(raw.meth, mSet=norm.meth)
	message("Data normalized.")
#o#
}

### Filter data ###

## Load ambiguous probes ##
exclude <- read.csv(non_specific_cg, quote='', header=T, stringsAsFactors=F)
exclude <- exclude$TargetID
exclude2 <- read.csv(non_specific_ch, quote='', header=T, stringsAsFactors=F)
exclude <- c(exclude, exclude2$TargetID)
#o#

## Remove blacklisted probes ##
if (blacklist != '') {
# Load blacklist
	remove <- scan(blacklist, what='character')
	exclude <- c(exclude, remove)
#o#
}

### Remove probes with interrogated CpGs mapping SNPs ###
if (removeEuropeanSNPs) {
# Process SNPs
	cpg.snps <- read.csv(polymorphic, stringsAsFactors=F)
	af <- 1/ncol(norm.meth)
	european.snps <- unique(cpg.snps$PROBE[cpg.snps[,'EUR_AF'] > af])
	european.snps <- european.snps[!is.na(european.snps)]
	exclude <- c(exclude, european.snps)
#o#
}

# Exclude probes
exclude <- unique(exclude[!is.na(exclude)])
exclude <- exclude[exclude %in% array.annot$Name]

# Mask probes
filtered.norm.meth <- norm.meth
assayDataElement(filtered.norm.meth, 'Meth')[exclude, ] <- NA
assayDataElement(filtered.norm.meth, 'Unmeth')[exclude, ] <- NA

### Remove measures with low detection ###
lowQ <- detectionP(raw.meth) > 0.01
assayDataElement(filtered.norm.meth, 'Meth')[lowQ] <- NA
assayDataElement(filtered.norm.meth, 'Unmeth')[lowQ] <- NA

### Extract genotyping probes ###
genotype.betas <- getSnpBeta(raw.meth)
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
processed.mval <- getM(filtered.norm.meth)
colnames(processed.mval) <- targets[colnames(processed.mval), 'Sample_Name']
save(filtered.norm.meth, file=file.path(wd, 'filtered_normalized_meth.RData'))
gz <- gzfile(file.path(wd, 'filtered_normalized_betas.gz'), 'w', compression=9)
write.table(processed.betas, gz, sep="\t", quote=F, row.names=T)
close(gz)
gz <- gzfile(file.path(wd, 'filtered_normalized_M.gz'), 'w', compression=9)
write.table(processed.mval, gz, sep="\t", quote=F, row.names=T)
close(gz)
message("Data ready in the output folder.")

pdata <- pData(filtered.norm.meth)
rownames(pdata) <- pdata$Sample_Name
write.csv(pdata, file.path(wd, 'extended_sample_sheet.csv'), quote=F, row.names=F)
#o#
