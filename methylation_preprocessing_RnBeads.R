#!/usr/bin/env Rscript-3.2.0

set.seed(seed)
library(RnBeads)
options(fftempdir=scratch)

rnb.options(import.table.separator=',', identifiers.column='Sample_Name', min.group.size=1)
raw.meth <- rnb.execute.import(list(idat_dir, sample.annotation), data.type='idat.dir')

#### Fetching array annotation ####
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
array.annot <- as.data.frame(IlluminaHumanMethylation450kanno.ilmn12.hg19@data$Other)
rnb.annot <- annotation(raw.meth)
array.annot <- cbind(rnb.annot, array.annot[row.names(rnb.annot), ])

## Output annotation ##
write.table(array.annot, file.path(wd, '450k_annotation.txt'), sep="\t", quote=F, row.names=T)

#### Filter data ####

### Remove ambiguous probes ###
filtered.raw.meth <- raw.meth
exclude <- scan(file.path(pipeline_dir, 'data', 'non-specific-probes-Illumina450k.csv'), what='character')

### Remove blacklisted probes ###
if (blacklist != '') {
	remove <- scan(blacklist, what='character')
	exclude <- c(exclude, remove)
}

exclude <- unique(exclude[!is.na(exclude)])
excluded <- row.names(meth(raw.meth, row.names=T)) %in% exclude
filtered.raw.meth <- remove.sites(raw.meth, excluded)
filtered.raw.betas <- meth(filtered.raw.meth, row.names=T)

### Remove outliers ###
if (! isEmpty(outliers)) filtered.raw.meth <- remove.samples(filtered.raw.meth, outliers)

### Output raw tables ###
save(filtered.raw.betas, file=file.path(wd, 'filtered_raw_betas.RData'))
write.table(filtered.raw.betas, file.path(wd, 'filtered_raw_betas.txt'), sep="\t", quote=F, row.names=T)

#### Normalize data ####
norm.meth <- rnb.execute.normalization(filtered.raw.meth, method=normalization, bgcorr.method=backgroundCorrection)

### Remove genotyping probes ###
genotype.betas <- meth(norm.meth, row.names=T)
remove.probe <- grepl('^rs', row.names(genotype.betas))
genotype.betas <- genotype.betas[remove.probe,]
filtered.norm.meth <- remove.sites(norm.meth, remove.probe)

## Output betas of genotyping probes ##
write.table(genotype.betas, file.path(wd, 'genotyping_betas.txt'), sep="\t", quote=F, row.names=T)

### Remove probes with interrogated CpGs mapping SNPs ###
if (removeEuropeanSNPs) {
	cpg.snps <- read.csv(file.path(pipeline_dir, 'data', 'polymorphic-CpGs-SNPs-Illumina450k.csv'), stringsAsFactors=F)
	af <- 1/ncol(meth(raw.meth))
	european.snps <- unique(cpg.snps$PROBE[cpg.snps$EUR_AF > af])
	european.snps <- european.snps[!is.na(european.snps)]
	remove.probe <- row.names(meth(filtered.norm.meth, row.names=T)) %in% european.snps
	filtered.norm.meth <- remove.sites(filtered.norm.meth, remove.probe)
}

filtered.betas <- meth(filtered.norm.meth, row.names=T)
transformed.betas <- filtered.betas
save(transformed.betas, file=file.path(wd, 'filtered_normalized_betas.RData'))
write.table(transformed.betas, file.path(wd, 'filtered_normalized_betas.txt'), sep="\t", quote=F, row.names=T)

### Correct batch efects ###
if (batchCorrection) {
	ph.dat <- pheno(filtered.norm.meth)
	ph.dat <- ph.dat[,colSums(! is.na(ph.dat)) != 0]
	ph.dat$Sentrix_ID <- factor(ph.dat$Sentrix_ID)
	transformed.betas <- transformed.betas[! apply(is.na(transformed.betas), 1, any),]
	mod <- model.matrix(~1, data=ph.dat)
	for (var in batch.vars) {
		transformed.betas <- ComBat(transformed.betas, batch=ph.dat[,var], mod=mod)
	}
	save(transformed.betas, file=file.path(wd, 'transformed_betas.RData'))
	write.table(transformed.betas, file.path(wd, 'batch_normalized_betas.txt'), sep="\t", quote=F, row.names=T)
}

#	ga.all <- GlobalAncova(corrected.betas, formula.full=~Sentrix_ID, formula.red=~1, model.dat=ph.dat, method='permutation', perm=10000)

#if (ga.all$test.result <= 0.05) {
#	illumina.cols <- c('Sample_Name', 'Sample_Well', 'Sample_Plate', 'Pool_ID', 'Sentrix_ID', 'Sentrix_Position')
#	interest.vars <- colnames(ph.dat)[! colnames(ph.dat) %in% c('barcode', illumina.cols)]
#	int.formula <- ~ 1
#	require(sva)
#	if (! isEmpty(interest.vars)) {
#		int.formula <- paste(interest.vars, collapse=' + ')
#		int.formula <- as.formula(paste0('~ ', int.formula))
#	}

#	mod <- model.matrix(int.formula, data=sva.vars)
#	mod0 <- model.matrix(~ Sentrix_ID + Sentrix_Position, data=sva.vars)
#	n.sv <- num.sv(filtered.betas[! apply(is.na(filtered.betas), 1, any),], mod, method='leek')
#	sva.betas <- sva(filtered.betas[! apply(is.na(filtered.betas), 1, any),], mod, mod0)
#	mod0 <- model.matrix(~ Sentrix_ID, data=sva.vars)
#}
#}

