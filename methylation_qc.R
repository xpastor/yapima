## Required variables: ##
## pipeline_dir
## qcdir
## raw.betas
## processed.betas
## pdata
## genotype.betas

source(file.path(pipeline_dir, 'functions.R'))

library(GenomicRanges)
library(cluster)
library(parallel)
library(ggplot2)
library(ComplexHeatmap)

dev.width <- .get_dev_width(raw.betas, name='Density')

#### Data QC ####

### Median intensities plot ###
qc.raw <- getQC(preprocessRaw(rgset))
#row.names(qc.raw) <- row.names(pdata)[match(row.names(qc.raw), pdata$Basename)]
qc.norm <- DataFrame(mMed=apply(log2(getMeth(norm.meth)), 2, median, na.rm=T), uMed=apply(log2(getUnmeth(norm.meth)), 2, median, na.rm=T))
#row.names(qc.norm) <- row.names(pdata)[match(row.names(qc.norm), pdata$Basename)]
plotQC <- function(qc, badSampleCutoff = 10.5, main=NULL)
{
	meds <- rowMeans(as.matrix(qc))
	bad <- meds < badSampleCutoff
	plot(qc$mMed, qc$uMed, xlim=c(8,14), ylim=c(8,14), xaxt = "n", yaxt = "n", xlab = "Meth median intensity (log2)", ylab = "Unmeth median intensity (log2)", col = ifelse(bad, "red", "black"), main=main)
	axis(side = 1, at = c(9, 11, 13))
	axis(side = 2, at = c(9, 11, 13))
	abline(badSampleCutoff * 2, -1, lty = 2)
	if (sum(bad) > 0)
		text(qc$mMed[bad], qc$uMed[bad] - 0.25, labels = row.names(qc)[bad], col = "red")
	legend("topleft", legend = c("good", "bad, with sample index"), pch = 1, col = c("black", "red"), bty = "n")
	invisible(NULL)
}
pdf(file.path(qcdir, 'median_intensities.pdf'))
plotQC(qc.raw, main='Raw intensities')
plotQC(qc.norm, main='Background corrected intensities')
dev.off()

### Bisulfite conversion plots ###
pdf(file.path(qcdir, 'bisulfite_conversion.pdf'), height=dev.width)
controlStripPlot(rgset, controls='BISULFITE CONVERSION I', sampNames=pData(rgset)$Sample_Name)
controlStripPlot(rgset, controls='BISULFITE CONVERSION II', sampNames=pData(rgset)$Sample_Name)
dev.off()

### Beta distribution plot ###
message('Plotting density plots of Beta values...')
pdf(file.path(qcdir, 'beta_distribution.pdf'))
plot(density(raw.betas[,1], na.rm=T), ylim=c(0,6), main='Raw samples', bty='n', lwd=0.1, xlab='Beta')
for(i in 2:ncol(raw.betas)) {
	lines(density(raw.betas[,i], na.rm=T), lwd=0.1)
}

plot(density(processed.betas[,1], na.rm=T), ylim=c(0,6), main='Normalized samples', bty='n', lwd=0.1, xlab='Beta')
for(i in 2:ncol(processed.betas)) {
	lines(density(processed.betas[,i], na.rm=T), lwd=0.1)
}
dev.off()
message('Finished.')
#gc()

### Beta density heatmap ###
message('Plotting density heatmaps of Beta values...')
heatmap.params <- list(cluster_rows=F, col=c('lightyellow', 'black'), name='Density')
pdf(file.path(qcdir, 'beta_distribution_heatmap.pdf'), width=dev.width)
densityHeatmap(raw.betas, ylab='Beta value', title='Raw beta distribution', range=c(0,1), cluster_columns=T, show_column_dend=T, column_names_gp=gpar(fontsize=7))
densityHeatmap(processed.betas, ylab='Beta value', title='Normalized beta distribution', range=c(0,1), cluster_columns=T, show_column_dend=T, column_names_gp=gpar(fontsize=7))
dev.off()
message('Finished.')
gc()

### PCA analysis ###
#pdata2 <- pdata[,colSums(! is.na(pdata)) != 0]
pdata2 <- pdata
pdata2$Slide <- as.character(pdata2$Slide)
raw.betas.narm <- raw.betas[! apply(is.na(raw.betas), 1, any),]
processed.betas.narm <- processed.betas[! apply(is.na(processed.betas), 1, any),]
raw.pca <- prcomp(t(raw.betas.narm))
pca <- prcomp(t(processed.betas.narm))
pdata2 <- pdata2[row.names(pca$x),]
gc()

## Batch variables ##
message('PCA plots of batch variables...')
pdata2$ArrayRow <- gsub('C..', '', pdata2$Array)
pdata2$ArrayColumn <- gsub('R..', '', pdata2$Array)

#plot.vars <- unique(c('Slide', 'ArrayRow', 'ArrayColumn', 'predictedSex', batch.vars))
plot.vars <- unique(c('Slide', 'ArrayRow', 'ArrayColumn', batch.vars))
if (usePredictedSex) {
	plot.vars <- c(plot.vars, 'predictedSex')
}
n.comp <- min(ncol(pca$x), 6)
n.grobs <- ceiling(n.comp/2)
width.unit <- 13/20
width <- n.grobs*2*4 + width.unit
if (n.comp > 2) {
	for (batch.var in plot.vars) {
		message(batch.var)
		groups <- pdata2[,batch.var]
		names(groups) <- row.names(pdata2)
		ggsave(filename=file.path(qcdir, paste0('raw_batch_PCA_', batch.var, '.png')), plot=plot.pca(raw.pca, groups, batch.var), width=width)
		ggsave(filename=file.path(qcdir, paste0('processed_batch_PCA_', batch.var, '.png')), plot=plot.pca(pca, groups, batch.var), width=width)
	}
}
message('Finished.')
     
## Interest variables ##
message('PCA plots of variables of interest...')
#interest.vars <- variablesOfInterest(pdata, batch.vars)
if (! isEmpty(interest.vars) & n.comp > 2) {
#filtered.betas.narm <- filtered.betas[! apply(is.na(filtered.betas), 1, any),]
	for (interest.var in interest.vars) {
		message(interest.var)
		groups <- pdata2[,interest.var]
		names(groups) <- row.names(pdata2)
		ggsave(filename=file.path(qcdir, paste0('raw_PCA_', interest.var, '.png')), plot=plot.pca(raw.pca, groups, interest.var), width=width)
		ggsave(filename=file.path(qcdir, paste0('processed_PCA_', interest.var, '.png')), plot=plot.pca(pca, groups, interest.var), width=width)
	}
}
message('Finished.')

### Correlation between samples ###
## All probes ##
#pdf(file.path(qcdir, 'samples_correlation.pdf'))
message('Plotting sample correlations...')
sample.cor <- cor(processed.betas, use='na.or.complete')
plot.vars <- interest.vars
if (usePredictedSex) plot.vars <- unique(c(plot.vars, 'predictedSex'))
categorical <- colnames(pdata)[sapply(pdata, class) %in% c('character', 'factor')]
plot.vars <- plot.vars[plot.vars %in% categorical]
dev.width <- .get_dev_width(sample.cor, name='Correlation', annotation_names=plot.vars)
dev.height <- .get_dev_width(sample.cor, name='AAA')
pdf(file.path(qcdir, 'samples_correlation.pdf'), width=dev.width, height=dev.height)
Heatmap2(mat=sample.cor, name='Correlation', column_annotation=pdata[,plot.vars,drop=F], row_annotation=pdata[,plot.vars,drop=F], column_title='All probes', col=c('lightyellow', 'black'))
dev.off()
message('Finished.')

### Sample genotyping ###
message('Genotyping samples...')
snps <- row.names(genotype.betas)
library(biomaRt)
#mart <- useMart('snp')
#snp.db <- useMart('snp', dataset='hsapiens_snp')
snp.db <- useMart('ENSEMBL_MART_SNP', dataset='hsapiens_snp', host="www.ensembl.org")
snp.genotype <- getBM(c('chr_name', 'chrom_start', 'chrom_end', 'refsnp_id', 'allele'), filters='snp_filter', values=snps, mart=snp.db)
rownames(snp.genotype) <- snp.genotype$refsnp_id
snp.unmeth <- gsub('[GC/]', '', snp.genotype$allele)
snp.meth <- ifelse(snp.unmeth=='A', 'G', 'C')
ref <- gsub('/..*', '', snp.genotype$allele)
genotypes.df <- data.frame(ref=ref, hipo=paste0(snp.unmeth, snp.unmeth), hemi=paste0(snp.meth,snp.unmeth), hyper=paste0(snp.meth,snp.meth), row.names=snp.genotype$refsnp_id, stringsAsFactors=F)

genotype.sample <- function(snps.betas) {
    snps <- row.names(snps.betas)
	library(fpc)
	library(cluster)
	asw <- numeric(3)
	for (k in 2:3) asw[[k]] <- pam(snps.betas,k) $ silinfo $ avg.width
	k.best <- which.max(asw)
	
	hc <- hclust(dist(snps.betas))
	haplotype <- cutree(hc, k.best)
	genotype.mean <- aggregate(snps.betas, list(haplotype), mean)
	hyper <- genotype.mean$Group.1[which.max(genotype.mean$x)]
	hipo <- genotype.mean$Group.1[which.min(genotype.mean$x)]
	haplotype <- ifelse(haplotype==hipo, 'hipo', ifelse(haplotype==hyper, 'hyper', 'hemi'))
}

smpl.genotypes <- apply(genotype.betas, 2, genotype.sample)
snp.idx <- match(row.names(smpl.genotypes), row.names(genotypes.df))
genotype.idx <- apply(smpl.genotypes, 2, match, colnames(genotypes.df))
genotypes <- matrix(genotypes.df[snp.idx + (genotype.idx - 1)*nrow(genotypes.df)], nrow=nrow(smpl.genotypes), ncol=ncol(smpl.genotypes), dimnames=list(row.names(smpl.genotypes), colnames(smpl.genotypes)))
genotypes <- data.frame(chrom=snp.genotype[rownames(genotypes), 'chr_name'], chromStart=snp.genotype[rownames(genotypes), 'chrom_start'], chromEnd=snp.genotype[rownames(genotypes), 'chrom_end'], name=row.names(genotypes), score=rep('.', nrow(genotypes)), strand=rep('.', nrow(genotypes)), genotypes, stringsAsFactors=F)
colnames(genotypes)[1] <- paste0('#', colnames(genotypes)[1])
write.table(genotypes, file.path(qcdir, 'genotypes.bed'), sep="\t", row.names=F, quote=F)
message('Finished.')

rm('genotype.betas', 'genotypes', 'genotypes.df', 'pca', 'polymorphic', 'processed.betas.narm', 'raw.betas', 'raw.betas.narm', 'raw.pca', 'rgset', 'smpl.genotypes', 'snp.genotype', 'snp.idx', 'snp.meth', 'snp.unmeth', 'snps')
gc()
