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

#### Data QC ####

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

### Beta density heatmap ###
message('Plotting density heatmaps of Beta values...')
heatmap.params <- list(cluster_rows=F, col=c('lightyellow', 'black'), name='Density')
dev.width <- .get_dev_width(raw.betas, name='Density')
pdf(file.path(qcdir, 'beta_distribution_heatmap.pdf'), width=dev.width)
h <- apply(raw.betas, 2, function(x) hist(x, breaks=seq(0,1,0.001), plot=F)$counts)
h <- h[nrow(h):1,]
do.call(Heatmap2, c(list(mat=h, column_title='Raw beta distribution'), heatmap.params))
#densityHeatmap(raw.betas, title='Raw beta distribution', range=c(0,1), column_names_gp=gpar(fontsize=7))

h <- apply(processed.betas, 2, function(x) hist(x, breaks=seq(0,1,0.001), plot=F)$counts)
h <- h[nrow(h):1,]
do.call(Heatmap2, c(list(mat=h, column_title='Normalized beta distribution'), heatmap.params))
#densityHeatmap(processed.betas, title='Normalized beta distribution', range=c(0,1), column_names_gp=gpar(fontsize=7))
dev.off()
message('Finished.')

### PCA analysis ###
#pdata2 <- pdata[,colSums(! is.na(pdata)) != 0]
pdata2 <- pdata
pdata2$Slide <- as.character(pdata2$Slide)
raw.betas.narm <- raw.betas[! apply(is.na(raw.betas), 1, any),]
processed.betas.narm <- processed.betas[! apply(is.na(processed.betas), 1, any),]
raw.pca <- prcomp(t(raw.betas.narm))
pca <- prcomp(t(processed.betas.narm))
pdata2 <- pdata2[row.names(pca$x),]

## Batch variables ##
message('PCA plots of batch variables...')
pdata2$ArrayRow <- gsub('C..', '', pdata2$Array)
pdata2$ArrayColumn <- gsub('R..', '', pdata2$Array)

plot.vars <- unique(c('Slide', 'ArrayRow', 'ArrayColumn', 'predictedSex', batch.vars))
n.comp <- min(ncol(pca$x), 6)
n.grobs <- ceiling(n.comp/2)
width.unit <- 13/20
width <- n.grobs*2*4 + width.unit
if (n.comp > 2) {
	for (batch.var in plot.vars) {
		groups <- pdata2[,batch.var]
		names(groups) <- row.names(pdata2)
		if (class(groups) %in% c('character', 'factor')) {
			ggsave(filename=file.path(qcdir, paste0('raw_batch_PCA_', batch.var, '.png')), plot=plot.pca(raw.pca, groups, batch.var), width=width)
			ggsave(filename=file.path(qcdir, paste0('processed_batch_PCA_', batch.var, '.png')), plot=plot.pca(pca, groups, batch.var), width=width)
		}
	}
}
message('Finished.')
     
## Interest variables ##
message('PCA plots of variables of interest...')
#interest.vars <- variablesOfInterest(pdata, batch.vars)
if (! isEmpty(interest.vars) & n.comp > 2) {
#filtered.betas.narm <- filtered.betas[! apply(is.na(filtered.betas), 1, any),]
	for (interest.var in interest.vars) {
		groups <- pdata2[,interest.var]
		names(groups) <- row.names(pdata2)
		if (class(groups) %in% c('character', 'factor')) {
			ggsave(filename=file.path(qcdir, paste0('raw_PCA_', interest.var, '.png')), plot=plot.pca(raw.pca, groups, interest.var), width=width)
			ggsave(filename=file.path(qcdir, paste0('processed_PCA_', interest.var, '.png')), plot=plot.pca(pca, groups, interest.var), width=width)
		}
	}
}
message('Finished.')

### Correlation between samples ###
## All probes ##
#pdf(file.path(qcdir, 'samples_correlation.pdf'))
message('Plotting sample correlations...')
sample.cor <- cor(processed.betas, use='na.or.complete')
plot.vars <- c(interest.vars, 'predictedSex')
dev.width <- .get_dev_width(sample.cor, name='Correlation', annotation_names=plot.vars)
dev.height <- .get_dev_width(sample.cor, name='AAA')
pdf(file.path(qcdir, 'samples_correlation.pdf'), width=dev.width, height=dev.height)
Heatmap2(sample.cor, name='Correlation', column_annotation=pdata[,plot.vars,drop=F], row_annotation=pdata[,plot.vars,drop=F], column_title='All probes', col=c('lightyellow', 'black'))
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
