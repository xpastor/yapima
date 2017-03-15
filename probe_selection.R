## Required variables: ##
## pipeline_dir
## processed.betas
## ncores
## pdata
## batch.vars
## interest.vars

#### Probe selection ####
source(file.path(pipeline_dir, 'functions.R'))

library(GenomicRanges)
library(cluster)
library(parallel)

### Optimal probe set ###
#clust.betas <- corrected.betas
## most variable probes ##
betas.pass <- processed.betas[annot.bed[rownames(processed.betas), 'score'] == 0,]
betas.sd <- apply(betas.pass, 1, sd)
betas.sorted <- betas.pass[order(betas.sd, decreasing=T),]

## bootstrap clustering ##
library(pvclust)
steps <- seq(1000, 12000, 1000)
cl <- makeCluster(ncores, type='FORK')
clust.obj <- lapply(steps, bootstrapClustering, betas.sorted, transpose=F, nboot=10000, ncores=ncores, cl=cl)
stopCluster(cl)
save(clust.obj, file=file.path(wd, 'pvclust.RData'))
scores <- score.clusters(clust.obj)
#o#
#interest.vars <- variablesOfInterest(pdata, batch.vars)

library(ComplexHeatmap)
plot.vars <- interest.vars
if (usePredictedSex) plot.vars <- c(plot.vars, 'predictedSex')
width.dev <- .get_dev_width(betas.sorted, annotation_names=plot.vars)
height.dev <- .get_dev_width(betas.sorted, name='AAA')

pdf(file.path(qcdir, 'pvclust_clusters.pdf'), width=width.dev)
for (i in 1:length(clust.obj)) {
	plot(clust.obj[[i]]$cluster, main=paste(clust.obj[[i]]$n, ' probes, score=', round(scores[i], 3), sep=''), cex=0.6)
	hc <- clust.obj[[i]]$cluster$hclust
	Heatmap2(head(betas.sorted, clust.obj[[i]]$n), name="Beta\nvalues", column_annotation=pdata[,plot.vars,drop=F], show_row_dend=F, cluster_columns=hc, show_row_names=F, heatmap_legend_param=list(at=seq(0,1,length.out=6)))
}
dev.off()

# Output selected probes
clust <- clust.obj[[which.max(scores)]]
nprobes <- clust$n
clust.betas <- head(betas.sorted, nprobes)
write.table(clust.betas, file.path(wd, paste0('betas_top_', nprobes, '_variable_probes.txt')), sep="\t", quote=F, row.names=T)
#o#

sample.cor.top <- cor(clust.betas, use='na.or.complete')
pdf(file.path(qcdir, paste0('sample_correlation_top_', nprobes,'_variable_probes.pdf')), width=width.dev, height=height.dev)
Heatmap2(sample.cor.top, name="Correlation", column_annotation=pdata[,plot.vars,drop=F], row_annotation=pdata[,plot.vars,drop=F], col=c('lightyellow', 'black'))
dev.off()

### PCA analysis on optimal probe set ###
#pca <- prcomp(t(clust.betas))
#message("Only variables without 'NA' values will be plotted.")
#ph.dat <- ph.dat[,colSums(! is.na(ph.dat)) != 0]
#filtered.betas.narm <- filtered.betas[! apply(is.na(filtered.betas), 1, any),]
pca <- prcomp(t(clust.betas))
pdata2 <- pdata[row.names(pca$x),]

library(ggplot2)
for (interest.var in plot.vars) {
	groups <- pdata2[,interest.var]
	names(groups) <- row.names(pdata2)
	if (class(groups) %in% c('character', 'factor')) {
		ggsave(file=file.path(qcdir, paste0('top', nprobes, '_variable_probes_PCA_', interest.var, '.png')), plot.pca(pca, groups, interest.var), width=20)
	}
}
