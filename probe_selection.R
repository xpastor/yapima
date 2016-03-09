## Required variables: ##
## pipeline_dir
## processed.betas
## ncores
## pdata
## batch.vars

#### Probe selection ####
source(file.path(pipeline_dir, 'qc_functions.R'))

library(GenomicRanges)
library(cluster)
library(parallel)

### Optimal probe set ###
#clust.betas <- corrected.betas
## most variable probes ##
betas.sd <- apply(processed.betas, 1, sd)
betas.sorted <- processed.betas[order(betas.sd, decreasing=T),]

## bootstrap clustering ##
library(pvclust)
steps <- seq(1000, 12000, 1000)
cl <- makeCluster(ncores, type='FORK')
clust.obj <- lapply(steps, bootstrapClustering, betas.sorted, transpose=F, nboot=10000, ncores=ncores, cl=cl)
stopCluster(cl)
save(clust.obj, file=file.path(wd, 'pvclust.RData'))
scores <- score.clusters(clust.obj)
pdf(file.path(qcdir, 'pvclust_clusters.pdf'), width=30)
for (i in 1:length(clust.obj)) {
	plot(clust.obj[[i]]$cluster, main=paste(clust.obj[[i]]$nprobes, ' probes, score=', round(scores[i], 3), sep=''))
}
dev.off()

clust <- clust.obj[[which.max(scores)]]
nprobes <- clust$nprobes
clust.betas <- head(betas.sorted, nprobes)
write.table(clust.betas, file.path(wd, paste0('betas_top_', nprobes, '_variable_probes.txt')), sep="\t", quote=F, row.names=T)
sample.cor.top <- cor(head(betas.sorted, nprobes), use='na.or.complete')
pheatmap(sample.cor.top, show_rownames=T, show_colnames=T, fontsize=6, filename=file.path(qcdir, 'sample_correlation_top_probes.pdf'), main=paste(nprobes, 'most variable probes'))
#o#

### PCA analysis on optimal probe set ###
#pca <- prcomp(t(clust.betas))
#message("Only variables without 'NA' values will be plotted.")
#ph.dat <- ph.dat[,colSums(! is.na(ph.dat)) != 0]
interest.vars <- variablesOfInterest(pdata, batch.vars)
if (! isEmpty(interest.vars)) {
#filtered.betas.narm <- filtered.betas[! apply(is.na(filtered.betas), 1, any),]
	pca <- prcomp(t(clust.betas))
	pdata2 <- pdata[row.names(pca$x),]

	for (interest.var in interest.vars) {
		groups <- pdata2[,interst.var]
		names(groups) <- row.names(pdata2)
		if (class(groups) %in% c('character', 'factor')) {
			ggsave(file=file.path(qcdir, paste0('top', nprobes, '_variable_probes_PCA_', interest.var, '.png')), plot.pca(pca, groups, interest.var), width=20)
		}
	}
}
