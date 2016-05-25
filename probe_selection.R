## Required variables: ##
## pipeline_dir
## processed.betas
## ncores
## pdata
## batch.vars
## interest.vars

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
#o#
#interest.vars <- variablesOfInterest(pdata, batch.vars)

library(ComplexHeatmap)
char.height <- convertHeight(grobHeight(textGrob("A", gp=gpar(fontsize = 7))), 'inch', valueOnly=T)
width.mat <- ncol(betas.sorted) * char.height
width.dev <- width.mat + convertWidth(unit(1, 'cm'), 'inch', valueOnly=T) + 10
if (! isEmpty(interest.vars)) {
	longest.title <- interest.vars[which.max(nchar(interest.vars))]
	width.title <- convertWidth(grobHeight(textGrob("A", gp=gpar(fontsize = 10))), 'inch', valueOnly=T) * (nchar(longest.title))
	width.dev <- width.dev + width.title
}

pdf(file.path(qcdir, 'pvclust_clusters.pdf'), width=width.dev)
for (i in 1:length(clust.obj)) {
	plot(clust.obj[[i]]$cluster, main=paste(clust.obj[[i]]$n, ' probes, score=', round(scores[i], 3), sep=''), cex=0.6)
	hc <- clust.obj[[i]]$cluster$hclust
	if (!isEmpty(interest.vars)) {
		ha <- HeatmapAnnotation(df=pdata[,interest.vars], gp=gpar(col='black'))
		Heatmap(head(betas.sorted, clust.obj[[i]]$n), name="Beta\nvalues", top_annotation=ha, show_row_hclust=F, cluster_columns=hc, show_row_names=F, column_names_gp=gpar(fontsize=7))
		#for(ann in colnames(heat.annot2)) {
		#   decorate_annotation(ann, {grid.text(ann, unit(1, 'npc') + unit(2, 'mm'), just='left')})
		#}
	} else {
		Heatmap(head(betas.sorted, clust.obj[[i]]$n), name="Beta\nvalues", show_row_hclust=F, cluster_columns=hc, show_row_names=F, column_names_gp=gpar(fontsize=7))
	}
}
dev.off()

# Output selected probes
clust <- clust.obj[[which.max(scores)]]
nprobes <- clust$n
clust.betas <- head(betas.sorted, nprobes)
write.table(clust.betas, file.path(wd, paste0('betas_top_', nprobes, '_variable_probes.txt')), sep="\t", quote=F, row.names=T)
#o#

sample.cor.top <- cor(head(betas.sorted, nprobes), use='na.or.complete')
pdf(file.path(qcdir, paste0('sample_correlation_top_', nprobes,'_variable_probes.pdf')), width=width.dev)
if (!isEmpty(interest.vars)) {
	col_annot <- HeatmapAnnotation(df=pdata[,interest.vars], gp=gpar(col='black'))
	row_annot <- rowAnnotation(df=pdata[,interest.vars], show_legend=F, gp=gpar(col='black'))
	draw(row_annot + Heatmap(sample.cor.top, name="Correlation", top_annotation=col_annot, column_names_gp=gpar(fontsize=7), row_names_gp=gpar(fontsize=7)), row_dend_side='left')
	#for(ann in colnames(heat.annot2)) {
	#   decorate_annotation(ann, {grid.text(ann, unit(1, 'npc') + unit(2, 'mm'), just='left')})
	#}
} else {
	Heatmap(sample.cor.top, name="Correlation", column_names_gp=gpar(fontsize=7), row_names_gp=gpar(fontsize=7))
}
dev.off()

### PCA analysis on optimal probe set ###
#pca <- prcomp(t(clust.betas))
#message("Only variables without 'NA' values will be plotted.")
#ph.dat <- ph.dat[,colSums(! is.na(ph.dat)) != 0]
if (! isEmpty(interest.vars)) {
#filtered.betas.narm <- filtered.betas[! apply(is.na(filtered.betas), 1, any),]
	pca <- prcomp(t(clust.betas))
	pdata2 <- pdata[row.names(pca$x),]

	for (interest.var in interest.vars) {
		groups <- pdata2[,interest.var]
		names(groups) <- row.names(pdata2)
		if (class(groups) %in% c('character', 'factor')) {
			ggsave(file=file.path(qcdir, paste0('top', nprobes, '_variable_probes_PCA_', interest.var, '.png')), plot.pca(pca, groups, interest.var), width=20)
		}
	}
}
