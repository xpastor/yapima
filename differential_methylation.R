## Required variables: ##
## pdata
## batch.vars
## processed.mval
## array.annot

#### Differential methylation analysis ####
message('Starting differential methylation analysis...')

annot.dm <- array.annot[,-match(c('chr', 'pos', 'Name', 'strand'), colnames(array.annot))]
#o#
library(ggplot2)

### Variable selection for analysis ###
#int.formula <- paste(interest.vars, collapse=' + ')
#int.formula <- as.formula(paste0('~ ', int.formula))
#design <- model.matrix(int.formula, data=pdata[,interest.vars,drop=F])
#diff.expr.mval <- processed.mval[,rownames(design)]

# Limma analysis
library(limma)
for (comparison in interest.vars) {
	design <- model.matrix(as.formula(paste('~', comparison)), data=pdata[,interest.vars,drop=F])
	diff.expr.mval <- processed.mval[,rownames(design)]
	message(paste0("Differential methylation analysis on '", comparison, "'..."))
	fit <- lmFit(diff.expr.mval, design)
	fit <- eBayes(fit)
	top <- NULL
	if (length(levels(pdata[,comparison])) == 2) {
		top <- topTable(fit, coef=2, number=Inf, sort.by='none')
	} else if (length(levels(pdata[,comparison])) > 2) {
		top <- topTableF(fit, number=Inf, sort.by='none')
	}
	colnames(top) <- paste(colnames(design)[2], colnames(top), sep='.')
	top <- cbind(annot.bed[row.names(top),], top, annot.dm[row.names(top),])
	filename <- paste(comparison, 'differentialMethylation.bed.gz', sep='_')
	gz <- gzfile(file.path(wd, filename), 'w', compression=9)
	write.table(top, gz, sep="\t", row.names=F, quote=F)
	close(gz)
#o#
	adjP <- grep('adj.P.Val', colnames(top))
	sig <- row.names(top)[top[,adjP] <= 0.05 & !is.na(top[,adjP])]
	sig.meth <- processed.mval[sig,] 
	sig.meth <- sig.meth[! apply(is.na(sig.meth), 1, any),]
	pca <- prcomp(t(sig.meth))
	groups <- pdata[,comparison]
	names(groups) <- row.names(pdata)
	ggsave(file=file.path(wd, paste0('differentialMethylation_PCA_', comparison, '.png')), plot.pca(pca, groups, comparison), width=20)
	scaled.mval <- apply(processed.mval[sig,], 1, scale, center=T)
	rownames(scaled.mval) <- colnames(processed.mval)
	hc <- hclust(dist(scaled.mval))
	hc.width <- .get_dev_width(processed.mval, name='', annotation_names=comparison)
	pdf(file.path(wd, paste0('differentialMethylation_cluster_', comparison, '.pdf')), height=4, width=hc.width)
	Heatmap2(matrix(nrow=0, ncol=ncol(processed.mval), dimnames=list(NULL, colnames(processed.mval))), column_annotation=pdata[,comparison,drop=F], cluster_columns=hc, column_dend_height=unit(5, 'cm'))
	dev.off()
}

message('Differential methylation analysis finished.')
