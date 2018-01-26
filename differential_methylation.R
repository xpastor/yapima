## Required variables: ##
## pdata
## batch.vars
## processed.mval
## array.annot

source(file.path(pipeline_dir, 'extractCoords.R'))
#### Differential methylation analysis ####
message('Starting differential methylation analysis...')

dmp.dir <- file.path(wd, 'DMP')
dir.create(dmp.dir, recursive=T)
dmr.dir <- file.path(wd, 'DMR')
dir.create(dmr.dir, recursive=T)

annot.dm <- as.data.frame(mcols(array.annot.gr[,-match(c('Name'), colnames(mcols(array.annot.gr)))]))
rownames(annot.dm) <- names(array.annot.gr)
annot.bed <- annot.bed[order(annot.bed[,1],annot.bed[,2]),]
#o#
library(ggplot2)

### Variable selection for analysis ###
#int.formula <- paste(interest.vars, collapse=' + ')
#int.formula <- as.formula(paste0('~ ', int.formula))
#design <- model.matrix(int.formula, data=pdata[,interest.vars,drop=F])
#diff.expr.mval <- processed.mval[,rownames(design)]

# Limma analysis
library(limma)

#design.vars <- c('Slide', batch.vars, interest.vars)
design.vars <- c(batch.vars, interest.vars)
my.formula <- paste(design.vars, collapse='+')
my.formula <- as.formula(paste('~', my.formula))
design <- model.matrix(my.formula, data=pdata[,design.vars, drop=F])
n.e. <- nonEstimable(design)
design <- design[, !colnames(design) %in% n.e.]
#no.replicates <- colnames(design)[apply(design, 2, function(x) any(table(x) == 1))]
categorical.vars <- interest.vars[sapply(pdata[, interest.vars], class) %in% c('character', 'factor')]
categorical.vars <- categorical.vars[!is.na(categorical.vars)]
no.replicates <- colnames(pdata[,categorical.vars, drop=F])[apply(pdata[,categorical.vars, drop=F], 2, function(x) any(table(x) == 1))]
if (length(no.replicates) > 0) {
	stop(paste0('The following variables contain groups without replicates: ', paste(no.replicates, collapse=', '), ". Replace them with 'NA' or reallocate them in other groups."))
}

coefs <- colnames(design)[-1]
diff.meth.mval <- processed.mval[row.names(annot.bed)[annot.bed$score == 0], rownames(design)]
fit <- lmFit(diff.meth.mval, design)
fit <- eBayes(fit)

betafit <- lmFit(processed.betas[row.names(diff.meth.mval),rownames(design)], design)
betafit <- eBayes(betafit)

## Prepare data for DMR analysis
library(DMRcate)
R.utils::reassignInPackage('extractCoords', 'DMRcate', extractCoords)

array.annot.gr <- sort(array.annot.gr)
annot.dm <- annot.dm[rownames(annot.bed),]
gc()
cate.array <- c(IlluminaHumanMethylation450k='450K', IlluminaHumanMethylationEPIC='EPIC')

for (comparison in interest.vars) {
	cat(paste0(comparison,"\n"))
	var.coefs <- grep(paste0('^', comparison), colnames(design), value=T)
	# DMP analysis
	top <- topTable(fit, coef=var.coefs, number=Inf, sort.by='none')
	coef.name <- gsub(paste0('^', comparison), paste0(comparison, '.'), var.coefs)
	if (length(var.coefs) == 1) {
		colnames(top) <- gsub('logFC', paste0(coef.name, '.logFC'), colnames(top))
	} else {
		colnames(top) <- gsub(paste0('^', comparison), paste0(comparison, '.'), colnames(top))
	}
	top <- cbind(annot.bed[row.names(top),], top, annot.dm[row.names(top),])
	colnames(top) <- gsub('X.chrom', '#chrom', colnames(top))
	filename <- paste0('DMP_', comparison, '.bed.gz')
	gz <- gzfile(file.path(dmp.dir, filename), 'w', compression=9)
	write.table(top, gz, sep="\t", row.names=F, quote=F)
	close(gz)
	adjP <- grep('adj.P.Val', colnames(top))
	sig <- row.names(top)[top[,adjP] <= 0.05 & !is.na(top[,adjP])]
	if (!isEmpty(sig)) {
#o#
		sig.meth <- processed.mval[sig,] 
		sig.meth <- sig.meth[! apply(is.na(sig.meth), 1, any),]
		pca <- prcomp(t(sig.meth))
		groups <- pdata[,comparison]
		names(groups) <- row.names(pdata)
		ggsave(file=file.path(dmp.dir, paste0('DMP_PCA_', comparison, '.png')), plot.pca(pca, groups, comparison), width=20)
		scaled.mval <- apply(processed.mval[sig,], 1, scale, center=T)
		rownames(scaled.mval) <- colnames(processed.mval)
		hc <- hclust(dist(scaled.mval))
		hc.width <- .get_dev_width(processed.mval, name='', annotation_names=comparison)
		pdf(file.path(dmp.dir, paste0('DMP_cluster_', comparison, '.pdf')), height=4, width=hc.width)
		Heatmap2(matrix(nrow=0, ncol=ncol(processed.mval), dimnames=list(NULL, colnames(processed.mval))), column_annotation=pdata[,comparison,drop=F], cluster_columns=hc, column_dend_height=unit(5, 'cm'))
		dev.off()

	#	rm('scaled.mval', 'sig', 'sig.meth')
		gc()

	# DMR analysis
#		chr <- as.character(seqnames(array.annot.gr))
#		pos <- end(array.annot.gr)
		top <- top[!is.na(top$adj.P.Val),]
		betatop <- topTable(betafit, coef=var.coefs, number=Inf)
		betafc <- if(length(var.coefs) == 1) {betatop[rownames(top),'logFC']} else {0}
		stat <- if(length(var.coefs) == 1) {top$t} else {sqrt(top$F)}
		cpg.annot <- data.frame(ID=rownames(top), stat=stat, CHR=top[, '#chrom'], pos=top[,'chromEnd'], betafc=betafc, indfdr=top$adj.P.Val, is.sig=top$adj.P.Val < 0.05)
		class(cpg.annot) <- 'annot'
		dmr <- dmrcate(cpg.annot, mc.cores=1)
		dmr.gr <- extractRanges(dmr, genome='hg19')
		dmr.gr <- dmr.gr[dmr.gr$no.cpgs > 1,]
		dmr.gr <- sort(dmr.gr)
		dmr.bed <- gr2bed(dmr.gr)
		dmr.bed <- cbind(dmr.bed, as.data.frame(mcols(dmr.gr)))
		colnames(dmr.bed)[1] <- '#chrom'
		if (length(var.coefs) > 1) {
			coef.name <- comparison
		}
		write.table(dmr.bed, file.path(dmr.dir, paste0('DMR_', coef.name, '.bed')), sep='\t', quote=F, row.names=F)
	#o#
		dmr.sig <- dmr.gr[order(dmr.gr$Stouffer)]
#		dmr.sig <- dmr.gr[dmr.gr$minfdr <= 0.01,]
#		dmr.score <- rank(dmr.sig$minfdr) + 2*rank(-abs(dmr.sig$meanbetafc)) + 2*rank(-dmr.sig$no.cpgs)
#		dmr.sig <- dmr.sig[order(dmr.score)]
		L <- cumsum(dmr.sig$no.cpgs)
		n.dmr <- sum(L <= 500)
		dmr.sig <- head(dmr.sig, n.dmr)
		overlap <- findOverlaps(array.annot.gr, dmr.sig)
		cpg.dmr <- array.annot.gr[queryHits(overlap),]
#		cpg.dmr$differentialDMR <- as.character(dmr.sig[subjectHits(overlap)])
		cpg.dmr <- cpg.dmr[annot.bed[names(cpg.dmr),'score']==0]
		
		pdf(file.path(dmr.dir, paste0('DMR_heatmap_', comparison, '.pdf')), height=10, width=hc.width)
		Heatmap2(processed.betas[names(cpg.dmr),], column_annotation=pdata[,comparison,drop=F], row_annotation=data.frame(Chromosome=as.character(seqnames(cpg.dmr)), row.names=names(cpg.dmr), stringsAsFactors=F), cluster_rows=F, show_row_names=F, name='Beta', row_dend_side='left')
		dev.off()
#		rm('cpg.dmr', 'dmr', 'dmr.bed', 'dmr.gr', 'dmr.sig', 'top.dmr')
	}
	gc()
}
message('Differential methylation analysis finished.')
