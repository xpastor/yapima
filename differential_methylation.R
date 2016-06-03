## Required variables: ##
## pdata
## batch.vars
## processed.mval
## array.annot

message('Starting differential methylation analysis...')

library(ggplot2)

#### Differential methylation analysis ####
library(limma)

### Variable selection for analysis ###
int.formula <- paste(interest.vars, collapse=' + ')
int.formula <- as.formula(paste0('~ ', int.formula))
design <- model.matrix(int.formula, data=pdata[,interest.vars,drop=F])
diff.expr.mval <- processed.mval[,rownames(design)]
#o#
if (! surrogateCorrection) {
	# Limma analysis
	fit <- lmFit(diff.expr.mval, design)
	fit <- eBayes(fit)
	coefs <- colnames(fit$coefficients)[-1]
	for (coef in 1:length(interest.vars)) {
		message(paste0("Differential methylation analysis on '", interest.vars[coef], "'..."))
		top <- topTable(fit, coef=coefs[1], number=Inf)
		colnames(top) <- paste(coefs[coef], colnames(top), sep='.')
		top <- cbind(top, array.annot[row.names(top),])
		filename <- paste(interest.vars[coef], 'diffMeth.gz', sep='_')
		gz <- gzfile(file.path(wd, filename), 'w', compression=9)
		write.table(top, gz, sep="\t", row.names=F, quote=F)
		close(gz)
		adjP <- grep('adj.P.Val', colnames(top))
	#o#
		sig <- row.names(top)[top[,adjP] <= 0.05 & !is.na(top[,adjP])]
#		sig.meth <- processed.mval[sig,] 
#		sig.meth <- sig.meth[! apply(is.na(sig.meth), 1, any),]
#		pca <- prcomp(t(sig.meth))
#		groups <- pdata[,interest.vars[coef]]
#		names(groups) <- row.names(pdata)
#		ggsave(file=file.path(wd, paste0('differentialMethylation_PCA_', coefs[coef], '.png')), plot.pca(pca, groups, interest.vars[coef]), width=20)
#		hc <- hclust(dist(t(processed.mval[sig,])))
	}
} else {
	### Selection of negative control probes ###
	library(ruv)
	library(missMethyl)
	test.mval <- diff.expr.mval[! apply(is.na(diff.expr.mval), 1, any),]
	inc <- getINCs(raw.meth) # 'missMethyl' package
	m.inc <- rbind(test.mval, inc)
	
	for (coef in 1:length(interest.vars)) {
		message(paste0("Differential methylation analysis on '", interest.vars[coef], "'..."))
		## Probe selection ##
		ctl <- rownames(m.inc) %in% rownames(inc)
		fit <- RUVfit(m.inc, design=design, coef=coef, ctl=ctl)
		fit <- RUVadj(fit)
		top <- topRUV(fit, num=Inf)
		ctl <- rownames(test.mval) %in% row.names(top)[top$p.ebayes.BH > 0.5]
		
		### Differential methylation analysis correcting for surrogate variables ###
		fit <- RUVfit(test.mval, design=design, coef=coef, ctl=ctl)
		fit <- RUVadj(fit)
		top <- topRUV(fit, number=Inf)
		hc <- hclust(dist(t(diff.expr.mval[sig,])))
		colnames(top) <- paste(interest.vars[coef], colnames(top), sep='.')
		top <- cbind(top, array.annot[row.names(top),])
		filename <- paste(interest.vars[coef], 'diffMeth.gz', sep='_')
		gz <- gzfile(file.path(wd, filename), 'w', compression=9)
		write.table(top, gz, sep="\t", row.names=F, quote=F)
		close(gz)
	#o#
		#fitvar <- varFit(processed.mval, design=design, coef=seq(length(interest.vars)))
		#topDV <- topVar(fitvar, coef=coef, num=nrow(top))
		sig <- row.names(top)[top$p.ebayes.BH <= 0.05]
#		sig.meth <- processed.mval[sig,]
#		sig.meth <- sig.meth[! apply(is.na(sig.meth), 1, any),]
#		pca <- prcomp(t(sig.meth))
#		groups <- pdata[,interest.vars[coef]]
#		names(groups) <- row.names(pdata)
#		ggsave(file=file.path(wd, paste0('differentialMethylation_PCA_', interest.vars[coef], '.png')), plot.pca(pca, groups, interest.vars[coef]), width=20)
	}
	sig.meth <- processed.mval[sig,] 
	sig.meth <- sig.meth[! apply(is.na(sig.meth), 1, any),]
	pca <- prcomp(t(sig.meth))
	groups <- pdata[,interest.vars[coef]]
	names(groups) <- row.names(pdata)
	ggsave(file=file.path(wd, paste0('differentialMethylation_PCA_', interest.vars[coef], '.png')), plot.pca(pca, groups, interest.vars[coef]), width=20)
}

message('Differential methylation analysis finished.')
