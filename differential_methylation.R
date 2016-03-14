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
interest.vars <- variablesOfInterest(pdata, batch.vars)
if (isEmpty(interest.vars)) {
	message("There's no factor eligible for a differential methylation analysis and it will be skipped.")
} else {
	int.formula <- paste(interest.vars, collapse=' + ')
	int.formula <- as.formula(paste0('~ ', int.formula))
	design <- model.matrix(int.formula, data=pdata[,interest.vars,drop=F])
#o#
	if (! surrogateCorrection) {
		# Limma analysis
		fit <- lmFit(processed.mval, design)
		fit <- eBayes(fit)
		coefs <- colnames(fit$coefficients)[-1]
		for (coef in 1:length(interest.vars)) {
			message(paste0("Differential methylation analysis on '", interest.vars[coef], "'..."))
			top <- topTable(fit, coef=coefs[1], number=Inf)
			sig <- row.names(top)[top$adj.P.Val <= 0.05 & !is.na(top$adj.P.Val)]
			colnames(top) <- paste(coefs[coef], colnames(top), sep='.')
			top <- cbind(top, array.annot[row.names(top),])
			filename <- paste(coefs[coef], 'diffMeth.gz', sep='_')
			gz <- gzfile(file.path(wd, filename), 'w', compression=9)
			write.table(top, gz, sep="\t", row.names=F, quote=F)
			close(gz)
		#o#
			sig.meth <- processed.mval[sig,] 
			sig.meth <- sig.meth[! apply(is.na(sig.meth), 1, any),]
			pca <- prcomp(t(sig.meth))
			ggsave(file=file.path(wd, paste0('differentialMethylation_PCA_', coefs[coef], '.png')), plot.pca(pca, groups, interest.vars[coef]), width=20)
#			hc <- hclust(dist(t(processed.mval[sig,])))
		}
	} else {
		### Selection of negative control probes ###
		library(ruv)
		library(missMethyl)
		test.mval <- processed.mval[! apply(is.na(processed.mval), 1, any),]
		inc <- getINCs(filtered.raw.meth) # 'missMethyl' package
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
#			hc <- hclust(dist(t(processed.mval[sig,])))
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
			sig.meth <- processed.mval[sig,]
			pca <- prcomp(t(processed.mval[sig,]))
			groups <- pdata[,interest.vars[coef]]
			names(groups) <- row.names(pdata)
			ggsave(file=file.path(wd, paste0('differentialMethylation_PCA_', interest.vars[coef], '.png')), plot.pca(pca, groups, interest.vars[coef]), width=20)
		}
	}
}

message('Differential methylation analysis finished.')
