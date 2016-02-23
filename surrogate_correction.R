## Required R objects ##
## pdata
## batch.vars
## filtered.raw.meth
## varianceProportion
## processed.mval

library(ruv)
library(missMethyl)
library(limma)
library(ggplot2)

### Variable selection for analysis ###
interest.vars <- variablesOfInterest(pdata, batch.vars)
int.formula <- paste(interest.vars, collapse=' + ')
pdata2 <- pdata
if (isEmpty(interest.vars)) {
	int.formula <- ~ 1
	pdata2 <- pdata[,interest.vars,drop=F]
} else {
	int.formula <- as.formula(paste0('~ ', int.formula))
}
design <- model.matrix(int.formula, data=pdata2)

#### Batch correction and analysis ####
### Selection of negative control probes ###
## k determination for RUVrinv ##
message('Determining negative control probes...')
inc <- getINCs(filtered.raw.meth) # 'missMethyl' package
m.inc <- rbind(processed.mval, inc)
ctl <- rownames(m.inc) %in% rownames(inc)
pca <- prcomp(t(m.inc))
k <- sum(summary(pca)$importance[2,] >= varianceProportion)

### Probe selection ###
fit <- RUVrinv(t(m.inc), design[,-1,drop=F], ctl=ctl, k=k) # 'ruv' package
fit <- missMethyl:::.toMArrayLM(fit) # 'missMethyl' package
fit <- RUVadj(fit) # 'missMethyl' package
ctl <- rownames(processed.mval) %in% row.names(m.inc)[!apply(fit$p.ebayes.BH <= 0.5, 1, any)]
message('Negative control probes ready.')

### Determining the coefficients for correction of surrogate variables ###
message('Correcting the data for surrogate variables...')
fit <- RUVrinv(t(processed.mval), design[,-1,drop=F], ctl=ctl, k=k) # 'ruv' package
#fit <- missMethyl:::.toMArrayLM(fit) # 'missMethyl' package
#fit <- RUVadj(fit) # 'missMethyl' package

### Transormation of M-values and beta-values ###
surrogated.mval <- removeBatchEffect(processed.mval, covariates=fit$W, design=design) # 'limma' package
surrogated.beta <- 2^surrogated.mval/(1+2^surrogated.mval)
message('Correction on surrogate variables finished.')

## Plot PCA of corrected values ##
message('Plotting PCA of corrected data...')
if (! isEmpty(interest.vars)) {
	pca <- prcomp(t(surrogated.beta))
	pdata2 <- pdata[row.names(pca$x),]
    for (interest.var in interest.vars) {
		groups <- pdata2[,interest.var]
		names(groups) <- row.names(pdata2)
		if (class(groups) %in% c('character', 'factor')) {
			ggsave(file=file.path(qcdir, paste0('surrogate_PCA_', interest.var, '.png')), plot.pca(pca, groups, interest.var), width=20)
		}
	}
}
message('PCA plots ready.')

### Ouptut transformed data ###
message('Saving data corrected for surrogate variables...')
save(surrogated.mval, surrogated.beta, file=file.path(wd, 'surrogated_data.RData'))
mval.out <- surrogated.mval
colnames(mval.out) <- paste0('Mval.', colnames(mval.out))
mval.out <- cbind(data.frame(mval.out, stringsAsFactors=F), array.annot[rownames(mval.out),])
write.table(mval.out, file.path(wd, 'Mvalue_surrogateCorrection.txt'), sep="\t", row.names=F, quote=F)
beta.out <- surrogated.beta
colnames(beta.out) <- paste0('Beta.', colnames(beta.out))
beta.out <- cbind(data.frame(beta.out, stringsAsFactors=F), array.annot[rownames(beta.out),])
write.table(beta.out, file.path(wd, 'Beta_surrogateCorrection.txt'), sep="\t", row.names=F, quote=F)
message('Data saved.')
