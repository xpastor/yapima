library(sva)
library(limma)
library(ggplot2)

edata <- processed.mval
edata.narm <- edata[apply(edata, 1, function(x)!any(is.na(x))),]

mod <- model.matrix(~ Sample + predictedSex, data=pdata)
mod0 <- model.matrix(~1,data=pdata)

n.sv <- num.sv(edata.narm, mod, method='leek')
svobj <- sva(edata.narm, mod, mod0, n.sv=n.sv)
pValues <- f.pvalue(edata.narm, mod, mod0)
qValues <- p.adjust(pValues, method="BH")

modSv <- cbind(mod, svobj$sv)
mod0Sv <- cbind(mod0, svobj$sv)
pValuesSv <- f.pvalue(edata.narm, modSv, mod0Sv)
qValuesSv <- p.adjust(pValuesSv, method="BH")

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

surrogated.mval <- removeBatchEffect(edata.narm, covariatese=svobj$sv, design=design)
surrogated.beta <- 2^surrogated.mval/(1+2^surrogated.mval)

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

