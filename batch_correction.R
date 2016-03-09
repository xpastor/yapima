library(sva)

### Correct batch efects ###
pdata.narm <- pdata[,colSums(! is.na(pdata)) != 0]
if ('Sentrix_ID' %in% batch.vars) {
	colnames(pdata.narm)[grep('Slide', colnames(pdata.narm))] <- 'Sentrix_ID'
}
if ('Sentrix_Position' %in% batch.vars) {
	colnames(pdata.narm)[grep('Array', colnames(pdata.narm))] <- 'Sentrix_Position'
}
processed.mval <- processed.mval[! apply(is.na(processed.mval), 1, any),]
mod <- model.matrix(~1, data=pdata.narm)
for (my.var in batch.vars) {
	message(paste0("Correcting for '", my.var, "' effect..."))
	batch <- pdata.narm[colnames(processed.mval), my.var]
	names(batch) <- colnames(processed.mval)
	processed.mval <- ComBat(processed.mval, batch=batch, mod=mod)
	message(paste0("Correction for '", my.var, "' effect finished."))
}
processed.beta <- 2^processed.mval/(1+2^processed.mval)
message("Saving the results...")
save(processed.mval, processed.beta, file=file.path(wd, 'batch_corrected.RData'))
write.table(processed.mval, file.path(wd, 'batch_normalized_M.txt'), sep="\t", quote=F, row.names=T)
write.table(processed.beta, file.path(wd, 'batch_normalized_betas.txt'), sep="\t", quote=F, row.names=T)
message("Batch corrected data already saved.")
