## Required variables: ##
## pipeline_dir
## qcdir
## exclude
## raw.betas
## norm.meth
## processed.betas
## pdata
## genotype.betas

#source(file.path(pipeline_dir, 'functions.R'))

library(GenomicRanges)
#library(cluster)
#library(parallel)

### CNV analysis ###
message('Running CNV analysis...')
library(conumee)
library(CopyNeutralIMA)
cnv.dir <- file.path(qcdir, 'CNV_report')
dir.create(cnv.dir, recursive=T)
#load(file.path(wd, 'filtered_normalized_meth.RData'))
cnv.intensity <- getMeth(norm.meth) + getUnmeth(norm.meth)
colnames(cnv.intensity) <- paste(targets[colnames(cnv.intensity), 'Sample_Name'], 'intensity', sep='.')
RGsetCtrl <- getCopyNeutralRGSet(array.type)
MsetCtrl <- preprocessENmix(RGsetCtrl, bgParaEst='est', dyeCorr='RELIC', exQCsample=F, exQCcpg=F, nCores=ncores)

exclude <- unique(exclude)
exclude.gr <- array.annot.gr[exclude]

cnv <- CNV.load(as.data.frame(cnv.intensity), names=sub('\\.intensity$', '', colnames(cnv.intensity)))
cnv.controls <- CNV.load(MsetCtrl)
conumee.array <- c(IlluminaHumanMethylation450k='450k', IlluminaHumanMethylationEPIC='EPIC')
cnv.annot <- CNV.create_anno(exclude_regions=exclude.gr, chrXY=T, array_type=conumee.array[array.type])
for (pid in names(cnv)) {
	cat(pid)
    cnv.analysis <- CNV.fit(cnv[pid], cnv.controls, cnv.annot)
	cnv.analysis <- CNV.bin(cnv.analysis)
	cnv.analysis <- CNV.detail(cnv.analysis)
	cnv.analysis <- CNV.segment(cnv.analysis)
	cnv.res <- CNV.write(cnv.analysis, what='segments')
	cnv.probes <- CNV.write(cnv.analysis, what='probes')
#	dir.create(file.path(qcdir, 'CNV_report', pid))
#	write.table(cnv.res, file.path(qcdir, 'CNV_report', pid, paste(pid, 'CNV_report.txt', sep='_')), sep="\t", quote=F, row.names=F)
	write.table(cnv.res, file.path(cnv.dir, paste(pid, 'CNV_report.txt', sep='_')), sep="\t", quote=F, row.names=F)
	pdf(file.path(cnv.dir, paste(pid, 'CNV_report.pdf', sep='_')), width=20, height=5)
#	png(file.path(qcdir, 'CNV_report', pid, paste(pid, 'whole_genome.png', sep='_')), width=2000, height=500)
	CNV.genomeplot(cnv.analysis)
#	dev.off()
	for (chr in row.names(cnv.analysis@anno@genome)) {
#		png(file.path(qcdir, 'CNV_report', pid, paste(pid, '_', chr, '.png', sep='')), width=2000, height=500)
		CNV.genomeplot(cnv.analysis, chr=chr)
#		dev.off()
	}
	dev.off()
}
message('Finished.')
