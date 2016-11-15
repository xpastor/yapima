## Required variables: ##
## pipeline_dir
## qcdir
## exclude
## raw.betas
## filtered.norm.meth
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
cnv.intensity <- getMeth(filtered.norm.meth) + getUnmeth(filtered.norm.meth)
colnames(cnv.intensity) <- paste(targets[colnames(cnv.intensity), 'Sample_Name'], 'intensity', sep='.')
#annotation(filtered.norm.meth)$array == 'IlluminaHumanMethylation450k'
library(CopyNumber450kData)
data(RGcontrolSetEx)
#annotation(filtered.norm.meth)$array == 'IlluminaHumanMethylationEPIC'
controls.norm <- preprocessENmix(RGcontrolSetEx)
controls.norm <- preprocessSWAN(RGcontrolSetEx, mSet=controls.norm)

exclude.gr <- array.annot.gr[exclude]

cnv <- CNV.load(as.data.frame(cnv.intensity), names=sub('\\.intensity$', '', colnames(cnv.intensity)))
cnv.controls <- CNV.load(controls.norm)
cnv.annot <- CNV.create_anno(exclude_regions=exclude.gr, chrXY=T)
dir.create(file.path(qcdir, 'CNV_report'), recursive=T)
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
	write.table(cnv.res, file.path(qcdir, 'CNV_report', paste(pid, 'CNV_report.txt', sep='_')), sep="\t", quote=F, row.names=F)
	pdf(file.path(qcdir, 'CNV_report', paste(pid, 'CNV_report.pdf', sep='_')), width=20, height=5)
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
