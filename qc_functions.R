#!/usr/bin/env Rscript-3.2.0

### Function for extracting interest variables ###
variablesOfInterest <- function(pdata, batch.vars=NULL)
{
	illumina.vars <- c('Sample_Name', 'Sample_Well', 'Sample_Plate', 'Sample_Group', 'Pool_ID', 'Basename', 'filenames', 'Sentrix_ID', 'Sentrix_Position')
	batch.vars <- c(batch.vars, 'Slide', 'Array')
	not.relevant <- c('predictedSex')
	interest.vars  <- colnames(pdata)[!colnames(pdata) %in% c(illumina.vars, batch.vars, not.relevant)]
	return(interest.vars)
}

### Functions for PCA plotting ###
.pca.grob <- function(pca, groups, comp1=1, comp2=2, legend=F)
{
	require(ggplot2)
	require(gridExtra)
	xvals <- pca$x[,comp1]
	yvals <- pca$x[,comp2]
	df <- data.frame(xvals=xvals, yvals=yvals, groups=groups)
	xtitle <- paste0('PC', comp1, ' (', round(summary(pca)$importance[2,comp1]*100, digits=1), '% of variance)')
	ytitle <- paste0('PC', comp2, ' (', round(summary(pca)$importance[2,comp2]*100, digits=1), '% of variance)')

	theme_pca <- theme_bw() + theme(axis.text=element_blank(), axis.ticks=element_blank(), legend.position='none', panel.grid=element_blank(), panel.border=element_rect(colour='black'))
	pca.plot <- ggplot(df) + geom_point(mapping=aes(x=xvals, y=yvals, colour=groups), size=4) + theme_pca + theme(plot.margin=unit(c(0,0,1,1), 'lines')) + xlab(xtitle) + ylab(ytitle)
	if (legend) {
		return(legend.plot <- pca.plot + theme(legend.position='right', legend.key=element_blank()) + labs(colour=''))
	} else {
		theme_density <- theme_pca + theme(panel.border=element_blank(), axis.title.x=element_blank())
		xdensity <- ggplot(df) + geom_density(aes(x=xvals, colour=groups), na.rm=T, adjust=2) + geom_density(aes(x=xvals), na.rm=T, adjust=2, linetype='dotted') + theme_density + theme(plot.margin=unit(c(2,0,0,1), 'lines')) + ylab('')
		ydensity <- ggplot(df) + geom_density(aes(x=yvals, colour=groups), na.rm=T, adjust=2) + geom_density(aes(x=yvals), na.rm=T, adjust=2, linetype='dotted') + theme_density + theme(plot.margin=unit(c(0,1,1,0), 'lines')) + coord_flip() + xlab('')

		return(arrangeGrob(xdensity, rectGrob(gp=gpar(lty='blank')), pca.plot, ydensity, ncol=2, nrow=2, widths=c(3,1), heights=c(1,3)))
	}
}

.legend.grob <- function(a.gplot) {
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
	return(tmp$grobs[[leg]])
}

plot.pca <- function(pca, groups, main=NULL) {
	require(gridExtra)
	grid.newpage()
	legend.plot <- .pca.grob(pca, groups, 1, 2, legend=T)
	main <- textGrob(main, gp=gpar(fontsize=20), just='top')
	return(arrangeGrob(.pca.grob(pca, groups, 1, 2), .pca.grob(pca, groups, 3, 4), .pca.grob(pca, groups, 5, 6), .legend.grob(legend.plot), widths=c(4, 4, 4, 1), ncol=4, main=main))
}

### Functions for bootstrap Clustering ###
bootstrapClustering <- function(n, mat, transpose=F, nboot=100, ncores=1, cl=cl)
{
	require(pvclust)
	print(n)
	top.mat <- head(mat, nprobes)
	if (transpose) {top.mat <- t(top.mat)}
	clust <- parPvclust(cl, top.mat, method.dist='euclidean', use.cor='na.or.complete', nboot=nboot)
	return(list(n=n, cluster=clust))
}

score.cluster <- function(clust, num.edges=nrow(clust$cluster$edges))
{
	if (is.null(num.edges)) {
		num.edges <- round(nrow(clust$cluster$edges)*0.25)
	}
	height <- clust$cluster$hclust$height
	pval <- clust$cluster$edges$au
	pval <- pval[order(height, decreasing=T)]
	pval <- head(pval, num.edges)
	height <- sort(height, decreasing=T)
	height <- head(height, num.edges)
	return(sum(pval * height / max(height)))
}

score.clusters <- function(clustList, num.edges=NULL)
{
	sapply(clustList, score.cluster, num.edges)
}

## Function to determine consensus for kmeans clustering ###

get.consensus.kmeans <- function(mat, k=2, ntrials=5, iterations=10, random.seed=F)
{
	seeds <- seq(ntrials)
	kmeans.trials <- matrix(nrow=nrow(mat), ncol=ntrials, dimnames=list(row.names(mat), NULL))
	for (i in seq(ntrials)) {
		if (! random.seed) set.seed(seeds[i])
		x <- kmeans(mat, centers=k, iter.max=iterations)$cluster
		kmeans.trials[names(x),i] <- x
	}
	d <- dist(kmeans.trials)
	hc <- hclust(d)
	consensus <- factor(cutree(hc, k))
	res <- list(k=k, trials=kmeans.trials, consensus=consensus)
}

get.range.kmeans <- function(mat, maxK=5, ...)
{
	res <- list()
	for (k in 2:maxK) {
		kmeans.res <- get.consensus.kmeans(mat, ..., k=k)
		res[[k-1]] <- kmeans.res
	}
	return(res)
}

get.best.kmeans <- function(mat, maxK=2, ntrials=5, iterations=10)
{
	## Silhouette plot ##
	require(fpc)
	asw <- numeric(20)
	for (k in 2:20) asw[[k]] <- pam(top.betas,k) $ silinfo $ avg.width
	best.k <- which.max(asw)

	## Affinity propagation clustering ##
	require(apcluster)
	apclus <- apcluster(negDistMat(r=2), top.betas)
	best.k <- c(best.k, length(apclus@clusters))
#	heatmap(apclus)
#	plot(apclus, top.betas)
	
	## Gap statistic ##
	require(cluster)
	gap <- clusGap(top.betas, kmeans, 10, B=500, verbose=interactive())
	best.k <- c(best.k, which.max(gap$Tab[,'gap']))
	maxSE()
#	plot(1:10, gap$Tab[,'gap'], type='b', xlab='Number of Clusters', ylab='Gap')

	k <- 
	get.consensus.kmeans(mat, k=k, ntrials=ntrials, iterations=iterations)
}

get.range.kmeans <- function(mat, maxK=5, ...)
{
	res <- list()
	for (k in 2:maxK) {
		kmeans.res <- get.consensus.kmeans(top.betas, ..., k=k)
		res[[k-1]] <- kmeans.res
	}
	return(res)
}
