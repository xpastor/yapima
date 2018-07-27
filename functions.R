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
#o#

### Function to convert GenomicRanges objects into BED-like data frames
gr2bed <- function(gr, name=NULL, score=NULL, additional=NULL)
{
	library(GenomicRanges)
	#ifNull <- rep('.', length(gr))
	if (is.null(name)) {
		name <- as.character(gr)
	}
	if (is.null(score)) {
		score <- '.'
	}
	ifNull <- as.character(gr)
	data.frame(chrom=as.character(seqnames(gr)), chromStart=start(gr)-1, chromEnd=end(gr), name=name, score=score, strand=strand(gr), stringsAsFactors=F)
}
#o#

### Functions for PCA plotting ###
.pca.grob <- function(pca, groups, comp1=1, comp2=2, legend=F)
{
	library(ggplot2)
	library(gridExtra)
	xvals <- pca$x[,comp1]
	yvals <- pca$x[,comp2]
	df <- data.frame(xvals=xvals, yvals=yvals, groups=groups)
	xtitle <- paste0('PC', comp1, ' (', round(summary(pca)$importance[2,comp1]*100, digits=1), '% of variance)')
	ytitle <- paste0('PC', comp2, ' (', round(summary(pca)$importance[2,comp2]*100, digits=1), '% of variance)')

#	theme_pca <- theme_bw() + theme(axis.text=element_blank(), axis.ticks=element_blank(), legend.position='none', panel.grid=element_blank(), panel.border=element_rect(colour='black'))
#	pca.plot <- ggplot(df) + geom_point(mapping=aes(x=xvals, y=yvals, colour=groups), size=4) + theme_pca + theme(plot.margin=unit(c(0,0,1,1), 'lines')) + xlab(xtitle) + ylab(ytitle)
	theme_pca <- theme_bw() + theme(legend.position='none', panel.grid=element_blank(), panel.border=element_rect(colour='black'))
	pca.plot <- ggplot(df) + geom_point(mapping=aes(x=xvals, y=yvals, colour=groups), size=4) + theme_pca + theme(axis.text=element_blank(), axis.ticks=element_blank(), plot.margin=unit(c(0,0,1,1), 'lines')) + xlab(xtitle) + ylab(ytitle)
	if (legend) {
		return(legend.plot <- pca.plot + theme(legend.position='right', legend.key=element_blank()) + labs(colour=''))
	} else {
#		theme_density <- theme_pca + theme(panel.border=element_blank(), axis.title.x=element_blank())
#		xdensity <- ggplot(df) + geom_density(aes(x=xvals, colour=groups), na.rm=T, adjust=2) + geom_density(aes(x=xvals), na.rm=T, adjust=2, linetype='dotted') + theme_density + theme(plot.margin=unit(c(2,0,0,1), 'lines')) + ylab('')
#		ydensity <- ggplot(df) + geom_density(aes(x=yvals, colour=groups), na.rm=T, adjust=2) + geom_density(aes(x=yvals), na.rm=T, adjust=2, linetype='dotted') + theme_density + theme(plot.margin=unit(c(0,1,1,0), 'lines')) + coord_flip() + xlab('')
		theme_aux_plot <- theme_pca + theme(panel.border=element_blank(), axis.title.x=element_blank())
		groups_n <- table(groups)
		groups_d <- names(groups_n)[groups_n >= 2]
		df2 <- df[df$groups %in% groups_d,]
		axis_plot <- ggplot(df2) + theme_aux_plot + theme(axis.text=element_blank(), axis.ticks=element_blank())
#		yplot <- ggplot(df2)
		if (class(groups) %in% c('factor', 'character')) {
			g <- ggplot_build(pca.plot)
			g <- g$data[[1]]
			g$group <- groups
			g <- unique(g[,c('colour', 'group')])
			color_map <- g$colour
			names(color_map) <- g$group
			xplot <- axis_plot + geom_density(aes(x=xvals, colour=groups), na.rm=T, adjust=2) + geom_density(aes(x=xvals), data=df, na.rm=T, adjust=2, linetype='dotted') + theme(plot.margin=unit(c(2,0,0,1), 'lines')) + ylab('')
			yplot <- axis_plot + geom_density(aes(x=yvals, colour=groups), na.rm=T, adjust=2) + geom_density(aes(x=yvals), data=df, na.rm=T, adjust=2, linetype='dotted') + theme(plot.margin=unit(c(0,1,1,0), 'lines')) + coord_flip() + xlab('')
		} else {
			xplot <- xplot + geom_point(aes(x=xvals, y=groups), na.rm=T) + geom_smooth(aes(x=xvals, y=groups), na.rm=T, se=FALSE, method='loess', span=1) + theme_aux_plot + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), plot.margin=unit(c(2,0,0,0), 'lines')) + ylab('')
			yplot <- yplot + geom_point(aes(x=yvals, y=groups), na.rm=T) + geom_smooth(aes(x=yvals, y=groups), na.rm=T, se=FALSE, method='loess', span=1) + theme_aux_plot + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), plot.margin=unit(c(0,1,1,0), 'lines')) + coord_flip() + xlab('')
		}
		return(arrangeGrob(xplot, rectGrob(gp=gpar(lty='blank')), pca.plot, yplot, ncol=2, nrow=2, widths=c(3,1), heights=c(1,3)))
	}
}

.legend.grob <- function(a.gplot) {
	tmp <- ggplot_gtable(ggplot_build(a.gplot))
	leg <- which(sapply(tmp$grobs, function(x) x$name) == 'guide-box')
	return(tmp$grobs[[leg]])
}

plot.pca <- function(pca, groups, main=NULL) {
	n.comp <- min(ncol(pca$x), 6)
	if (n.comp > 1) {
		library(grid)
		library(gridExtra)
		grid.newpage()
		legend.plot <- .pca.grob(pca, groups, 1, 2, legend=T)
		main <- textGrob(main, gp=gpar(fontsize=20), just='top')
		n.grobs <- n.comp%/%2
		rep.comp <- n.comp < 6 & n.comp%%2 == 1
		list.grobs <- list()
		for (i in seq(n.grobs)) {
			list.grobs[[i]] <- .pca.grob(pca, groups, 2*i-1, 2*i)
		}
		if (rep.comp) {
			list.grobs[[length(list.grobs)+1]] <- .pca.grob(pca, groups, n.comp-1, n.comp)
		}
		widths <- rep(4, length(list.grobs))
		list.grobs[[length(list.grobs) + 1]] <-.legend.grob(legend.plot)
		widths <- c(widths, 1)
#		return(arrangeGrob(.pca.grob(pca, groups, 1, 2), .pca.grob(pca, groups, 3, 4), .pca.grob(pca, groups, 5, 6), .legend.grob(legend.plot), widths=c(4, 4, 4, 1), ncol=4, top=main))
		return(arrangeGrob(grobs=list.grobs, widths=widths, ncol=length(widths), top=main))
	}
}

### Function to produce Heatmaps with ComplexHeatmap ###
Heatmap2 <- function(mat, ..., column_annotation=NULL, row_annotation=NULL, column_names_gp=gpar(fontsize=7), row_names_gp=column_names_gp, row_dend_side='right', row_names_side='left', heatmap_legend_param=list(color_bar='continuous', legend_height=unit(3, 'cm')))
{
	library(ComplexHeatmap)
	heatmap.params <- list(matrix=mat, column_names_gp=column_names_gp, row_names_gp=row_names_gp, row_dend_side=row_dend_side, row_names_side=row_names_side, heatmap_legend_param=heatmap_legend_param, ...)
	annotation.width <- unit(2, 'mm')
	if (! is.null(column_annotation)) {
		ha_cols <- .annotation_colors(column_annotation)
		ha <- HeatmapAnnotation(df=column_annotation, col=ha_cols, gp=gpar(col='black'), na_col='white')
		heatmap.params$top_annotation <- ha
#		annotation.width <- grobWidth(textGrob(colnames(column_annotation), gp=gpar(fontsize=12))) + unit(1, 'mm') - grobWidth(textGrob(rownames(mat), gp=gpar(fontsize=7)))
		annotation.width <- grobWidth(textGrob(colnames(column_annotation), gp=gpar(fontsize=12))) + unit(1, 'mm') - grobWidth(textGrob(rownames(mat), gp=gpar(fontsize=7)))
	}
	hm <- do.call(Heatmap, args=heatmap.params)
	rows_hc <- ifelse(is.null(heatmap.params$cluster_rows), T, heatmap.params$cluster_rows)
	show_rows_hc <- ifelse(is.null(heatmap.params$show_row_dend), T, heatmap.params$show_row_dend)
	if (convertWidth(annotation.width, 'mm', valueOnly=T) < 0) annotation.width <- unit(2, 'mm')
	padding <- unit.c(unit(2, 'mm'), annotation.width, unit(2,'mm'), unit(2, 'mm'))
	if (is.null(row_annotation)) {
		draw(hm, padding=padding)
	} else {
		row_annot_params <- list(df=row_annotation, gp=column_names_gp, show_legend=!identical(row_annotation, column_annotation))
		if (identical(column_annotation, row_annotation)) {
			col <- list()
			for (anno in names(ha@anno_list)) {
				if (class(row_annotation[,anno]) %in% c('character', 'factor')) {
					col[[anno]] <- ha@anno_list[[anno]]@color_mapping@colors
				} else {
					col[[anno]] <- .annotation_colors(row_annotation[,anno, drop=F])[[anno]]
				}
			}
			row_annot_params <- c(row_annot_params, col=list(col), na_col='white')
		} else {
			col <- .annotation_colors(row_annotation)
			row_annot_params <- c(row_annot_params, col=list(col))
		}
		row_annot <- do.call(rowAnnotation, row_annot_params)
		draw(hm + row_annot, row_dend_side=row_dend_side, padding=padding)
	}
	if (! is.null(column_annotation)) {
		for(ann in colnames(column_annotation)) {
			if(rows_hc & show_rows_hc) {
					decorate_annotation(ann, {grid.text(ann, unit(0, 'npc') - unit(2, 'mm'), just='right', gp=gpar(fontsize=10))})
			} else {
				  	decorate_annotation(ann, {grid.text(ann, unit(1, 'npc') + unit(2, 'mm'), just='left', gp=gpar(fontsize=10))})
			}
		}
	}
	gc()
}

.get_dev_width <- function(mat, name='matrix_0', annotation_names=NULL, fontsize=7)
{
	library(grid)
	char.height <- convertHeight(grobHeight(textGrob("A", gp=gpar(fontsize = fontsize))), 'inch', valueOnly=T)
	width.mat <- ncol(mat) * char.height
	nchar.name <- max(nchar(unlist(strsplit(name, '\n'))))
	width.name <- convertWidth(grobHeight(textGrob('A', gp=gpar(fontsize=10))), 'inch', valueOnly=T) * nchar.name
	if (width.name < 0) width.name <- 0
	width.dev <- width.mat + convertWidth(unit(1, 'cm'), 'inch', valueOnly=T) + width.name + 4
	if (! is.null(annotation_names)) {
	    length.title <- max(nchar(annotation_names))
		width.title <- convertWidth(grobHeight(textGrob("A", gp=gpar(fontsize = 10))), 'inch', valueOnly=T) * length.title
		width.dev <- width.dev + width.title
	}
	return(width.dev)
}

.annotation_colors <- function(df)
{
#	library(circlize)
	cols <- c('darkgrey', 'black', 'red', 'yellow', 'blue', 'orange', 'cyan', 'magenta', 'darkgreen', 'khaki')
	annot_cols <- list()
	groups <- NULL
	categorical <- colnames(df)[sapply(df, class) %in% c('character', 'factor')]
	for (annot in categorical) {
		groups <- NULL
		x <- df[,annot]
		if (is.character(x)) {
			groups <- unique(x)
		} else if (is.factor(x)) {
			groups <- levels(x)
		}
		groups <- groups[!is.na(groups)]
		if (length(groups) <= length(cols)) {
			group_cols <- cols[seq(length(groups))]
		} else {
			group_cols <- rainbow(length(groups))
		}
		names(group_cols) <- groups
		annot_cols[[annot]] <- group_cols
	}
	other <- colnames(df)[!colnames(df) %in% categorical]
	cols <- cols[-1]
	numeric_cols <- seq_along(other)
	names(numeric_cols) <- other
#	numeric_cols <- lapply(numeric_cols, function(x) colorRampPalette(c('white', x)))
#	numeric_cols <- lapply(numeric_cols, function(x) colorRamp(c('white', x)))
#	numeric_cols <- lapply(numeric_cols, function(x) colorRamp2(c(0,10), c('white', x)))
	annot_cols <- c(annot_cols, numeric_cols)
	for (i in seq_along(other)) {
		my.col <- cols[(i-1) %% length(cols) + 1]
		my.col.func <- colorRampPalette(c('white', my.col))

		annot_cols[[other[i]]] <- .col_mapping(c('white', my.col), range(df[,other[i]], na.rm=T))
	}
	annot_cols <- annot_cols[colnames(df)]
	return(annot_cols)
}

.col_mapping <- function(colors, breaks)
{
	if (length(breaks) != 2) {
		stop('The vector \'ends\' must have 2 values.')
	}
	ends <- sort(breaks)
	attr <- list(breaks=breaks, colors=colors)
	fun <- function(x) {
		width.breaks <- breaks[2] - breaks[1]
		x[x < breaks[1]] <- breaks[1]
		x[x > breaks[2]] <- breaks[2]
		x <- (x-breaks[1])/width.breaks
		my.col <- colorRamp(colors, alpha=T)(x)
		my.col <- apply(my.col, 1, function(rgb) paste(as.hexmode(as.integer(rgb)), collapse=''))
		my.col <- paste0('#', toupper(my.col))
		my.col[is.na(x)] <- 'FFFFFFFF'
		return(my.col)
	}
	attributes(fun) <- attr
	return(fun)
}
