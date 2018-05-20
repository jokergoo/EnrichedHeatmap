
# == title
# Class for a single heatmap
#
# == details
# The `EnrichedHeatmap-class` is inherited from `ComplexHeatmap::Heatmap-class`.
#
# == methods
# The `EnrichedHeatmap-class` provides following methods:
#
# - `EnrichedHeatmap`: constructor method.
# - `draw,EnrichedHeatmap-method`: draw a single heatmap.
#
# == seealso
# `EnrichedHeatmapList-class`
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
EnrichedHeatmap = setClass("EnrichedHeatmap",
	slots = getClass("Heatmap")@slots,
	contains = "Heatmap")

# == title
# Enriched scores
#
# == param
# -mat a normalized matrix from `normalizeToMatrix`
#
# == details
# The function calculates how the signal is enriched in the target by weighting
# the distance to the target.
#
# For a numeric vector, assume the vector is denoted as combination of three sub-vectors
# ``c(x1, x2, x3)`` with length ``n1``, ``n2`` and ``n3``, 
# where ``x1`` are data points in upstream windows, ``x2`` are data points in target windows and 
# ``x3`` are data points in downstream windows, the enriched score is calcualted as 
#
# sum(x_1i* i/n1) + sum(x_3j* (n3 - j + 1)/n3) + sum(x_2k * abs(n2/2 - abs(k - n2/2)))
#
# where the first two terms are the distance to the start or end position of the target
# by weighting the distance to the position that if it is closer to the start or end position
# of the target, it has higher weight. The second term weight the distance to the center point
# of the target and similar, if it is closer to the center position, it has higher weight.
#
# == value
# A numeric vector
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
enriched_score = function(mat) {

	calc_score = function(x1, x2, x3) {
		n1 = length(x1)
		n2 = length(x2)
		n3 = length(x3)

		x1_index = seq_along(x1)
		x2_index = seq_along(x2)
		x3_index = seq_along(x3)

		l1 = !is.na(x1)
		x1 = x1[l1]
		x1_index = x1_index[l1]

		l2 = !is.na(x2)
		x2 = x2[l2]
		x2_index = x2_index[l2]

		l3 = !is.na(x3)
		x3 = x3[l3]
		x3_index = x3_index[l3]

		if(length(n1) && length(n2)) {
			sum(x1 * x1_index/n1) + 
				sum(x2 * abs(n2/2 - abs(x2_index - n2/2))) + 
				sum(x3 * rev(x3_index)/n3)
		} else if(!length(n1) && length(n2)) {
			sum(x2 * abs(n2/2 - abs(x2_index - n2/2))) + 
				sum(x3 * rev(x3_index)/n3)
		} else if(length(n1) && !length(n2)) {
			sum(x1 * x1_index/n1) + 
				sum(x2 * abs(n2/2 - abs(x2_index - n2/2)))
		} else {
			sum(x2 * abs(n2/2 - abs(x2_index - n2/2)))
		}
	}

	upstream_index = attr(mat, "upstream_index")
	downstream_index = attr(mat, "downstream_index")
	target_index = attr(mat, "target_index")

	score = apply(mat, 1, function(x) {
				x1 = x[upstream_index]
				x2 = x[target_index]
				x3 = x[downstream_index]
				
				calc_score(x1, x2, x3)
			})
	return(score)
}


# == title
# Constructor method for EnrichedHeatmap class
# 
# == param
# -mat a matrix which is returned by `normalizeToMatrix`
# -col color settings. If the signals are categorical, color should be a vector with category levels as names.
# -top_annotation a specific annotation which is always put on top of the enriched heatmap and is constructed by `anno_enriched`
# -top_annotation_height the height of the top annotation
# -row_order row order. Default rows are ordered by enriched scores calculated from `enriched_score`
# -pos_line whether draw vertical lines which represent the positions of ``target``
# -pos_line_gp graphic parameters for the position lines
# -axis_name names for axis which is below the heatmap. If the targets are single points, ``axis_name`` is a vector
#         of length three which corresponds to upstream, target itself and downstream. If the
#         targets are regions with width larger than 1, ``axis_name`` should be a vector of length four which 
#        corresponds to upstream, start of targets, end of targets and downstream.
# -axis_name_rot rotation for axis names
# -axis_name_gp graphic parameters for axis names
# -border whether show border of the heatmap
# -cluster_rows clustering on rows are turned off by default
# -show_row_dend whether show dendrograms on rows if apply hierarchical clustering on rows
# -row_dend_reorder weight for reordering the row dendrogram. It is reordered by enriched scores by default.
# -show_row_names whether show row names
# -heatmap_legend_param a list of settings for heatmap legends. ``at`` and ``labels`` can not be set here.
# -... pass to `ComplexHeatmap::Heatmap`
#
# == details
# `EnrichedHeatmap-class` is inherited from `ComplexHeatmap::Heatmap-class`. Following parameters are 
# set with pre-defined values:
#
# -``cluster_columns`` enforced to be ``FALSE``
# -``show_column_names`` enforced to be ``FALSE``
# -``bottom_annotation`` enforced to be ``NULL`` 
# -``column_title_side`` enforced to be ``top``
#
# A `EnrichedHeatmap-class` object is also a `ComplexHeatmap::Heatmap-class` object, thus, most of the 
# arguments in `ComplexHeatmap::Heatmap` are usable in `EnrichedHeatmap` such as
# to apply clustering on rows, or to split rows by data frame or k-means clustering. Users can also 
# add more than one heatmaps by ``+`` operator. For a detailed demonstration, please go to the vignette.
#
# == value
# An `EnrichedHeatmap-class` object which is inherited from `ComplexHeatmap::Heatmap-class`.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# load(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap"))
# mat3 = normalizeToMatrix(meth, cgi, value_column = "meth", mean_mode = "absolute",
#     extend = 5000, w = 50, smooth = TRUE)
# EnrichedHeatmap(mat3, name = "methylation", column_title = "methylation near CGI")
# EnrichedHeatmap(mat3, name = "meth1") + EnrichedHeatmap(mat3, name = "meth2")
# # for more examples, please go to the vignette
EnrichedHeatmap = function(mat, col, top_annotation = HeatmapAnnotation(enriched = anno_enriched()),
	top_annotation_height = unit(2, "cm"),
	row_order = order(enriched_score(mat), decreasing = TRUE), pos_line = TRUE, 
	pos_line_gp = gpar(lty = 2), axis_name = NULL, axis_name_rot = 0, 
	axis_name_gp = gpar(fontsize = 10), border = TRUE, cluster_rows = FALSE, 
	row_dend_reorder = -enriched_score(mat),
	show_row_dend = FALSE, show_row_names = FALSE, 
	heatmap_legend_param = list(), ...) {

	upstream_index = attr(mat, "upstream_index")
	downstream_index = attr(mat, "downstream_index")
	target_index = attr(mat, "target_index")

	if(is.null(upstream_index) || is.null(downstream_index)) {
		stop("`mat` should be generated by `normalizeToMatrix()`.")
	}
	
	n1 = length(upstream_index)
	n2 = length(target_index)
	n3 = length(downstream_index)
	n = n1 + n2 + n3

	extend = attr(mat, "extend")
	if(is.null(axis_name)) {
		if(n1 && n2 && n3) {
			axis_name = c(paste0("-", extend[1]), "start", "end", extend[2])
		} else if(n1 && !n2 && n3) {
			axis_name = c(paste0("-", extend[1]), "start", extend[2])
		} else if(!n1 && n2 && n3) {
			axis_name = c("start", "end", extend[2])
		} else if(n1 && n2 && !n3) {
			axis_name = c(paste0("-", extend[1]), "start", "end")
		} else if(!n1 && n2 && !n3) {
			axis_name = c("start", "end")
		} else if(n1 && !n2 && !n3) {
			axis_name = c(paste0("-", extend[1]), "start")
		} else if(!n1 && !n2 && n3) {
			axis_name = c("end", extend[2])
		}
	}

	axis_name_rot = axis_name_rot %% 360
	if(axis_name_rot > 90 && axis_name_rot < 270) axis_name_rot = (axis_name_rot + 180) %% 360

	axis_fun = function() {
		grid.lines(c(0.5/n, (n-0.5)/n), c(1, 1))
		if(n1 && n2 && n3) {
			grid.segments(c(0.5/n, (n1+0.5)/n, (n1+n2-0.5)/n, (n-0.5)/n), 
				          unit(1, "npc") - unit(c(1, 1, 1, 1), "mm"), 
				          c(0.5/n, (n1+0.5)/n, (n1+n2-0.5)/n, (n-0.5)/n), 
				          c(1, 1, 1, 1))
			if(axis_name_rot == 0) {
				grid.text(axis_name,
					      c(0.5/n, (n1+0.5)/n, (n1+n2-0.5)/n, (n-0.5)/n),
					      unit(1, "npc") - unit(c(2, 2, 2, 2), "mm"), gp = axis_name_gp,
					      hjust = c(0, 0.5, 0.5, 1), vjust = 1)
			} else {
				if(axis_name_rot > 0 & axis_name_rot <= 90) {
					hjust = c(1, 1, 1, 1)
					vjust = 0.5
				} else {
					hjust = c(0, 0, 0, 0)
					vjust = 0.5
				}
				grid.text(axis_name,
					      c(0.5/n, (n1+0.5)/n, (n1+n2-0.5)/n, (n-0.5)/n),
					      unit(1, "npc") - unit(c(2, 2, 2, 2), "mm"), gp = axis_name_gp, rot = axis_name_rot,
					      hjust = hjust, vjust = vjust)
			}
		} else if(n1 && !n2 && n3) {
			grid.segments(c(0.5/n, (n1+0.5)/n, (n-0.5)/n), 
				          unit(1, "npc") - unit(c(1, 1, 1), "mm"), 
				          c(0.5/n, (n1+0.5)/n, (n-0.5)/n), 
				          c(1, 1, 1))
			if(axis_name_rot == 0) {
				grid.text(axis_name,
					      c(0.5/n, (n1+0.5)/n, (n-0.5)/n),
					      unit(1, "npc") - unit(c(2, 2, 2), "mm"), gp = axis_name_gp,
					      hjust = c(0, 0.5, 1), vjust = 1)
			} else {
				if(axis_name_rot > 0 & axis_name_rot <= 90) {
					hjust = c(1, 1, 1)
					vjust = 0.5
				} else {
					hjust = c(0, 0, 0)
					vjust = 0.5
				}
				grid.text(axis_name,
					      c(0.5/n, (n1+0.5)/n, (n-0.5)/n),
					      unit(1, "npc") - unit(c(2, 2, 2), "mm"), gp = axis_name_gp, rot = axis_name_rot,
					      hjust = hjust, vjust = vjust)
			}
		} else if(!n1 && n2 && n3) {
			grid.segments(c(0.5/n, (n1+n2-0.5)/n, (n-0.5)/n), 
				          unit(1, "npc") - unit(c(1, 1, 1), "mm"), 
				          c(0.5/n, (n1+n2-0.5)/n, (n-0.5)/n), 
				          c(1, 1, 1))
			if(axis_name_rot == 0) {
				grid.text(axis_name,
					      c(0.5/n, (n1+n2-0.5)/n, (n-0.5)/n),
					      unit(1, "npc") - unit(c(2, 2, 2), "mm"), gp = axis_name_gp,
					      hjust = c(0, 0.5, 1), vjust = 1)
			} else {
				if(axis_name_rot > 0 & axis_name_rot <= 90) {
					hjust = c(1, 1, 1)
					vjust = 0.5
				} else {
					hjust = c(0, 0, 0)
					vjust = 0.5
				}
				grid.text(axis_name,
					      c(0.5/n, (n1+n2-0.5)/n, (n-0.5)/n),
					      unit(1, "npc") - unit(c(2, 2, 2), "mm"), gp = axis_name_gp, rot = axis_name_rot,
					      hjust = hjust, vjust = vjust)
			}
		} else if(n1 && n2 && !n3) {
			grid.segments(c(0.5/n, (n1+0.5)/n, (n1+n2-0.5)/n), 
				          unit(1, "npc") - unit(c(1, 1, 1), "mm"), 
				          c(0.5/n, (n1+0.5)/n, (n1+n2-0.5)/n), 
				          c(1, 1, 1))
			if(axis_name_rot == 0) {
				grid.text(axis_name,
					      c(0.5/n, (n1+0.5)/n, (n1+n2-0.5)/n),
					      unit(1, "npc") - unit(c(2, 2, 2), "mm"), gp = axis_name_gp,
					      hjust = c(0, 0.5, 1), vjust = 1)
			} else {
				if(axis_name_rot > 0 & axis_name_rot <= 90) {
					hjust = c(1, 1, 1)
					vjust = 0.5
				} else {
					hjust = c(0, 0, 0)
					vjust = 0.5
				}
				grid.text(axis_name,
					      c(0.5/n, (n1+0.5)/n, (n1+n2-0.5)/n),
					      unit(1, "npc") - unit(c(2, 2, 2), "mm"), gp = axis_name_gp, rot = axis_name_rot,
					      hjust = hjust, vjust = vjust)
			}
		} else {
			grid.segments(c(0.5/n, (n-0.5)/n), 
				          unit(1, "npc") - unit(c(1, 1), "mm"), 
				          c(0.5/n, (n-0.5)/n), 
				          c(1, 1))
			if(axis_name_rot == 0) {
				grid.text(axis_name,
					      c(0.5/n, (n-0.5)/n),
					      unit(1, "npc") - unit(c(2, 2), "mm"), gp = axis_name_gp,
					      hjust = c(0, 1), vjust = 1)
			} else {
				if(axis_name_rot > 0 & axis_name_rot <= 90) {
					hjust = c(1, 1)
					vjust = 0.5
				} else {
					hjust = c(0, 0)
					vjust = 0.5
				}
				grid.text(axis_name,
					      c(0.5/n, (n-0.5)/n),
					      unit(1, "npc") - unit(c(2, 2), "mm"), gp = axis_name_gp, rot = axis_name_rot,
					      hjust = hjust, vjust = vjust)
			}
		}
		minor_ticks = calc_minor_ticks(mat)
		if(length(minor_ticks)) {
	        grid.segments(minor_ticks, unit(1, "npc") - unit(0.5, "mm"), minor_ticks, 1)
	    }
	}
	
	if(axis_name_rot == 0) {
		axis_height = grobHeight(textGrob("a")) + unit(4, "mm")
	} else {
	 	axis_height = max(grobWidth(textGrob(axis_name, gp = gpar(axis_name_gp))))*abs(sin(axis_name_rot/180*pi)) + unit(4, "mm")
	}

	class(mat) = NULL

	signal_is_categorical = attr(mat, "signal_is_categorical")
	if(is.null(signal_is_categorical)) signal_is_categorical = FALSE
	if(signal_is_categorical) {
		signal_level = attr(mat, "signal_level")
		n_level = length(signal_level)
		if(missing(col)) {
			default_col = colorRamp2(c(1, (1+n_level)/2, n_level), c("blue", "white", "red"))(1:n_level)
			col = structure(default_col, names = 1:n_level)
			col = c("0" = "white", col)
		} else {
			if(is.null(names(col))) {
				if(length(col) != n_level) {
					stop("If `col` has no name, it must have same length as number of levels in signals.")
				} else {
					col = structure(col, names = 1:n_level)
				}
			} else {
				if(!identical(sort(names(col)), sort(signal_level))) {
					stop("Names of `col` should be same as the levels in signals.")
				} else {
					col = structure(col[signal_level], names = 1:n_level)
				}
			}
			col = c("0" = "white", col)
		}
		if(is.null(heatmap_legend_param$at)) {
			heatmap_legend_param$at = 1:n_level
		}
		if(is.null(heatmap_legend_param$labels)) {
			heatmap_legend_param$labels = signal_level
		}
	}
	
	ht = Heatmap(mat, col, row_order = row_order, cluster_columns = FALSE, cluster_rows = cluster_rows,
			show_row_names = show_row_names, show_column_names = FALSE, bottom_annotation = NULL, 
			column_title_side = "top", show_row_dend = show_row_dend, row_dend_reorder = row_dend_reorder,
			top_annotation = top_annotation, top_annotation_height = top_annotation_height, 
			heatmap_legend_param = heatmap_legend_param, ...)

	# additional parameters specific for `EnrichedHeatmap` class
	ht@heatmap_param$pos_line = pos_line
	ht@heatmap_param$pos_line_gp = pos_line_gp
	ht@heatmap_param$axis_fun = axis_fun
	ht@heatmap_param$axis_height = axis_height
	ht@heatmap_param$border = border

	return(changeClassName(ht, "EnrichedHeatmap"))
}

# == title
# Draw the single heatmap with default parameters
#
# == param
# -object an `EnrichedHeatmap-class` object.
#
# == details
# Actually it calls `draw,EnrichedHeatmap-method`, but only with default parameters. If users want to customize the heatmap,
# they can pass parameters directly to `draw,EnrichedHeatmap-method`.
#
# == value
# An `EnrichedHeatmapList-class` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# # see documentation of EnrichedHeatmap
# NULL
setMethod(f = "show",
	signature = "EnrichedHeatmap",
	definition = function(object) {

	# `draw` method is inherited from `Heatmap` class
	draw(object)

})

# == title
# Draw a single heatmap
#
# == param
# -object an `EnrichedHeatmap-class` object.
# -internal only used internally.
# -... pass to `ComplexHeatmap::draw`,HeatmapList-method.
#
# == detail
# The function creates an `EnrichedHeatmapList-class` object which only contains a single heatmap
# and call `draw,EnrichedHeatmapList-method` to make the final heatmap.
#
# == value
# An `EnrichedHeatmapList-class` object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# # see documentation of EnrichedHeatmap
# NULL
setMethod(f = "draw",
	signature = "EnrichedHeatmap",
	definition = function(object, internal = FALSE, ...) {
	
	if(internal) {
		object = changeClassName(object, "Heatmap")
		draw(object, internal = internal)
	} else {
		ht_list = new("HeatmapList")
	    ht_list = ht_list + object
	    draw(ht_list, ...)
	}
})


# == title
# Annotation function to show the enrichment
#
# == param
# -gp graphic parameters. There are two non-standard parameters: ``neg_col`` and ``pos_col``. 
#     If these two parameters are defined, the positive signals and negatie signals are visualized separatedly.
#     The graphic parameters can be set as vectors when the heatmap or heatmap list is split into several row clusters.
# -pos_line whether to draw vertical lines which represent positions of ``target``
# -pos_line_gp graphic parameters for the position lines
# -yaxis whether show yaxis
# -ylim ranges on y-axis, by default it is inferred from the data
# -value the method to summarize signals from columns of the noramlized matrix
# -yaxis_side side of y-axis
# -yaxis_facing facing of the axis ticks and labels. It can be set to avoid overlapping text when multiple
#               heatmaps are plotted together
# -yaxis_gp graphic parameters for y-axis
# -show_error whether show error regions which are one standard error to the mean value. Color of error
#            area is same as the corresponding lines with 75 percent transparency.
#
# == details
# This annotation functions shows mean values (or depends on the method set in ``value`` argument) of columns in the normalized matrix
# which summarises the enrichment of the signals to the targets.
#
# If rows are splitted, the enriched lines are calculated for each row cluster and there will also be multiple lines in this annotation viewport.
#
# It should only be placed as column annotation of the enriched heatmap.
#
# == values
# A column annotation function which can be set to ``top_annotation`` argument in `EnrichedHeatmap`.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# load(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap"))
# tss = promoters(genes, upstream = 0, downstream = 1)
# mat1 = normalizeToMatrix(H3K4me3, tss, value_column = "coverage", 
#     extend = 5000, mean_mode = "w0", w = 50, keep = c(0, 0.99))
# EnrichedHeatmap(mat1, col = c("white", "red"), name = "H3K4me3",
#     top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
#     km = 3, row_title_rot = 0)
#
anno_enriched = function(gp = gpar(col = "red"), pos_line = TRUE, pos_line_gp = gpar(lty = 2),
	yaxis = TRUE, ylim = NULL, value = c("mean", "sum", "abs_mean", "abs_sum"), 
	yaxis_side = "right", yaxis_facing = ifelse(yaxis_side == "right", "right", "left"), 
	yaxis_gp = gpar(fontsize = 8), show_error = FALSE) {

	# in case of lazy loading
	gp = gp
	pos_line = pos_line
	pos_line_gp = pos_line_gp
	yaxis = yaxis
	ylim = ylim
	yaxis_side = yaxis_side
	yaxis_facing = yaxis_facing
	yaxis_gp = yaxis_gp
	show_error = show_error

	by_sign = FALSE
	if(!is.null(gp$neg_col) && !is.null(gp$pos_col)) {
		by_sign = TRUE
		show_error = FALSE
	} else if(!is.null(gp$neg_col) && is.null(gp$pos_col)) {
		stop("Since you defined `neg_col` in `gp`, you should also define `pos_col`.")
	} else if(is.null(gp$neg_col) && !is.null(gp$pos_col)) {
		stop("Since you defined `pos_col` in `gp`, you should also define `neg_col`.")
	}

	value = match.arg(value)[1]
	function(index) {

		ht = get("object", envir = parent.frame(n = 5))
		mat = ht@matrix
		
		if(by_sign) {
			mat_pos = mat
			mat_pos[mat_pos < 0] = 0
			mat_neg = mat
			mat_neg[mat_neg > 0] = 0
			mat_pos = abs(mat_pos)
			mat_neg = abs(mat_neg)
		}

		upstream_index = attr(mat, "upstream_index")
		downstream_index = attr(mat, "downstream_index")
		target_index = attr(mat, "target_index")
		signal_is_categorical = attr(mat, "signal_is_categorical")
		signal_level = attr(mat, "signal_level")
		if(is.null(signal_is_categorical)) signal_is_categorical = FALSE
	
		n1 = length(upstream_index)
		n2 = length(target_index)
		n3 = length(downstream_index)
		n = n1 + n2 + n3

		if(signal_is_categorical) {
			mat_list = lapply(seq_along(signal_level), function(i) {
				mat0 = mat
				mat0[mat0 != i] = 0
				mat0[mat0 > 0] = 1
				mat0
			})
			yl = lapply(ht@row_order_list, function(i) {
				if(value == "sum" || value == "abs_sum") {
					y = sapply(mat_list, function(m) {
						colSums(m[i, , drop = FALSE], na.rm = TRUE)
					})
					
				} else if(value == "mean" || value == "abs_mean") {
					y = sapply(mat_list, function(m) {
						colMeans(m[i, , drop = FALSE], na.rm = TRUE)
					})
				}
				y
			})
			show_error = FALSE
		} else if(by_sign) {
			if(value == "sum" || value == "abs_sum") {
				y_pos = sapply(ht@row_order_list, function(i) {
					colSums(mat_pos[i, , drop = FALSE], na.rm = TRUE)
				})
				y_neg = sapply(ht@row_order_list, function(i) {
					colSums(mat_neg[i, , drop = FALSE], na.rm = TRUE)
				})
				show_error = FALSE
			} else if(value == "mean" || value == "abs_mean") {
				y_pos = sapply(ht@row_order_list, function(i) {
					colMeans(mat_pos[i, , drop = FALSE], na.rm = TRUE)
				})
				y_neg = sapply(ht@row_order_list, function(i) {
					colMeans(mat_neg[i, , drop = FALSE], na.rm = TRUE)
				})
			}
		} else {
			if(value == "sum") {
				y = sapply(ht@row_order_list, function(i) {
					colSums(mat[i, , drop = FALSE], na.rm = TRUE)
				})
				show_error = FALSE
			} else if(value == "abs_sum") {
				y = sapply(ht@row_order_list, function(i) {
					colSums(abs(mat[i, , drop = FALSE]), na.rm = TRUE)
				})
				show_error = FALSE
			} else if(value == "mean") {
				y = sapply(ht@row_order_list, function(i) {
					colMeans(mat[i, , drop = FALSE], na.rm = TRUE)
				})
			} else if(value == "abs_mean") {
				y = sapply(ht@row_order_list, function(i) {
					colMeans(abs(mat[i, , drop = FALSE]), na.rm = TRUE)
				})
			}
		}

		if(show_error) {
			y_se = sapply(ht@row_order_list, function(i) {
				colSds(mat[i, , drop = FALSE], na.rm = TRUE)/sqrt(length(i))
			})
			if(is.null(ylim)) {
				ylim = range(c(y+y_se, y-y_se), na.rm = TRUE)
				ylim[2] = ylim[2] + (ylim[2] - ylim[1]) * 0.05
			}
		} else {
			if(is.null(ylim)) {
				if(signal_is_categorical) {
					ylim = range(unlist(yl), na.rm = TRUE)
				} else {
					if(by_sign) {
						ylim = range(c(y_pos, y_neg), na.rm = TRUE)
					} else {
						ylim = range(y, na.rm = TRUE)
					}
				}
				ylim[2] = ylim[2] + (ylim[2] - ylim[1]) * 0.05
			}
		}

		pushViewport(viewport(xscale = c(0, n), yscale = ylim))
		grid.rect(gp = gpar(col = "black", fill = NA))
		if(signal_is_categorical) {
			gp = recycle_gp(gp, length(ht@row_order_list))
			mat_col = ht@matrix_color_mapping@colors
			for(k in seq_along(yl)) {
				y = yl[[k]]
				for(i in seq_len(ncol(y))) {
					gp2 = subset_gp(gp, k)
					gp2$col = mat_col[as.character(i)]
					grid.lines(seq_len(n)-0.5, y[,i], default.units = "native", gp = gp2)
				}
			}
		} else {
			gp = recycle_gp(gp, length(ht@row_order_list))
			if(by_sign) {
				for(i in seq_len(ncol(y_pos))) {
					gp2 = subset_gp(gp, i); gp2$col = gp2$pos_col
					grid.lines(seq_len(n)-0.5, y_pos[,i], default.units = "native", gp = gp2)
					gp2 = subset_gp(gp, i); gp2$col = gp2$neg_col
					grid.lines(seq_len(n)-0.5, y_neg[,i], default.units = "native", gp = gp2)
				}
			} else {
				for(i in seq_len(ncol(y))) {
					if(show_error) {
						line_col = col2rgb(subset_gp(gp, i)$col, alpha = TRUE)[, 1]
						line_col[4] = floor(line_col[4]*0.25)
						grid.polygon(c(seq_len(n)-0.5, rev(seq_len(n)-0.5)), c(y[,i]+y_se[,i], rev(y[,i]-y_se[,i])), 
							default.units = "native", gp = gpar(col = NA, fill = rgb(line_col[1], line_col[2], line_col[3], line_col[4], maxColorValue = 255)))
					}
					grid.lines(seq_len(n)-0.5, y[,i], default.units = "native", gp = subset_gp(gp, i))
				}
			}
		}

		if(pos_line) {
		    if(n1 && n2 && n3) {
                grid.lines(rep((n1+0.5)/n, 2), c(0, 1), gp = pos_line_gp)
                grid.lines(rep((n1+n2-0.5)/n, 2), c(0, 1), gp = pos_line_gp)
            } else if(n1 && !n2 && n3) {
                grid.lines(rep((n1+0.5)/n, 2), c(0, 1), gp = pos_line_gp)
            } else if(!n1 && n2 && n3) {
                grid.lines(rep((n1+n2-0.5)/n, 2), c(0, 1), gp = pos_line_gp)
            } else if(n1 && n2 && !n3) {
                grid.lines(rep((n1-0.5)/n, 2), c(0, 1), gp = pos_line_gp)
            }
		}
		if(yaxis) {
			le1 = grid.pretty(ylim)
			le2 = pretty(ylim, n = 3)
			if(abs(length(le1) - 5) < abs(length(le2) - 5)) {
				le = le1
			} else {
				le = le2
			}
			breaks = le[le >= ylim[1] & le <= ylim[2]]
			n_break = length(breaks)
			if(yaxis_side == "right") {
				if(yaxis_facing == "right") {
					grid.segments(unit(1, "npc"), breaks, unit(1, "npc") + unit(1, "mm"), breaks, default.units = "native", gp = yaxis_gp)
					grid.text(breaks, unit(1, "npc") + unit(2, "mm"), breaks, default.units = "native", gp = yaxis_gp, just = "left")
				} else {
					grid.segments(unit(1, "npc"), breaks, unit(1, "npc") - unit(1, "mm"), breaks, default.units = "native", gp = yaxis_gp)
					grid.text(breaks, unit(1, "npc") - unit(2, "mm"), breaks, default.units = "native", gp = yaxis_gp, just = "right")
				}
			} else {
				if(yaxis_facing == "right") {
					grid.segments(unit(0, "npc"), breaks, unit(0, "npc") + unit(1, "mm"), breaks, default.units = "native", gp = yaxis_gp)
					grid.text(breaks, unit(0, "npc") + unit(2, "mm"), breaks, default.units = "native", gp = yaxis_gp, just = "left")
				} else {
					grid.segments(unit(0, "npc"), breaks, unit(0, "npc") - unit(1, "mm"), breaks, default.units = "native", gp = yaxis_gp)
					grid.text(breaks, unit(0, "npc") - unit(2, "mm"), breaks, default.units = "native", gp = yaxis_gp, just = "right")
				}
			}
		}

		# minor_ticks = calc_minor_ticks(mat)
		# if(length(minor_ticks)) {
	 #        grid.segments(minor_ticks, unit(0.5, "mm"), minor_ticks, 0)
	 #    }

	    upViewport()
	}
}

recycle_gp = function(gp, n = 1) {
	g = lapply(gp, function(x) {
		if(length(x) == 1 && n > 1) {
			rep(x, n)
		} else {
			x
		}
	})
	class(g) = "gpar"
	return(g)
}

subset_gp = function(gp, i = 1) {
	g = lapply(gp, function(x) {
		x[i]
	})
	class(g) = "gpar"
	return(g)
}
