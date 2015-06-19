
CentralizedHeatmap = setClass("CentralizedHeatmap",
	slots = getClass("Heatmap")@slots,
	contains = "Heatmap")

CentralizedHeatmap = function(mat, pos_line = TRUE, pos_line_gp = gpar(lty = 2), 
	axis_name = NULL, axis_name_rot = 0, axis_name_gp = gpar(), ...) {

	upstream_index = attr(mat, "upstream_index")
	downstream_index = attr(mat, "downstream_index")
	body_index = attr(mat, "body_index")

	if(is.null(upstream_index) || is.null(downstream_index)) {
		stop("wrong `mat`.")
	}

	score = apply(mat, 1, function(x) {
			x1 = x[upstream_index]
			x2 = x[body_index]
			x3 = x[downstream_index]
			n1 = length(x1)
			n2 = length(x2)
			n3 = length(x3)

			if(n2) {
				sum(x1 * seq_along(x1)/n1) + sum(x2 * abs(n2/2 - abs(seq_along(x2) - n2/2))) + sum(x3 * rev(seq_along(x3))/n3)
			} else {
				sum(x1 * seq_along(x1)/n1) + sum(x3 * rev(seq_along(x3))/n3)
			}
		})

	od = order(score, decreasing = TRUE)

	n1 = length(upstream_index)
	n2 = length(body_index)
	n3 = length(downstream_index)
	n = n1 + n2 + n3

	extend = attr(mat, "extend")
	if(is.null(axis_name)) {
		if(n2) {
			axis_name = c(paste0("-", extend[1]), "start", "end", extend[2])
		} else {
			axis_name = c(paste0("-", extend[1]), "start", extend[2])
		}
	}

	axis_name_rot = axis_name_rot %% 360
	if(axis_name_rot > 90 && axis_name_rot < 270) axis_name_rot = (axis_name_rot + 180) %% 360

	axis_fun = function() {
		
		grid.lines(c(0, 1), c(1, 1))
		if(n2) {
			grid.segments(c(0, n1/n, (n1+n2)/n, 1), 
				          unit(1, "npc") - unit(c(1, 1, 1, 1), "mm"), 
				          c(0, n1/n, (n1+n2)/n, 1), 
				          c(1, 1, 1, 1))
			if(axis_name_rot == 0) {
				grid.text(axis_name,
					      c(0, n1/n, (n1+n2)/n, 1),
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
					      c(0, n1/n, (n1+n2)/n, 1),
					      unit(1, "npc") - unit(c(2, 2, 2, 2), "mm"), gp = axis_name_gp, rot = axis_name_rot,
					      hjust = hjust, vjust = vjust)
			}
		} else {
			grid.segments(c(0, n1/n, 1), 
				          unit(1, "npc") - unit(c(1, 1, 1), "mm"), 
				          c(0, n1/n, 1), 
				          c(1, 1, 1))
			if(axis_name_rot == 0) {
				grid.text(axis_name,
					      c(0, n1/n, 1),
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
					      c(0, n1/n, 1),
					      unit(1, "npc") - unit(c(2, 2, 2), "mm"), gp = axis_name_gp, rot = axis_name_rot,
					      hjust = hjust, vjust = vjust)
			}
		}
	}
	
	if(axis_name_rot == 0) {
		axis_height = grobHeight(textGrob("a")) + unit(4, "mm")
	} else {
	 	axis_height = max(grobWidth(textGrob(axis_name, gp = gpar(axis_name_gp))))*abs(sin(axis_name_rot/180*pi)) + unit(4, "mm")
	}

	ht = Heatmap(mat, row_order = od, cluster_columns = FALSE, cluster_rows = FALSE,
			show_row_names = FALSE, show_column_names = FALSE, ...)

	# additional parameters specific for `CentralizedHeatmap` class
	ht@heatmap_param$pos_line = pos_line
	ht@heatmap_param$pos_line_gp = pos_line_gp
	ht@heatmap_param$axis_fun = axis_fun
	ht@heatmap_param$axis_height = axis_height

	return(changeClassName(ht, "CentralizedHeatmap"))
}

setMethod(f = "show",
	signature = "CentralizedHeatmap",
	definition = function(object) {

	# `draw` method is inherited from `Heatmap` class
	draw(object)

})

setMethod(f = "draw",
	signature = "CentralizedHeatmap",
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
