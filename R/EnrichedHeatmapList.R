
# == title
# Class for a list of heatmaps
#
# == details
# The structure of `EnrichedHeatmapList-class` is the same as
# `ComplexHeatmap::HeatmapList-class` and the class is inherited from `ComplexHeatmap::HeatmapList-class`.
#
# == methods
# The `EnrichedHeatmapList-class` provides following methods:
#
# - `draw,EnrichedHeatmapList-method`: draw a list of heatmaps.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
EnrichedHeatmapList = setClass("EnrichedHeatmapList",
	slots = getClass("HeatmapList")@slots,
	contains = "HeatmapList")

# == title
# Constructor method for EnrichedHeatmapList class
#
# == param
# -... arguments
#
# == details
# There is no public constructor method for the `EnrichedHeatmapList-class`.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# # no example
# NULL
EnrichedHeatmapList = function(...) {
    new("EnrichedHeatmapList", ...)
}

# == title
# Add heatmaps or row annotations to a heatmap list
#
# == param
# -x an `EnrichedHeatmap-class` object, a `Heatmap-class` object, a `HeatmapAnnotation-class` object or a `HeatmapList-class` object.
# -y an `EnrichedHeatmap-class` object, a `Heatmap-class` object, a `HeatmapAnnotation-class` object or a `HeatmapList-class` object.
#
# == detail
# It overwrites ``+.AdditiveUnit`` in the ComplexHeatmap package.
#
# == value
# A `ComplexHeatmap::HeatmapList-class` object or an `EnrichedHeatmapList-class` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# # users should not use it directly
# NULL
"+.AdditiveUnit" = function(x, y) {
    if(inherits(x, "EnrichedHeatmap") || 
       inherits(x, "EnrichedHeatmapList") ||
       inherits(y, "EnrichedHeatmap") ||
       inherits(y, "EnrichedHeatmapList")) {
    	
    	# should return a `EnrichedHeatmapList` object
        if(is.null(x)) {
            ht_list = new("HeatmapList")
            ht_list = add_heatmap(ht_list, y)
        } else if(is.null(y)) {
            ht_list = new("HeatmapList")
            ht_list = add_heatmap(ht_list, x)
        } else {
    	   ht_list = add_heatmap(x, y)
        }
        changeClassName(ht_list, "EnrichedHeatmapList")
    } else {
    	ComplexHeatmap::`+.AdditiveUnit`(x, y)
    }
}

# == title
# Draw a list of heatmaps with default parameters
#
# == param
# -object an `EnrichedHeatmapList-class` object.
#
# == details
# Actually it calls `draw,EnrichedHeatmapList-method`, but only with default parameters. If users want to customize the heatmap,
# they can pass parameters directly to `draw,EnrichedHeatmapList-method`.
#
# == value
# An `EnrichedHeatmapList-class` object.
#
# == example
# # see documentation of EnrichedHeatmap
# NULL
setMethod(f = "show",
	signature = "EnrichedHeatmapList",
	definition = function(object) {

	# re-define `draw` for HeatmapList
	draw(object)

})

# == title
# Draw a list of heatmaps
#
# == param
# -object an `EnrichedHeatmapList-class` object
# -padding padding of the plot. Elements correspond to bottom, left, top, right paddings.
# -... pass to `ComplexHeatmap::make_layout,HeatmapList-method`
# -newpage whether to create a new page
#
# == detail
# It calls `ComplexHeatmap::draw,HeatmapList-method` to make the plot but with some adjustment
# specificly for enriched heatmaps.
#
# == value
# An `EnrichedHeatmapList` object
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# # see documentation of EnrichedHeatmap
# NULL
setMethod(f = "draw",
    signature = "EnrichedHeatmapList",
    definition = function(object, padding = unit(c(2, 2, 2, 2), "mm"), ..., newpage= TRUE) {

    if(newpage) {
        grid.newpage()
    }

    object = make_layout(object, ...)

    # which heatmap are EnrichedHeatmap
	Enriched_heatmap_index = which(sapply(object@ht_list, function(ht) {
		inherits(ht, "EnrichedHeatmap")
	}))

	max_axis_height = max(do.call("unit.c", lapply(object@ht_list[Enriched_heatmap_index], function(ht) {
		ht@heatmap_param$axis_height
	})))
    
    normal_heatmap_index = which(sapply(object@ht_list, function(ht) {
		inherits(ht, "Heatmap") & !inherits(ht, "EnrichedHeatmap")
	}))
    
    if(length(normal_heatmap_index) == 0) {
    	# if all heatmaps are Enriched heatmaps, the just put the
    	# axis in the bottom colnames component
    	for(i in Enriched_heatmap_index) {
    		ht = object@ht_list[[i]]

    		# change the size of bottom_annotation
    		ht@layout$layout_column_names_bottom_height = max_axis_height
            ht@layout$layout_index = rbind(ht@layout$layout_index, c(6, 4))
            ht@layout$graphic_fun_list = c(ht@layout$graphic_fun_list, function(ht) ht@heatmap_param$axis_fun())

            object@ht_list[[i]] = ht
    	}
    } else {
    	i = normal_heatmap_index[1]
    	# bottom_height are all fixed length
    	bottom_height = sum(component_height(object@ht_list[[i]], 6:9)) # column_names, annotation, dendrogram and title

    	# assume nothing is allowed to plotted below the Enriched heatmap
    	if(compare_height(bottom_height, max_axis_height) == -1) {

    		for(i in normal_heatmap_index) {
    			ht = object@ht_list[[i]]
    			ht@layout$layout_column_title_bottom_height = ht@layout$layout_column_title_bottom_height + bottom_height - max_axis_height
    			object@ht_list[[i]] = ht
    		}

    	}

        # axis is in the bottom colnames component
        for(i in Enriched_heatmap_index) {
            ht = object@ht_list[[i]]
            ht@layout$layout_index = rbind(ht@layout$layout_index, c(6, 4))
            ht@layout$graphic_fun_list = c(ht@layout$graphic_fun_list, function(ht) {
                pushViewport(viewport(name = paste0(ht@name, "_axis"), y = unit(1, "npc"), 
                    height = ht@heatmap_param$axis_height, just = "top"))
                ht@heatmap_param$axis_fun()
                upViewport()
            })
            object@ht_list[[i]] = ht
        }
    }

    object = selectMethod("draw", "HeatmapList")(object, padding = padding, newpage = FALSE, ...)
    changeClassName(object, "EnrichedHeatmapList")

    # add borders and pos lines
    for(i in Enriched_heatmap_index) {
        ht = object@ht_list[[i]]
        heatmap_name = ht@name

        upstream_index = attr(ht@matrix, "upstream_index")
        downstream_index = attr(ht@matrix, "downstream_index")
        target_index = attr(ht@matrix, "target_index")
        n1 = length(upstream_index)
        n2 = length(target_index)
        n3 = length(downstream_index)
        n = n1 + n2 + n3

        for(k in seq_along(ht@row_order_list)) {
            if(ht@heatmap_param$border) {
                decorate_heatmap_body(heatmap_name, {
                    grid.rect(gp = gpar(col = "black", fill = NA))
                }, slice = k)
            }
            if(ht@heatmap_param$pos_line) {
                decorate_heatmap_body(heatmap_name, {
                    grid.lines(rep((n1-0.5)/n, 2), c(0, 1), gp = ht@heatmap_param$pos_line_gp)
                    if(n2) grid.lines(rep((n1+n2-0.5)/n, 2), c(0, 1), gp = ht@heatmap_param$pos_line_gp)
                }, slice = k)
            }
        }
    }

    return(invisible(object))
})

