
# == title
# Class for a list of heatmaps
#
# == details
# The structure of `CentralizedHeatmapList-class` is the same as
# `HeatmapList-class` and the class is inherited from `HeatmapList-class`.
# Then it allows to support all features that are provided by `HeatmapList-class`.
# Also some functions for `HeatmapList-class` are overwitten to provide adjustment
# specifically for centralized heatmaps.
#
# == methods
# The `CentralizedHeatmapList-class` provides following methods:
#
# - `Heatmap`: constructor method.
# - `draw,Heatmap-method`: draw a single heatmap.
# - `add_heatmap,Heatmap-method` append heatmaps and row annotations to a list of heatmaps.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
CentralizedHeatmapList = setClass("CentralizedHeatmapList",
	slots = getClass("HeatmapList")@slots,
	contains = "HeatmapList")

# == title
# Constructor method for CentralizedHeatmapList class
#
# == param
# -... arguments
#
# == details
# There is no public constructor method for the `CentralizedHeatmapList-class`.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
CentralizedHeatmapList = function(...) {
    new("CentralizedHeatmapList", ...)
}

# == title
# 
"+.AdditiveUnit" = function(x, y) {
    if(inherits(x, "CentralizedHeatmap") || 
       inherits(x, "CentralizedHeatmapList") ||
       inherits(y, "CentralizedHeatmap") ||
       inherits(y, "CentralizedHeatmapList")) {
    	
    	# should return a `CentralizedHeatmapList` object
    	ht_list = add_heatmap(x, y)
        changeClassName(ht_list, "CentralizedHeatmapList")
    } else {
    	ComplexHeatmap::`+.AdditiveUnit`(x, y)
    }
}

setMethod(f = "show",
	signature = "CentralizedHeatmapList",
	definition = function(object) {

	# re-define `draw` for HeatmapList
	draw(object)

})

# overwrite the `draw` method for `HeatmapList` class
setMethod(f = "draw",
    signature = "CentralizedHeatmapList",
    definition = function(object, padding = unit(c(2, 2, 2, 2), "mm"), ..., newpage= TRUE) {

    if(newpage) {
        grid.newpage()
    }

    object = make_layout(object, ...)

    # which heatmap are CentralizedHeatmap
	centralized_heatmap_index = which(sapply(object@ht_list, function(ht) {
		inherits(ht, "CentralizedHeatmap")
	}))

	max_axis_height = max(do.call("unit.c", lapply(object@ht_list[centralized_heatmap_index], function(ht) {
		ht@heatmap_param$axis_height
	})))
    
    normal_heatmap_index = which(sapply(object@ht_list, function(ht) {
		inherits(ht, "Heatmap") & !inherits(ht, "CentralizedHeatmap")
	}))
    
    if(length(normal_heatmap_index) == 0) {
    	# if all heatmaps are centralized heatmaps, the just put the
    	# axis in the bottom annotation component
    	for(i in centralized_heatmap_index) {
    		ht = object@ht_list[[i]]

    		# change the size of bottom_annotation
    		ht@layout$layout_column_names_bottom_height = max_axis_height
            ht@layout$layout_index = rbind(ht@layout$layout_index, c(6, 4))
            ht@layout$graphic_fun_list = c(ht@layout$graphic_fun_list, function(object) ht@heatmap_param$axis_fun())

            object@ht_list[[i]] = ht
    	}
    } else {
    	i = normal_heatmap_index[1]
    	# bottom_height are all fixed length
    	bottom_height = component_height(object@ht_list[[i]], 6:9) # column_names, annotation, dendrogram and title

    	# assume nothing is allowed to plotted below the centralized heatmap
    	if(compare_unit(bottom_height, max_axis_height) == -1) {

    		for(i in normal_heatmap_index) {
    			ht = object@ht_list[[i]]
    			ht@layout$layout_column_title_bottom_height = ht@layout$layout_column_title_bottom_height + bottom_height - max_axis_height
    			object@ht_list[[i]] = ht
    		}

    		# add in column rownames component
    		for(i in centralized_heatmap_index) {
    			ht = object@ht_list[[i]]
    			ht@layout$layout_index = rbind(ht@layout$layout_index, c(6, 4))
            	ht@layout$graphic_fun_list = c(ht@layout$graphic_fun_list, function(object) {
	    			pushViewport(viewport(name = paste0(ht@name, "_axis"), y = unit(1, "npc"), 
	    				height = ht@heatmap_param$axis_height, just = "top"))
	    			ht@heatmap_param$axis_fun()
	    			upViewport()
	    		})
	    		object@ht_list[[i]] = ht
    		}
    	}
    }

    if(length(padding) == 1) {
        padding = rep(padding, 4)
    } else if(length(padding) == 2) {
        padding = rep(padding, 2)
    } else if(length(padding) != 4) {
        stop("`padding` can only have length of 1, 2, 4")
    }

    layout = grid.layout(nrow = 7, ncol = 7, widths = component_width(object, 1:7), heights = component_height(object, 1:7))
    pushViewport(viewport(layout = layout, name = "global", width = unit(1, "npc") - padding[2] - padding[4],
        height = unit(1, "npc") - padding[1] - padding[3]))
    ht_layout_index = object@layout$layout_index
    ht_graphic_fun_list = object@layout$graphic_fun_list
    
    for(j in seq_len(nrow(ht_layout_index))) {
        pushViewport(viewport(layout.pos.row = ht_layout_index[j, 1], layout.pos.col = ht_layout_index[j, 2]))
        ht_graphic_fun_list[[j]](object)
        upViewport()
    }

    # add decorations (axis, and lines)

    upViewport()

    for(i in centralized_heatmap_index) {
    	
    	ht = object@ht_list[[i]]
    	heatmap_name = ht@name
    	upstream_index = attr(ht@matrix, "upstream_index")
		downstream_index = attr(ht@matrix, "downstream_index")
		body_index = attr(ht@matrix, "body_index")
		n1 = length(upstream_index)
		n2 = length(body_index)
		n3 = length(downstream_index)
		n = n1 + n2 + n3

		if(ht@heatmap_param$pos_line) {
	    	for(j in seq_along(ht@row_order_list)) {
	    		seekViewport(paste0(heatmap_name, "_heatmap_body_", j))
	    		grid.lines(rep(n1/n, 2), c(0, 1), gp = ht@heatmap_param$pos_line_gp)
	    		if(n2) grid.lines(rep((n1+n2)/n, 2), c(0, 1), gp = ht@heatmap_param$pos_line_gp)
	    	}
	    }
    }
})

