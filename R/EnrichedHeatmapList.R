
# == title
# Extarct Enrichment Annotation Graphics as a Separated Plot
#
# == param
# -ht_list The heatmap list returned by `draw,HeatmapList-method`.
# -which The index of enriched heatmap in the heatmap list. The value can be an integer index or a character index (the name of the heatmap).
# -newpage Whether call `grid::grid.newpage` to create a new page?
# -padding Padding of the plot.
#
# == details
# The extracted plot is exactly the same as that on the enriched heatmap. 
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
extract_anno_enriched = function(ht_list, which = NULL, newpage = TRUE, padding = NULL) {

    if(!any(sapply(ht_list@ht_list, is_enriched_heatmap))) {
        stop_wrap("`ht_list` should contain at least one enriched heatmap.")
    }

    if(!ht_list@layout$initialized) {
        stop_wrap("`ht_list` should be returned by `draw()` function.")
    }

    if(newpage) grid.newpage()
    if(is.null(which)) which = which(sapply(ht_list@ht_list, is_enriched_heatmap))[1]
    object = ht_list@ht_list[[which]]
    if(!is_enriched_heatmap(object)) {
        stop_wrap(paste0("heamtap ", which, " is not an enriched heatmap."))
    }

    column_title = object@column_title
    if(length(column_title) == 0) column_title = object@name
    title_height = 2*grobHeight(textGrob(column_title))
    axis_height = object@heatmap_param$axis_height
    
    ha = ht_list@ht_list[[which]]@top_annotation
    anno = ha@anno_list[[1]]@fun
    anno@height = unit(1, "npc")
    left_ext = anno@extended[2]
    right_ext = anno@extended[4]

    if(is.null(padding)) padding = unit.c(left_ext, right_ext) + unit(2, "mm")
    
    # viewprot for title
    pushViewport(viewport(y = 1, x = padding[1], 
        height = title_height, 
        width = unit(1, "npc") - padding[1] - padding[2], just = c("left", "top")))
    grid.text(column_title)
    upViewport()

    # viewport the enriched lines
    pushViewport(viewport(y = axis_height, x = padding[1], 
        height = unit(1, "npc") - axis_height - title_height, 
        width = unit(1, "npc") - padding[1] - padding[2], 
        just = c("left", "bottom")))
    x = calc_minor_ticks(object@matrix)
    if(length(x)) {
        grid.segments(x, 0, x, 1, gp = gpar(col = "#CCCCCC", lty = 2))
    }
    
    f1 = function() draw(anno, seq_len(ncol(object@matrix)))
    f2 = function() f1()
    f3 = function() f2()
    f4 = function() f3()
    f5 = function() f4()
    f6 = function() f5()
    f7 = function() f6()
    f7()

    object@heatmap_param$axis_fun()
    upViewport()
}


