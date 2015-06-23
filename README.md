# Make Enriched Heatmaps

Enriched heatmap is a special type of heatmap which visualizes the enrichment of genomic signals on specific target regions. It is broadly used to visualize e.g. how histone marks are enriched to specific sites.

There are several tools that can make such heatmap (e.g. [ngs.plot](https://github.com/shenlab-sinai/ngsplot) or [deepTools](https://github.com/fidelram/deepTools)). Here we implement Enriched heatmap by [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap) package. Since this type of heatmap is just a normal heatmap but with some fixed settings, with the functionality of ComplexHeatmap, it would be much easier to customize the heatmap as well as concatenating a list of heatmaps to show correspondance between differnet data sources.

Examples:

```{r}
EnrichedHeatmap(mat1, col = c("white", "red"), name = "H3K4me3", km = 3, width = 1,
    top_annotation = HeatmapAnnotation(lines = anno_enriched()), 
    top_annotation_height = unit(2, "cm"), row_title_rot = 0,
    column_title = "H3K4me3") + 
EnrichedHeatmap(mat2, name = "methylation", width = 1,
    column_title = "Methylation") +
Heatmap(log2(rpkm+1), col = c("white", "orange"), name = "log2(rpkm+1)", 
    show_row_names = FALSE, width = unit(5, "mm"))
```

![EnrichedHeatmap](https://cloud.githubusercontent.com/assets/449218/8305941/ee147c18-19b3-11e5-8732-863d17e53c60.png)
