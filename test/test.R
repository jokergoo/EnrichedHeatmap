
library(ComplexHeatmap)
library(GenomicFeatures)
library(GenomicRanges)

tss = promoters(genes, upstream = 0, downstream = 1)

x1 = normalizeToMatrix(H3K4me3, tss, value_column = "coverage", extend = 5000, mean_mode = "w0", w = 100, smooth = TRUE)
x2 = normalizeToMatrix(meth, tss, value_column = "meth", extend = 5000, w = 100, empty_value = 0.5, smooth = TRUE)

CentralizedHeatmap(x1, name = "H3K4me3", column_title = "histone modification") + 
CentralizedHeatmap(x2, name = "methylation", column_title = "methylation")

CentralizedHeatmap(x1, name = "H3K4me3", column_title = "histone modification", 
	split = sample(letters[1:2], length(tss), replace = TRUE))

x = normalizeToMatrix(meth, cgi, value_column = "meth", extend = 5000, w = 100, empty_value = 0.5)

CentralizedHeatmap(x2, name = "methylation", column_title = "methylation", axis_name_rot = -90)

