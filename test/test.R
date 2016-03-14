
library(ComplexHeatmap)
library(GenomicFeatures)
library(GenomicRanges)
library(locfit)
library(GetoptLong)

load(paste0(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap")))

tss = promoters(genes, upstream = 0, downstream = 1)

x1 = normalizeToMatrix(H3K4me3, tss, value_column = "coverage", extend = 5000, mean_mode = "w0", w = 100)
x2 = normalizeToMatrix(meth, tss, value_column = "meth", mean_mode = "absolute", extend = 5000, w = 100, smooth = TRUE)

EnrichedHeatmap(x1, name = "H3K4me3", column_title = "histone modification") + 
EnrichedHeatmap(x2, name = "methylation", column_title = "methylation")

EnrichedHeatmap(x1, name = "H3K4me3", column_title = "histone modification", 
	split = sample(letters[1:2], length(tss), replace = TRUE))

x = normalizeToMatrix(meth, cgi, value_column = "meth", mean_mode = "absolute", extend = 5000, w = 50, smooth = TRUE)
EnrichedHeatmap(x, name = "methylation")

x = normalizeToMatrix(meth, cgi, value_column = "meth", mean_mode = "absolute", extend = c(0, 5000), smooth = TRUE)
EnrichedHeatmap(x, name = "methylation")


x = normalizeToMatrix(meth, cgi, value_column = "meth", mean_mode = "absolute", extend = c(5000, 0), smooth = TRUE)
EnrichedHeatmap(x, name = "methylation")

x = normalizeToMatrix(meth, cgi, value_column = "meth", mean_mode = "absolute", extend = 0, smooth = TRUE)
EnrichedHeatmap(x, name = "methylation")

x = normalizeToMatrix(meth[1:100], cgi, value_column = "meth", mean_mode = "absolute", extend = 0, smooth = TRUE)


kmeans(mat1)$clusters
EnrichedHeatmap(mat1, col = c("white", "red"), name = "H3K4me3", 
	split = kmeans(mat1, centers =2)$cluster)
