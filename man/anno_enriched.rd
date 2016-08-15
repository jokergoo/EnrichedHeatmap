\name{anno_enriched}
\alias{anno_enriched}
\title{
Annotation function to show the enrichment
}
\description{
Annotation function to show the enrichment
}
\usage{
anno_enriched(gp = gpar(col = "red"), pos_line = TRUE, pos_line_gp = gpar(lty = 2),
    yaxis = TRUE, ylim = NULL, value = c("mean", "sum", "abs_mean", "abs_sum"), yaxis_side = "right",
    yaxis_gp = gpar(fontsize = 8), show_error = FALSE)
}
\arguments{

  \item{gp}{graphic parameters}
  \item{pos_line}{whether draw vertical lines which represent the position of \code{target}}
  \item{pos_line_gp}{graphic parameters}
  \item{yaxis}{whether show yaxis}
  \item{ylim}{ranges on y-axis}
  \item{value}{what type of value corresponds to the y-axis}
  \item{yaxis_side}{side of y-axis}
  \item{yaxis_gp}{graphic parameters for yaxis}
  \item{show_error}{whether show error regions which are +-1 se to the mean value. Color of error area is same as the corresponding lines with 75 percent transparency.}

}
\details{
This annotation functions shows mean values of columns in the normalized matrix
which represents the enrichment of the signals to the targets.

If rows are splitted, there will also be multiple lines in this annotation.

It should only be placed as column annotation of the Enriched Heatmap.
}
\value{
A column annotation function which can be set to \code{top_annotation} argument in \code{\link{EnrichedHeatmap}}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
load(paste0(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap")))
tss = promoters(genes, upstream = 0, downstream = 1)
mat1 = normalizeToMatrix(H3K4me3, tss, value_column = "coverage", 
    extend = 5000, mean_mode = "w0", w = 50, trim = c(0, 0.01))
EnrichedHeatmap(mat1, col = c("white", "red"), name = "H3K4me3",
    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
    top_annotation_height = unit(2, "cm"),
    km = 3, row_title_rot = 0)
}
