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
    yaxis_facing = ifelse(yaxis_side == "right", "right", "left"),
    yaxis_gp = gpar(fontsize = 8), show_error = FALSE)
}
\arguments{

  \item{gp}{graphic parameters. There are two non-standard parameters: \code{neg_col} and \code{pos_col}.  If these two parameters are defined, the positive signals and negatie signals are visualized separatedly. The graphic parameters can be set as vectors when the heatmap or heatmap list is split into several row clusters.}
  \item{pos_line}{whether to draw vertical lines which represent positions of \code{target}}
  \item{pos_line_gp}{graphic parameters for the position lines}
  \item{yaxis}{whether show yaxis}
  \item{ylim}{ranges on y-axis, by default it is inferred from the data}
  \item{value}{the method to summarize signals from columns of the noramlized matrix}
  \item{yaxis_side}{side of y-axis}
  \item{yaxis_facing}{facing of the axis ticks and labels. It can be set to avoid overlapping text when multiple heatmaps are plotted together}
  \item{yaxis_gp}{graphic parameters for y-axis}
  \item{show_error}{whether show error regions which are one standard error to the mean value. Color of error area is same as the corresponding lines with 75 percent transparency.}

}
\details{
This annotation functions shows mean values (or depends on the method set in \code{value} argument) of columns in the normalized matrix
which represents the enrichment of the signals to the targets.

If rows are splitted, the enriched lines are calculated for each row cluster and there will also be multiple lines in this annotation viewport.

It should only be placed as column annotation of the enriched heatmap.
}
\value{
A column annotation function which can be set to \code{top_annotation} argument in \code{\link{EnrichedHeatmap}}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
load(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap"))
tss = promoters(genes, upstream = 0, downstream = 1)
mat1 = normalizeToMatrix(H3K4me3, tss, value_column = "coverage", 
    extend = 5000, mean_mode = "w0", w = 50, keep = c(0, 0.99))
EnrichedHeatmap(mat1, col = c("white", "red"), name = "H3K4me3",
    top_annotation = HeatmapAnnotation(lines = anno_enriched(gp = gpar(col = 2:4))), 
    km = 3, row_title_rot = 0)
}
