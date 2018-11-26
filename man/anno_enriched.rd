\name{anno_enriched}
\alias{anno_enriched}
\title{
Annotation Function to Show the Enrichment
}
\description{
Annotation Function to Show the Enrichment
}
\usage{
anno_enriched(gp = gpar(col = "red"), pos_line = NULL, pos_line_gp = NULL,
    ylim = NULL, value = c("mean", "sum", "abs_mean", "abs_sum"),
    yaxis = TRUE, axis = yaxis, axis_param = list(side = "right"),
    show_error = FALSE, height = unit(2, "cm"), ...)
}
\arguments{

  \item{gp}{Graphic parameters. There are two non-standard parameters: \code{neg_col} and \code{pos_col}.  If these two parameters are defined, the positive signals and negatie signals are visualized separatedly. The graphic parameters can be set as vectors when the heatmap or heatmap list is split into several row clusters.}
  \item{pos_line}{Whether draw vertical lines which represent positions of \code{target}?}
  \item{pos_line_gp}{Graphic parameters for the position lines.}
  \item{ylim}{Ranges on y-axis. By default it is inferred from the data.}
  \item{value}{The method to summarize signals from columns of the normalized matrix.}
  \item{yaxis}{Deprecated, use \code{axis} instead.}
  \item{axis}{Whether show axis?}
  \item{axis_param}{parameters for controlling axis. See \code{\link[ComplexHeatmap]{default_axis_param}} for all possible settings and default parameters.}
  \item{show_error}{Whether show error regions which are one standard error to the mean value? Color of error area is same as the corresponding lines with 75 percent transparency.}
  \item{height}{Height of the annotation.}
  \item{...}{Other arguments.}

}
\details{
This annotation functions shows mean values (or depends on the method set in \code{value} argument) of columns in the normalized matrix
which summarises the enrichment of the signals to the targets.

If rows are splitted, the enriched lines are calculated for each row cluster and there will also be multiple lines in this annotation viewport.

It should only be placed as column annotation of the enriched heatmap.
}
\value{
A column annotation function which should be set to \code{top_annotation} argument in \code{\link{EnrichedHeatmap}}.
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
