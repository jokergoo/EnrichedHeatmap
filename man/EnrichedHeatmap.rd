\name{EnrichedHeatmap}
\alias{EnrichedHeatmap}
\title{
Constructor method for EnrichedHeatmap class
}
\description{
Constructor method for EnrichedHeatmap class
}
\usage{
EnrichedHeatmap(mat, top_annotation = HeatmapAnnotation(enriched = anno_enriched()),
    top_annotation_height = unit(2, "cm"),
    row_order = order(enriched_score(mat), decreasing = TRUE), pos_line = TRUE,
    pos_line_gp = gpar(lty = 2), axis_name = NULL, axis_name_rot = 0,
    axis_name_gp = gpar(fontsize = 10), border = TRUE, cluster_rows = FALSE,
    row_dend_reorder = -enriched_score(mat),
    show_row_dend = FALSE, show_row_names = FALSE, ...)
}
\arguments{

  \item{mat}{a matrix which is returned by \code{\link{normalizeToMatrix}}}
  \item{top_annotation}{a specific annotation which is always put on top of the enriched heatmap and is constructed by \code{\link{anno_enriched}}}
  \item{top_annotation_height}{the height of the top annotation}
  \item{row_order}{row order. Default rows are ordered by enriched scores calculated from \code{\link{enriched_score}}}
  \item{pos_line}{whether draw vertical lines which represent the positions of \code{target}}
  \item{pos_line_gp}{graphic parameters for the position lines}
  \item{axis_name}{names for axis which is below the heatmap. If the targets are single points, \code{axis_name} is a vector of length three which corresponds to upstream, target itself and downstream. If the targets are regions with width larger than 1, \code{axis_name} should be a vector of length four which  corresponds to upstream, start of targets, end of targets and downstream.}
  \item{axis_name_rot}{rotation for axis names}
  \item{axis_name_gp}{graphic parameters for axis names}
  \item{border}{whether show border of the heatmap}
  \item{cluster_rows}{clustering on rows are turned off by default}
  \item{show_row_dend}{whether show dendrograms on rows if apply hierarchical clustering on rows}
  \item{row_dend_reorder}{weight for reordering the row dendrogram. It is reordered by enriched scores by default.}
  \item{show_row_names}{whether show row names}
  \item{...}{pass to \code{\link[ComplexHeatmap]{Heatmap}}}

}
\details{
\code{\link{EnrichedHeatmap-class}} is inherited from \code{\link[ComplexHeatmap]{Heatmap-class}}. Following parameters are 
set with pre-defined values:

\describe{
  \item{\code{cluster_columns}}{enforced to be \code{FALSE}}
  \item{\code{show_column_names}}{enforced to be \code{FALSE}}
  \item{\code{bottom_annotation}}{enforced to be \code{NULL} }
  \item{\code{column_title_side}}{enforced to be \code{top}}
}

A \code{\link{EnrichedHeatmap-class}} object is also a \code{\link[ComplexHeatmap]{Heatmap-class}} object, thus, most of the 
arguments in \code{\link[ComplexHeatmap]{Heatmap}} are usable in \code{\link{EnrichedHeatmap}} such as
to apply clustering on rows, or to split rows by data frame or k-means clustering. Users can also 
add more than one heatmaps by \code{+} operator. For a detailed demonstration, please go to the vignette.
}
\value{
An \code{\link{EnrichedHeatmap-class}} object which is inherited from \code{\link[ComplexHeatmap]{Heatmap-class}}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
load(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap"))
mat3 = normalizeToMatrix(meth, cgi, value_column = "meth", mean_mode = "absolute",
    extend = 5000, w = 50, smooth = TRUE)
EnrichedHeatmap(mat3, name = "methylation", column_title = "methylation near CGI")
EnrichedHeatmap(mat3, name = "meth1") + EnrichedHeatmap(mat3, name = "meth2")
# for more examples, please go to the vignette
}
