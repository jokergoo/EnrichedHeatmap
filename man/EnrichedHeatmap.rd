\name{EnrichedHeatmap}
\alias{EnrichedHeatmap}
\title{
Constructor Method for the Enriched Heatmap
}
\description{
Constructor Method for the Enriched Heatmap
}
\usage{
EnrichedHeatmap(mat,
    col,
    top_annotation = HeatmapAnnotation(enriched = anno_enriched()),
    row_order = order(enriched_score(mat), decreasing = TRUE),
    pos_line = TRUE,
    pos_line_gp = gpar(lty = 2),
    axis_name = NULL,
    axis_name_rot = 0,
    axis_name_gp = gpar(fontsize = 10),
    border = TRUE,
    cluster_rows = FALSE,
    row_dend_reorder = -enriched_score(mat),
    show_row_dend = FALSE,
    show_row_names = FALSE,
    heatmap_legend_param = list(),
    ...)
}
\arguments{

  \item{mat}{A matrix which is returned by \code{\link{normalizeToMatrix}}.}
  \item{col}{Color settings. If the signals are categorical, color should be a vector with category levels as names.}
  \item{top_annotation}{A special annotation which is always put on top of the enriched heatmap and is constructed by \code{\link{anno_enriched}}.}
  \item{row_order}{Row order. Default rows are ordered by enriched scores calculated from \code{\link{enriched_score}}.}
  \item{pos_line}{Whether draw vertical lines which represent the positions of \code{target}?}
  \item{pos_line_gp}{Graphic parameters for the position lines.}
  \item{axis_name}{Names for axis which is below the heatmap. If the targets are single points, \code{axis_name} is a vector of length three which corresponds to upstream, target itself and downstream. If the targets are regions with width larger than 1, \code{axis_name} should be a vector of length four which  corresponds to upstream, start of targets, end of targets and downstream.}
  \item{axis_name_rot}{Rotation for axis names.}
  \item{axis_name_gp}{Graphic parameters for axis names.}
  \item{border}{Whether show the border of the heatmap?}
  \item{cluster_rows}{Clustering on rows are turned off by default.}
  \item{show_row_dend}{Whether show dendrograms on rows if hierarchical clustering is applied on rows?}
  \item{row_dend_reorder}{Weight for reordering the row dendrogram. It is reordered by enriched scores by default.}
  \item{show_row_names}{Whether show row names?}
  \item{heatmap_legend_param}{A list of settings for heatmap legends. \code{at} and \code{labels} can not be set here.}
  \item{...}{Other arguments passed to \code{\link[ComplexHeatmap]{Heatmap}}.}

}
\details{
The enriched heatmap is essentially a normal heatmap but with several special settings. Following parameters are 
set with pre-defined values:

\describe{
  \item{\code{cluster_columns}}{enforced to be \code{FALSE}}
  \item{\code{show_column_names}}{enforced to be \code{FALSE}}
  \item{\code{bottom_annotation}}{enforced to be \code{NULL} }
}

\code{\link{EnrichedHeatmap}} calls \code{\link[ComplexHeatmap]{Heatmap}}, thus, most of the 
arguments in \code{\link[ComplexHeatmap]{Heatmap}} are usable in \code{\link{EnrichedHeatmap}} such as
to apply clustering on rows, or to split rows by a data frame or k-means clustering. Users can also 
add more than one heatmaps by \code{+} operator. Enriched heatmaps and normal heatmaps can be
concatenated mixed.

For detailed demonstration, please go to the vignette.
}
\value{
A \code{\link{Heatmap-class}} object.
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
