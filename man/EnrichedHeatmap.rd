\name{EnrichedHeatmap}
\alias{EnrichedHeatmap}
\title{
Constructor method for EnrichedHeatmap class
}
\description{
Constructor method for EnrichedHeatmap class
}
\usage{
<<<<<<< HEAD
EnrichedHeatmap(mat, top_annotation = HeatmapAnnotation(enriched = anno_enriched()),
    top_annotation_height = unit(2, "cm"),
    score_fun = enriched_score, row_order = NULL, pos_line = TRUE,
    pos_line_gp = gpar(lty = 2), axis_name = NULL, axis_name_rot = 0,
    axis_name_gp = gpar(fontsize = 10), border = TRUE, cluster_rows = FALSE,
    show_row_dend = FALSE, show_row_names = FALSE, ...)
=======
EnrichedHeatmap(mat, score_fun = enriched_score, row_order = NULL, pos_line = TRUE,
    pos_line_gp = gpar(lty = 2), axis_name = NULL, axis_name_rot = NULL,
    axis_name_gp = gpar(fontsize = 10), border = TRUE, cluster_rows = FALSE,
    show_row_dend = FALSE, ...)
>>>>>>> bioc/master
}
\arguments{

  \item{mat}{a matrix which is returned by \code{\link{normalizeToMatrix}}}
<<<<<<< HEAD
  \item{top_annotation}{a specific annotation which is always put on top of the enriched heatmap and is constructed by \code{\link{anno_enriched}}}
  \item{top_annotation_height}{the height of the top annotation}
  \item{score_fun}{a function which calculates enriched scores for rows in \code{mat}. This function can be self-defined, refer to \code{\link{enriched_score}} to find out how to design it. Note if row clustering is turned on, this argument is ignored.}
  \item{row_order}{row order. If it is specified, \code{score_fun} is ignored.}
  \item{pos_line}{whether draw vertical lines which represent the positions of \code{target}}
  \item{pos_line_gp}{graphic parameters for the position lines}
=======
  \item{score_fun}{a function which calculates enriched scores for rows in \code{mat}. This function can be self-defined, take a look at \code{\link{enriched_score}} to find out how to design it. Note if row clustering is turned on, this argument is ignored.}
  \item{row_order}{row order. If it is specified, \code{score_fun} is ignored.}
  \item{pos_line}{whether draw vertical lines which represent the position of \code{target}}
  \item{pos_line_gp}{graphic parameters for lines}
>>>>>>> bioc/master
  \item{axis_name}{names for axis which is below the heatmap. If the targets are single points, \code{axis_name} is a vector of length three which corresponds to upstream, target itself and downstream. If the targets are regions with width larger than 1, \code{axis_name} should be a vector of length four which  corresponds to upstream, start of targets, end of targets and downstream.}
  \item{axis_name_rot}{rotation for axis names}
  \item{axis_name_gp}{graphic parameters for axis names}
  \item{border}{whether show border of the heatmap}
  \item{cluster_rows}{clustering on rows are turned off by default}
<<<<<<< HEAD
  \item{show_row_dend}{whether show dendrograms on rows if apply hierarchical clustering on rows}
  \item{show_row_names}{whether show row names}
=======
  \item{show_row_dend}{whether show dendrograms on rows}
>>>>>>> bioc/master
  \item{...}{pass to \code{\link[ComplexHeatmap]{Heatmap}}}

}
\details{
\code{\link{EnrichedHeatmap-class}} is inherited from \code{\link[ComplexHeatmap]{Heatmap-class}}. Following parameters are 
set with pre-defined values:

\describe{
<<<<<<< HEAD
  \item{\code{row_order}}{the rows are sorted by the enriched score which is calculated by \code{score_fun}. The sorting is applied decreasingly.}
  \item{\code{cluster_columns}}{enforced to be \code{FALSE}}
=======
  \item{\code{row_order}}{the rows are sorted by the enriched score which is calcualted by \code{score_fun}. The sorting is applied decreasingly.}
  \item{\code{cluster_columns}}{enforced to be \code{FALSE}}
  \item{\code{show_row_names}}{enforced to be \code{FALSE}}
>>>>>>> bioc/master
  \item{\code{show_column_names}}{enforced to be \code{FALSE}}
  \item{\code{bottom_annotation}}{enforced to be \code{NULL} }
  \item{\code{column_title_side}}{enforced to be \code{top}}
}

<<<<<<< HEAD
A \code{\link{EnrichedHeatmap-class}} object is also a \code{\link[ComplexHeatmap]{Heatmap-class}} object, thus, most of the 
arguments in \code{\link[ComplexHeatmap]{Heatmap}} are usable in \code{\link{EnrichedHeatmap}} such as
to apply clustering on rows, or to split rows by data frame or k-means clustering. Users can also 
add more than one heatmaps by \code{+} operator. For a detailed demonstration, please go to the vignette.
=======
With above pre-defined values, no graphics will be drawn below the heatmap, then the space
below the heatmap can be used to add a new graph which contains the axis. A (or two) line which corresponds to 
the position of \code{target} will be added to the heatmap body as well.

Same as the \code{\link[ComplexHeatmap]{Heatmap-class}}, users can make more controls on the heatmap such as
apply clustering on rows, or split rows by data frame or k-means clustering. Users can also 
add more than one heatmaps by \code{+} operator.

For a detailed demonstration, please go to the vignette.
>>>>>>> bioc/master
}
\value{
An \code{\link{EnrichedHeatmap-class}} object which is inherited from \code{\link[ComplexHeatmap]{Heatmap-class}}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
<<<<<<< HEAD
load(system.file("extdata", "chr21_test_data.RData", package = "EnrichedHeatmap"))
mat3 = normalizeToMatrix(meth, cgi, value_column = "meth", mean_mode = "absolute",
    extend = 5000, w = 50, background = 0.5)
=======
load(paste0(system.file("extdata", "chr21_test_data.RData", 
    package = "EnrichedHeatmap")))
mat3 = normalizeToMatrix(meth, cgi, value_column = "meth", mean_mode = "absolute",
    extend = 5000, w = 50, empty_value = 0.5)
>>>>>>> bioc/master
EnrichedHeatmap(mat3, name = "methylation", column_title = "methylation near CGI")
EnrichedHeatmap(mat3, name = "meth1") + EnrichedHeatmap(mat3, name = "meth2")
# for more examples, please go to the vignette
}
