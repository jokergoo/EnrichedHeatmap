\name{EnrichedHeatmap}
\alias{EnrichedHeatmap}
\title{
Constructor method for EnrichedHeatmap class

}
\description{
Constructor method for EnrichedHeatmap class

}
\usage{
EnrichedHeatmap(mat, score_fun = enriched_score, pos_line = TRUE,
    pos_line_gp = gpar(lty = 2), axis_name = NULL, axis_name_rot = 0,
    axis_name_gp = gpar(fontsize = 10), border = TRUE, cluster_rows = FALSE, ...)}
\arguments{

  \item{mat}{a matrix which is returned by \code{\link{normalizeToMatrix}}}
  \item{score_fun}{score function which calcualte enriched scores for rows in \code{mat}}
  \item{pos_line}{whether draw vertical lines which represent the position of \code{target}}
  \item{pos_line_gp}{graphical parameters for lines}
  \item{axis_name}{names for axis}
  \item{axis_name_rot}{rotation for axis names}
  \item{axis_name_gp}{graphical parameters for axis names}
  \item{border}{whether show border of the heatmap}
  \item{cluster_rows}{clustering on rows are turned off by default}
  \item{...}{pass to \code{\link[ComplexHeatmap]{Heatmap}}}
}
\details{
\code{EnrichedHeatmap-class} is inherited from \code{Heatmap-class}. Following parameters are 
set with pre-defined values:

\describe{
  \item{row_order}{the rows are sorted by the enriched score which is calcualted by \code{score_fun} argument.The sorting is applied decreasingly.}
  \item{cluster_columns}{enforced to be \code{FALSE}}
  \item{show_row_names}{enforced to be \code{FALSE}}
  \item{show_column_names}{enforced to be \code{FALSE}}
  \item{bottom_annotation}{enforced to be \code{NULL} }
  \item{column_title_side}{enforced to be \code{top}}
}

With above pre-defined values, no graphics will be drawn below the heatmap, then the space
below the heatmap is used to add a new graph which contains the axis. A (or two) line which corresponds to 
the position of \code{target} will be added to the heatmap body as well.

Same as the \code{\link[ComplexHeatmap]{Heatmap-class}}, users can make controls on the heatmap such as
apply clustering on rows (this will ignore the \code{row_order}), split rows by data frame or k-means clustering.

}
\value{
A \code{\link{EnrichedHeatmap-class}} object which is inherited from \code{\link[ComplexHeatmap]{Heatmap-class}}

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
