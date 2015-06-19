\name{CentralizedHeatmap}
\alias{CentralizedHeatmap}
\title{
Constructor method for CentralizedHeatmap class

}
\description{
Constructor method for CentralizedHeatmap class

}
\usage{
CentralizedHeatmap(mat, score_fun = score_fun, pos_line = TRUE, pos_line_gp = gpar(lty = 2),
    axis_name = NULL, axis_name_rot = 0, axis_name_gp = gpar(), ...)}
\arguments{

  \item{mat}{a matrix which is returned by \code{\link{normalizeToMatrix}}}
  \item{score_fun}{score function which calcualte scores by rows in \code{mat}}
  \item{pos_line}{whether draw vertical lines which represent the position of center}
  \item{pos_line_gp}{graphical parameters for lines}
  \item{axis_name}{names for axis}
  \item{axis_name_rot}{rotation for axis names}
  \item{axis_name_gp}{graphical parameters for axis names}
  \item{...}{pass to \code{\link[ComplexHeatmap]{Heatmap}}}
}
\details{
xxx

}
\value{
a \code{\link{CentralizedHeatmap-class}} object which is inherited from \code{\link[ComplexHeatmap]{HeatmapList-class}}

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
