\name{draw-CentralizedHeatmap-method}
\alias{draw,CentralizedHeatmap-method}
\title{
Draw a single heatmap

}
\description{
Draw a single heatmap

}
\usage{
\S4method{draw}{CentralizedHeatmap}(object, internal = FALSE, ...)}
\arguments{

  \item{object}{a \code{\link{CentralizedHeatmap-class}} object.}
  \item{internal}{only used inside the calling of \code{\link{draw,CentralizedHeatmapList-method}}. Only heatmap without legends will be drawn.}
  \item{...}{pass to \code{\link[ComplexHeatmap]{draw,HeatmapList-method}}.}
}
\details{
The function creates a \code{\link{CentralizedHeatmapList-class}} object which only contains a single heatmap
and call \code{\link{draw,CentralizedHeatmapList-method}} to make the final heatmap.

}
\value{
This function returns no value.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
