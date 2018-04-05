\name{draw-EnrichedHeatmap-method}
\alias{draw,EnrichedHeatmap-method}
\title{
Draw a single heatmap
}
\description{
Draw a single heatmap
}
\usage{
\S4method{draw}{EnrichedHeatmap}(object, internal = FALSE, ...)
}
\arguments{

  \item{object}{an \code{\link{EnrichedHeatmap-class}} object.}
  \item{internal}{only used internally.}
  \item{...}{pass to \code{\link[ComplexHeatmap]{draw}},HeatmapList-method.}

}
\details{
The function creates an \code{\link{EnrichedHeatmapList-class}} object which only contains a single heatmap
and call \code{\link{draw,EnrichedHeatmapList-method}} to make the final heatmap.
}
\value{
An \code{\link{EnrichedHeatmapList-class}} object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# see documentation of EnrichedHeatmap
NULL
}
