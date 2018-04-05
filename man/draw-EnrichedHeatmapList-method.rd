\name{draw-EnrichedHeatmapList-method}
\alias{draw,EnrichedHeatmapList-method}
\title{
Draw a list of heatmaps
}
\description{
Draw a list of heatmaps
}
\usage{
\S4method{draw}{EnrichedHeatmapList}(object, padding = unit(c(2, 2, 2, 2), "mm"),
    newpage= TRUE, ...)
}
\arguments{

  \item{object}{an \code{\link{EnrichedHeatmapList-class}} object}
  \item{padding}{padding of the plot. The four values correspond to bottom, left, top, right paddings.}
  \item{newpage}{whether to create a new page}
  \item{...}{pass to \code{\link[ComplexHeatmap]{make_layout}}, HeatmapList-method or \code{\link[ComplexHeatmap]{draw}}, HeatmapList-method}

}
\details{
It calls \code{\link[ComplexHeatmap]{draw}}, HeatmapList-method to make the plot but with some adjustment
specificly for enriched heatmaps.
}
\value{
An \code{\link{EnrichedHeatmapList}} object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# see documentation of EnrichedHeatmap
NULL
}
