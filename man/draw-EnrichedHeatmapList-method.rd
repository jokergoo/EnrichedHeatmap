\name{draw-EnrichedHeatmapList-method}
\alias{draw,EnrichedHeatmapList-method}
\title{
Draw a list of heatmaps

}
\description{
Draw a list of heatmaps

}
\usage{
\S4method{draw}{EnrichedHeatmapList}(object, padding = unit(c(2, 2, 2, 2), "mm"), ..., newpage= TRUE)}
\arguments{

  \item{object}{a \code{\link{EnrichedHeatmapList-class}} object}
  \item{padding}{padding of the plot. Elements correspond to bottom, left, top, right paddings.}
  \item{...}{pass to \code{\link[ComplexHeatmap]{make_layout,HeatmapList-method}}}
  \item{newpage}{whether to create a new page}
}
\details{
It calls \code{\link[ComplexHeatmap]{draw,HeatmapList-method}} to make the plot but with some adjustment.

}
\value{
No value is returned.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
\examples{
# see documentation of EnrichedHeatmap
NULL

}
