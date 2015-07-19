\name{show-EnrichedHeatmapList-method}
\alias{show,EnrichedHeatmapList-method}
\title{
Draw a list of heatmaps with default parameters

}
\description{
Draw a list of heatmaps with default parameters

}
\usage{
\S4method{show}{EnrichedHeatmapList}(object)}
\arguments{

  \item{object}{an \code{\link{EnrichedHeatmapList-class}} object.}
}
\details{
Actually it calls \code{\link{draw,EnrichedHeatmapList-method}}, but only with default parameters. If users want to customize the heatmap,
they can pass parameters directly to \code{\link{draw,EnrichedHeatmapList-method}}.

}
\value{
This function returns no value.

}
\examples{
# see documentation of EnrichedHeatmap
NULL

}
