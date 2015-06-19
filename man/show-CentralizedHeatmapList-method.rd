\name{show-CentralizedHeatmapList-method}
\alias{show,CentralizedHeatmapList-method}
\title{
Draw a list of heatmaps with default parameters

}
\description{
Draw a list of heatmaps with default parameters

}
\usage{
\S4method{show}{CentralizedHeatmapList}(object)}
\arguments{

  \item{object}{a \code{\link{CentralizedHeatmapList-class}} object.}
}
\details{
Actually it calls \code{\link{draw,CentralizedHeatmapList-method}}, but only with default parameters. If users want to customize the heatmap,
they can pass parameters directly to \code{\link{draw,CentralizedHeatmapList-method}}.

}
\value{
This function returns no value.

}
