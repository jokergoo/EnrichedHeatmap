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

  \item{object}{a \code{\link{EnrichedHeatmap-class}} object.}
  \item{internal}{only used inside the calling of \code{\link{draw,EnrichedHeatmapList-method}}.}
  \item{...}{pass to \code{\link[ComplexHeatmap]{draw,HeatmapList-method}}.}

}
\details{
The function creates a \code{\link{EnrichedHeatmapList-class}} object which only contains a single heatmap and call \code{\link{draw,EnrichedHeatmapList-method}} to make the final heatmap.  


}
\value{
This function returns no value.  


}
\author{
Zuguang Gu <z.gu@dkfz.de>  


}
\section{Example}{
# see documentation of \code{\link{EnrichedHeatmap}} NULL 


}
