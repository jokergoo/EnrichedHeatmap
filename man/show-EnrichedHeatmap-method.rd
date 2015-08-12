\name{show-EnrichedHeatmap-method}
\alias{show,EnrichedHeatmap-method}
\title{
Draw the single heatmap with default parameters
}
\description{
Draw the single heatmap with default parameters
}
\usage{
\S4method{show}{EnrichedHeatmap}(object)
}
\arguments{

  \item{object}{an \code{\link{EnrichedHeatmap-class}} object.}

}
\details{
Actually it calls \code{\link{draw,EnrichedHeatmap-method}}, but only with default parameters. If users want to customize the heatmap,
they can pass parameters directly to \code{\link{draw,EnrichedHeatmap-method}}.
}
\value{
This function returns no value.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# see documentation of EnrichedHeatmap
NULL
}
