\name{show-CentralizedHeatmap-method}
\alias{show,CentralizedHeatmap-method}
\title{
Draw the single heatmap with default parameters

}
\description{
Draw the single heatmap with default parameters

}
\usage{
\S4method{show}{CentralizedHeatmap}(object)}
\arguments{

  \item{object}{a \code{\link{CentralizedHeatmap-class}} object.}
}
\details{
Actually it calls \code{\link{draw,CentralizedHeatmap-method}}, but only with default parameters. If users want to customize the heatmap,
they can pass parameters directly to \code{\link{draw,CentralizedHeatmap-method}}.

}
\value{
This function returns no value.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
