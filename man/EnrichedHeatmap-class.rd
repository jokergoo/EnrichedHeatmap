\name{EnrichedHeatmap-class}
\docType{class}
\alias{EnrichedHeatmap-class}
\title{
Class for a single heatmap
}
\description{
Class for a single heatmap
}
\details{
The structure of \code{\link{EnrichedHeatmap-class}} is the same as
\code{\link[ComplexHeatmap]{HeatmapList-class}} and the class is inherited from \code{\link[ComplexHeatmap]{Heatmap-class}}.

The \code{\link{EnrichedHeatmap-class}} pre-defines some parameters for \code{\link[ComplexHeatmap]{Heatmap-class}} such as
the order of rows and supressing column clustering. Also there are several
new parameters that are attached in the object.
}
\section{Methods}{
The \code{\link{EnrichedHeatmap-class}} provides following methods:

\itemize{
  \item \code{\link{EnrichedHeatmap}}: constructor method.
  \item \code{\link{draw,EnrichedHeatmap-method}}: draw a single heatmap.
}}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
