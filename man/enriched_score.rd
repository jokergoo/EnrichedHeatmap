\name{enriched_score}
\alias{enriched_score}
\title{
Enriched scores
}
\description{
Enriched scores
}
\usage{
enriched_score(x1, x2, x3)
}
\arguments{

  \item{x1}{a vector corresponding to values in upstream windows}
  \item{x2}{a vector corresponding to values in target windows}
  \item{x3}{a vector corresponding to values in downstream windows}

}
\details{
The function calculates how the signal is enriched in the target by weighting
the distance to the target.

For a numeric vector, assume the vector is denoted as combination of three sub-vectors
\code{c(x1, x2, x3)} with length \code{n1}, \code{n2} and \code{n3}, 
where \code{x1} are data points in upstream windows, \code{x2} are data points in target windows and 
\code{x3} are data points in downstream windows, the enriched score is calcualted as

sum(x_1i* i/n1) + sum(x_3j* (n3 - j + 1)/n3) + sum(x_2k * abs(n2/2 - abs(k - n2/2)))

where the first two terms are the distance to the start or end position of the target
by weighting the distance to the position that if it is closer to the start or end position
of the target, it has higher weight. The second term weight the distance to the center point
of the target and similar, if it is closer to the center position, it has higher weight.
}
\seealso{
This \code{\link{enriched_score}} is the default scoring function for \code{score_fun} argument in \code{\link{EnrichedHeatmap}}
function. It is also an example function for implementing customized scoreing function.
Basically, to be a score function which calculates enriched score, it should accept three arguments
which are the values in upstream windows, the target windows and downstream windows 
The user-defined function should return a single value. Rows are sorted decreasingly by the enriched scores.
}
\value{
A numeric value.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
enriched_score(c(1, 2, 3), c(1, 2, 1), c(3, 2, 1))
enriched_score(c(3, 2, 1), c(2, 1, 2), c(1, 2, 3))
}
