\name{enriched_score}
\alias{enriched_score}
\title{
Enriched Scores
}
\description{
Enriched Scores
}
\usage{
enriched_score(mat)
}
\arguments{

  \item{mat}{A normalized matrix from \code{\link{normalizeToMatrix}}.}

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
\value{
A numeric vector.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
