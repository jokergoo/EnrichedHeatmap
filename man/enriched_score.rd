\name{enriched_score}
\alias{enriched_score}
\title{
Enriched scores for matrix rows

}
\description{
Enriched scores for matrix rows

}
\usage{
enriched_score(x1, x2, x3)}
\arguments{

  \item{x1}{a vector corresponding to values in upstream windows}
  \item{x2}{a vector corresponding to values in target windows}
  \item{x3}{a vector corresponding to values in downstream windows}
}
\details{
After obtaining a matrix from \code{\link{normalizeToMatrix}}, enrichment scores
are calculated to order rows in the matrix.

}
\value{
a numeric value.

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
\examples{
NULL

}
