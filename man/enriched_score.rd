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
The function calculates how the signal is enriched in the targets.
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
