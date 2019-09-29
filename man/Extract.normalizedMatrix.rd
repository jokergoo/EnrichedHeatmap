\name{[.normalizedMatrix}
\alias{[.normalizedMatrix}
\alias{Extract.normalizedMatrix}
\title{
Subset normalized matrix by rows
}
\description{
Subset normalized matrix by rows
}
\usage{
\method{[}{normalizedMatrix}(x, i, j, drop = FALSE)
}
\arguments{

  \item{x}{the normalized matrix returned by \code{\link{normalizeToMatrix}}}
  \item{i}{row index}
  \item{j}{column index}
  \item{drop}{whether drop the dimension}

}
\value{
A \code{normalizedMatrix} class object.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
