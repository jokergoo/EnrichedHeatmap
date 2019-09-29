\name{discretize}
\alias{discretize}
\title{
Discretize a Continuous Matrix to a Discrete Matrix
}
\description{
Discretize a Continuous Matrix to a Discrete Matrix
}
\usage{
discretize(mat, rule, right_closed = FALSE)
}
\arguments{

  \item{mat}{A normalize matrix from \code{\link{normalizeToMatrix}}.}
  \item{rule}{A list of intervals which provide mapping between continuous values to discrete values. Note the order of intervals determines the order of corresponding discrete levels.}
  \item{right_closed}{Is the interval right closed?}

}
\details{
Assuming we have a normalized matrix with both positive values and negative values, we only 
want to see the enrichment of the windows/regions showing significant positive values and 
negative values and we are only interested in the direction of the values while not the value itself,
then we can define the rule as:

  \preformatted{
    rule = list(
        "positive" = c(0.5, Inf),
        "negative" = c(-Inf, -0.5)
    )  }

And we can convert the continuous matrix to a discrete matrix and visualize it:

  \preformatted{
    mat2 = discretize(mat, rule)
    EnrichedHeatmap(mat2, col = c("positive" = "red", "negative" = "green"))  }

Another example is to discretize the signals to discrete levels according to the intensities:

  \preformatted{
    rule = list(
        "very_high" = c(100, Inf),
        "high" = c(50, 100),
        "intermediate" = c(25, 50),
        "low" = c(1e-6, 25)
    )  }
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
