\name{dist_by_closeness}
\alias{dist_by_closeness}
\title{
Distance by closeness
}
\description{
Distance by closeness
}
\usage{
dist_by_closeness(mat)
}
\arguments{

  \item{mat}{a numeric matrix where the distance is calculated by rows}

}
\details{
For two rows in the matrix, assume x_1, x_2, ..., x_n1 are the column index of none-zero values in row 1
and y_1, y_2, ... y_n2 are the column index for non-zero values in row 2, 
the distance between the two rows based on the closeness is calculated as:

d_closeness = sum_i sum_j(|x_i - y_j|) / (n_1*n_2)
}
\value{
A \code{\link[stats]{dist}} object
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
x1 = c(0, 0, 0, 0, 1, 1, 1, 0, 0, 0)
x2 = c(0, 0, 0, 1, 1, 1, 0, 0, 0, 0)
x3 = c(1, 0, 0, 0, 1, 1, 0, 0, 0, 0)
m = rbind(x1, x2, x3)
dist(m)
dist_by_closeness(m)
}
