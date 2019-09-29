\name{default_smooth_fun}
\alias{default_smooth_fun}
\title{
Default Smoothing function
}
\description{
Default Smoothing function
}
\usage{
default_smooth_fun(x)
}
\arguments{

  \item{x}{Input numeric vector.}

}
\details{
The smoothing function is applied to every row in the normalized matrix. For this default smoothing function,
\code{\link[locfit]{locfit}} is first tried on the vector. If there is error, \code{\link[stats]{loess}} smoothing is tried afterwards.
If both smoothing are failed, there will be an error.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
