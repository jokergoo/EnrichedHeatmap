\name{default_smooth_fun}
\alias{default_smooth_fun}
\title{
Default smooth function
}
\description{
Default smooth function
}
\usage{
default_smooth_fun(x)
}
\arguments{

  \item{x}{input numeric vector}

}
\details{
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
