\name{default_smooth_fun}
\alias{default_smooth_fun}
\title{
<<<<<<< HEAD
Default smoothing function
}
\description{
Default smoothing function
=======
Default smooth function
}
\description{
Default smooth function
>>>>>>> bioc/master
}
\usage{
default_smooth_fun(x)
}
\arguments{

  \item{x}{input numeric vector}

}
\details{
<<<<<<< HEAD
The smoothing function is applied to every row in the normalized matrix. For this default smoothing function,
=======
The smooth function is applied to every row in the normalized matrix. For this default smooth function,
>>>>>>> bioc/master
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
