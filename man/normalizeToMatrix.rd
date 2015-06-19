\name{normalizeToMatrix}
\alias{normalizeToMatrix}
\title{
Normalize regions to a matrix

}
\description{
Normalize regions to a matrix

}
\usage{
normalizeToMatrix(gr, center, extend = 5000, w = extend/50, value_column = NULL,
    empty_value = 0, mean_mode = c("absolute", "weighted", "w0"), show_body = any(width(center) > 1),
    body_ratio = 0.1, smooth = FALSE, span = 0.75)}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{center}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{extend}{extension to the upstream and downstream of \code{center}. It can be a vector of length one or two.}
  \item{w}{window size for splitting upstream and downstream}
  \item{value_column}{index for column in \code{gr} that map to colors. If the value is \code{NULL}, an internal columnwhich contains 1 will be attached.}
  \item{empty_value}{values for windows that don't overlap with \code{gr}}
  \item{mean_mode}{when a window overlaps with more than one regions in \code{gr}, how to calculate the mean values in this window.}
  \item{show_body}{whether show \code{center}}
  \item{body_ratio}{the ratio of width of \code{center} in the heatmap}
  \item{smooth}{whether apply smoothing in every row in the matrix. The smoothing is applied by \code{\link[stats]{loess}}}
  \item{span}{degree of smoothing, pass to \code{\link[stats]{loess}}.}
}
\details{
Following illustrates different settings for \code{mean_mode}:

  \preformatted{
       4      5      2     values
    ++++++   +++   +++++   gr
      ================     window (16bp)

    absolute: (4 + 5 + 2)/3
    weighted: (4*4 + 5*3 + 2*3)/(4 + 3 + 3)
    w0:       (4*4 + 5*3 + 2*3)/16
  }

}
\value{
A matrix with following additional attributes:

\describe{
  \item{upstream_index}{column index corresponding to upstream}
  \item{body_index}{column index corresponding to body}
  \item{downstream_index}{column index corresponding to downstream}
  \item{extend}{extension on upstream and downstream}
}

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
\examples{
gr = GRanges(seqnames = "chr1", 
	  ranges = IRanges(start = c(1, 4, 7, 11, 14, 17, 21, 24, 27),
                     end = c(2, 5, 8, 12, 15, 18, 22, 25, 28)))
center = GRanges(seqnames = "chr1", ranges = IRanges(start = 10, end = 20))
normalizeToMatrix(gr, center, extend = 10, w = 2)

}
