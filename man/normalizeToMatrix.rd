\name{normalizeToMatrix}
\alias{normalizeToMatrix}
\title{
Normalize associations between genomic signals and target regions into a matrix
}
\description{
Normalize associations between genomic signals and target regions into a matrix
}
\usage{
normalizeToMatrix(signal, target, extend = 5000, w = extend/50, value_column = NULL,
    mapping_column = NULL, empty_value = 0, mean_mode = c("absolute", "weighted", "w0"),
    include_target = any(width(target) > 1), target_ratio = 0.1, smooth = FALSE,
    s = 1, trim = 0.01, ...)
}
\arguments{

  \item{signal}{a \code{\link[GenomicRanges]{GRanges}} object which is the genomic signals.}
  \item{target}{a \code{\link[GenomicRanges]{GRanges}} object.}
  \item{extend}{extended base pairs to the upstream and downstream of \code{target}. It can be a vector of length one or two. If it is length one, it means extension to the upstream and downstream are the same.}
  \item{w}{window size for splitting upstream and downstream, and probably \code{target} itself.}
  \item{value_column}{column index in \code{signal} that will be mapped to colors. If it is \code{NULL}, an internal column which all contains 1 will be attached.}
  \item{mapping_column}{mapping column to restrict overlapping between \code{signal} and \code{target}. By default it tries to look for all regions in \code{signal} that overlap with every target.}
  \item{empty_value}{values for small windows that don't overlap with \code{signal}. }
  \item{mean_mode}{when a window is not perfectly overlapped to \code{signal}, how to correspond  the values to this window. See 'Details' section for a detailed explanation.}
  \item{include_target}{whether include \code{target} in the heatmap. If the width of all regions in \code{target} is 1, \code{include_target} is enforced to \code{FALSE}.}
  \item{target_ratio}{the ratio of width of \code{target} part compared to the full heatmap}
  \item{smooth}{whether apply smoothing on rows in the matrix. The smoothing is applied by \code{\link[locfit]{locfit}}. Please note the data range will change, you need to adjust values in the new matrix afterward.}
  \item{s}{\code{\link[GenomicRanges]{findOverlaps}} sometimes uses a lot of memory. \code{target} is splitted into \code{s} parts and each part is processed serialized (note it will be slow!).}
  \item{trim}{percent of extreme values to remove, currently it is disabled.}
  \item{...}{pass to \code{\link[locfit]{locfit}}}

}
\details{
In order to visualize associations between \code{signal} and \code{target}, the data is transformed into a matrix
and visualized as a heatmap afterward.

Upstream and downstream also with the target body are splitted into a list of small windows and overlap
to \code{signal}. Since regions in \code{signal} and small windows do not always 100 percent overlap, averaging should be applied.

Following illustrates different settings for \code{mean_mode}:

  \preformatted{
       4      5      2     values in signal
    ++++++   +++   +++++   signal
      ================     window (16bp)
        4     3     3      overlap

    absolute: (4 + 5 + 2)/3
    weighted: (4*4 + 5*3 + 2*3)/(4 + 3 + 3)
    w0:       (4*4 + 5*3 + 2*3)/16  }
}
\value{
A matrix with following additional attributes:

\describe{
  \item{upstream_index}{column index corresponding to upstream of \code{target}}
  \item{target_index}{column index corresponding to \code{target}}
  \item{downstream_index}{column index corresponding to downstream of \code{target}}
  \item{extend}{extension on upstream and downstream}
  \item{smooth}{whether smoothing was applied on the matrix}
}

The matrix is wrapped into a simple \code{normalizeToMatrix} class.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
signal = GRanges(seqnames = "chr1", 
	  ranges = IRanges(start = c(1, 4, 7, 11, 14, 17, 21, 24, 27),
                     end = c(2, 5, 8, 12, 15, 18, 22, 25, 28)),
    score = c(1, 2, 3, 1, 2, 3, 1, 2, 3))
target = GRanges(seqnames = "chr1", ranges = IRanges(start = 10, end = 20))
normalizeToMatrix(signal, target, extend = 10, w = 2)
normalizeToMatrix(signal, target, extend = 10, w = 2, include_target = TRUE)
normalizeToMatrix(signal, target, extend = 10, w = 2, value_column = "score")
}
