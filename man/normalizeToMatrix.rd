\name{normalizeToMatrix}
\alias{normalizeToMatrix}
\title{
Normalize associations between genomic signals and target regions into a matrix
}
\description{
Normalize associations between genomic signals and target regions into a matrix
}
\usage{
normalizeToMatrix(signal, target, extend = 5000, w = max(extend)/50,
    value_column = NULL, mapping_column = NULL, background = ifelse(smooth, NA, 0),
    mean_mode = c("absolute", "weighted", "w0", "coverage"), include_target = any(width(target) > 1),
    target_ratio = min(c(0.6, mean(width(target))/(sum(extend)) + mean(width(target)))),
    k = min(c(20, min(width(target)))), smooth = FALSE, smooth_fun = default_smooth_fun,
    keep = c(0, 1))
}
\arguments{

  \item{signal}{a \code{\link[GenomicRanges]{GRanges}} object.}
  \item{target}{a \code{\link[GenomicRanges]{GRanges}} object.}
  \item{extend}{extended base pairs to the upstream and downstream of \code{target}. It can be a vector of length one or two. Length one means same extension to the upstream and downstream.}
  \item{w}{window size for splitting upstream and downstream.}
  \item{value_column}{column index in \code{signal} that is mapped to colors. If it is not set, it assumes values for all all signal regiosn are 1.}
  \item{mapping_column}{mapping column to restrict overlapping between \code{signal} and \code{target}. By default it tries to look for all regions in \code{signal} that overlap with every target.}
  \item{background}{values for windows that don't overlap with \code{signal}. }
  \item{mean_mode}{when a window is not perfectly overlapped to \code{signal}, how to summarize  values to the window. See 'Details' section for a detailed explanation.}
  \item{include_target}{whether include \code{target} in the heatmap. If the width of all regions in \code{target} is 1, \code{include_target} is enforced to \code{FALSE}.}
  \item{target_ratio}{the ratio of \code{target} in the full heatmap. If the value is 1, \code{extend} will be reset to 0.}
  \item{k}{number of windows only when \code{target_ratio = 1} or \code{extend == 0}, otherwise ignored.}
  \item{smooth}{whether apply smoothing on rows in the matrix. }
  \item{smooth_fun}{the smoothing function that is applied to each row in the matrix. This self-defined function accepts a numeric vector (may contains \code{NA} values) and returns a vector with same length. If the smoothing is failed, the function should call \code{\link[base]{stop}} to throw errors so that \code{\link{normalizeToMatrix}} can catch how many rows are failed in smoothing.  See the default \code{\link{default_smooth_fun}} for example.}
  \item{keep}{values in the normalized matrix to keep. The value is a vector of two percent values. Values less than the first percentile is replaces with the first pencentile and values larger than the second percentile is replaced with the second percentile.}

}
\details{
In order to visualize associations between \code{signal} and \code{target}, the data is transformed into a matrix
and visualized as a heatmap by \code{\link{EnrichedHeatmap}} afterwards.

Upstream and downstream also with the target body are splitted into a list of small windows and overlap
to \code{signal}. Since regions in \code{signal} and small windows do not always 100 percent overlap, there are four different average modes:

Following illustrates different settings for \code{mean_mode} (note there is one signal region overlapping with other signals):

  \preformatted{
      40      50     20     values in signal
    ++++++   +++    +++++   signal
           30               values in signal
         ++++++             signal
      =================     window (17bp), there are 4bp not overlapping to any signal region.
        4  6  3      3      overlap

    absolute: (40 + 30 + 50 + 20)/4
    weighted: (40*4 + 30*6 + 50*3 + 20*3)/(4 + 6 + 3 + 3)
    w0:       (40*4 + 30*6 + 50*3 + 20*3)/(4 + 6 + 3 + 3 + 4)
    coverage: (40*4 + 30*6 + 50*3 + 20*3)/17  }
}
\value{
A matrix with following additional attributes:

\describe{
  \item{\code{upstream_index}}{column index corresponding to upstream of \code{target}}
  \item{\code{target_index}}{column index corresponding to \code{target}}
  \item{\code{downstream_index}}{column index corresponding to downstream of \code{target}}
  \item{\code{extend}}{extension on upstream and downstream}
  \item{\code{smooth}}{whether smoothing was applied on the matrix}
  \item{\code{failed_rows}}{index of rows which are failed for smoothing}
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
