\name{normalizeToMatrix}
\alias{normalizeToMatrix}
\title{
Normalize Associations between Genomic Signals and Target Regions into a Matrix
}
\description{
Normalize Associations between Genomic Signals and Target Regions into a Matrix
}
\usage{
normalizeToMatrix(signal, target, extend = 5000, w = max(extend)/100,
    value_column = NULL, mapping_column = NULL, background = ifelse(smooth, NA, 0), 
    empty_value = NULL, mean_mode = c("absolute", "weighted", "w0", "coverage"), 
    include_target = any(width(target) > 1),
    target_ratio = min(c(0.4, mean(width(target))/(sum(extend) + mean(width(target))))),
    k = min(c(20, min(width(target)))), smooth = FALSE, smooth_fun = default_smooth_fun,
    keep = c(0, 1), limit = NULL, trim = NULL, flip_upstream = FALSE, verbose = TRUE)
}
\arguments{

  \item{signal}{A \code{\link[GenomicRanges]{GRanges-class}} object.}
  \item{target}{A \code{\link[GenomicRanges]{GRanges-class}} object.}
  \item{extend}{Extended base pairs to the upstream and/or downstream of \code{target}. It can be a vector of length one or two. Length one means same extension to the upstream and downstream.}
  \item{w}{Window size for splitting upstream and downstream, measured in base pairs}
  \item{value_column}{Column index in \code{signal} that is mapped to colors. If it is not set, it assumes values for all signal regions are 1.}
  \item{mapping_column}{Mapping column to restrict overlapping between \code{signal} and \code{target}. By default it tries to look for all regions in \code{signal} that overlap with every target.}
  \item{background}{Values for windows that don't overlap with \code{signal}. }
  \item{empty_value}{Deprecated, please use \code{background} instead.}
  \item{mean_mode}{When a window is not perfectly overlapped to \code{signal}, how to summarize  values to the window. See 'Details' section for a detailed explanation.}
  \item{include_target}{Whether include \code{target} in the heatmap? If the width of all regions in \code{target} is 1, \code{include_target} is enforced to \code{FALSE}.}
  \item{target_ratio}{The ratio of \code{target} columns in the normalized matrix. If the value is 1, \code{extend} will be reset to 0.}
  \item{k}{Number of windows only when \code{target_ratio = 1} or \code{extend == 0}, otherwise ignored.}
  \item{smooth}{Whether apply smoothing on rows in the matrix?}
  \item{smooth_fun}{The smoothing function that is applied to each row in the matrix. This self-defined function accepts a numeric vector (may contain \code{NA} values) and returns a vector with same length. If the smoothing is failed, the function should call \code{\link[base]{stop}} to throw errors so that \code{\link{normalizeToMatrix}} can catch how many rows are failed in smoothing.  See the default \code{\link{default_smooth_fun}} for example.}
  \item{keep}{Percentiles in the normalized matrix to keep. The value is a vector of two percent values. Values less than the first percentile is replaces with the first pencentile and values larger than the second percentile is replaced with the second percentile.}
  \item{limit}{Similar as \code{keep}, but it provides boundary for absolute values. The value should be a vector of length two.}
  \item{trim}{Deprecated, please use \code{keep} instead.}
  \item{flip_upstream}{Sometimes whether the signals are on the upstream or the downstream of the targets are not important and users only want to show the relative distance to targets. If the value is set to \code{TRUE}, the upstream part in the normalized matrix is flipped and added to the downstream part The flipping is only allowed when the targets are single-point targets or the targets are excluded in the normalized matrix (by setting \code{include_target = FALSE}). If the extension for the upstream and downstream is not equal, the smaller extension is used for the final matrix.}
  \item{verbose}{Whether to print help messages.}

}
\details{
In order to visualize associations between \code{signal} and \code{target}, the data is transformed into a matrix
and visualized as a heatmap by \code{\link{EnrichedHeatmap}} afterwards.

Upstream and downstream also with the target body are splitted into a list of small windows and overlap
to \code{signal}. Since regions in \code{signal} and small windows do not always 100 percent overlap, there are four different averaging modes:

Following illustrates different settings for \code{mean_mode} (note there is one signal region overlapping with other signals):

  \preformatted{
      40      50     20     values in signal regions
    ++++++   +++    +++++   signal regions
           30               values in signal region
         ++++++             signal region
      =================     a window (17bp), there are 4bp not overlapping to any signal regions.
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
  \item{\code{failed_rows}}{index of rows which are failed after smoothing}
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
