\name{normalizeToMatrix}
\alias{normalizeToMatrix}
\title{
Normalize associations between genomic regions and target regions into a matrix

}
\description{
Normalize associations between genomic regions and target regions into a matrix

}
\usage{
normalizeToMatrix(gr, target, extend = 5000, w = extend/50, value_column = NULL, mapping_column = NULL,
    empty_value = 0, mean_mode = c("absolute", "weighted", "w0"), include_target = any(width(target) > 1),
    target_ratio = 0.1, smooth = FALSE, span = 0.5, s = 1, trim = 0.01)}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{target}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{extend}{extended base pairs to the upstream and downstream of \code{target}. It can be a vector of length one or two.}
  \item{w}{window size for splitting upstream and downstream in \code{target}.}
  \item{value_column}{index for column in \code{gr} that will be mapped to colors. If it is \code{NULL}, an internal columnwhich contains 1 will be attached.}
  \item{mapping_column}{mapping column to restrict overlapping between \code{gr} and \code{target}}
  \item{empty_value}{values for windows that don't overlap with \code{gr}}
  \item{mean_mode}{when a window is not perfectkt matched to one region in \code{gr}, how to calculate the mean values in this window. See 'Details' section for a detailed explanation.}
  \item{include_target}{whether include \code{target} in the heatmap. If the width of all regions in \code{target} is 1, \code{include_target}is enforced to \code{FALSE}.}
  \item{target_ratio}{the ratio of width of \code{target} compared to 'upstream + target + downstream' in the heatmap}
  \item{smooth}{whether apply smoothing in every row in the matrix. The smoothing is applied by \code{\link[stats]{loess}}. Pleasenote the data range will change, you need to adjust values in the new matrix afterwards.}
  \item{span}{degree of smoothing, pass to \code{\link[stats]{loess}}.}
  \item{s}{\code{\link[GenomicRanges]{findOverlaps}} sometimes uses a lot of memory. \code{target} is splitted into \code{s} parts and eachpart is processed serialized.}
  \item{trim}{percent of extreme values to remove}
}
\details{
In order to visualize associations between \code{gr} and \code{target}, the data is transformed into a matrix
and visualized as a heatmap.

Following illustrates different settings for \code{mean_mode}:

  \preformatted{
       4      5      2     values in gr
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
  \item{upstream_index}{column index corresponding to upstream of \code{target}}
  \item{target_index}{column index corresponding to \code{target}}
  \item{downstream_index}{column index corresponding to downstream of \code{target}}
  \item{extend}{extension on upstream and downstream}
  \item{smooth}{whether smoothing was applied on the matrix}
}

}
\author{
Zuguang Gu <z.gu@dkfz.de>

}
\examples{
gr = GRanges(seqnames = "chr1", 
	  ranges = IRanges(start = c(1, 4, 7, 11, 14, 17, 21, 24, 27),
                     end = c(2, 5, 8, 12, 15, 18, 22, 25, 28)),
    score = c(1, 2, 3, 1, 2, 3, 1, 2, 3))
target = GRanges(seqnames = "chr1", ranges = IRanges(start = 10, end = 20))
normalizeToMatrix(gr, target, extend = 10, w = 2)
normalizeToMatrix(gr, target, extend = 10, w = 2, include_target = TRUE)
normalizeToMatrix(gr, target, extend = 10, w = 2, value_column = "score")

}
