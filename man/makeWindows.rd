\name{makeWindows}
\alias{makeWindows}
\title{
Split regions into windows  


}
\description{
Split regions into windows  


}
\usage{
makeWindows(gr, w = NULL, k = NULL, direction = c("normal", "reverse"),
    short.keep = FALSE)
}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object.}
  \item{w}{window size, a value larger than 1 means the number of base pairs and a value between 0 and 1 is the percent to the current region.}
  \item{k}{number of partitions for each region. If it is set, all other arguments are ignored.}
  \item{direction}{where to start the splitting. See 'Details' section.}
  \item{short.keep}{if the the region can not be splitted equally under the window size,  whether to keep the windows that are smaller than the window size. See 'Details' section.}

}
\details{
Following illustrates the meaning of \code{direction} and \code{short.keep}:  

  \preformatted{
    ----------  a region
    aaabbbccc   direction = "normal",  short.keep = FALSE
    aaabbbcccd  direction = "normal",  short.keep = TRUE
     aaabbbccc  direction = "reverse", short.keep = FALSE
    abbbcccddd  direction = "reverse", short.keep = TRUE
  }

There is an additional column \code{.row} attached which contains the correspondance between small windows and original regions in \code{gr}  


}
\value{
A \code{\link[GenomicRanges]{GRanges}} object.  


}
\author{
Zuguang gu <z.gu@dkfz.de>  


}
\section{Example}{
gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 21), end = c(10, 20, 30))) makeWindows(gr, w = 2) makeWindows(gr, w = 0.2) makeWindows(gr, w = 3) makeWindows(gr, w = 3, direction = "reverse") makeWindows(gr, w = 3, short.keep = TRUE) makeWindows(gr, w = 3, direction = "reverse", short.keep = TRUE) makeWindows(gr, w = 12) makeWindows(gr, w = 12, short.keep = TRUE)  

makeWindows(gr, k = 2) makeWindows(gr, k = 3)  

gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 31), end = c(10, 30, 70))) makeWindows(gr, w = 2) makeWindows(gr, w = 0.2)  


}
\examples{
gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 21), end = c(10, 20, 30)))
makeWindows(gr, w = 2)
makeWindows(gr, w = 0.2)
makeWindows(gr, w = 3)
makeWindows(gr, w = 3, direction = "reverse")
makeWindows(gr, w = 3, short.keep = TRUE)
makeWindows(gr, w = 3, direction = "reverse", short.keep = TRUE)
makeWindows(gr, w = 12)
makeWindows(gr, w = 12, short.keep = TRUE)

makeWindows(gr, k = 2)
makeWindows(gr, k = 3)

gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 31), end = c(10, 30, 70)))
makeWindows(gr, w = 2)
makeWindows(gr, w = 0.2)
}
