\name{makeWindows}
\alias{makeWindows}
\title{
Split regions into windows
}
\description{
Split regions into windows
}
\usage{
makeWindows(query, w = NULL, k = NULL, direction = c("normal", "reverse"),
    short.keep = FALSE)
}
\arguments{

  \item{query}{a \code{\link[GenomicRanges]{GRanges}} object.}
  \item{w}{window size, a value larger than 1 means the number of base pairs and a value between 0 and 1 is the percent to the current region.}
  \item{k}{number of partitions for each region. If it is set, all other arguments are ignored.}
  \item{direction}{where to start the splitting. See 'Details' section.}
  \item{short.keep}{if the the region can not be split equally under the window size,  the argument controls whether to keep the windows that are smaller than the window size. See 'Details' section.}

}
\details{
Following illustrates the meaning of \code{direction} and \code{short.keep}:

  \preformatted{
    -->-->-->-  one region, split by 3bp window (">" means the direction of the sequence)
    aaabbbccc   direction = "normal",  short.keep = FALSE
    aaabbbcccd  direction = "normal",  short.keep = TRUE
     aaabbbccc  direction = "reverse", short.keep = FALSE
    abbbcccddd  direction = "reverse", short.keep = TRUE  }
}
\value{
A \code{\link[GenomicRanges]{GRanges}} object with two additional columns attached:

\itemize{
  \item \code{.i_query} which contains the correspondance between small windows and original regions in \code{query}
  \item \code{.i_window} which contains the index of the small window on the current region.
}
}
\author{
Zuguang gu <z.gu@dkfz.de>
}
\examples{
query = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 21), end = c(10, 20, 30)))
makeWindows(query, w = 2)
makeWindows(query, w = 0.2)
makeWindows(query, w = 3)
makeWindows(query, w = 3, direction = "reverse")
makeWindows(query, w = 3, short.keep = TRUE)
makeWindows(query, w = 3, direction = "reverse", short.keep = TRUE)
makeWindows(query, w = 12)
makeWindows(query, w = 12, short.keep = TRUE)
makeWindows(query, k = 2)
makeWindows(query, k = 3)
query = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 31), end = c(10, 30, 70)))
makeWindows(query, w = 2)
makeWindows(query, w = 0.2)
}
