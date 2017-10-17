\name{copyAttr}
\alias{copyAttr}
\title{
Copy attributes to another object
}
\description{
Copy attributes to another object
}
\usage{
copyAttr(x, y)
}
\arguments{

  \item{x}{object 1}
  \item{y}{object 2}

}
\details{
The \code{\link{normalizeToMatrix}} object is actually a matrix but with more additional attributes attached.
When manipulating such matrix, there are some circumstances that the attributes are lost.
This function is used to copy these specific attributes when dealing with the matrix.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
gr = GRanges(seqnames = c("chr5", "chr5"),
	ranges = IRanges(start = c(98, 98),
	                 end = c(104, 104)))
target = GRanges(seqnames = "chr5",
	ranges = IRanges(start = 100, 
		             end = 100))
mat1 = normalizeToMatrix(gr, target, extend = 6, w = 1)
# attributes removed and you cannot use it for EnrichedHeatmap()
mat2 = mat1[]
# copy attributes to mat2 and now mat3 can be used for EnrichedHeatmap()
mat3 = copyAttr(mat1, mat2)
}
