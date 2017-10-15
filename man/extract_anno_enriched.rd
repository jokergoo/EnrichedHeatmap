\name{extract_anno_enriched}
\alias{extract_anno_enriched}
\title{
Extarct enrichment annotation graphic as a separate plot
}
\description{
Extarct enrichment annotation graphic as a separate plot
}
\usage{
extract_anno_enriched(ht_list, which = NULL, newpage = TRUE)
}
\arguments{

  \item{ht_list}{the heatmap list returned by \code{\link{draw,EnrichedHeatmapList-method}}}
  \item{which}{the index of enriched heamtap in the heatmap list. The value can be an integer index or a character index}
  \item{newpage}{whether call \code{\link[grid]{grid.newpage}} to create a new page}

}
\details{
The extracted plot is exactly the same as the one on the heatmap.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL

}
