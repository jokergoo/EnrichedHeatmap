\name{extract_anno_enriched}
\alias{extract_anno_enriched}
\title{
Extarct Enrichment Annotation Graphics as a Separated Plot
}
\description{
Extarct Enrichment Annotation Graphics as a Separated Plot
}
\usage{
extract_anno_enriched(ht_list, which = NULL, newpage = TRUE, padding = NULL)
}
\arguments{

  \item{ht_list}{The heatmap list returned by \code{\link{draw,HeatmapList-method}}.}
  \item{which}{The index of enriched heatmap in the heatmap list. The value can be an integer index or a character index (the name of the heatmap).}
  \item{newpage}{Whether call \code{\link[grid]{grid.newpage}} to create a new page?}
  \item{padding}{Padding of the plot.}

}
\details{
The extracted plot is exactly the same as that on the enriched heatmap.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
# There is no example
NULL
}
