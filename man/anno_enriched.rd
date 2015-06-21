\name{anno_enriched}
\alias{anno_enriched}
\title{
Annotation function to show the enrichment

}
\description{
Annotation function to show the enrichment

}
\usage{
anno_enriched(mat, gp = gpar(col = "red"), pos_line = TRUE, pos_line_gp = gpar(lty = 2),
    yaxis = TRUE, ylim = NULL, yaxis_side = "right", yaxis_gp = gpar(fontsize = 8), show_error = FALSE)}
\arguments{

  \item{mat}{tha matrix returned by \code{\link{normalizeToMatrix}}}
  \item{gp}{graphical parameters for the line}
  \item{pos_line}{whether draw vertical lines which represent the position of \code{target}}
  \item{pos_line_gp}{graphical parameters for lines}
  \item{yaxis}{whether show yaxis}
  \item{ylim}{ranges on y-axis}
  \item{yaxis_side}{side of y-axis}
  \item{yaxis_gp}{graphical parameters for yaxis}
  \item{show_error}{whether show error regions which are 1 sd to the mean value}
}
\details{
It should only be placed as column annotation of the Enriched Heatmap.

}
