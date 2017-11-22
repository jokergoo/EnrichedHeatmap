\name{getSignalsFromList}
\alias{getSignalsFromList}
\title{
Get signals from a list
}
\description{
Get signals from a list
}
\usage{
getSignalsFromList(lt, fun = function(x) mean(x, na.rm = TRUE))
}
\arguments{

  \item{lt}{a list of normalized matrices which are returned by \code{\link{normalizeToMatrix}}. Matrices in the list should be generated with same settings (e.g. they should use same target regions, same extension to targets and same number of windows).}
  \item{fun}{a user-defined function to summarize signals.}

}
\details{
Let's assume you have a list of histone modification signals for different samples and you want
to visualize the mean pattern across samples. You can first normalize histone mark signals for each sample and then
calculate means values across all samples. In following example code, \code{hm_gr_list} is a list of \code{GRanges} objects
which contain positions of histone modifications, \code{tss} is a \code{GRanges} object containing positions of gene TSS.

  \preformatted{
    mat_list = NULL
    for(i in seq_along(hm_gr_list)) \{
        mat_list[[i]] = normalizeToMatrix(hm_gr_list[[i]], tss, value_column = "density")
    \}  }

If we compress the list of matrices as a three-dimension array where the first dimension corresponds to genes,
the second dimension corresponds to windows and the third dimension corresponds to samples, the mean signal
across all sample can be calculated on the third dimension. Here \code{\link{getSignalsFromList}} simplifies this job.

Applying \code{getSignalsFromList()} to \code{mat_list}, it gives a new normalized matrix which contains mean signals across all samples and can
be directly used in \code{EnrichedHeatmap()}.

  \preformatted{
    mat_mean = getSignalsFromList(mat_list)
    EnrichedHeatmap(mat_mean)  }

The correlation between histone modification and gene expression can
also be calculated on the third dimension of the array. In the user-defined function \code{fun}, \code{x} is the vector for gene i
and window j in the array, and \code{i} is the index of current gene.

  \preformatted{
    mat_corr = getSignalsFromList(mat_list, 
        fun = function(x, i) cor(x, expr[i, ], method = "spearman"))  }

Then \code{mat_corr} here can be used to visualize how gene expression is correlated to histone modification around TSS.

  \preformatted{
    EnrichedHeatmap(mat_corr)  }
}
\value{
A \code{\link{normalizeToMatrix}} object which can be directly used for \code{\link{EnrichedHeatmap}}.
}
\author{
Zuguang Gu <z.gu@dkfz.de>
}
\examples{
NULL
}
