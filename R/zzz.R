.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")

  	msg = paste0("========================================
", pkgname, " version ", version, "
Bioconductor page: http://bioconductor.org/packages/EnrichedHeatmap/
Github page: https://github.com/jokergoo/EnrichedHeatmap
Documentation: http://bioconductor.org/packages/EnrichedHeatmap/

If you use it in published research, please cite:
Gu, Z. EnrichedHeatmap: an R/Bioconductor package for comprehensive 
visualization of genomic signal associations. BMC Genomics 2018.
========================================
")	
    # suppressPackageStartupMessages(require(ComplexHeatmap))
    packageStartupMessage(msg)
}
