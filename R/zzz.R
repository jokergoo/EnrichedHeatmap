.onAttach = function(libname, pkgname) {
    version = packageDescription(pkgname, fields = "Version")

  	msg = paste0("========================================
", pkgname, " version ", version, "
Bioconductor page: http://bioconductor.org/packages/EnrichedHeatmap/
Github page: https://github.com/jokergoo/EnrichedHeatmap
Documentation: http://bioconductor.org/packages/EnrichedHeatmap/
========================================
")	
    # suppressPackageStartupMessages(require(ComplexHeatmap))
    packageStartupMessage(msg)
}
