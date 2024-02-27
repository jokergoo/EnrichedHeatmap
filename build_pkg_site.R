if(file.exists(".Rbuildignore")) {
	ln = readLines(".Rbuildignore")
	if(!any(ln == "^_pkgdown\\.yml$")) {
		ln = c(ln, "^_pkgdown\\.yml$")
	}
	if(!any(ln == "^docs$")) {
		ln = c(ln, "^docs$")
	}
	if(!any(ln == "^pkgdown$")) {
		ln = c(ln, "^pkgdown$")
	}
	if(!any(ln == "build_pkg_site.R")) {
		ln = c(ln, "build_pkg_site.R")
	}
	writeLines(ln, ".Rbuildignore")
} else {
	writeLines("
^_pkgdown\\.yml$
^docs$
^pkgdown$", ".Rbuildignore")
}

pkgname = read.dcf("DESCRIPTION")[1, "Package"]

vig_files = list.files(path = "vignettes", full.names = TRUE)
vig_files = vig_files[basename(vig_files) != paste0(pkgname, ".Rmd")]

ln = readLines(".Rbuildignore")

for(v in vig_files) {
	if(!any(ln == v)) {
		ln = c(ln, v)
	}
}
writeLines(ln, ".Rbuildignore")


options(rmarkdown.html_vignette.check_title = FALSE)
pkgdown::build_site(run_dont_run = TRUE)
