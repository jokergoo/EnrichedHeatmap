 
stop_wrap = ComplexHeatmap:::stop_wrap
warning_wrap = ComplexHeatmap:::warning_wrap


recycle_gp = function(gp, n = 1) {
	g = lapply(gp, function(x) {
		if(length(x) == 1 && n > 1) {
			rep(x, n)
		} else {
			x
		}
	})
	class(g) = "gpar"
	return(g)
}

subset_gp = function(gp, i = 1) {
	g = lapply(gp, function(x) {
		x[i]
	})
	class(g) = "gpar"
	return(g)
}
