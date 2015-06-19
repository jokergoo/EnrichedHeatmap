

# == title
# Normalize regions to a matrix
#
# == param
# -gr a `GenomicRanges::GRanges` object
# -center a `GenomicRanges::GRanges` object
# -extend extension to the upstream and downstream of ``center``. It can be a vector of length one or two.
# -w window size for splitting upstream and downstream
# -value_column index for column in ``gr`` that map to colors. If the value is ``NULL``, an internal column
#         which contains 1 will be attached.
# -empty_value values for windows that don't overlap with ``gr``
# -mean_mode when a window overlaps with more than one regions in ``gr``, how to calculate 
#       the mean values in this window.
# -show_body  whether show ``center``
# -body_ratio  the ratio of width of ``center`` in the heatmap
# -smooth whether apply smoothing in every row in the matrix. The smoothing is applied by `stats::loess`
# -span degree of smoothing, pass to `stats::loess`.
#
# == details
# Following illustrates different settings for ``mean_mode``:
#
#        4      5      2     values
#     ++++++   +++   +++++   gr
#       ================     window (16bp)
#
#     absolute: (4 + 5 + 2)/3
#     weighted: (4*4 + 5*3 + 2*3)/(4 + 3 + 3)
#     w0:       (4*4 + 5*3 + 2*3)/16
#
# == value
# A matrix with following additional attributes:
#
# -upstream_index column index corresponding to upstream
# -body_index column index corresponding to body
# -downstream_index column index corresponding to downstream
# -extend extension on upstream and downstream
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr = GRanges(seqnames = "chr1", 
# 	  ranges = IRanges(start = c(1, 4, 7, 11, 14, 17, 21, 24, 27),
#                      end = c(2, 5, 8, 12, 15, 18, 22, 25, 28)))
# center = GRanges(seqnames = "chr1", ranges = IRanges(start = 10, end = 20))
# normalizeToMatrix(gr, center, extend = 10, w = 2)
normalizeToMatrix = function(gr, center, extend = 5000, w = extend/50, value_column = NULL,
    empty_value = 0, mean_mode = c("absolute", "weighted", "w0"), show_body = any(width(center) > 1),
    body_ratio = 0.1, smooth = FALSE, span = 0.75) {
  
	if(any(extend %% w > 0)) {
		stop("`extend` should be divisible by `w`.\n")
	}
  
	if(length(extend) == 1) extend = c(extend, extend)
  	
  	if(all(width(center) <= 1)) {
  		# do not need to separate upstream and downstream
  		suppressWarnings(both <- promoters(center, upstream = extend[1], downstream = extend[2]))
		strand(both) = "*"

		mat_both = makeMatrix(gr, both, w = w, value_column = value_column, empty_value = empty_value, mean_mode = mean_mode)
		i = round(extend[1]/(extend[1] + extend[2]) * ncol(mat_both))  # assume
		if(i < 2 | ncol(mat_both) - i < 2) {
			stop("Maybe `w` is too large or one of `extend` is too small.")
		}
		mat_upstream = mat_both[, 1:i, drop = FALSE]
		mat_downstream = mat_both[, (i+1):ncol(mat_both), drop = FALSE]
	  
  	} else {
		# extend and normalize in upstream 
		suppressWarnings(upstream <- promoters(center, upstream = extend[1], downstream = 0))
		strand(upstream) = "*"
	  
		mat_upstream = makeMatrix(gr, upstream, w = w, value_column = value_column, empty_value = empty_value, mean_mode = mean_mode)
	  
		# extend and normalize in downstream
		if(all(width(center) == 1) && !show_body) {
			e = ifelse(strand(center) == "-", start(center), end(center))
		} else {
			e = ifelse(strand(center) == "-", start(center) - 1, end(center) + 1)
		}
		end_center = GRanges(seqnames = seqnames(center),
	                         ranges = IRanges(start = e, end = e),
	                         strand = strand(center))
		suppressWarnings(downstream <- promoters(end_center, upstream = 0, downstream = extend[2]))
		strand(downstream) = "*"
	  
		mat_downstream = makeMatrix(gr, downstream, w = w, value_column = value_column, empty_value = empty_value,
	                              mean_mode = mean_mode)
	}

	if(show_body) {
		k = (ncol(mat_upstream) + ncol(mat_downstream)) * body_ratio/(1-body_ratio)
		mat_body = makeMatrix(gr, center, k = k, value_column = value_column, empty_value = empty_value, mean_mode = mean_mode)
	} else {
		mat_body = matrix(0, nrow = length(center), ncol = 0)
	}

	if(show_body) {
		mat = cbind(mat_upstream, mat_downstream)
		# apply smoothing on rows in mat
		if(smooth) mat = t(apply(mat, 1, function(x) loess(x ~ seq_along(x), span = span)$fitted))

		attr(mat, "upstream_index") = seq_len(ncol(mat_upstream))
		attr(mat, "body_index") = seq_len(ncol(mat_body)) + ncol(mat_upstream)
		attr(mat, "downstream_index") = seq_len(ncol(mat_downstream)) + ncol(mat_upstream)
		attr(mat, "extend") = extend
	} else {
  		mat = cbind(mat_upstream, mat_body, mat_downstream)
	  	# apply smoothing on rows in mat
		if(smooth) mat = t(apply(mat, 1, function(x) loess(x ~ seq_along(x), span = span)$fitted))
		
		attr(mat, "upstream_index") = seq_len(ncol(mat_upstream))
		attr(mat, "body_index") = numeric(0)
		attr(mat, "downstream_index") = seq_len(ncol(mat_downstream)) + ncol(mat_upstream) + ncol(mat_body)		
		attr(mat, "extend") = extend
	}

	

	return(mat)
}

# 
# -gr input regions
# -reference the upstream part or body part
# -window absolute size (100) or relative size(0.1)
# -value_column
# -mean_mode how to calculate mean value in a window
#
# == example
# gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 4, 7), end = c(2, 5, 8)))
# reference = GRanges(seqnames = "chr1", ranges =IRanges(start = 1, end = 10))
# makeMatrix(gr, reference, w = 2)
#
makeMatrix = function(gr, reference, w = NULL, k = NULL, value_column = NULL, empty_value = 0,
                      mean_mode = c("absolute", "weighted", "w0"), 
                      direction = c("normal", "reverse")) {
  
	if(is.null(value_column)) {
		mcols(gr) = NULL
		mcols(gr)$.value = rep(1, length(gr))
		value_column = ".value"
	}
  
	# split `reference` into small windows
	reference_windows = makeWindows(reference, w = w, k = k, direction = direction)
  
	# overlap `gr` to `reference_windows`
	mtch = findOverlaps(gr, reference_windows)
	mtch = as.matrix(mtch)

	# add a `value` column in `reference_windows` which is the mean value for intersected gr
	# in `reference_window`
	m_gr = gr[ mtch[, 1] ]
	m_reference_windows = reference_windows[ mtch[, 2] ]
	v = mcols(m_gr)[[value_column]]
  
	mean_mode = match.arg(mean_mode)[1]
  
	if(mean_mode == "w0") {
		mintersect = pintersect(m_gr, m_reference_windows)
		p = width(mintersect)/width(m_reference_windows)
		x = tapply(p*v, mtch[, 2], sum)
	} else if(mean_mode == "absolute") {
		x = tapply(v, mtch[, 2], mean)
	} else {
		mintersect = pintersect(m_gr, m_reference_windows)
		w = width(mintersect)
		x = tapply(w*v, mtch[, 2], sum) / tapply(w, mtch[, 2], sum)
	}
  
	v2 = rep(empty_value, length(reference_windows))
	v2[ as.numeric(names(x)) ] = x
  
	reference_windows$.value = v2
  
	# transform into a matrix
	tb = table(reference_windows$.orow)
	reference_strand = strand(reference)
	column_index = mapply(as.numeric(names(tb)), tb, FUN = function(i, n) {
		if(as.vector(reference_strand[i] == "-")) {
			rev(seq_len(n))
		} else {
			seq_len(n)
		}
	})
  
	# is column_index has the same length for all regions in reference?
	# if extension of upstream are the same or split body into k pieces,
	# then all column index has the same length
	# if it is not the same, throw error!
	if(!is.matrix(column_index)) {
		stop("numbers of columns are not the same.")
	}
  
	mat = matrix(empty_value, nrow = length(reference), ncol = dim(column_index)[1])
	mat[ reference_windows$.orow + (as.vector(column_index) - 1)* nrow(mat) ] = reference_windows$.value
  
	return(mat)
}

# == title
# Split regions into windows
#
# == param
# -gr a `GenomicRanges::GRanges` object. Regions in the object will be splitted into windows
# -w window size, a value larger than 1 means the number of base pairs and a value between 0 and 1
#    is the percent to the current region.
# -k number of partitions for each region. If it is set, all other arguments are ignored.
# -direction where to start the splitting. See 'Details' section.
# -short.keep if the the region can not be splitted equally under the window size, 
#             whether to keep the windows that are smaller than the window size. See 'Details' section.
#
# == details
# Following illustrates the meaning of ``direction`` and ``short.keep``:
#
#     ----------  a region
#     aaabbbccc   direction = "normal",  short.keep = FALSE
#     aaabbbcccd  direction = "normal",  short.keep = TRUE
#      aaabbbccc  direction = "reverse", short.keep = FALSE
#     abbbcccddd  direction = "reverse", short.keep = TRUE
#     
# There is an additional column ``.orow`` attached which contains the correspondance between small windows
# and original regions in ``gr``
#
# == value
# A `GenomicRanges::GRanges` object.
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
# == example
# gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 21), end = c(10, 20, 30)))
# makeWindows(gr, w = 2)
# makeWindows(gr, w = 0.2)
# makeWindows(gr, w = 3)
# makeWindows(gr, w = 3, direction = "reverse")
# makeWindows(gr, w = 3, short.keep = TRUE)
# makeWindows(gr, w = 3, direction = "reverse", short.keep = TRUE)
# makeWindows(gr, w = 12)
# makeWindows(gr, w = 12, short.keep = TRUE)
#
# makeWindows(gr, k = 2)
# makeWindows(gr, k = 3)
# 
# gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 31), end = c(10, 30, 70)))
# makeWindows(gr, w = 2)
# makeWindows(gr, w = 0.2)
#
makeWindows = function(gr, w = NULL, k = NULL, direction = c("normal", "reverse"), 
	short.keep = FALSE) {

	direction = match.arg(direction)[1]
  
	if(is.null(w) & is.null(k)) {
		stop("You should define either `w` or `k`.")
	}
	ostart = start(gr)
	oend = end(gr)
  
	if(!is.null(k)) {
		pos = mapply(ostart, oend, FUN = function(s, e) {
			x = seq(s, e, length = k+1)
			y = round(x[-1])
			x = round(x[-length(x)])
			y[-k] = ifelse(y[-k] > x[-k], y[-k] - 1, y[-k])
			return(list(start = x, end = y))
		})
	} else {

		if(direction == "normal") {
			if(w >= 1) {
				w = as.integer(w)
				pos = mapply(ostart, oend, FUN = function(s, e) {
					x = seq(s, e, by = w)
					y = x + w - 1
					y = ifelse(y > e, e, y)
					if(!short.keep) {
						l = (y-x+1) == w
						x = x[l]
						y = y[l]
					}
					return(list(start = x, end = y))
				})
			} else if(w > 0 & w < 1) {
				pos = mapply(ostart, oend, FUN = function(s, e) {
					w = as.integer(round(e - s + 1)*w)
					x = seq(s, e, by = w)
					y = x + w - 1
					y = ifelse(y > e, e, y)
					if(!short.keep) {
						l = (y-x+1) == w
						x = x[l]
						y = y[l]
					}
					return(list(start = x, end = y))
				})
			} else {
				stop("`w` is wrong.")
			}
		} else {
			if(w >= 1) {
				w = as.integer(w)
				pos = mapply(ostart, oend, FUN = function(s, e) {
					y = seq(e, s, by = -w)
					x = y - w + 1
					x = ifelse(x < s, s, x)
					if(!short.keep) {
						l = (y-x+1) == w
						x = x[l]
						y = y[l]
					}
					return(list(start = rev(x), end = rev(y)))
				})
			} else if(w > 0 & w < 1) {
				pos = mapply(ostart, oend, FUN = function(s, e) {
					w = as.integer(round(e - s + 1)*w)
					y = seq(e, s, by = -w)
					x = y - w + 1
					x = ifelse(x < s, s, x)
					if(!short.keep) {
						l = (y-x+1) == w
						x = x[l]
						y = y[l]
					}
					return(list(start = rev(x), end = rev(y)))
				})
			} else {
				stop("`w` is wrong.")
			}
		}
	}
  
	# check start and end

	start = unlist(pos[1, ])
	end = unlist(pos[2, ])
	orow = rep(seq_len(ncol(pos)), times = sapply(pos[1, ], length))
	chr = seqnames(gr)[orow]
	strand = strand(gr)[orow]

	gr = GRanges(seqnames = chr,
		         ranges = IRanges(start = start,
		         	              end = end),
		         strand = strand,
		         .orow = orow)
	return(gr)

}

