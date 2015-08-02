

# == title
# Normalize associations between genomic signals and target regions into a matrix
#
# == param
# -signal a `GenomicRanges::GRanges` object which is the genomic signals.
# -target a `GenomicRanges::GRanges` object.
# -extend extended base pairs to the upstream and downstream of ``target``. It can be a vector of length one or two.
# -w window size for splitting upstream and downstream in ``target``.
# -value_column index for column in ``signal`` that will be mapped to colors. If it is ``NULL``, an internal column
#         which contains 1 will be attached.
# -mapping_column mapping column to restrict overlapping between ``signal`` and ``target``. By default it tries to look for
#           all regions in ``signal`` that overlap with every target.
# -empty_value values for windows that don't overlap with ``signal``. 
# -mean_mode when a window is not perfectly matched to ``signal``, how to calculate 
#       the mean values in this window. See 'Details' section for a detailed explanation.
# -include_target  whether include ``target`` in the heatmap. If the width of all regions in ``target`` is 1, ``include_target``
#               is enforced to ``FALSE``.
# -target_ratio  the ratio of width of ``target`` compared to 'upstream + target + downstream' in the heatmap
# -smooth whether apply smoothing in rows in the matrix. The smoothing is applied by `stats::loess`. Please
#         note the data range will change, you need to adjust values in the new matrix afterwards.
# -span degree of smoothing, pass to `stats::loess`.
# -s `GenomicRanges::findOverlaps` sometimes uses a lot of memory. ``target`` is splitted into ``s`` parts and each
#     part is processed serialized (it will be slow!).
# -trim percent of extreme values to remove
#
# == details
# In order to visualize associations between ``signal`` and ``target``, the data is transformed into a matrix
# and visualized as a heatmap.
#
# Upstream and downstream also with the target body are splitted into a list of small windows and overlap
# to ``signal``. Since regions in ``signal`` and small windows do not always 100 percent overlap, averaging should be applied.
# 
# Following illustrates different settings for ``mean_mode``:
#
#        4      5      2     values in signal
#     ++++++   +++   +++++   signal
#       ================     window (16bp)
#         4     3     3      overlap
#
#     absolute: (4 + 5 + 2)/3
#     weighted: (4*4 + 5*3 + 2*3)/(4 + 3 + 3)
#     w0:       (4*4 + 5*3 + 2*3)/16
#
# == value
# A matrix with following additional attributes:
#
# -upstream_index column index corresponding to upstream of ``target``
# -target_index column index corresponding to ``target``
# -downstream_index column index corresponding to downstream of ``target``
# -extend extension on upstream and downstream
# -smooth whether smoothing was applied on the matrix
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# signal = GRanges(seqnames = "chr1", 
# 	  ranges = IRanges(start = c(1, 4, 7, 11, 14, 17, 21, 24, 27),
#                      end = c(2, 5, 8, 12, 15, 18, 22, 25, 28)),
#     score = c(1, 2, 3, 1, 2, 3, 1, 2, 3))
# target = GRanges(seqnames = "chr1", ranges = IRanges(start = 10, end = 20))
# normalizeToMatrix(signal, target, extend = 10, w = 2)
# normalizeToMatrix(signal, target, extend = 10, w = 2, include_target = TRUE)
# normalizeToMatrix(signal, target, extend = 10, w = 2, value_column = "score")
#
normalizeToMatrix = function(signal, target, extend = 5000, w = extend/50, value_column = NULL, 
	mapping_column = NULL, empty_value = 0, mean_mode = c("absolute", "weighted", "w0"), 
	include_target = any(width(target) > 1), target_ratio = 0.1, smooth = FALSE, 
	span = 0.5, s = 1, trim = 0.01) {

	signal_name = as.character(substitute(signal))
	target_name = as.character(substitute(target))

	if(s > 1) {
		n = length(target)
		if(s > n) s = n
		x = seq(1, n, by = s)
		if(x < n) x = c(x, n)
		start_index = x[-length(x)]
		end_index = x[-1] - 1
		end_index[length(end_index)] = x[length(x)]

		lt = lapply(seq_along(start_index), function(i) {
			normalizeToMatrix(signal, target[ start_index[i]:end_index[i] ], extend = extend, w = w, value_column = value_column, mapping_column = mapping_column,
				empty_value = empty_value, mean_mode = mean_mode, include_target = include_target,
				target_ratio = target_ratio, smooth = smooth, span = span, trim = 0)
		})

		upstream_index = attr(lt[[1]], "upstream_index")
		target_index = attr(lt[[1]], "target_index")
		downstream_index = attr(lt[[1]], "downstream_index")
		extend = attr(lt[[1]], "extend")
		smooth = attr(lt[[1]], "smooth")

		mat = do.call("rbind", lt)

		attr(mat, "upstream_index") = upstream_index
		attr(mat, "target_index") = target_index
		attr(mat, "downstream_index") = downstream_index
		attr(mat, "extend") = extend
		attr(mat, "smooth") = smooth
		attr(mat, "signal_name") = signal_name
		attr(mat, "target_name") = target_name

		if(trim > 0) {
	  		q1 = quantile(mat, trim/2, na.rm = TRUE)
	  		q2 = quantile(mat, 1 - trim/2, na.rm = TRUE)
	  		mat[mat <= q1] = q1
	  		mat[mat >= q2] = q2
	  	}
	  	class(mat) = c("normalizeToMatrix", "matrix")
		return(mat)
		
	}

	if(all(width(target) == 1)) {
		include_target = FALSE
	}
  
	if(any(extend %% w > 0)) {
		stop("`extend` should be divisible by `w`.\n")
	}
  
	if(length(extend) == 1) extend = c(extend, extend)
  	
  	if(all(width(target) <= 1)) {
  		# do not need to separate upstream and downstream
  		suppressWarnings(both <- promoters(target, upstream = extend[1], downstream = extend[2]))

		mat_both = makeMatrix(signal, both, w = w, value_column = value_column, mapping_column = mapping_column, empty_value = empty_value, mean_mode = mean_mode)
		i = round(extend[1]/(extend[1] + extend[2]) * ncol(mat_both))  # assume
		if(i < 2 | ncol(mat_both) - i < 2) {
			stop("Maybe `w` is too large or one of `extend` is too small.")
		}
		mat_upstream = mat_both[, 1:i, drop = FALSE]
		mat_downstream = mat_both[, (i+1):ncol(mat_both), drop = FALSE]
	  
  	} else {
		# extend and normalize in upstream 
		suppressWarnings(upstream <- promoters(target, upstream = extend[1], downstream = 0))
	  
		mat_upstream = makeMatrix(signal, upstream, w = w, value_column = value_column, mapping_column = mapping_column, empty_value = empty_value, mean_mode = mean_mode)
	  
		# extend and normalize in downstream
		if(all(width(target) == 1) && !include_target) {
			e = ifelse(strand(target) == "-", start(target), end(target))
		} else {
			e = ifelse(strand(target) == "-", start(target) - 1, end(target) + 1)
		}
		end_target = GRanges(seqnames = seqnames(target),
	                         ranges = IRanges(start = e, end = e),
	                         strand = strand(target))
		suppressWarnings(downstream <- promoters(end_target, upstream = 0, downstream = extend[2]))
		names(downstream) = names(target)
	  
		mat_downstream = makeMatrix(signal, downstream, w = w, value_column = value_column, mapping_column = mapping_column, empty_value = empty_value,
	                              mean_mode = mean_mode)
	}

	if(include_target) {
		k = (ncol(mat_upstream) + ncol(mat_downstream)) * target_ratio/(1-target_ratio)
		mat_target = makeMatrix(signal, target, k = k, value_column = value_column, mapping_column = mapping_column, empty_value = empty_value, mean_mode = mean_mode)
	} else {
		mat_target = matrix(0, nrow = length(target), ncol = 0)
	}

  	mat = cbind(mat_upstream, mat_target, mat_downstream)
  	# apply smoothing on rows in mat
	if(smooth) mat = t(apply(mat, 1, function(x) loess(x ~ seq_along(x), span = span)$fitted))

	
	upstream_index = seq_len(ncol(mat_upstream))
	target_index = seq_len(ncol(mat_target)) + ncol(mat_upstream)	
	downstream_index = seq_len(ncol(mat_downstream)) + ncol(mat_upstream) + ncol(mat_target)

	attr(mat, "upstream_index") = upstream_index
	attr(mat, "target_index") = target_index
	attr(mat, "downstream_index") = downstream_index
	attr(mat, "extend") = extend
	attr(mat, "smooth") = smooth
	attr(mat, "signal_name") = signal_name
	attr(mat, "target_name") = target_name

  	rownames(mat) = names(target)
  	if(ncol(mat_target)) {
  		colnames(mat) = c(paste0("u", seq_along(upstream_index)), paste0("t", seq_along(target_index)), paste0("d", seq_along(downstream_index)))
  	} else {
  		colnames(mat) = c(paste0("u", seq_along(upstream_index)), paste0("d", seq_along(downstream_index)))
  	}

  	if(trim > 0) {
  		q1 = quantile(mat, trim/2, na.rm = TRUE)
  		q2 = quantile(mat, 1 - trim/2, na.rm = TRUE)
  		mat[mat <= q1] = q1
  		mat[mat >= q2] = q2
  	}
	class(mat) = c("normalizeToMatrix", "matrix")
	return(mat)
}

# 
# -gr input regions
# -target the upstream part or body part
# -window absolute size (100) or relative size(0.1)
# -value_column
# -mean_mode how to calculate mean value in a window
#
# == example
# gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 4, 7), end = c(2, 5, 8)))
# target = GRanges(seqnames = "chr1", ranges =IRanges(start = 1, end = 10))
# makeMatrix(gr, target, w = 2)
#
makeMatrix = function(gr, target, w = NULL, k = NULL, value_column = NULL, mapping_column = mapping_column, empty_value = 0,
    mean_mode = c("absolute", "weighted", "w0"), direction = c("normal", "reverse")) {
  
	if(is.null(value_column)) {
		gr$..value = rep(1, length(gr))
		value_column = "..value"
	}
  
	# split `target` into small windows
	target_windows = makeWindows(target, w = w, k = k, direction = direction)
 	strand(target_windows) = "*"
 	strand(gr) = "*"
 	
	# overlap `gr` to `target_windows`
	mtch = findOverlaps(gr, target_windows)
	mtch = as.matrix(mtch)
	
	# add a `value` column in `target_window` which is the mean value for intersected gr
	# in `target_window`
	m_gr = gr[ mtch[, 1] ]
	m_target_windows = target_windows[ mtch[, 2] ]

	if(!is.null(mapping_column)) {
		
		mapping = mcols(m_gr)[[mapping_column]]
		if(is.numeric(mapping)) {
			l = mapping == m_target_windows$.row
		} else {
			if(is.null(names(target))) {
				stop("`mapping_column` in `gr` is mapped to the names of `target`, which means `target` should have names.")
			} else {
				l = mapping == names(target)[m_target_windows$.row]
			}
		}

		m_gr = m_gr[l]
		m_target_windows = m_target_windows[l]
		mtch = mtch[l , , drop = FALSE]
	}

	v = mcols(m_gr)[[value_column]]
  
	mean_mode = match.arg(mean_mode)[1]
  
	if(mean_mode == "w0") {
		mintersect = pintersect(m_gr, m_target_windows)
		p = width(mintersect)/width(m_target_windows)
		x = tapply(p*v, mtch[, 2], sum)
	} else if(mean_mode == "absolute") {
		x = tapply(v, mtch[, 2], mean)
	} else {
		mintersect = pintersect(m_gr, m_target_windows)
		w = width(mintersect)
		x = tapply(w*v, mtch[, 2], sum) / tapply(w, mtch[, 2], sum)
	}
  
	v2 = rep(empty_value, length(target_windows))
	v2[ as.numeric(names(x)) ] = x
  
	target_windows$..value = v2

	# transform into a matrix
	tb = table(target_windows$.row)
	target_strand = strand(target)
	column_index = mapply(as.numeric(names(tb)), tb, FUN = function(i, n) {
		if(as.vector(target_strand[i] == "-")) {
			rev(seq_len(n))
		} else {
			seq_len(n)
		}
	})
  
	# is column_index has the same length for all regions in target?
	# if extension of upstream are the same or split body into k pieces,
	# then all column index has the same length
	# if it is not the same, throw error!
	if(!is.matrix(column_index)) {
		stop("numbers of columns are not the same.")
	}
  
	mat = matrix(empty_value, nrow = length(target), ncol = dim(column_index)[1])
	mat[ target_windows$.row + (as.vector(column_index) - 1)* nrow(mat) ] = target_windows$..value

	# findOverlaps may use a lot of memory
	rm(list = setdiff(ls(), "mat"))
	gc(verbose = FALSE)

	return(mat)
}

# == title
# Split regions into windows
#
# == param
# -gr a `GenomicRanges::GRanges` object.
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
#     ----------  a region, split by 3bp window
#     aaabbbccc   direction = "normal",  short.keep = FALSE
#     aaabbbcccd  direction = "normal",  short.keep = TRUE
#      aaabbbccc  direction = "reverse", short.keep = FALSE
#     abbbcccddd  direction = "reverse", short.keep = TRUE
#     
# There is an additional column ``.row`` attached which contains the correspondance between small windows
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
# makeWindows(gr, k = 2)
# makeWindows(gr, k = 3)
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
	ncol = unlist(lapply(pos[1, ], seq_along))  # which window from left to right
	chr = seqnames(gr)[orow]
	strand = strand(gr)[orow]

	gr = GRanges(seqnames = chr,
		         ranges = IRanges(start = start,
		         	              end = end),
		         strand = strand,
		         .row = orow,
		         .column = ncol)
	return(gr)

}

# == title
# Subset normalized matrix by rows
#
# == param
# -x the normalized matrix returned by `normalizeToMatrix`
# -i row index
# -j column index, disabled
# -drop disabled
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
"[.normalizeToMatrix" = function(x, i, j, drop = FALSE) {
	
	attr = attributes(x)
	attributes(x) = NULL
	for(bb in intersect(names(attr), c("dim", "dimnames"))) {
		attr(x, bb) = attr[[bb]]
	}
	if(!missing(i) && !missing(j)) {
		return(x[i, j, drop = FALSE])
	}
	if(nargs() == 2) {
		return(x[i])
	}
	if(nargs() == 3 && missing(i)) {
		return(x[, j])
	}
	if(missing(i)) {
		return(x[i, j, drop = drop])
	}
	x = x[i, , drop = FALSE]
	for(bb in setdiff(names(attr), c("dim", "dimnames"))) {
		attr(x, bb) = attr[[bb]]
	}
	return(x)
}

# == title
# Print normalized matrix
#
# == param
# -x the normalized matrix returned by `normalizeToMatrix`
# -... other arguments
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
print.normalizeToMatrix = function(x, ...) {
	upstream_index = attr(x, "upstream_index")
	target_index = attr(x, "target_index")
	downstream_index = attr(x, "downstream_index")
	extend = attr(x, "extend")
	smooth = attr(x, "smooth")
	signal_name = attr(x, "signal_name")
	target_name = attr(x, "target_name")

	cat("Normalize ", signal_name, " to ", target_name, ":\n", sep = "")
	cat("  Upstream ", extend[1], "bp\n", sep = "")
	cat("  Downstream ", extend[2], "bp\n", sep = "")
	if(length(target_index) == 0) {
		cat("  Not show target regions\n", sep = "")
	} else {
		cat("  Show target regions.\n", sep = "")
	}
	cat("  ", nrow(x), " signal regions\n", sep = "")
}

# == title
# Copy attributes of a normalized matrix to another
#
# == param
# -x x
# -y y
#
copyAttr = function(x, y) {
	if(!identical(dim(x), dim(y))) {
		stop("x and y should have same dimension.\n")
	}
	attr = attributes(x)
	for(bb in names(attr)) {
		attr(y, bb) = attr[[bb]]
	}
	attr(y, "signal_name") = "\b"
	return(y)
}

