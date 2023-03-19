

# == title
# Normalize Associations between Genomic Signals and Target Regions into a Matrix
#
# == param
# -signal A `GenomicRanges::GRanges-class` object.
# -target A `GenomicRanges::GRanges-class` object.
# -extend Extended base pairs to the upstream and/or downstream of ``target``. It can be a vector of length one or two.
#         Length one means same extension to the upstream and downstream.
# -w Window size for splitting upstream and downstream, measured in base pairs
# -value_column Column index in ``signal`` that is mapped to colors. If it is not set, it assumes values for all signal regions are 1.
# -mapping_column Mapping column to restrict overlapping between ``signal`` and ``target``. By default it tries to look for
#           all regions in ``signal`` that overlap with every target.
# -background Values for windows that don't overlap with ``signal``. 
# -empty_value Deprecated, please use ``background`` instead.
# -mean_mode When a window is not perfectly overlapped to ``signal``, how to summarize 
#        values to the window. See 'Details' section for a detailed explanation.
# -include_target  Whether include ``target`` in the heatmap? If the width of all regions in ``target`` is 1, ``include_target``
#               is enforced to ``FALSE``.
# -target_ratio The ratio of ``target`` columns in the normalized matrix. If the value is 1, ``extend`` will be reset to 0.
# -k Number of windows only when ``target_ratio = 1`` or ``extend == 0``, otherwise ignored.
# -smooth Whether apply smoothing on rows in the matrix?
# -smooth_fun The smoothing function that is applied to each row in the matrix. This self-defined function accepts a numeric
#    vector (may contain ``NA`` values) and returns a vector with same length. If the smoothing is failed, the function
#    should call `base::stop` to throw errors so that `normalizeToMatrix` can catch how many rows are failed in smoothing. 
#    See the default `default_smooth_fun` for example.
# -keep Percentiles in the normalized matrix to keep. The value is a vector of two percent values. Values less than the first
#       percentile is replaces with the first pencentile and values larger than the second percentile is replaced with the
#       second percentile.
# -limit Similar as ``keep``, but it provides boundary for absolute values. The value should be a vector of length two.
# -trim Deprecated, please use ``keep`` instead.
# -flip_upstream Sometimes whether the signals are on the upstream or the downstream of the targets
#      are not important and users only want to show the relative distance to targets. If the value is set
#      to ``TRUE``, the upstream part in the normalized matrix is flipped and added to the downstream part
#      The flipping is only allowed when the targets are single-point targets or the targets are excluded
#      in the normalized matrix (by setting ``include_target = FALSE``). If the extension for the upstream
#      and downstream is not equal, the smaller extension is used for the final matrix.
# -verbose Whether to print help messages.
#
# == details
# In order to visualize associations between ``signal`` and ``target``, the data is transformed into a matrix
# and visualized as a heatmap by `EnrichedHeatmap` afterwards.
#
# Upstream and downstream also with the target body are splitted into a list of small windows and overlap
# to ``signal``. Since regions in ``signal`` and small windows do not always 100 percent overlap, there are four different averaging modes:
# 
# Following illustrates different settings for ``mean_mode`` (note there is one signal region overlapping with other signals):
#
#       40      50     20     values in signal regions
#     ++++++   +++    +++++   signal regions
#            30               values in signal region
#          ++++++             signal region
#       =================     a window (17bp), there are 4bp not overlapping to any signal regions.
#         4  6  3      3      overlap
#
#     absolute: (40 + 30 + 50 + 20)/4
#     weighted: (40*4 + 30*6 + 50*3 + 20*3)/(4 + 6 + 3 + 3)
#     w0:       (40*4 + 30*6 + 50*3 + 20*3)/(4 + 6 + 3 + 3 + 4)
#     coverage: (40*4 + 30*6 + 50*3 + 20*3)/17
#
# == value
# A matrix with following additional attributes:
#
# -``upstream_index`` column index corresponding to upstream of ``target``
# -``target_index`` column index corresponding to ``target``
# -``downstream_index`` column index corresponding to downstream of ``target``
# -``extend`` extension on upstream and downstream
# -``smooth`` whether smoothing was applied on the matrix
# -``failed_rows`` index of rows which are failed after smoothing
#
# The matrix is wrapped into a simple ``normalizeToMatrix`` class.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# signal = GRanges(seqnames = "chr1", 
#     ranges = IRanges(start = c(1, 4, 7, 11, 14, 17, 21, 24, 27),
#                      end = c(2, 5, 8, 12, 15, 18, 22, 25, 28)),
#     score = c(1, 2, 3, 1, 2, 3, 1, 2, 3))
# target = GRanges(seqnames = "chr1", ranges = IRanges(start = 10, end = 20))
# normalizeToMatrix(signal, target, extend = 10, w = 2)
# normalizeToMatrix(signal, target, extend = 10, w = 2, include_target = TRUE)
# normalizeToMatrix(signal, target, extend = 10, w = 2, value_column = "score")
#
normalizeToMatrix = function(signal, target, extend = 5000, w = max(extend)/100, 
	value_column = NULL, mapping_column = NULL, background = ifelse(smooth, NA, 0), 
	empty_value = NULL, mean_mode = c("absolute", "weighted", "w0", "coverage"), 
	include_target = any(width(target) > 1), 
	target_ratio = min(c(0.4, mean(width(target))/(sum(extend) + mean(width(target))))), 
	k = min(c(20, min(width(target)))), smooth = FALSE, smooth_fun = default_smooth_fun,
	keep = c(0, 1), limit = NULL, trim = NULL, flip_upstream = FALSE, verbose = TRUE) {

	signal_name = deparse(substitute(signal))
	target_name = deparse(substitute(target))
	
	if(length(extend) == 1) extend = c(extend, extend)

	## if value is categorical data
	if(!is.null(value_column)) {
		x = mcols(signal)[, value_column]
		if(is.factor(x) || is.character(x)) {

			if(flip_upstream) {
				stop_wrap("The categorical signal is not allowed to flip the upstream signals.")
			}

			mat = normalize_with_category_variable(signal, target, x, extend = extend, w = w, mapping_column = mapping_column, background = 0,
				empty_value = 0, target_ratio = target_ratio, k = k, smooth = smooth, smooth_fun = smooth_fun, keep = keep, limit = limit, trim = trim)
			attributes(mat)$signal_name = signal_name
			attributes(mat)$target_name = target_name
			return(mat)
		}
	}

	if(abs(target_ratio - 1) < 1e-6 || abs(target_ratio) >= 1) {
		if(!all(extend == 0)) warning_wrap("Rest `extend` to 0 when `target_ratio` is larger than or euqal to 1.")
		extend = c(0, 0)
	} else if(all(extend == 0)) {
		warning_wrap("Reset `target_ratio` to 1 when `extend` is 0.")
		target_ratio = 1
	}
	if(abs(target_ratio) > 1) target_ratio = 1

	target_is_single_point = all(width(target) <= 1)

	if(target_is_single_point) {
		if(include_target) {
			warning_wrap("Width of `target` are all 1, `include_target` is set to `FALSE`.")
		}
		include_target = FALSE
	}
  
	if(extend[1] > 0) {
		if(extend[1] %% w > 0) {
			warning_wrap("Length of upstream extension is not completely divisible by `w`.")
			extend[1] = extend[1] - extend[1] %% w
		}
	}
	if(extend[2] > 0) {
		if(extend[2] %% w > 0) {
			warning_wrap("Length of downstream extension is not completely divisible by `w`.")
			extend[2] = extend[2] - extend[2] %% w
		}
	}

	if(!missing(empty_value) && missing(background)) {
		if(!is.null(empty_value)) {
			background = empty_value
		}
	}

	if(!missing(trim) && missing(keep)) {
		if(!is.null(trim)) {
			if(length(trim) == 1) trim = rep(trim, 2)
			keep = c(trim[1], 1 - trim[2])
		}
	}

	.seq = function(start, end, by = 1) {
		if(end < start) {
			return(integer(0))
		} else {
			seq(start, end, by = by)
		}
	}
  	
  	if(target_is_single_point) {
  		# do not need to separate upstream and downstream
  		# and it makes the boundary between upstream and downstream smoothing
  		suppressWarnings(both <- promoters(target, upstream = extend[1], downstream = extend[2] + 1))

		mat_both = makeMatrix(signal, both, w = w, value_column = value_column, mapping_column = mapping_column, 
			background = background, mean_mode = mean_mode)
		i = round(extend[1]/(extend[1] + extend[2]) * ncol(mat_both))  # assume
		# if(i < 2 | ncol(mat_both) - i < 2) {
		# 	stop("Maybe `w` is too large or one of `extend` is too small.")
		# }
		mat_upstream = mat_both[, .seq(1, i), drop = FALSE]
		mat_downstream = mat_both[, .seq(i+1, ncol(mat_both)), drop = FALSE]

  	} else {
		# extend and normalize in upstream 
		if(extend[1] <= 0) {
			mat_upstream = matrix(0, nrow = length(target), ncol = 0)
		} else {
			suppressWarnings(upstream <- promoters(target, upstream = extend[1], downstream = 0))
			mat_upstream = makeMatrix(signal, upstream, w = w, value_column = value_column, mapping_column = mapping_column, 
				background = background, mean_mode = mean_mode)
		}

		# extend and normalize in downstream
		e = ifelse(strand(target) == "-", start(target) - 1, end(target) + 1)
		end_target = GRanges(seqnames = seqnames(target),
	                         ranges = IRanges(start = e, end = e),
	                         strand = strand(target))
		if(extend[2] <= 0) {
			mat_downstream = matrix(0, nrow = length(target), ncol = 0)
		} else {
			suppressWarnings(downstream <- promoters(end_target, upstream = 0, downstream = extend[2] ))
			names(downstream) = names(target)
		  
			mat_downstream = makeMatrix(signal, downstream, w = w, value_column = value_column, mapping_column = mapping_column, 
				background = background, mean_mode = mean_mode)
		}
	}

	if(include_target) {
		if(!all(extend == 0)) {
			k = round((ncol(mat_upstream) + ncol(mat_downstream)) * target_ratio/(1-target_ratio))
			if(k < 1) k = 1
		} 
		mat_target = makeMatrix(signal, target, k = k, value_column = value_column, mapping_column = mapping_column, background = background, mean_mode = mean_mode)
	} else {
		mat_target = matrix(0, nrow = length(target), ncol = 0)
	}

  	mat = cbind(mat_upstream, mat_target, mat_downstream)

  	# guess whether the signal is methylation
  	if(is.null(limit) && smooth) {
  		xx = mat[!is.na(mat)]
  		if(all(xx >= 0 & xx <= 1)) {
  			if(verbose) {
  				message_wrap("All signal values are within [0, 1], so we assume it is methylation signal. Automatically set limit [0, 1] to the smoothed values. If this is not the case, set argument `limit = NA` in the function to remove the limits. Set `verbose = FALSE` to turn off this message.")
  			}
  			limit = c(0, 1)
  		}
  	}
  	if(identical(limit, NA)) {
  		limit = NULL
  	}

  	# apply smoothing on rows in mat
  	failed_rows = NULL
  	all_positive = all(mat >= 0, na.rm = TRUE)
  	
	if(smooth) {
		i_row = 0
		ow = options("warn")[[1]]
		mat = t(apply(mat, 1, function(x) {
			i_row <<- i_row + 1

			oe = try(x <- suppressWarnings(smooth_fun(x)), silent = TRUE)
			if(inherits(oe, "try-error")) {
				failed_rows <<- c(failed_rows, i_row)
			}
			if(all_positive) x[x < 0] = 0
			return(x)
		}))
		options(warn = ow)
		
		if(!is.null(failed_rows)) {
			if(length(failed_rows) == 1) {
				msg = paste(strwrap(paste0("Smoothing is failed for one row because there are very few signals overlapped to it. Please use `failed_rows(mat)` to get the index of the failed row and consider to remove it.\n")), collapse = "\n")
			} else {
				msg = paste(strwrap(paste0("Smoothing are failed for ", length(failed_rows), " rows because there are very few signals overlapped to them. Please use `failed_rows(mat)` to get the index of failed rows and consider to remove them.\n")), collapse = "\n")
			}
			msg = paste0("\n", msg, "\n")
			warning_wrap(msg)
		}
	}
	
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
	attr(mat, "target_is_single_point") = target_is_single_point
	attr(mat, "failed_rows") = failed_rows
	attr(mat, "background") = background
	attr(mat, "signal_is_categorical") = FALSE

	.paste0 = function(a, b) {
		if(length(a) == 0 || length(b) == 0) {
			return(NULL)
		} else {
			paste0(a, b)
		}
	}

	# dimension names are mainly for debugging
  	rownames(mat) = names(target)
  	if(ncol(mat_target)) {
  		colnames(mat) = c(.paste0("u", seq_along(upstream_index)), .paste0("t", seq_along(target_index)), .paste0("d", seq_along(downstream_index)))
  	} else {
  		colnames(mat) = c(.paste0("u", seq_along(upstream_index)), .paste0("d", seq_along(downstream_index)))
  	}
	q1 = quantile(mat, keep[1], na.rm = TRUE)
	q2 = quantile(mat, keep[2], na.rm = TRUE)
	mat[mat <= q1] = q1
	mat[mat >= q2] = q2
	if(!is.null(limit)) {
		mat[mat < limit[1]] = limit[1]
		mat[mat > limit[2]] = limit[2]
	}
  	
  	class(mat) = c("normalizedMatrix", "matrix")

  	if(flip_upstream) {
  		mat = flip_upstream(mat)
  	}

	return(mat)
}

normalize_with_category_variable = function(signal, target, category, ...) {
	if(!is.factor(category)) {
		category = factor(category)
	}

	level = levels(category)
	n_level = length(level)
	mat_list = lapply(as.list(split(signal, category)), function(gr) {
		normalizeToMatrix(gr, target, value_column = NULL, mean_mode = "coverage", ...)
	})[level]

	arr = array(, dim = c(dim(mat_list[[1]]), n_level))
	for(i in 1:n_level) {
		arr[, , i] = mat_list[[i]]
	}

	mat = apply(arr, c(1, 2), function(x) {
		if(all(x == 0)) {
			return(0)
		} else {
			which.max(x)
		}
	})
	mat = copyAttr(mat_list[[1]], mat)
	attributes(mat)$signal_is_categorical = TRUE
	attributes(mat)$signal_level = level

	return(mat)
}

flip_upstream = function(mat) {
	if(attributes(mat)$signal_is_categorical) {
		stop_wrap("A normalized matrix with categorical signals is not allowed to flip the upstream signals.")
	}

	upstream_index = attr(mat, "upstream_index")
	target_index = attr(mat, "target_index")
	downstream_index = attr(mat, "downstream_index")
	extend = attr(mat, "extend")

	if(length(target_index) != 0) {
		stop_wrap("Upstream can only be flipped when the targets are single points or excluded.")
	}

	k = min(length(upstream_index), length(downstream_index))
	m = mat[, upstream_index[seq(length(upstream_index), length(upstream_index) - k + 1)]] + 
	    mat[, downstream_index[1:k]]
	attr(m, "upstream_index") = integer(0)
	attr(m, "downstream_index") = 1:k
	attr(m, "extend") = c(0, min(extend))
	attr(m, "target_is_single_point") = attr(mat, "target_is_single_point")
	attr(m, "target_index") = integer(0)
	attr(m, "smooth") = attr(mat, "smooth")
	attr(m, "signal_name") = attr(mat, "signal_name")
	attr(m, "target_name") = attr(mat, "target_name")
	attr(m, "failed_rows") = attr(mat, "failed_rows")
	attr(m, "background") = attr(mat, "background")
	colnames(m) = paste0("d", 1:k)
	attr(m, "upstream_flipped") = TRUE

	class(m) = c("normalizedMatrix", "matrix")
	return(m)
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
makeMatrix = function(gr, target, w = NULL, k = NULL, value_column = NULL, mapping_column = NULL, background = 0,
    mean_mode = c("absolute", "weighted", "w0", "coverage"), direction = c("normal", "reverse")) {
  
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

	# subset `m_gr` and `m_target_windows` based on `mapping_column`
	if(!is.null(mapping_column)) {
		
		mapping = mcols(m_gr)[[mapping_column]]
		if(is.numeric(mapping)) {
			l = mapping == m_target_windows$.i_query
		} else {
			if(is.null(names(target))) {
				stop("`mapping_column` in `gr` is mapped to the names of `target`, which means `target` should have names.")
			} else {
				l = mapping == names(target)[m_target_windows$.i_query]
			}
		}

		m_gr = m_gr[l]
		m_target_windows = m_target_windows[l]
		mtch = mtch[l , , drop = FALSE]
	}

	# the value associated with `gr`
	v = mcols(m_gr)[[value_column]]

	mean_mode = match.arg(mean_mode)[1]

	if(length(mtch)) {
		if(mean_mode == "w0") {
			mintersect = pintersect(m_gr, m_target_windows)
			w = width(mintersect)
			target_windows_list = split(ranges(m_gr), mtch[, 2])
			target_windows2 = target_windows[as.numeric(names(target_windows_list))]
			cov = coverage(target_windows_list, shift = -start(target_windows2), width = width(target_windows2))
			#non_intersect_width = sapply(cov, function(x) sum(x == 0))
			non_intersect_width = sapply(cov@listData, function(x) {ind = x@values == 0;sum(x@lengths[ind])})
			x = tapply(w*v, mtch[, 2], sum, na.rm = TRUE) / (tapply(w, mtch[, 2], sum, na.rm = TRUE) + non_intersect_width)
		} else if(mean_mode == "coverage") {
			mintersect = pintersect(m_gr, m_target_windows)
			p = width(mintersect)/width(m_target_windows)
			x = tapply(p*v, mtch[, 2], sum, na.rm = TRUE)
		} else if(mean_mode == "absolute") {
			x = tapply(v, mtch[, 2], mean, na.rm = TRUE)
		} else {
			mintersect = pintersect(m_gr, m_target_windows)
			w = width(mintersect)
			x = tapply(w*v, mtch[, 2], sum, na.rm = TRUE) / tapply(w, mtch[, 2], sum, na.rm = TRUE)
		}
		v2 = rep(background, length(target_windows))
		v2[ as.numeric(names(x)) ] = x
	} else {
		v2 = rep(background, length(target_windows))
	}

	target_windows$..value = v2

	# transform into a matrix
	tb = table(target_windows$.i_query)
	target_strand = strand(target)
	column_index = mapply(as.numeric(names(tb)), tb, FUN = function(i, n) {
		if(as.vector(target_strand[i] == "-")) {
			rev(seq_len(n))
		} else {
			seq_len(n)
		}
	}, SIMPLIFY = FALSE)
	column_index = do.call("cbind", column_index)

	# is column_index has the same length for all regions in target?
	# if extension of upstream are the same or split body into k pieces,
	# then all column index has the same length
	# if it is not the same, throw error!
	if(!is.matrix(column_index)) {
		stop("numbers of columns are not the same.")
	}
  
	mat = matrix(background, nrow = length(target), ncol = dim(column_index)[1])
	mat[ target_windows$.i_query + (as.vector(column_index) - 1)* nrow(mat) ] = target_windows$..value

	# findOverlaps may use a lot of memory
	rm(list = setdiff(ls(), "mat"))
	gc(verbose = FALSE)

	return(mat)
}

# == title
# Split Regions into Windows
#
# == param
# -query A `GenomicRanges::GRanges-class` object.
# -w Window size. A value larger than 1 means the number of base pairs and a value between 0 and 1
#    is the percent to the current region.
# -k Number of partitions for each region. If it is set, all other arguments are ignored.
# -direction Where to start the splitting? See 'Details' section.
# -short.keep If the the region can not be split equally under the window size, 
#             the argument controls whether to keep the windows that are smaller than the window size. See 'Details' section.
#
# == details
# Following illustrates the meaning of ``direction`` and ``short.keep``:
#
#     -->-->-->-  one region, split by 3bp window (">" represents the direction of the sequence)
#     aaabbbccc   direction = "normal",  short.keep = FALSE
#     aaabbbcccd  direction = "normal",  short.keep = TRUE
#      aaabbbccc  direction = "reverse", short.keep = FALSE
#     abbbcccddd  direction = "reverse", short.keep = TRUE
#     
#
# == value
# A `GenomicRanges::GRanges-class` object with two additional columns attached:
# 
# - ``.i_query`` which contains the correspondance between small windows and original regions in ``query``
# - ``.i_window`` which contains the index of the small window on the current region.
#
# == author
# Zuguang gu <z.gu@dkfz.de>
#
# == example
# query = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 21), end = c(10, 20, 30)))
# makeWindows(query, w = 2)
# makeWindows(query, w = 0.5)
# makeWindows(query, w = 3)
# makeWindows(query, w = 3, direction = "reverse")
# makeWindows(query, w = 3, short.keep = TRUE)
# makeWindows(query, w = 3, direction = "reverse", short.keep = TRUE)
# makeWindows(query, w = 12)
# makeWindows(query, w = 12, short.keep = TRUE)
# makeWindows(query, k = 2)
# makeWindows(query, k = 3)
# query = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 31), end = c(10, 30, 70)))
# makeWindows(query, w = 2)
# makeWindows(query, w = 0.2)
#
makeWindows = function(query, w = NULL, k = NULL, direction = c("normal", "reverse"), 
	short.keep = FALSE) {

	direction = match.arg(direction)[1]
  
	if(is.null(w) & is.null(k)) {
		stop("You should define either `w` or `k`.")
	}
	ostart = start(query)
	oend = end(query)
  
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
	i_query = rep(seq_len(ncol(pos)), times = sapply(pos[1, ], length))
	i_window = unlist(lapply(pos[1, ], seq_along))  # which window from left to right
	chr = seqnames(query)[i_query]
	strand = strand(query)[i_query]

	gr = GRanges(seqnames = chr,
		         ranges = IRanges(start = start,
		         	              end = end),
		         strand = strand,
		         .i_query = i_query,
		         .i_window = i_window)
	return(gr)

}

# == title
# Subset normalized matrix by rows
#
# == param
# -x the normalized matrix returned by `normalizeToMatrix`
# -i row index
# -j column index
# -drop whether drop the dimension
#
# == value
# A ``normalizedMatrix`` class object.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
"[.normalizedMatrix" = function(x, i, j, drop = FALSE) {
	
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
# Bind Matrix by Rows
#
# == param
# -... Matrices
# -deparse.level Not used.
#
# == value
# A ``normalizedMatrix`` class object.
#
# == author
# z.gu@dkfz.de
#
rbind.normalizedMatrix = function(..., deparse.level = 1) {
	mat_list = list(...)
	mat_list = mat_list[sapply(mat_list, function(x) !is.null(x))]
	mat_list2 = lapply(mat_list, function(x) {
		attributes(x) = attributes(x)["dim"]
		x
	})
	rbind_matrix = selectMethod("rbind", signature = "matrix")
	mat = do.call("rbind_matrix", mat_list2)
	mat = copyAttr(mat_list[[1]], mat)
	return(mat)
}

# == title
# Print the Normalized Matrix
#
# == param
# -x The normalized matrix returned by `normalizeToMatrix`.
# -... Other arguments.
#
# == value
# No value is returned.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
print.normalizedMatrix = function(x, ...) {

	upstream_index = attr(x, "upstream_index")
	target_index = attr(x, "target_index")
	downstream_index = attr(x, "downstream_index")
	extend = attr(x, "extend")
	smooth = attr(x, "smooth")
	signal_name = attr(x, "signal_name")
	target_name = attr(x, "target_name")
	target_is_single_point = attr(x, "target_is_single_point")
	if(is.null(target_is_single_point)) {
		target_is_single_point = FALSE
	}
	upstream_flipped = attr(x, "upstream_flipped")
	if(is.null(upstream_flipped)) upstream_flipped = FALSE

	op = qq.options(READ.ONLY = FALSE)
    on.exit(qq.options(op))
    qq.options(code.pattern = "@\\{CODE\\}")

	qqcat("Normalize @{signal_name} to @{target_name}:\n")
	if(upstream_flipped) {
		qqcat("  Extension @{extend[2]} bp (@{length(downstream_index)} window@{ifelse(length(upstream_index) > 1, 's', '')})\n")
		qqcat("    upstream is flipped to downstream.\n")
	} else {
		qqcat("  Upstream @{extend[1]} bp (@{length(upstream_index)} window@{ifelse(length(upstream_index) > 1, 's', '')})\n")
		qqcat("  Downstream @{extend[2]} bp (@{length(downstream_index)} window@{ifelse(length(upstream_index) > 1, 's', '')})\n")	
	}
	if(target_is_single_point) {
			qqcat("  Include target regions (width = 1)\n")
	} else {
		if(length(target_index) == 0) {
			qqcat("  Not include target regions\n")
		} else {
			qqcat("  Include target regions (@{length(target_index)} window@{ifelse(length(target_index) > 1, 's', '')})\n")
		}
	}
	qqcat("  @{nrow(x)} target region@{ifelse(nrow(x) > 1, 's', '')}\n")
	signal_is_categorical = attr(x, "signal_is_categorical")
	if(is.null(signal_is_categorical)) signal_is_categorical = FALSE
	if(signal_is_categorical) {
		qqcat("  signal is categorical (@{length(attr(x, 'signal_level'))} levels)\n")
	}
}

# == title
# Copy Attributes to Another Object
#
# == param
# -x Object 1.
# -y Object 2.
#
# == details
# The `normalizeToMatrix` object is actually a matrix but with more additional attributes attached.
# When manipulating such matrix, there are some circumstances that the attributes are lost.
# This function is used to copy these specific attributes when dealing with the matrix.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# gr = GRanges(seqnames = c("chr5", "chr5"),
# 	ranges = IRanges(start = c(98, 98),
# 	                 end = c(104, 104)))
# target = GRanges(seqnames = "chr5",
# 	ranges = IRanges(start = 100, 
# 		             end = 100))
# mat1 = normalizeToMatrix(gr, target, extend = 6, w = 1)
# # attributes removed and you cannot use it for EnrichedHeatmap()
# mat2 = mat1[]
# # copy attributes to mat2 and now mat3 can be used for EnrichedHeatmap()
# mat3 = copyAttr(mat1, mat2)
copyAttr = function(x, y) {
	if(!identical(ncol(x), ncol(y))) {
		stop("x and y should have same number of columns.\n")
	}
	attr = attributes(x)
	for(bb in setdiff(names(attr), c("dim"))) {
		if(bb == "dimnames") {
			attr(y, bb)[[2]] = attr[[bb]][[2]]  # set same column names
		} else {
			attr(y, bb) = attr[[bb]]
		}
	}
	attr(y, "signal_name") = "\b"
	return(y)
}

# == title
# Get Signals from a List
#
# == param
# -lt A list of normalized matrices which are returned by `normalizeToMatrix`. Matrices in the list should be generated with same settings (e.g. they
#     should use same target regions, same extension to targets and same number of windows).
# -fun A user-defined function to summarize signals.
#
# == details
# Let's assume you have a list of histone modification signals for different samples and you want
# to visualize the mean pattern across samples. You can first normalize histone mark signals for each sample and then
# calculate means values across all samples. In following example code, ``hm_gr_list`` is a list of ``GRanges`` objects
# which contain positions of histone modifications, ``tss`` is a ``GRanges`` object containing positions of gene TSS.
#
#     mat_list = NULL
#     for(i in seq_along(hm_gr_list)) {
#         mat_list[[i]] = normalizeToMatrix(hm_gr_list[[i]], tss, value_column = "density")
#     }
#
# If we compress the list of matrices as a three-dimension array where the first dimension corresponds to genes,
# the second dimension corresponds to windows and the third dimension corresponds to samples, the mean signal
# across all sample can be calculated on the third dimension. Here `getSignalsFromList` simplifies this job.
#
# Applying ``getSignalsFromList()`` to ``mat_list``, it gives a new normalized matrix which contains mean signals across all samples and can
# be directly used in ``EnrichedHeatmap()``.
#
#     mat_mean = getSignalsFromList(mat_list)
#     EnrichedHeatmap(mat_mean)
#
# The correlation between histone modification and gene expression can
# also be calculated on the third dimension of the array. In the user-defined function ``fun``, ``x`` is the vector for gene i
# and window j in the array, and ``i`` is the index of current gene.
# 
#     mat_corr = getSignalsFromList(mat_list, 
#         fun = function(x, i) cor(x, expr[i, ], method = "spearman"))
#
# Then ``mat_corr`` here can be used to visualize how gene expression is correlated to histone modification around TSS.
#
#     EnrichedHeatmap(mat_corr)
#
#
# == value
# A `normalizeToMatrix` object which can be directly used for `EnrichedHeatmap`.
# 
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# NULL
getSignalsFromList = function(lt, fun = function(x) mean(x, na.rm = TRUE)) {

	if(!inherits(lt, "list")) {
		stop_wrap("`lt` should be a list of objects which are returned by `normalizeToMatrix()`.")
	}

	if(!all(sapply(lt, inherits, "normalizedMatrix"))) {
		stop_wrap("`lt` should be a list of objects which are returned by `normalizeToMatrix()`.")
	}

	n = length(lt)
	if(n > 1) {
		for(i in seq_len(n-1)) {
			attr1 = attributes(lt[[i]])[ c("upstream_index", "target_index", "downstream_index", "extend") ]
			attr2 = attributes(lt[[i+1]])[ c("upstream_index", "target_index", "downstream_index", "extend") ]
			if(!identical(attr1, attr2)) {
				stop_wrap("Objects in `lt` should have same settings.")
			}
		}
	}

	flag = 0
	for(i in seq_len(n)) {
		tm = lt[[i]]
	    if(!flag) {
	    	arr = array(dim = c(dim(tm), length(lt)))
	    	flag = 1
	    }
	    arr[, , i] = tm
	}

	if(identical(fun, mean)) {
		fun = function(x) mean(x, na.rm = TRUE)
	} else if(identical(fun, median)) {
		fun = function(x) median(x, na.rm = TRUE)
	} else if(identical(fun, max)) {
		fun = function(x) max(x, na.rm = TRUE)
	} else if(identical(fun, min)) {
		fun = function(x) min(x, na.rm = TRUE)
	} 

	if(length(as.list(fun)) == 2) {
		m = apply(arr[, , ,drop = FALSE], c(1, 2), fun)
	} else if(length(as.list(fun)) == 3) {
		m = matrix(nrow = nrow(lt[[1]]), ncol = ncol(lt[[1]]))
		for(i in seq_len(nrow(m))) {
			for(j in seq_len(ncol(m))) {
				m[i, j] = fun(arr[i, j, ], i)
			}
		}
	} else {
		stop_wrap("`fun` can only have one or two arguments.")
	}
	m = copyAttr(lt[[1]], m)
	return(m)
}

# == title
# Default Smoothing function
#
# == param
# -x Input numeric vector.
#
# == details
# The smoothing function is applied to every row in the normalized matrix. For this default smoothing function,
# `locfit::locfit` is first tried on the vector. If there is error, `stats::loess` smoothing is tried afterwards.
# If both smoothing are failed, there will be an error.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
default_smooth_fun = function(x) {
	l = !is.na(x)
	if(sum(l) >= 2) {
		oe1 = try(x <- suppressWarnings(predict(locfit(x[l] ~ lp(seq_along(x)[l], nn = 0.1, h = 0.8)), seq_along(x))), silent = TRUE)
		if(inherits(oe1, "try-error")) {
			oe2 = try(x <-  suppressWarnings(predict(loess(x[l] ~ seq_along(x)[l], control = loess.control(surface = "direct")), seq_along(x))))

			if(inherits(oe2, "try-error")) {
				stop_wrap("error when doing locfit or loess smoothing")
			} else {
				return(x)
			}
		} else {
			return(x)
		}
	} else {
		stop_wrap("Too few data points.")
	}
	return(x)
}


# == title
# Discretize a Continuous Matrix to a Discrete Matrix
#
# == param
# -mat A normalize matrix from `normalizeToMatrix`.
# -rule A list of intervals which provide mapping between continuous values to discrete values.
#       Note the order of intervals determines the order of corresponding discrete levels.
# -right_closed Is the interval right closed?
#
# == details
# Assuming we have a normalized matrix with both positive values and negative values, we only 
# want to see the enrichment of the windows/regions showing significant positive values and 
# negative values and we are only interested in the direction of the values while not the value itself,
# then we can define the rule as:
#
#     rule = list(
#         "positive" = c(0.5, Inf),
#         "negative" = c(-Inf, -0.5)
#     )
#
# And we can convert the continuous matrix to a discrete matrix and visualize it:
#
#     mat2 = discretize(mat, rule)
#     EnrichedHeatmap(mat2, col = c("positive" = "red", "negative" = "green"))
#
# Another example is to discretize the signals to discrete levels according to the intensities:
#     
#     rule = list(
#         "very_high" = c(100, Inf),
#         "high" = c(50, 100),
#         "intermediate" = c(25, 50),
#         "low" = c(1e-6, 25)
#     )
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
discretize = function(mat, rule, right_closed = FALSE) {
	if(!inherits(mat, "normalizedMatrix")) {
		stop_wrap("`mat` should be generated by `normalizeToMatrix().")
	}

	if(!inherits(rule, "list")) {
		stop_wrap("`rule` should be defined as a list of intervals.")
	}

	if(!all(sapply(rule, length) == 2)) {
		stop_wrap("All intervals in `rule` should have length of 2.")
	}

	mat2 = matrix(0, nrow = nrow(mat), ncol = ncol(mat))
	for(i in seq_along(rule)) {
		interval = rule[[i]]
		if(right_closed) {
			l = mat > interval[1] & mat <= interval[2]
		} else {
			l = mat >= interval[1] & mat < interval[2]
		}
		mat2[l] = i
	}

	mat2 = copyAttr(mat, mat2)
	attributes(mat2)$signal_is_categorical = TRUE
	attributes(mat2)$signal_level = names(rule)

	return(mat2)
}

# == title
# Convert a Normal Matrix to a normalizedMatrix Object
#
# == param
# -mat A matrix generated by other software.
# -k_upstream Number of windows in the upstream.
# -k_downstream Number of windows in the downstream.
# -k_target Number of windows in the target.
# -extend Extension to the target. The length should be 1 (if one of ``k_upstream`` or ``k_downstream`` is zero).
#  or 2 (if both of ``k_upstream`` and ``k_downstream`` are non-zero).
# -signal_name The name of signal regions. It is only used for printing the object.
# -target_name The name of the target names. It is only used for printing the object.
# -background The background value in the matrix.
# -smooth Whether apply smoothing on rows in the matrix. 
# -smooth_fun The smoothing function that is applied to each row in the matrix. This self-defined function accepts a numeric
#    vector (may contain ``NA`` values) and returns a vector with same length. If the smoothing is failed, the function
#    should call `base::stop` to throw errors so that `normalizeToMatrix` can catch how many rows are failed in smoothing. 
#    See the default `default_smooth_fun` for example.
# -keep Percentiles in the normalized matrix to keep. The value is a vector of two percent values. Values less than the first
#       percentile is replaces with the first pencentile and values larger than the second percentile is replaced with the
#       second percentile.
# -trim Deprecated, please use ``keep`` instead.
#
# == details
# If users use the matrix from other software, they can use this function to convert it to the ``normalizedMatrix`` object
# and visualize it afterwards.
#
# == value
# A ``normalizedMatrix`` object.
#
# == author
# z.gu@dkfz.de
#
as.normalizedMatrix = function(mat, k_upstream = 0, k_downstream = 0, k_target = 0,
	extend, signal_name = "signals", target_name = "targets",
	background = NA, smooth = FALSE, smooth_fun = default_smooth_fun, 
	keep = c(0, 1), trim = NULL) {

	if(k_upstream + k_target + k_downstream != ncol(mat)) {
		stop_wrap("sum of `k_upstream`, `k_target` and `k_downstream` should be equal to the col of `mat`.")
	}

	# apply smoothing on rows in mat
  	failed_rows = NULL
  	
	if(smooth) {
		mat[mat == background] = NA
		i_row = 0
		ow = options("warn")[[1]]
		mat = t(apply(mat, 1, function(x) {
			i_row <<- i_row + 1

			oe = try(x <- suppressWarnings(smooth_fun(x)), silent = TRUE)
			if(inherits(oe, "try-error")) {
				failed_rows <<- c(failed_rows, i_row)
			}
			return(x)
		}))
		options(warn = ow)
		
		if(!is.null(failed_rows)) {
			if(length(failed_rows) == 1) {
				msg = paste(strwrap(paste0("Smoothing is failed for one row because there are very few signals overlapped to it. Please use `attr(mat, 'failed_rows')` to get the index of the failed row and consider to remove it.\n")), collapse = "\n")
			} else {
				msg = paste(strwrap(paste0("Smoothing are failed for ", length(failed_rows), " rows because there are very few signals overlapped to them. Please use `attr(mat, 'failed_rows')` to get the index of failed rows and consider to remove them.\n")), collapse = "\n")
			}
			msg = paste0("\n", msg, "\n")
			warning_wrap(msg)
		}
		background = NA
	}

	if(k_upstream > 0) {
		upstream_index = seq_len(k_upstream)
	} else {
		upstream_index = integer(0)
	}
	if(k_target > 0) {
		target_index = seq(k_upstream + 1, k_upstream + k_target)
	} else {
		target_index = integer(0)
	}
	if(k_downstream > 0) {
		downstream_index = seq(k_upstream + k_target + 1, k_upstream + k_target + k_downstream)
	} else {
		downstream_index = integer(0)
	}

	if(length(extend) == 1) {
		if(length(upstream_index) == 0 && length(downstream_index) == 0) {
			extend = c(0, 0)
		} else if(length(upstream_index) == 0) {
			extend = c(0, extend)
		} else if(length(downstream_index) == 0) {
			extend = c(extend, 0)
		} else {
			extend = rep(extend, 2)
		}
	} else if(length(extend) > 2) {
		stop_wrap("length of `extend` should only be 1 or 2.")
	}
	
	attr(mat, "upstream_index") = upstream_index
	attr(mat, "target_index") = target_index
	attr(mat, "downstream_index") = downstream_index
	attr(mat, "extend") = extend
	attr(mat, "smooth") = smooth
	attr(mat, "signal_name") = signal_name
	attr(mat, "target_name") = target_name
	attr(mat, "target_is_single_point") = FALSE
	attr(mat, "failed_rows") = failed_rows
	attr(mat, "background") = background
	attr(mat, "signal_is_categorical") = FALSE

	.paste0 = function(a, b) {
		if(length(a) == 0 || length(b) == 0) {
			return(NULL)
		} else {
			paste0(a, b)
		}
	}

	# dimension names are mainly for debugging
  	if(!attr(mat, "target_is_single_point")) {
  		colnames(mat) = c(.paste0("u", seq_along(upstream_index)), .paste0("t", seq_along(target_index)), .paste0("d", seq_along(downstream_index)))
  	} else {
  		colnames(mat) = c(.paste0("u", seq_along(upstream_index)), .paste0("d", seq_along(downstream_index)))
  	}
	q1 = quantile(mat, keep[1], na.rm = TRUE)
	q2 = quantile(mat, keep[2], na.rm = TRUE)
	mat[mat <= q1] = q1
	mat[mat >= q2] = q2
  	
	class(mat) = c("normalizedMatrix", "matrix")
	return(mat)
}

# == title
# Indices of Rows Failed from Smoothing
#
# == param
# -m Matrix from `normalizeToMatrix`.
#
# == value
# A numeric vector or ``NULL``.
#
failed_rows = function(m) {
	attr(m, "failed_rows")
}
