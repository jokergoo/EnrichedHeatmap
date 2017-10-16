signal = GRanges(seqnames = "chr1", 
	  ranges = IRanges(start = c(1, 4, 7, 11, 14, 17, 21, 24, 27),
                     end = c(2, 5, 8, 12, 15, 18, 22, 25, 28)),
    score = c(1, 2, 3, 1, 2, 3, 1, 2, 3))
target = GRanges(seqnames = "chr1", ranges = IRanges(start = c(10, 25), end = c(20, 30)))
normalizeToMatrix(signal, target, extend = 10, w = 2)
normalizeToMatrix(signal, target, extend = 10, w = 2, include_target = TRUE)
normalizeToMatrix(signal, target, extend = 10, w = 2, value_column = "score")


normalizeToMatrix(signal, target, extend = 10, w = 2, target_ratio = 0.5)
normalizeToMatrix(signal, target, extend = 10, w = 2, target_ratio = 1)
normalizeToMatrix(signal, target, w = 2, extend = c(5, 5))
normalizeToMatrix(signal, target, w = 2, extend = c(0, 5))
normalizeToMatrix(signal, target, w = 2, extend = c(5, 0))
normalizeToMatrix(signal, target, w = 2, extend = c(0, 0))



