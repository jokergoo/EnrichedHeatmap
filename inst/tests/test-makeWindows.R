

test_that("test w", {
	gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 21), end = c(10, 20, 30)))

	x = makeWindows(gr, w = 2)[1:5]
	expect_equal(start(x), c(1, 3, 5, 7, 9))
	expect_equal(end(x), c(2, 4, 6, 8, 10))

	x = makeWindows(gr, w = 0.2)[1:5]
	expect_equal(start(x), c(1, 3, 5, 7, 9))
	expect_equal(end(x), c(2, 4, 6, 8, 10))

	x = makeWindows(gr, w = 3)[1:5]
	expect_equal(start(x), c(1, 4, 7, 11, 14))
	expect_equal(end(x), c(3, 6, 9, 13, 16))

	x = makeWindows(gr, w = 3, direction = "reverse")[1:5]
	expect_equal(start(x), c(2, 5, 8 ,12, 15))
	expect_equal(end(x), c(4, 7, 10, 14, 17))

	x = makeWindows(gr, w = 3, short.keep = TRUE)[1:5]
	expect_equal(start(x), c(1, 4, 7, 10, 11))
	expect_equal(end(x), c(3, 6, 9, 10, 13))

	x = makeWindows(gr, w = 3, direction = "reverse", short.keep = TRUE)[1:5]
	expect_equal(start(x), c(1, 2, 5, 8, 11))
	expect_equal(end(x), c(1, 4, 7, 10, 11))

	x = makeWindows(gr, w = 12)
	expect_equal(length(x), 0)

	gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 31), end = c(10, 30, 70)))

	x = makeWindows(gr, w = 2)
	expect_equal(length(x), 35)

	x = makeWindows(gr, w = 0.2)
	expect_equal(length(x), 15)

})

test_that("test k", {
	gr = GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 11, 21), end = c(10, 20, 30)))
	x = makeWindows(gr, k = 2)[1:5]
	expect_equal(start(x), c(1, 6, 11, 16, 21))
	expect_equal(end(x), c(5, 10, 15, 20, 25))

	x = makeWindows(gr, k = 3)[1:5]
	expect_equal(start(x), c(1, 4, 7, 11, 14))
	expect_equal(end(x), c(3, 6, 10, 13, 16))
})
