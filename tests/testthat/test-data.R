test_that("teaching datasets have expected dimensions", {
  data(depress, package = "ibist")
  data(fluorosis, package = "ibist")
  data(midbp, package = "ibist")
  data(mihdl, package = "ibist")
  data(nrs, package = "ibist")
  data(rds, package = "ibist")

  expect_equal(dim(depress), c(1730L, 9L))
  expect_equal(dim(fluorosis), c(24L, 5L))
  expect_equal(dim(midbp), c(14L, 4L))
  expect_equal(dim(mihdl), c(24L, 4L))
  expect_equal(dim(nrs), c(91795L, 6L))
  expect_equal(dim(rds), c(16L, 4L))
})

test_that("aggregated count datasets have valid counts and exposure", {
  data(fluorosis, package = "ibist")
  data(midbp, package = "ibist")
  data(mihdl, package = "ibist")
  data(rds, package = "ibist")

  expect_true(all(fluorosis$count >= 0))
  expect_true(all(midbp$cases >= 0))
  expect_true(all(midbp$pyears > 0))
  expect_true(all(mihdl$cases >= 0))
  expect_true(all(mihdl$pyears > 0))
  expect_true(all(rds$count >= 0))
})
