test_that("wang.paired.ci matches published ExactCIdiff examples", {
  ci1 <- wang.paired.ci(3, 1, 0, conf.level = 0.95,
                        precision = 0.0001)
  ci2 <- wang.paired.ci(2, 0, 2, conf.level = 0.95,
                        precision = 0.0001)
  ci3 <- wang.paired.ci(1, 1, 2, conf.level = 0.95,
                        precision = 0.0001)

  expect_equal(ci1$estimate, 0.75)
  expect_equal(ci1$ExactCI, c(-0.2494, 0.9937), tolerance = 1e-4)
  expect_equal(ci2$ExactCI, c(-0.8649, 0.8649), tolerance = 1e-4)
  expect_equal(ci3$ExactCI, c(-0.8832, 0.6430), tolerance = 1e-4)
})

test_that("wang.paired.reject agrees with individual interval inversion", {
  N <- 4
  tabs <- expand.grid(n11 = 0:N, n10 = 0:N, n01 = 0:N)
  tabs$n00 <- N - tabs$n11 - tabs$n10 - tabs$n01
  tabs <- tabs[tabs$n00 >= 0, , drop = FALSE]
  rownames(tabs) <- NULL

  individual <- vapply(
    seq_len(nrow(tabs)),
    function(i) {
      ci <- wang.paired.ci(
        n10 = tabs$n10[i],
        t = tabs$n11[i] + tabs$n00[i],
        n01 = tabs$n01[i],
        precision = 0.0001
      )$ExactCI
      ci[1] > 0 || ci[2] < 0
    },
    logical(1)
  )

  bulk <- wang.paired.reject(tabs, precision = 0.0001)
  expect_identical(bulk, individual)
})

test_that("wang.paired.ci validates inputs", {
  expect_error(wang.paired.ci(-1, 1, 0), "nonnegative")
  expect_error(wang.paired.ci(1, 1, 0, conf.level = 1), "conf.level")
  expect_error(wang.paired.ci(1, 1, 0, grid.one = 1), "grid")
  expect_error(wang.paired.reject(data.frame(x = 1)), "n11")
})
