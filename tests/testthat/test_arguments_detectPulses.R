library(NullTurnover)

context('ensure detectPulses_ functions have identical arguments')

test_that("argument lists are identical", {
  expect_identical(names(formals(detectPulses)), names(formals(detectPulsesChiSq)))
})