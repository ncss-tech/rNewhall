context("moistureCondition")

m <- matrix(rep(0, times = 64), ncol = 8, byrow = TRUE)


test_that("working as expected", {

  # standard 8x8 configuration

  # "dry"
  m[c(9, 17, 25)] <- 0
  expect_true(
    moistureCondition(m) == 1
  )

  # "dry/moist"
  m[c(9, 17, 25)] <- c(0, 1, 0)
  expect_true(
    moistureCondition(m) == 2
  )

  # "moist"
  m[c(9, 17, 25)] <- 1
  expect_true(
    moistureCondition(m) == 3
  )

})


test_that("expected errors", {

  # incorrect input length
  expect_error(moistureCondition(1))

})
