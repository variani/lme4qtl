context("Argumet relmat")

test_that("relmatLmer", {
  data(dat40)

  expect_warning(relmatLmer(trait1 ~ (1|ID), dat40, relmat = list(MyID = kin2)),
    "not all relmat ID variables")
})

test_that("relmatGlmer", {
  data(dat40)

  expect_warning(relmatGlmer(trait1bin ~ (1|ID), dat40, relmat = list(MyID = kin2),
      family = binomial),
    "not all relmat ID variables")
})
