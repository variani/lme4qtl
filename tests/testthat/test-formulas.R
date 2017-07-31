context("Formula")

test_that("basic formula", {
  f <- trait1 ~ -1 + AGE + (0 + SEX|ID) + (1|FAMID)

  out <- parseFormula(f)

  expect_true(length(out$labels) == 3) # term `-1` is not counted
  expect_true(all(c("ID", "FAMID") %in% out$ranef2))
})
