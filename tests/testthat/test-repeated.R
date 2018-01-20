context("repeated measurements")

test_that("relmatLmer/relmatGlmer", {
  stopifnot(require(pedigreemm))
  
  data(milk, package = "pedigreemm")

  milk <- within(milk, {
    id <- as.character(id)
    sdMilk <- milk / sd(milk) # outcome for `lmer`
    sdMilk_counts <- round(sdMilk, 0) # outcome for `glmer`
  })

  ids <- with(milk, id[sire %in% 1:3]) # for more ids: c(1:3, 319-321)
  milk_subset <- subset(milk, id %in% ids)

  milk_subset <- within(milk_subset, {
     herd <- droplevels(herd)
     herd_id <- paste0("herd", seq(1, length(herd)))
  })
  
  A <- pedigreemm::getA(pedCowsR)

  # lmer
  m1 <- try(relmatLmer(sdMilk ~ lact + log(dim) + (1|id) + (1|herd), 
    data = milk_subset, relmat = list(id = A)))

  m2 <- try(relmatGlmer(sdMilk_counts ~ lact + log(dim) + (1|id) + (1|herd),
    data = milk_subset, relmat = list(id = A), 
    family = poisson, nAGQ = 0))
  
  expect_true(class(m1) == "lmerMod")
  expect_true(class(m2) == "glmerMod")
})

