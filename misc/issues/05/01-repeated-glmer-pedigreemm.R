# Addressing issue #5 using data & models
# from http://www.animalsciencepublications.org/publications/jas/abstracts/88/2/497

### inc
library(pedigreemm)

#library(lme4qtl)
library(devtools)
load_all("~/git/variani/lme4qtl/")

### data
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

### genetic relationship matrix
A <- getA(pedCowsR)

### fit `*lmer` models
m1 <- relmatLmer(sdMilk ~ lact + log(dim) + (1|id) + (1|herd),
  data = milk_subset, relmat = list(id = A))

m2 <- try(relmatGlmer(sdMilk ~ lact + log(dim) + (1|id) + (1|herd),
  data = milk_subset, relmat = list(id = A), 
  family = "poisson"))

# > packageVersion("lme4qtl")
#[1] ‘0.1.9’

#> m2
#[1] "Error : nrow(Ztlist[[i]])%%ncol(Ztlist[[i]]) == 0 is not TRUE\n"
#attr(,"class")
#[1] "try-error"
#attr(,"condition")
#<simpleError: nrow(Ztlist[[i]])%%ncol(Ztlist[[i]]) == 0 is not TRUE>

