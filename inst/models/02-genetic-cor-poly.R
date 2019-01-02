### inc
library(lme4)

library(dplyr)
library(magrittr)

library(ggplot2)
library(cowplot)

### load simulated dataset `dat40`
data(dat40)

m1 <- relmatLmer(trait1 ~ (1|ID), dat40, relmat = list(ID = kin2))
h1 <- with(VarProp(m1), prop[grp == "ID"])

m2 <- relmatLmer(trait2 ~ (1|ID), dat40, relmat = list(ID = kin2))
h2 <- with(VarProp(m2), prop[grp == "ID"])

### sex-specificity model
m3 <- relmatLmer(trait2 ~ (0 + SEX|ID), dat40, relmat = list(ID = kin2))

# prepare data & fit bi-variate model
mat <- as.matrix(kin2)
mat_upper <- cbind(mat, mat)
mat_all <- rbind(mat_upper, mat_upper)
tkin2 <- Matrix(mat_all)

dat <- bind_rows(
    mutate(dat40, trait = trait1, tindex = 1, tname = "trait1"),
    mutate(dat40, trait = trait2, tindex = 2, tname = "trait2")) %>%
  mutate(ID2 = ID, RID = ID, TID = paste(tindex, ID, sep = "_"), TRID = TID)

rownames(tkin2) <- dat$TID
colnames(tkin2) <- dat$TID

m12 <- relmatLmer(trait ~ (0 + tname|ID), dat, relmat = list(ID = kin2))  

#m12 <- relmatLmer(trait ~ (0 + tname|TID), dat, relmat = list(TID = tkin2))  
#m12 <- relmatLmer(trait ~ (0 + tname|TID) + (0 + dummy(tname)|TRID), dat, relmat = list(TID = tkin2))

#m12 <- relmatLmer(trait ~ (1|ID) + (0 + tindex|ID2), dat, relmat = list(ID = kin2))
