### data
data(dat40, package = "lme4qtl")

head(dat40)
image(kin2[1:15, 1:15])

#### fit model 1
m1 <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))

### get a point estimate of heritability, the proportion of variance explained by (1|ID)
lme4::VarCorr(m1)
lme4qtl::VarProp(m1)

### try to profile and get CI...
prof1 <- profile(m1, which = "theta_", prof.scale = "varcov")
head(prof1) # many NAs

### fit model 2: dropping not significant terms AGE + SEX + (1|FAMID)
m2 <- relmatLmer(trait1 ~ (1|ID), dat40, relmat = list(ID = kin2))

prof2 <- profile(m2, which = "theta_", prof.scale = "varcov")
head(prof2)

prof2 <- profile(m2, which = "theta_", prof.scale = "sdcor")
head(prof2)

confint(prof2, level = 0.95)

prop2 <- try(varpropProf(prof2, prof.scale = "sdcor")) # error: system is computationally singular

### fit model 3
m3 <- relmatLmer(trait1 ~ (1|FAMID), dat40)

prof3 <- profile(m3, which = "theta_", prof.scale = "sdcor")
head(prof3)

confint(prof3, level = 0.95)

prop3 <- varpropProf(prof3, prof.scale = "sdcor")
