###
# Reproducing the errors reported in https://github.com/variani/lme4qtl/issues/7

data(dat40)

# model 1: two random effects
m1 <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))

p1 <- profile(m1, which = "theta_", prof.scale = "varcov")

# model 2: one random effect
m2 <- relmatLmer(trait1 ~ AGE + SEX + (1|ID), dat40, relmat = list(ID = kin2))

p2 <- profile(m2, which = "theta_", prof.scale = "varcov")

# model 3: the simplest model
m3 <- relmatLmer(trait1 ~ (1|ID), dat40, relmat = list(ID = kin2))

p3 <- profile(m3, which = "theta_")

# model 4: lme4 example
m4 <- lmer(Yield ~ 1|Batch, Dyestuff)

p4 <- profile(m4, which = "beta_")

# model 5: model 1 with only FAM effect (so it can be fitted by lme4)
m5 <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID), dat40)

p5 <- profile(m5, which = "theta_", prof.scale = "varcov")

pp5 <- varpropProf(p5)

ci5 <- confint(pp5, level = 0.95)

#> ci5
#                2.5 %    97.5 %
#.sigprop01 0.04146613 0.3965989
#.sigmaprop 0.76210966 0.9112842

