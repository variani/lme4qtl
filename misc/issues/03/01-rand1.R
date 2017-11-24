dataX <- read.table("lme4qtl_data.csv", sep=",", dec=".", header=TRUE, skipNul = TRUE)
library(ape)
tree <- read.nexus("tree_all.tre")
phylo_all <- vcv.phylo(tree, cor = T) 

library(coxme)
m_lmekin <- lmekin(variable ~ mass*habitat + (1|species), data=dataX,
                      varlist=phylo_all, method = "REML")
print(m_lmekin)

#  Log-likelihood = -111.105 
#  n= 1868 
#
# Model:  variable ~ mass * habitat + (1 | species) 
#
# Fixed coefficients
#                   Value  Std Error      z       p
# (Intercept)    1.0926508 0.22466203   4.86 1.2e-06
# mass           0.1483706 0.05843213   2.54 1.1e-02
# habitatU      -0.2663249 0.01286789 -20.70 0.0e+00
# mass:habitatU -0.1065064 0.02695353  -3.95 7.8e-05

library(lme4qtl)
m_lme4qtl <- relmatLmer(variable ~ mass*habitat + (1|species), data=dataX, 
                     relmat = list(species = phylo_all), REML = T)

summary(m_lme4qtl)[["logLik"]]

# 'log Lik.' -107.7412 (df=6)

summary(m_lme4qtl)[["coefficients"]]

#               Estimate Std. Error    t value
# (Intercept)    0.9814953 0.23427874   4.189434
# mass           0.3642359 0.15830427   2.300860
# habitatU      -0.2655202 0.01284164 -20.676496
# mass:habitatU -0.1103905 0.02686246  -4.109473


library(visreg)
# CI bands for mass are very wide even though the model indicates a signficant interaction:
visreg(m_lme4qtl, "mass", "habitat", overlay=F, partial=T, band=T, xlab="centered and scaled log10(mass)")

