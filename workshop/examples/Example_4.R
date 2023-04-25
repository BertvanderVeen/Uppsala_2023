library(gllvm)
data("microbialdata")

# Associated article: https://www.frontiersin.org/articles/10.3389/fmicb.2017.00012/full
# microbialdata is a list with two components: microbialdata$Y (counts) and microbialdata$Xenv (predictors)
# Models from the presentation are included and can be loaded as:
# load("ftNULL.RData")
# load("ftCN.RData")
# load("ftCN2.RData")
# load("ftCN3.RData")

# Prepare data
Ysoil <- microbialdata$Y
# Center & scale predictors
Xenv <- scale(microbialdata$Xenv[,1:3])
# Separate variable for the row-effects
sDesign <- data.frame(Site=microbialdata$Xenv$Site)
# Color vector for ordiplot later
ph <- microbialdata$Xenv$pH
rbPal <- colorRampPalette(c('mediumspringgreen', 'blue'))
Colorsph <- rbPal(20)[as.numeric(cut(ph, breaks = 20))]
pchr = NULL
pchr[microbialdata$Xenv$Region == "Kil"] = 1
pchr[microbialdata$Xenv$Region == "NyA"] = 2
pchr[microbialdata$Xenv$Region == "Aus"] = 3

########### Analysis ###########
# Note: often you might want to include a row-effect or offset
# to account for differences in library size.
# An offset requires knowing libray size a-priori, which we here do not
# But can be included as a matrix (same size as Y) in the "offset" argument of gllvm
# A row-effect per sample is not possibel (currently) in combination with
# a structured row-effect (as we are doing here)

# Fit uconstrained ordination or "JSDM"
# Large number of species, so it might take a while
# Turning sd.errors off will speed things up
# And "EVA" is faster than the default "VA" or "LA"
ftNULL <- gllvm(Ysoil, studyDesign = sDesign, family = "negative.binomial", row.eff = ~(1|Site), num.lv = 2, sd.errors = FALSE, method = "EVA")

# Ordination plot colored by pH
gllvm::ordiplot(ftNULL, main = "Ordination of sites, color: pH",
                symbols = TRUE, pch = pchr, s.colors = Colorsph)

# Let's fit a concurrent ordination and see if this relationship holds up
# The only thing we need to do it add predictors X, and change num.lv to num.lv.c
ftCN <- gllvm(Ysoil, X=Xenv, studyDesign = sDesign, family = "negative.binomial", row.eff = ~(1|Site), num.lv.c = 2, sd.errors = FALSE, method = "EVA")

# Erm, now we actually need standard errors to evaluate uncertainty around the effect
# Luckily, we can also get those afterwards:
# NOTE: THIS TAKES TIME (15min or so?)
sdErr <- se(ftCN) # This takes a long time for large datasets
ftCN$sd <- sdErr$sd
ftCN$Hess <- sdErr$Hess
confint(ftCN,"LvXcoef") # now we can look at the CIs

# Assess how much variation is explained by the predictors
1 - getResidualCov(ftCN)$trace/getResidualCov(ftNULL)$trace

# Or a partial r^2 that also incorporates uncertainty:
source("partR2.gllvm.R")
partR2(ftCN)

# We can add shrinkage. This can make ordination plots nicer.
# We can add shrinkage to "select" predictors
ftCN2 <- update(ftCN,randomB="P")

# OR, LV-based shrinkage
ftCN3 <- update(ftCN, randomB="LV")

# More information on using "randomB" will be available soon in the gllvm vignette.