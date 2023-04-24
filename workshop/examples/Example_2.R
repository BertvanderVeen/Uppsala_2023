library(vegan)
library(gllvm)
data("pyrifos")

# Associated article: https://setac.onlinelibrary.wiley.com/doi/10.1002/etc.5620180207
# Pyrifos is a data frame of abundances
# Models from the presentation are included and can be loaded as:
# load("pyrifos.unconstrained.RData")
# load("pyrifos.time.RData")

# Prepare the data
pyrifos<-(exp(pyrifos)-1)/10 # Backtransform: we don't want to use GLLVM on transformed data

# Create predictors
ditch <- gl(12, 1, length=132)
week <- gl(11, 12, labels=c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24))
dose <- factor(rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11))
X <- data.frame(ditch,week,dose)

########### Analysis ###########
# Now we can start the analysis
# Let's immediatley fit a model with 2LVs
# You can afterall not make a good ordination plot with fewer.
mod <- gllvm(y = pyrifos, num.lv = 2, family = "tweedie")

# The ordination plot
# Here, sites (rows) are colored by "week since treatment"
# ordiplot re-rotates the model solution by default. This can be turned off with rotate = FALSE
gllvm::ordiplot(mod,s.col=viridis::viridis(length(unique(X$week)))[as.numeric(as.factor(X$week))])

# Now we will include "week" as a predictor
# This should remove its effects from the ordination plot
mod.week<-gllvm(y = pyrifos, X = X, formula = ~as.numeric(week), num.lv = 2, family = "tweedie")
gllvm::ordiplot(mod.week,s.col=viridis::viridis(length(unique(X$week)))[as.numeric(as.factor(X$week))])

# It looks like some effect has indeed been removed
# Let's examine the year effect per species
coefplot(mod.week) # Grey = CI crosses zero

# We can also make a biplot (i.e., plot with two quantities)
gllvm::ordiplot(mod.week, biplot =T, predict.region=T, s.col = scales::alpha("gray", 0.5), col.ellips = scales::alpha("gray", 0.5))

# Compare models
anova(mod, mod.week) # This will warn that we are doing something stupid
AIC(mod, mod.time) # Safer to compare models with a large number of parameters differences