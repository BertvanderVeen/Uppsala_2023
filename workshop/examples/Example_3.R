library(gllvm)
data("Skabbholmen")

# Associated article: https://link.springer.com/article/10.1007/BF00038697
# Skabbholmen is a list with 3 components: Y (cover classes), X (transectID, Elevation, Year)
# and "species" which includes the full species names
# Models from the presentation are included and can be loaded as:
# load("Skabbholem.UO.RData")
# load("Skabbholmen.CO.RData")
# load("Skabbholmen.CN.RData")

# Prepare the data
# If we leave Year the same, we will get in trouble
# As the coefficients will be very large (intercept at year 0)
# Resulting in lack of convergence
Skabbholmen$X$Year <- Skabbholmen$X$Year - min(Skabbholmen$X$Year)
# Centering & standardising variables is -very- important to improve convergence
Skabbholmen$X$Elevation <- scale(Skabbholmen$X$Elevation)

########### Analysis ###########
# Unconstrained ordination with 2LVs and random row-effect (num.lv)
model.UO <- gllvm(y = Skabbholmen$Y, 
                  studyDesign = Skabbholmen$X, num.lv = 2, 
                  row.eff = ~(1|transectID), family = "ordinal", 
                  zeta.struc = "common", n.init = 3, trace = T)
# Constrained ordination with 2LVs and random row-effect (num.RR)
# NOTE: DO NOT INCLUDE ROW-EFFECT IN BOTH X AND "STUDYDESIGN" (known bug)
model.CO <- gllvm(y = Skabbholmen$Y, X = Skabbholmen$X, 
                  num.RR = 2, lv.formula = ~Elevation+Year, 
                  row.eff = ~(1|transectID), family = "ordinal", 
                  zeta.struc = "common", n.init = 3, trace = T)

# Concurrent ordination with 2LVs and random row-effect (num.lv.c)
model.CN <- gllvm(y = Skabbholmen$Y, X = Skabbholmen$X, 
                  num.lv.c = 2, lv.formula = ~Elevation+Year, 
                  row.eff = ~(1|transectID), family = "ordinal", 
                  zeta.struc = "common", n.init = 3, trace = T)

# Model comparison with LRT
# Can only compare -nested- models
# UO is nested in CN but not in CO
anova(model.UO, model.CN)
# and, CO is nested in CN
anova(model.CO, model.CN) #CN better in both instances

# Procrustes rotation
# This can be useful
# Not for model selection, but for assessing how similar ordinations are
# the "getLV" function retrieves the predicted for the ordination axes
# below, 1 means that two ordinations are very different, and 0 that they are same
# This also works for non-nested models
vegan::procrustes(getLV(model.UO), getLV(model.CO), symmetric=T) # Very different
vegan::procrustes(getLV(model.UO), getLV(model.CN), symmetric=T) # Somewhere in the middle
vegan::procrustes(getLV(model.CO), getLV(model.CN), symmetric=T) # Also different

# Examine results for CN
gllvm::ordiplot(model.CN, biplot=T) # can make arrows brighter by arrow.CI=F
coefplot(model.CN)

summary(model.CN, rotate = TRUE) # the "rotate" option let's you choose
# between examining the raw model coefficients (FALSE) or the ones
# rotated to what "ordiplot" shows (TRUE)
# we can also turn off the rotation in ordiplot instead (rotate = FALSE in ordiplot)