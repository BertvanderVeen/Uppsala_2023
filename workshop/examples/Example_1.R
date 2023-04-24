library(gllvm)
library(Momocs)
data(apodemus)

# Associated article: https://doi.org/10.1046/j.1365-2699.2003.00932.x
# Apodemus is a tibble that contains the responses as the first slot
# and information of individuals in the second slot
# Models from the presentation are included and can be loaded as a list load("mice.mods.RData")

# Step 1: plot the data
# First, we need to arrange the data in a 3d array
dat3d <- array(unlist(apodemus[[1]]),dim=c(nrow(apodemus[[1]][[1]]),2,length(apodemus[[1]])))

# Second to a matrix
dat<-sapply(seq(dim(dat3d)[3]), function(x) matrix(unlist(apply(dat3d[,,x],1,c,simplify = F)),ncol=1))

# Let's plot the data first, always good practice to look at data before analysis
# First we make an empty plot
plot(NA, type="n", xlim = c(5, 25), ylim = c(0, 40), ylab = "y", xlab = "x", main = "Data points") 
cols <- scales::alpha(c("blue", "black", "brown"), 0.5)

# Now iterate over the individuals and plot the coordinates (x,y)
# We color by species (3)
for(i in 1:ncol(dat)){
  points(x = dat3d[,1,i], y = dat3d[,2,i],  col = cols[as.numeric(unlist(apodemus[[2]][,3]))][i])
}

# Add a line for readability
abline(h = c(10, 30), lty = "dashed", col = "red")

########### Analysis ###########
# Now we can start the analysis
# First we will fit a GLLVM with 1 latent variable; starting simple
# You can fit more complex models by increasing the value for "num.lv"
# Below, n.init makes sure the model tries different sets of starting values
# And pick the "best" solution. trace = T keeps us up to date of the progress
library(gllvm)
model.one <- gllvm(y = dat, family = "gamma", num.lv = 1, starting.val = "random", n.init = 3, trace = T)
model.one

# Have a look at residuals
plot(model.one) # Well that does not look great!

# Our model didn't fit the data well
# Maybe it's too simple, so let's fit a model with 4 LVs
# This will take a little longer
model.four <- gllvm(y = dat, num.lv = 4, family = "gamma", starting.val = "random", n.init = 3, trace = T)
plot(model.four)

# Now we can compare the models. Here are some options:
AIC(model.one, model.four);AICc(model.one, model.four);BIC(model.one, model.four)

# model.four was best, let's have a look at the genotypic correlations:
correlations <- getResidualCor(model.four) # get estimated correlation matrix
# and plot the correlations
corrplot::corrplot(correlations, type="lower", order="AOE")