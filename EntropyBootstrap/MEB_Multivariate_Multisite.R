#==================================================================================================
#==================================================================================================
# MULTISITE MULTIVARIATE STOCHASTIC VARIABLE GENERATOR USING MAXIMUM ENTROPY BOOTSTRAP
# Roshan K. Srivastav
# Facility for Intelligent Decision Support
# Department of Civil Engineering
# University of Western Ontario
# July, 2014
# 
# Citation: Roshan K. Srivastav, Slobodan P. Simonovic (2014), Multi-site, multivariate weather 
# generator using maximum entropy bootstrap, Climate Dynamics, doi: 10.1007/s00382-014-2157-x
#
# Roshan K. Srivastav, Slobodan P. Simonovic (2014), An Analytical Procedure for Multi-Site, 
# Multi-Season Streamflow Generation using Maximum Entropy Bootstrapping, 
# Environmental Modeling and Software Journal, 59, 1-17.
# http://dx.doi.org/10.1016/j.envsoft.2014.05.005
#
#
# Current Address
# HydroSystems Research Laboratory
# Department of Civil and Environmental Engineering
# IIT Tirupati, Tirupati, India
# Email: roshan@iittp.ac.in
#==================================================================================================
#==================================================================================================


# Required Library
library(meboot)

#========================================================================
#========================================================================

# Read Data: Columns represent variables at different spatial locations
x<-read.csv("data.csv", row.names = NULL)

# Number of Replicates
reps <- 100

# Trimmed mean percentage
trim.per <- 10
#========================================================================
#========================================================================
# Apply maximum entropy bootstrap in combination with PCA

# Fit PCA
fit.pca <- princomp(x.scale, scores=TRUE)

# Extract Scores
z.pca <- fit.pca$scores

# Extract Loadings
l.pca <- fit.pca$loadings

# Generate Ensembles using Maximum Entropy Bootstrap using first PC 
z1 <- meboot(z.pca[,1], reps=reps, trim = 10)
z1.meb <- z1$ensemble

# Generate Replicates
replicates <- list()
for(i in 1:reps){
  z.i <- z.pca
  z.i[,1] <- z1.meb[,i]
  x.recon.z <- z.i %*% t(l.pca)
  x.recon.z <- sweep(sweep(x.recon.z, 2, attr(x.scale, 'scaled:scale'),'*'),2, attr(x.scale, 'scaled:center'), '+')
  replicates[[i]] <- x.recon.z
}


