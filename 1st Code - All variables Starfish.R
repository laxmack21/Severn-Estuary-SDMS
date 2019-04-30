library(raster)
# Package to handle raster-formatted spatial data
library(rasterVis)
# The rasterVis package complements the raster package, providing a set of methods for enhanced visualization and interaction
# Defines visualisation methods with 'levelplot'
library(dismo)
# Dismo has the SDM analyses for maxent and support vector machines used by R
library(rgeos)
# To define circles with a radius around the subsampled points
# geos is a geometry engine, need to install package to access these capabilities (such as defining circumfrances)
library(rJava)
library(rgdal)
# Provides access to projection/transformation operations from a different library
# Coordinate referancing system**
library(sp)
# Coordinate referancing system
library(ncdf4)
# Opens access to read and write on netCDF files
library(kernlab)
# Required for support vector machines
# installed and running BUT UNSURE of function
library(grDevices)
# For colouring maps
library(colorRamps)
#Allows easy construction of color palettes


#Loading data for project now
#Ensure WD is in correct place WILL BE IN NEW PLACE FOR EACH SPECIES
setwd("~/Documents/UoY/Dissertation/Common Starfish")
locs = read.csv("Common Starfish Severn.csv", header=T, sep = ",")

#loading severn files
#had to add the file location for R to access the severn files, is this right?
dry_always<-raster("Severn_unaltered copy/always_dry_masked.tif")
tidal_range<-raster("Severn_unaltered copy/tidal_range_masked.tif")
subtidal<-raster("Severn_unaltered copy/subtidal_masked.tif")
min_elev<-raster("Severn_unaltered copy/min_elev_masked.tif")
max_velocity<-raster("Severn_unaltered copy/max_vel_masked.tif")
max_elev<-raster("Severn_unaltered copy/max_elev_masked.tif")
mask_2<-raster("Severn_unaltered copy/mask2.tif")
intertidal<-raster("Severn_unaltered copy/intertidal_masked.tif")
depth<-raster("Severn_unaltered copy/bathy_masked.tif")
avg_velocity<-raster("Severn_unaltered copy/av_vel_masked.tif")
#ALL raster data is uploaded here
mask<-depth

#DO NOT HAVE 'distance_to_coast' comparison in our data set as in MaxEnt Code
#DO NOT HAVE 'lat and lon' tifs as in MaxEnt Code


# Extract depth values to table of species co-ordinates
locs_ext=extract(depth, locs[,c("X","Y")])
#this has created a VALUE of depth for each single point as dictated by x&y coordinates from species data
#now each species seen has a depth based on its coordinates in the depth raster file we are given!!

# Build a data frame of species occurrence data and depth data
locs = data.frame(locs, locs_ext)
# added locs_ext to the final column in locs file so now coordinates for species can be coupled with their depth in teh same file

# Remove points with NA values for depth, i.e. on land
locs = subset(locs, !is.na(locs_ext))
e = extent(depth)
#subset extracted all values and rows with 'na' from the locs_ext column
# WHAT DOES EXTENT DO?!
# without using the 'mask' technique above will this still remove all 'land' data above?
#what is "e"?? - is it simply giving the 'extent' of the data set in a min and max of x and y?


# Create sequences of X and Y values to define a grid
# this a 1x1 km grid
xgrid = seq(e@xmin,e@xmax,1000)
ygrid = seq(e@ymin,e@ymax,1000)
#"seq()" works by 'from', 'to', 'by incremental step'
#generated a sequence from xmin value to xmax value in "e" that increase by 1000

# Identify occurrence points within each grid cell, then draw one at random
subs = c()
for(i in 1:(length(xgrid)-1)) {
  for(j in 1:(length(ygrid)-1)) {
    gridsq = subset(locs, Y > ygrid[j] & Y < ygrid[j+1] & X > xgrid[i] & X < xgrid[i+1])
    if(dim(gridsq)[1]>0) {
      subs = rbind(subs, gridsq[sample(1:dim(gridsq)[1],1 ), ])
    }
  }
}
dim(locs);dim(subs) # Confirm that you have a smaller dataset than you started with (1st number)
#for is an argument that will loop a desired action on a given value in a vector
#length will get value the legth of vectors and factors in a defined object
##this a loop going through x values (every 1000m) and at each new x square, looping through all the y's related to that x (and so on for all the x values)
#gridsq is a complex way of saying the square is greater than the start of one x/y value and less than the next one after it
#rbind & cbind combine/create a matrix by rows (rbind) or columns (cbind) of the two seperate vector sets


# Assign correct co-ordinate reference system to subset
coordinates <- cbind(subs$X, subs$Y)
subs_df <- SpatialPointsDataFrame(coordinates, subs, proj4string=CRS("+proj=utm +zone=30 ellps=WGS84"))
#cbind of subs$X and subs$Y created a new data set/matrix called coordinates that only has coordinate data in it!

# we create 20,000 random "background points". There are other ways to do this, but start with this.
#NOTE
psa <- randomPoints(mask, 20000, ext=e)
#need to make sure all is up-to-date: previous error due to 'dismo' not being updated

# Stack raster layers into one variable
#NOTE (make it match your environmental variables from above)
env_uk<-stack(depth,max_elev,max_velocity,min_elev,avg_velocity,dry_always,intertidal,subtidal,tidal_range)
#NEED TO CHECK THAT IS ALL OF THEM

# Pull environmental data for the sumbsampled-presence points from the raster stack
presence_uk= extract(env_uk, subs_df[,c("X","Y")])
#Warning messages: transforming SpatialPoints to the CRS of the Raster?

# Pull environmental data for the pseudo-absence points from the raster stack
pseudo_uk = extract(env_uk, psa)

# Build some useful dataframes, with two columns of coordinates followed by the environmental variables. For the presence points:
presence_uk = data.frame(X=subs_df$X, Y=subs_df$Y, presence_uk)
#HOW IS THIS DIFFERENT TO ABOVE FUCNTION WITH "EXTRACT"?

# Convert psa from atomic vector matrix to data.frame
psapoints=data.frame(psa)
# Bind co-ordinates
coordinates <- cbind(psapoints$x, psapoints$y)
# Create spatial data frame of pseudo absences
psadf <- SpatialPointsDataFrame(coordinates, psapoints, proj4string=CRS("+proj=utm +zone=30 ellps=WGS84"))

# Build dataframe, with two columns of coordinates followed by the 5 environmental variables. For the pseudo-absences:
psadfx = psadf@coords
colnames(psadfx) = c("X","Y")
pseudo_uk = data.frame(cbind(psadfx,pseudo_uk))


# Vector of group assignments splitting the subsampled presence points data fram with environmental data into 5 groups
group_p = kfold(presence_uk, 5)
#kfold partitions a data set k times (in this case 5 times) for model testing purposes

# Repeat above step for pseudo-absence points
group_a = kfold(pseudo_uk, 5)

# create output required for the loop
evaluations = list(5)
models = list(5)


# where it says maxent - you may need to swap this for other functions if you're exploring different models
# Note that some model may need different inputs etc. Read the docs to figure this out.

# This is our k-fold test. You will want to spend a bit of time making predictions on each of the 5 sub-models
# created here to check you can make decent predictions even with missing data
for (test in 1:5) {
  
  # Then we use test and the kfold groupings to divide the presence and absence points:
  train_p = presence_uk[group_p!=test, c("X","Y")]
  train_a = pseudo_uk[group_a!=test, c("X","Y")]
  test_p = presence_uk[group_p==test, c("X","Y")]
  test_a = pseudo_uk[group_a==test, c("X","Y")]
  
  # Now, estimate a maxent model using the "training" points and the environmental data. This may take a few moments to run:
  models[test] = maxent(env_uk, p=train_p, a=train_a)
  
  # To validate the model, we use the appropriately named function.
  # Produces warning message about implicit list embedding being depreciated. May fail in future versions of R
  evaluations[test] = evaluate(test_p, test_a, models[[test]], env_uk)
}


# print out the AUC for the k-fold tests
# ideally should be > 0.75 for all
cat("K-FOLD AUC: ")
for (test in 1:5) {
  cat(paste0(evaluations[[test]]@auc,","))
}
cat("\n")

#IF ONE WANTED to visualise the 1st model NEED TO FINISH WITH THE OTHER 4 KFOLDS
## pred <- predict(models[[1]], env_uk)
## plot(pred)
#-could visualise all k-fold models to show they are the same (strong statistical result)
#BUT DONT NEED TO


# Assess Spatial Sorting Bias (SSB)
pres_train_me <- train_p
pres_test_me <- test_p
back_train_me <- train_a
back_test_me <- test_a
sb <- ssb(pres_test_me, back_test_me, pres_train_me)
sb[,1] / sb[,2]

#creates a model of spacial biasing to compare to given preditions 
# Adjust for SSB if present via distance based point-wise sampling
i <- pwdSample(pres_test_me, back_test_me, pres_train_me, n=1, tr=0.1)
pres_test_pwd_me <- pres_test_me[!is.na(i[,1]), ]
back_test_pwd_me <- back_test_me[na.omit(as.vector(i)), ]
sb2 <- ssb(pres_test_pwd_me, back_test_pwd_me, pres_train_me)
sb2[1]/ sb2[2]


#creates full model without any K fold statistics etc 
pres_points = presence_uk[c("X","Y")]
abs_points = pseudo_uk[c("X","Y")]
# create full maxent with all points
model <- maxent(env_uk, p=pres_points, a=abs_points)
#turn model into prediction that can be plotted into a raster
pred_PredFull <- predict(model, env_uk)
#to see model and obtain jpeg
plot(pred_PredFull)

#Gives AUC for full model (pred_PredFull)
evaluate_full <- evaluate(presence_uk[c("X","Y")], pseudo_uk[c("X","Y")], model, env_uk)
#see what AUC is by typing it in
evaluate_full

#see what the specific sensitivity values is of species
#use value givenas a base level or higher that one would expect to see species (compare to evaluate_full)
message(threshold(evaluate_full)$spec_sens)


#check response curves to see if they change the FULL MODEL
response(model)

#creates a raster file in the WD 
#will be useful when putting a file into qgis!
writeRaster(pred_PredFull, filename="pred1_me.tif", options="INTERLEAVE=BAND", overwrite=TRUE)

