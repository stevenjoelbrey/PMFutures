#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args=2003
}
year <- args

################################################################################
# assign_ecoregion_to_border_points.R
# This is run right after assign_ecoregion_to_FPAFOD.R. It assigns ecoregion to
# points that did not fall inside of polygons. 

# Learning how to calculate distance between points and polygons was aquired by
# this message thread:
# https://gis.stackexchange.com/questions/225102/calculate-distance-between-points-and-nearest-polygon-in-r
################################################################################

library(maps)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(geosphere)


dataDir   <- "Data/FPA_FOD/"
f <- paste0(dataDir, "FPA_FOD_",year,".RData")
load(f)

# We need PR and HI out of this dataset, as it is outside the scope of this work.
# barplot(prop.table(table(as.factor(FPA_FOD$STATE))))
PRMask <- FPA_FOD$STATE == "PR" | FPA_FOD$STATE == "HI"
FPA_FOD <- FPA_FOD[!PRMask, ]

# Load the region polygons
load("Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData")
SPDF@data$NA_L2CODE <- as.numeric(as.character(SPDF@data$NA_L2CODE))

# Get the centriods for fast subsetting
centriods <- gCentroid(SPDF, byid=TRUE)

# Get fire coordinates 
coordinates <- cbind(FPA_FOD$LONGITUDE, FPA_FOD$LATITUDE)
fireLocations <- SpatialPoints(coordinates)
fireLocations@proj4string <- SPDF@proj4string

# Where are they missing an assignment? 
ptsMissingIndex <- which(is.na(FPA_FOD$NA_L2CODE))
ptsMissing      <- fireLocations[ptsMissingIndex]

# How many assigments do we need to make?
nMissing <- length(ptsMissing)

print(paste("Trying calculation for",nMissing, "pts"))
t1 <- Sys.time()

#for (i in round(runif(100, min=1, max=nMissing))){
for (i in 1:nMissing){

  print(paste(i/nMissing*100))
  
  pt <- ptsMissing[i]
  
  # This is super fast, so get rid of polygons that are really far away from this
  # point. 
  centriodDist <- distHaversine(pt, centriods)
  cuttoffDistance <- sort(centriodDist)[10]
  SPDFMask <- centriodDist <= cuttoffDistance
  
  # Now this calculation can be done on many fewer polygons
  dist.mat <- geosphere::dist2Line(p = pt, line = SPDF[SPDFMask,])
  
  # Assign the closest polygon in the subset ecoregions polygons
  ecoregionAssignment <- SPDF[SPDFMask,]@data$NA_L2CODE[dist.mat[,4]]
  
  # Place assignment back into the original data 
  FPA_FOD$NA_L2CODE[ptsMissingIndex[i]] <- ecoregionAssignment
  
  # # Sanity check figure.
  # quartz()
  # plot(SPDF[SPDFMask,], col="blue") # all polygon options
  # plot(SPDF[SPDF@data$NA_L2CODE==ecoregionAssignment,],add=T, col="green") # chosen
  # plot(pt, add=T, pch=3, col="red")
  # title(i)

}

tf <- Sys.time()
print(paste("It took", tf - t1, "to complete calculation."))

saveName <- paste0(dataDir, "FPA_FOD_borderfix_", year,".RData")
save(FPA_FOD, file=saveName)



