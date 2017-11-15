
# https://gis.stackexchange.com/questions/225102/calculate-distance-between-points-and-nearest-polygon-in-r

# Let's set up some sample data:
library(sp)
library(spdep)
example(columbus)
plot(columbus)

# We'll now get some points. Click 20 times on the map, some inside and some outside polygons:
pts = locator(20,type="p")

# Convert to Spatial data type:
spts = SpatialPoints(pts)

# Now use rgeos to compute point-polygon distances, and take the minimum for each point:
library(rgeos)

apply(gDistance(spts, columbus,byid=TRUE),2,min)

# There you go. That's the minimum distance from each of the 20 points to any of 
# the polygons. The zeroes are for points inside one of the polygons.
# Note you should use data in a projected spatial coordinate system and not lat-long.

################################################################################
# For lat lon data, try this alternative approach
################################################################################
library(sp)
library(geosphere)

# some country polygons to try on
data(wrld_simpl, package = "maptools")
wrld_subset <- wrld_simpl[wrld_simpl@data$ISO2 %in% c("RO","HU","AT","DE","FR"),]

# Generate random points (in and out)
set.seed(2017)
pts <- sp::makegrid(wrld_subset, n = 5)

plot(wrld_subset)

# compute the shortest distance between points and polygons
# (from ?dist2Line): "returns matrix with distance and lon/lat of the nearest point" & 
# "the ID (index) of (one of) the nearest objects"; distance is in meters (default)
dist.mat <- geosphere::dist2Line(p = pts, line = wrld_subset)

# bind results with original points
pts.wit.dist <- cbind(pts, dist.mat)
pts.wit.dist[1:3,]
##      x1   x2 distance        lon      lat ID
## 1 -10.0 40.2 767133.6 -1.7808770 43.35992  2
## 2  -0.2 40.2 282022.2  0.1894124 42.71846  2
## 3   9.6 40.2 134383.0  9.1808320 41.36472  2


# Plot some results to get an idea

pts.sp <- sp::SpatialPoints(coords = pts[,c("x1","x2")], # order matters
                            proj4string = wrld_subset@proj4string)
plot(pts.sp, col="red")
plot(wrld_subset, add=TRUE)
# plot arrows to indicate the direction of the great-circle-distance
for (i in 1:nrow(pts.wit.dist)) {
  arrows(x0 = pts.wit.dist[i,1], 
         y0 = pts.wit.dist[i,2], 
         x1 = pts.wit.dist[i,4], 
         y1 = pts.wit.dist[i,5],
         length = 0.1,
         col = "green")
}

