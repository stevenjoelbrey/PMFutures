# mergeUSandCanadaPolygons.R

# This script gets US and canada shapefiles and merges them for easy access to 
# all political borders in North America. 

library(stringr)
library(maps)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(raster)


Canada <- getData('GADM', country='CAN', level=1)
Canada_df <- Canada@data
Canada_simple <- rgeos::gSimplify(Canada, tol = 0.1, topologyPreserve = TRUE)
Canada_simple_spdf <- SpatialPolygonsDataFrame(Canada_simple, Canada_df)

USA <- getData('GADM', country='USA', level=1)
USA_df <- USA@data
USA_simple <- rgeos::gSimplify(USA, tol = 0.1, topologyPreserve = TRUE)
USA_simple_spdf <- SpatialPolygonsDataFrame(USA_simple, USA_df)

# Merge the U.S. and Canada
#north_america <- raster::union(USA_simple_spdf, Canada_simple_spdf)
north_america <- rbind(USA_simple_spdf, Canada_simple_spdf)

quartz()
plot(north_america)
# nCanada <- length(Canada_simple_spdf)
# north_america <- USA_simple_spdf
# for (i in 1:nCanada){
#   print(i)
#   print(paste("Merging:", Canada_df$NAME_1[i]))
#   north_america <- raster::union(north_america, Canada_simple_spdf[i,])
#   
# }

save(north_america, file="Data/GIS/north_america.RData")