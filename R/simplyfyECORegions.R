# simplyfyECORegions.R

library(stringr)
library(maps)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)


# Load the level II ecoregions of the U.S. 
layername <- "NA_CEC_Eco_Level2"
layerDir  <- "Data/GIS/na_cec_eco_l2/"
data_projected <- readOGR(dsn=layerDir, layer=layername) 

# Change projection to longlat
assign("data_projected", spTransform(data_projected, CRS("+proj=longlat")))
polyData <- data_projected@data

IDs <- polyData$NA_L2CODE
uniqueIDs <- unique(IDs)

firstMatchIndex <- match(uniqueIDs, IDs)
data_subset <- polyData[firstMatchIndex,]

# Length and area descriptors no longer correct, since the match() command just
# grabbed us one row where the polygon matched. So get rid of these.
keeperColumns <- !(names(polyData) %in% c("Shape_Leng","Shape_Area"))
data_subset <- data_subset[,keeperColumns]
row.names(data_subset) <- uniqueIDs # will make matching polygons and data easy

# Ok now merge the polygons based on the IDs we want to join them by
newPolyObject <- unionSpatialPolygons(data_projected, IDs)
newPolyObject_simple <- gSimplify(newPolyObject, tol = 0.1, topologyPreserve = TRUE)

# NOTE: tol argument chose based on value that leaves no visible plotting gap
# NOTE: in region 6.2 and 10.1 plotted side by side. 

# plot(newPolyObject[1,], border="black", col="black")
# plot(newPolyObject_simple[1,], border="red")

# TODO: Now we want to simplyfy these polygons while keeping the same IDs

# Now create a spatial polygons dataframe with the subset data and merged polys
SPDF <- SpatialPolygonsDataFrame(newPolyObject_simple, data_subset)

m <- uniqueIDs == 6.2

plot(SPDF[m,], col="forestgreen")
map('world', add=T)
map("state", add=T)
plot(SPDF[uniqueIDs == 10.1,],add=T, col="grey")
plot(SPDF[uniqueIDs == 13.1,],add=T, col="red")


save(SPDF, file="Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData")

