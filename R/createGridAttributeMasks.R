# createGridAttributeMasks.R

# This script will be used to mask the 0.25 and 0.75 degree grids used in this
# analysis. Each grid will have a layer with the following descriptive attributes
# 1) Area, these will be from GFED4s grids
# 2) state/province
# 3) ecoregion
# 4) elevation 

library(stringr)
library(maps)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(raster)
library(ncdf4)

# Get the quarter degree grid that all data lives on (GFED4s grid). All attribute
# layers will live on the same lat lon grid as this. 
fgrid <- "/Volumes/Brey_external/GFED4s/GFED4.1s_burned_area_2004.nc"
nc <- nc_open(fgrid)
grid_area <- ncvar_get(nc, "grid_area")
longitude <- ncvar_get(nc, "longitude")
latitude  <- ncvar_get(nc, "latitude")
nc_close(nc)

# Load political borders of North America (sans Mexico)
load("Data/GIS/north_america.RData")

# Ecoregions
load("Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData")

# Load the levation data
felev <- "Data/GIS/elev.0.25-deg.nc"
nc <- nc_open(felev)
elevation <- ncvar_get(nc, "data")
elevationUnits <- ncatt_get(nc,"data","units")
elev_longitude <- ncvar_get(nc, "lon")
elev_latitude  <- ncvar_get(nc, "lat")
nc_close(nc)

# These data need to be shifted in the longitude direction to match the GFED
# grid.
elevation_new <- elevation
elevation_new[] <- NA
elevation_new[1:360,] <- elevation[361:720,]
elevation_new[361:720,] <- elevation[1:360,]

quartz()
image(elevation)

elev_latitude - latitude



# point over polygon

# if no assignment, what is closest? 