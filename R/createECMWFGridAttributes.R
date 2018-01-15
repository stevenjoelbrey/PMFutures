# createECMWFGridAttributes.R

# This script will be used to mask the 0.75 degree grids used in this
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
library(geosphere)
library(fields)

# # Get the elevation grid
fgrid <- "Data/GIS/elev.0.75-deg.nc"
nc <- nc_open(fgrid)
elevation <- ncvar_get(nc, "data")
longitude <- ncvar_get(nc, "lon")
latitude  <- ncvar_get(nc, "lat")
nc_close(nc)

flipper <- length(latitude):1
quartz()
image.plot(longitude, latitude[flipper], elevation[,flipper])
title("Elevation, not transformed")

# Make sure these dimensions match ecmwf file exactly
externalFile <- "/Volumes/Brey_external/era_interim_nc_daily/u10_2011.nc"
nc <- nc_open(externalFile)
longitude_era <- ncvar_get(nc, "longitude")
latitude_era  <- ncvar_get(nc, "latitude")
nc_close(nc)


dlat <- unique(latitude - latitude_era)
dlon <- unique(longitude - longitude_era) 
if(dlat != 0 | dlon != 0){
  stop("Elevation and ECMWF meteorology do not live on the same grid!")
}

################################################################################
# Calculate area of grid cells 
# http://badc.nerc.ac.uk//help/coordinates/cell-surf-area.html
################################################################################
lat_rad <- latitude * pi/180
lon_rad <- longitude * pi/180
dx <- diff(lon_rad)[1]
R <- 6371*1e3 # radius of earth m

nLat <- length(lat_rad)
area <- rep(NA, nLat)

for (i in 1:(nLat-1)){
  area[i] <- R^2 * dx * (sin(lat_rad[i]) - sin(lat_rad[i+1]))
}
area[nLat] <- area[1]

areaMax <- max(area)

# At equator, where each deg of lat and lon are 111km
mathAreaMax <- (0.75*111*1e3)^2

PercentDifferent <- abs(mathAreaMax - areaMax)/areaMax * 100
print(paste("% Difference in max grid cell area calculated from assumed 111km/deg at equator=", PercentDifferent))

# Assign to each latitude, does not change in longitude
grid_area <- elevation
grid_area[] <- NA
for (j in 1:length(lon_rad)){
  grid_area[j,] <- area
}
# alternatively
#grid_area <- t(replicate(length(lon_rad), area))

################################################################################
# Now overlap calculations for state and ecoregions
################################################################################
lon <- longitude - 180 # not shifting data, just the labeling of this dimension 
psudeoLon <- rep(NA, length(longitude))
psudeoLon[1:(sum(lon>=0))] <- lon[lon >=0]
psudeoLon[(sum(lon>=0)+1):length(lon)] <- lon[lon < 0]

lon <- psudeoLon
lat <- latitude

# point over polygon calculations. Now we need to assign ecoregions as well as
# political borders. 
# Speed this up by only looping through North America coordinates 
lon_index <- which(lon >= -168 & lon <= -52) # needs to be the shifted grid
lat_index <- which(latitude >= 25 & latitude <= 80)

quartz()
f <- length(lat_index):1
image.plot(longitude[lon_index], latitude[lat_index][f], elevation[lon_index, lat_index[f]])
title("This map needs to look like the united states")

# Load ecoregions "SPDF". Use this CRS for all spatial objects
load("Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData")
# Handle the annoying way codes are stored as factors
SPDF@data$NA_L2CODE <- as.numeric(as.character(SPDF@data$NA_L2CODE))
ecoregion_centriods <- gCentroid(SPDF, byid=TRUE)

# Load "north_america"
load("Data/GIS/north_america.RData")
north_america@proj4string <- SPDF@proj4string
north_america_centriods <- gCentroid(north_america, byid=TRUE)

# Make dummy grids that will be filled with overlap values! 
ecoregion <- grid_area
ecoregion[] <- 0
state <- ecoregion
state[] <- ""

nIterations <- length(lon_index) * length(lat_index)
count <- 0
print("Starting work on the main, large loop.")
t1 <- Sys.time()
for (xi in lon_index){
  
  gridLon <- lon[xi]
  
  for (yi in lat_index){
    
    gridLat <- lat[yi]
    
    # Make this grid location a spatial point opject
    pt <- SpatialPoints(cbind(gridLon, gridLat), proj4string = SPDF@proj4string)
    
    ############################################################################        
    # State overlap calculation
    ############################################################################   
    centriodDist <- distHaversine(pt, north_america_centriods)
    cuttoffDistance <- sort(centriodDist)[3] # States do not have wierd shapes so this small number should be fine
    StateMask <- centriodDist <= cuttoffDistance
    # Perform overlap calculation for close points
    state[xi, yi] <- over(pt, north_america[StateMask,])$NAME_1

    ############################################################################
    # Ecoregion overlap calculation 
    ############################################################################
    # Eliminate points that are far away
    centriodDist <- distHaversine(pt, ecoregion_centriods)
    cuttoffDistance <- sort(centriodDist)[7] # ecoregions have wierd shapes so keep many around
    SPDFMask <- centriodDist <= cuttoffDistance
    # Perform calculation on the remaining, close points
    ecoregion[xi,yi] <- over(pt, SPDF[SPDFMask,])$NA_L2CODE

    # Keep track of our progress
    count <- count + 1
    
    if(count %% 100 == 0){
      print(paste("Percent Complete:", count/nIterations*100))
    }
    
  }
  
}

t2 <- Sys.time()
dt <- t2 - t1
print("Since starting the loop:")
print(dt)

dataDir <- "Data/grid_attributes/"

# Sadly we cannot save entire state names in NC data. This will result in fill
# value getting written. So we need to make a unique numeric ID for each state
stateNames <- unique(as.character(state))
ORDER <- order(stateNames)
stateNamesOrdered <- stateNames[ORDER]
# create an array to store numeric state data
state_numeric <- matrix(0, nrow = dim(state)[1], ncol = dim(state)[2])

stateID <- 0
stateID_ <- numeric()
for (s in stateNamesOrdered){
  stateID_ <- append(stateID_, stateID)
  mask <- s == state
  state_numeric[mask] <- stateID
  stateID <- stateID + 1
}

# Save state dictionary as a csv
state_dict <- data.frame(stateID=stateID_, 
                         stateName=stateNamesOrdered)

write.csv(state_dict, file=paste0(dataDir, "state_dict_75.csv"), 
          row.names = FALSE)

# path and file name, set dname
ncpath <- dataDir
ncname <- "grid_attributes_75x75"  
ncfname <- paste(ncpath, ncname, ".nc", sep="")

# create and write the netCDF file -- ncdf4 version
# define dimensions
londim <- ncdim_def("longitude","degrees_east", as.double(longitude)) 
latdim <- ncdim_def("latitude","degrees_north", as.double(latitude)) 

# define variables
fillvalue <- 1e32
dlname <- "grid area meters squared"
grid_area.def <- ncvar_def("grid_area","m**2",list(londim,latdim),fillvalue,dlname,prec="float")

dlname <- "ecoregion"
ecoregion.def <- ncvar_def("ecoregion","none",list(londim,latdim),fillvalue,dlname,prec="float")

dlname <- "state"
state.def <- ncvar_def("state","none",list(londim,latdim),fillvalue,dlname, prec="float")

dlname <- "elevation"
elevation.def <- ncvar_def("elevation","m",list(londim,latdim),fillvalue,dlname,prec="float")


# create netCDF file and put arrays. All sperarate working. Together not. 
# Likes grid_area.df + elevation.def + ecoregion.def
ncout <- nc_create(ncfname, list(grid_area.def, 
                                 elevation.def, 
                                 ecoregion.def, 
                                 state.def), 
                   force_v4=TRUE)

# put variables
ncvar_put(ncout, grid_area.def, grid_area)
ncvar_put(ncout, ecoregion.def, ecoregion)
ncvar_put(ncout, state.def, state_numeric)
ncvar_put(ncout, elevation.def, elevation)


# put additional attributes into dimension and data variables
ncatt_put(ncout,"longitude","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"latitude","axis","Y")

# add global attributes
ncatt_put(ncout,0,"title","grid attributes for ecmwf grid. Created by createECMWFGridAttributeMasks.R")
ncatt_put(ncout,0,"institution","Colorado State University")
history <- paste("Steven J. Brey", date(), sep=", ")

# Get a summary of the created file:
ncout

# close the file, writing data to disk
nc_close(ncout)