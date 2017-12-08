# createGridAttributeMasks.R

# This script will be used to mask the 0.25 and 0.75 degree grids used in this
# analysis. Each grid will have a layer with the following descriptive attributes
# 1) Area, these will be from GFED4s grids
# 2) state/province
# 3) ecoregion
# 4) elevation 

# TODO: Make sure no land ecoregions are left as "0" meaning water label


library(stringr)
library(maps)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(raster)
library(ncdf4)
library(geosphere) # for calculating distance of point to polygon edge

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

# Load the elevation data
felev <- "Data/GIS/elev.0.25-deg.nc"
nc <- nc_open(felev)
elevation <- ncvar_get(nc, "data")
elevationUnits <- ncatt_get(nc,"data","units")
elev_longitude <- ncvar_get(nc, "lon")
elev_latitude  <- ncvar_get(nc, "lat")
nc_close(nc)

# These elevation data need to be shifted in the longitude direction to match 
# the GFED grid.
elevation_new <- elevation
elev_longitude_new <- longitude
elevation_new[] <- NA
elev_longitude_new[] <- NA
# Start shifting, one half of the world at a time
elevation_new[1:720,] <- elevation[721:1440,]
elevation_new[721:1440,] <- elevation[1:720,]
# Do the same for the longitude array, provides direct and easy check
elev_longitude_new[1:720] <- longitude[721:1440]
elev_longitude_new[721:1440] <- longitude[1:720]

# Make sure all values have been replaced. 
if(sum(is.na(elevation_new)) > 0 ){
  warning("There is an NA value present in the new elevation data")
  # where? 
  #elevation_new[is.na(elevation_new)] <- 1e6
  quartz()
  image(is.na(elevation_new))
  title("Locations where new elevation data is still NA")
}

# # Sanity check figures
# quartz()
# image(elevation)
# title("original elevation nc data")
# 
# quartz()
# image(elevation_new)
# title("shifted elevation nc data")
elevation <- elevation_new

unique(elev_latitude - latitude)
# I do not care about ozone. Anywhere less than zero will be set to 0
# elevation[elevation<0] <- 0 

# NOTE: setting this to NA makes it really easy to see land
# NOTE: I am setting all below 

# point over polygon calculations. Now we need to assign ecoregions as well as
# political borders. 
# Speed this up by only looping through North America coordinates 
lon_index <- which(longitude >= -168 & longitude <= -52)
lat_index <- which(latitude >= 25 & latitude <= 80)

# quartz()
# image(elevation[lon_index, lat_index])
# title("Spatial subset where we will be making overlap calculations")


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
  
  gridLon <- longitude[xi]
  
  for (yi in lat_index){
    
    gridLat <- latitude[yi]
    
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
    #state[xi, yi] <- over(pt, north_america)$NAME_1 # SLOW
    
    ############################################################################
    # Ecoregion overlap calculation 
    ############################################################################
    # Eliminate points that are far away
    centriodDist <- distHaversine(pt, ecoregion_centriods)
    cuttoffDistance <- sort(centriodDist)[25] # ecoregions have wierd shapes so keep many around
    SPDFMask <- centriodDist <= cuttoffDistance
    # Perform calculation on the remaining, close points
    ecoregion[xi,yi] <- over(pt, SPDF[SPDFMask,])$NA_L2CODE
    
    # # If a state has been assigned but an ecoregion has not, then that means
    # # we have an on the border type situation
    # if(is.na(ecoregion[xi,yi]) & !is.na(state[xi, yi]) ){
    #   print("No ecoregion selected, requires expensive calculations")
    #   dist.mat <- geosphere::dist2Line(p = pt, line = SPDF)
    #   ecoregion[xi,yi] <- SPDF[dist.mat[4],]$NA_L2CODE
    # }
    
    #ecoregion[xi,yi] <- over(pt, SPDF)$NA_L2CODE # VERY SLOW
    
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


################################################################################
# Write out the nc data
################################################################################

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
  print(s)
  stateID_ <- append(stateID_, stateID)
  mask <- s == state
  state_numeric[mask] <- stateID
  stateID <- stateID + 1
}

# Save state dictionary as a csv
state_dict <- data.frame(stateID=stateID_, 
                         stateName=stateNamesOrdered)

write.csv(state_dict, file=paste0(dataDir, "state_dict.csv"), 
          row.names = FALSE)

# path and file name, set dname
ncpath <- dataDir
ncname <- "grid_attributes_25x25"  
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

# Same dimensions as variables too
latitude.def <- ncvar_def("lat","m", list(latdim), fillvalue,"lat", prec="float")
longitude.def <- ncvar_def("lon","m", list(londim), fillvalue,"lon", prec="float")

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
# ncvar_put(ncout, latitude.def, latitude)
# ncvar_put(ncout, longitude.def, longitude)


# put additional attributes into dimension and data variables
ncatt_put(ncout,"longitude","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"latitude","axis","Y")

# add global attributes
ncatt_put(ncout,0,"title","grid attributes for GFED4s grid. Create by createGridAttributeMasks.R")
ncatt_put(ncout,0,"institution","Colorado State University")
history <- paste("Steven J. Brey", date(), sep=", ")

# Get a summary of the created file:
ncout

# close the file, writing data to disk
nc_close(ncout)



