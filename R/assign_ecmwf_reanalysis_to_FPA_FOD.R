#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  year1 <- 1992
  year2 <- 2015
}else{
  year1 <- args[1]
  year2 <- args[2]
}

print("-----------------------------------------------------------------------")
print("This script assigns ecmwf reanalysis data to FPA-FOD fires.")
print(args)
print("-----------------------------------------------------------------------")

# TODO: wind gusts

# The base of the FPA_FOD fires that are getting weather assigned here are from
# the script:
#   readFPAFODFireFeatures.R

# execute via command line, this takes about 3 hours, e.g. below  
# Rscript --vanilla R/assign_ecmwf_reanalysis_to_FPA_FOD.R 1992 2015
################################################################################
# assign_ecmwf_reanalysis_to_FPA_FOD.R

# This script is assigns ecmwf daily weather attributes to FPA-FOD, fire 
# ingition conditions. 

# Assgin:
  # RH
  # T
  # precip, TODO: how much precip in season so far
  # days_since_rain_2003_2016.nc
  # dew point 

library(stringr)
library(maps)
library(ncdf4)
library(fields)
library(geosphere)

# Load FPA-FOD data, the one with all attributes (ecoregions) assgined.
FPAFOD_file <- paste0("Data/FPA_FOD/FPA_FOD_", year1, "_", year2, ".RData")
load(FPAFOD_file)
nRow     <- dim(FPA_FOD)[1]
fireDate <- FPA_FOD$DISCOVERY_DATE
fireLat  <- FPA_FOD$LATITUDE
fireLon  <- FPA_FOD$LONGITUDE

print("unique Months in the data START_MONTH column")
print(sort(unique(FPA_FOD$START_MONTH)))

print("unique years in the data FIRE_YEAR column")
print(sort(unique(FPA_FOD$FIRE_YEAR)))

# NOTE: This only works since all fires are treated as points, and only exist
# NOTE: in the western hemisphere. 
fireLonAdjusted <- fireLon + 360

# Set path to ecmwf era-interim data 
ncDir <- "/Volumes/Brey_external/era_interim_nc_daily_merged/"

print("Loading lots of large nc data")

nc_file <- paste0(ncDir,"t2m_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)

# Handle ecmwf time with origin 1900-01-01 00:00:0.0
ecmwf_hours <- ncvar_get(nc, "time")
ecmwf_seconds <- ecmwf_hours * 60^2

# make time useful unit
t0 <- as.POSIXct("1900-01-01 00:00:0.0", tz="UTC")
ecmwfDate <- t0 + ecmwf_seconds

# We only want to load through 2013
tf <- which(ecmwfDate == as.POSIXct(paste0(year2, "-12-31"), tz="UTC"))

# Now actually load the data
ecmwf_latitude <- ncvar_get(nc, "latitude")
ecmwf_longitude <- ncvar_get(nc, "longitude")
nLat <- length(ecmwf_latitude)
nLon <- length(ecmwf_longitude)

# Figure out the lon and lat subsets. These hard coded numbers are known limits
# to the spatial extent of North America when in lat lon [0 360] space. 
lon1 <- which(ecmwf_longitude == 180)
lon2 <- which(ecmwf_longitude == 310.5)
lonCount <- (lon2 - lon1) + 1

lat1 <- which(ecmwf_latitude == 87)
lat2 <- which(ecmwf_latitude == 15)
latCount <- (lat2 - lat1) + 1

# Get the spatial dims again using the new lonCount and latCount variables 
ecmwf_latitude  <- ncvar_get(nc, "latitude", start=lat1, count=latCount)
ecmwf_longitude <- ncvar_get(nc, "longitude", start=lon1, count=lonCount)

t2m <- ncvar_get(nc, "t2m", start=c(lon1,lat1, 1), count=c(lonCount, latCount, tf))

nc_close(nc)

# To keep things as clear as possible, subset the time array so that they ALL
# match in terms of dimensions. 
ecmwfDate <- ecmwfDate[1:tf]
if(length(ecmwfDate) != dim(t2m)[3]){
  stop("The ecmwfDate date array and variable array do not match in length")
}

################################################################################
# Plot a single time slice of t2m to make sure that you are getting the land
# area you think you are...
quartz()
f <- length(ecmwf_latitude):1
image.plot(ecmwf_longitude, ecmwf_latitude[f], t2m[,f,180])
title(paste("This map should show T over North America for", ecmwfDate[180]))

# RH% 
nc_file <- paste0(ncDir,"rh2m_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
rh2m <- ncvar_get(nc, "rh2m", start=c(lon1,lat1, 1), count=c(lonCount, latCount, tf))
nc_close(nc)

# days since rain 
nc_file <- paste0(ncDir,"days_since_rain_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
daysSinceRain <- ncvar_get(nc, "days_since_rain", start=c(lon1,lat1, 1), count=c(lonCount, latCount, tf))
nc_close(nc)

# total precip
nc_file <- paste0(ncDir,"tp_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
tp <- ncvar_get(nc, "tp", start=c(lon1,lat1, 1), count=c(lonCount, latCount, tf))
nc_close(nc)

# dew point
nc_file <- paste0(ncDir,"d2m_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
d2m <- ncvar_get(nc, "d2m", start=c(lon1,lat1, 1), count=c(lonCount, latCount, tf))
nc_close(nc)

# Evaporation
nc_file <- paste0(ncDir,"e_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
evap <- ncvar_get(nc, "e", start=c(lon1,lat1, 1), count=c(lonCount, latCount, tf))
nc_close(nc)

# Wind speed
nc_file <- paste0(ncDir,"u10_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
u10 <- ncvar_get(nc, "u10", start=c(lon1,lat1, 1), count=c(lonCount, latCount, tf))
nc_close(nc)

nc_file <- paste0(ncDir,"v10_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
v10 <- ncvar_get(nc, "v10", start=c(lon1,lat1, 1), count=c(lonCount, latCount, tf))
nc_close(nc)

windSpeed <- sqrt(u10^2 + v10^2)
# save workspace memory 
rm(v10, u10)

################################################################################
# Plot a time slice of wind. 
quartz()
f <- length(ecmwf_latitude):1
image.plot(ecmwf_longitude, ecmwf_latitude[f], windSpeed[,f,180])
title(paste("This map should show windspeed over North America for", ecmwfDate[180]))

# Get the nc grid attributes. Might as well get eleation from a higher resolution
# dataset 
# NOTE: The elevation grid being loaded is on a lon -180:180 grid!!! 
nc_file <- "Data/grid_attributes/grid_attributes_25x25.nc"
nc <- nc_open(nc_file)
elev  <- ncvar_get(nc, "elevation")
elev_lat <- ncvar_get(nc, "latitude")
elev_lon <- ncvar_get(nc, "longitude")
nc_close(nc)

quartz()
f <- length(elev_lat):1
image.plot(elev_lon, elev_lat[f], elev[,f])
title(paste("This is the levation to be assigned. GRID IS DIFFERENT"))

# Make sure the fireLon and elev_lon can be paired by ploting together. 
# Plot the grids and data together to make sure the grids are the same
# Sanity check the placement of everything by visualized ecmwf and fire locations
quartz(width=8, height=5)
elev_flipper <- length(elev_lat):1
image.plot(elev_lon, elev_lat[elev_flipper], elev[,elev_flipper])
points(fireLon, fireLat, pch=".", col="black")
title("The fire locations (black dots, should be over thge U.S. )")

# Make sure the fireLonAdjusted and ecmwf_longitude can be paired by ploting 
# together. 
quartz(width=8, height=5)
flipper                <- length(ecmwf_latitude):1
ecmwf_latitude_flipped <- ecmwf_latitude[flipper]
t2m_flipped            <- t2m[,flipper, 180] # and time sliced 

image.plot(ecmwf_longitude, ecmwf_latitude_flipped, t2m_flipped)
points(fireLonAdjusted, fireLat, pch=".", col="black")
title("The fire locations (black dots, should be over the U.S. )")

# NOTE: If the previous two images show fires and data agreeing of spatial coords
# NOTE: then we are good to go with these long coord pairs. 

# Get the assignment loop working, first just for temperature
v <- rep(NA, nRow)
t2m_assigned           <- v
rh2m_assigned          <- v
daysSinceRain_assigned <- v
tp_assigned            <- v
d2m_assigned           <- v
elev_assigned          <- v

# Long term predictive measures
tm2_lastMonth  <- v
rh2m_lastMonth <- v
tp_lastMonth   <- v
evap_lastMonth <- v

# After fire start days metric. Namely I am going to look at precip and wind
windSpeed_MonthAfter <- v
tp_MonthAfter        <- v
t2m_monthAfter       <- v

# Loop through every single fire and assign its past current and future weather
nECMWFTime <- length(ecmwfDate)
for (i in 1:nRow){ 
  
  # Find the fire match in space and time in the reanalysis data 
  xi <- which.min(abs(fireLonAdjusted[i] - ecmwf_longitude))
  yi <- which.min(abs(fireLat[i] - ecmwf_latitude))
  ti <- which(fireDate[i] == ecmwfDate)
  
  # Check the distance on the assigned grid points
  dist_meters <- distHaversine(c(fireLonAdjusted[i]-180, fireLat[i]), 
                               c(ecmwf_longitude[xi]-180, ecmwf_latitude[yi]))
  dist_km <- dist_meters/1000.0
  if(dist_km > 60){
    
    # Plot the error
    quartz()
    f <- length(ecmwf_latitude):1
    image.plot(ecmwf_longitude, ecmwf_latitude[f], t2m[,f,ti])
    points(fireLonAdjusted[i], fireLat[i], pch=3, cex=2, col="black")
    points(ecmwf_longitude[xi], ecmwf_latitude[yi])
    title("Fire Point and lat point are too far apart.")
    # Stop the code
    stop(paste("At row:", i, "distance=", dist_km, "km. Exceeds 60 km tolerance."))

  }
  
  # Assign each environmental variable that lives on the adjusted lon grid
  t2m_assigned[i] <- t2m[xi, yi, ti]
  rh2m_assigned[i] <- rh2m[xi, yi, ti]
  daysSinceRain_assigned[i] <- daysSinceRain[xi, yi, ti]
  tp_assigned[i] <- tp[xi, yi, ti]
  d2m_assigned[i] <- d2m[xi, yi, ti]
  
  # Assign environmental variables with a timescale greater than say 
  pastIndicies <- (ti-29):ti # length == 30
  if(pastIndicies[1] > 0){
    
    tm2_lastMonth[i] <- mean(t2m[xi, yi, pastIndicies])
    tp_lastMonth[i]  <- sum(tp[xi, yi, pastIndicies])
    evap_lastMonth[i]<- sum(evap[xi, yi, pastIndicies])
    rh2m_lastMonth[i]<- mean(rh2m[xi, yi, pastIndicies])
    
  }

  # After fire start date indicies 
  futureIndicies <- ti:(ti+29) # length == 30
  if(futureIndicies[30] <= nECMWFTime){
    
    windSpeed_MonthAfter[i] <- mean(windSpeed[xi, yi, futureIndicies])
    tp_MonthAfter[i]        <- sum(tp[xi, yi, futureIndicies])
    t2m_monthAfter[i]       <- mean(t2m[xi, yi, futureIndicies])
    
  }
  
  # Not ecmwf, but elevation too, no time dimension
  # NOTE: non adjusted lon
  xxi <- which.min(abs(fireLon[i] - elev_lon))
  yyi <- which.min(abs(fireLat[i] - elev_lat))
  
  # Check the distance on the assigned grid points
  dist_meters <- distHaversine(c(fireLon[i], fireLat[i]), 
                               c(elev_lon[xxi], elev_lat[yyi]))
  dist_km <- dist_meters/1000.0

  # Different grid different tolerance (27.82987 between centers). length of 
  # hypotenuse is then 19.67869 km 
  if(dist_km > 20){
    stop(paste("Elevation assignent for i =", i, "is not within error tolerance."))
  }
  
  elev_assigned[i] <- elev[xxi, yyi]
  
  # Output progress to the screen
  if(i %% 1000 == 0){
    print(paste("Percent Complete: ", i/nRow*100))
  }
  
}

# Assign each of the new environmental arrays to FPA_FOD
FPA_FOD$t2m             <- t2m_assigned
FPA_FOD$rh2m            <- rh2m_assigned
FPA_FOD$days_since_rain <- daysSinceRain_assigned
FPA_FOD$tp              <- tp_assigned
FPA_FOD$d2m             <- d2m_assigned

# Past environment
FPA_FOD$tm2_lastMonth  <- tm2_lastMonth
FPA_FOD$tp_lastMonth   <- tp_lastMonth
FPA_FOD$rh2m_lastMonth <- rh2m_lastMonth
FPA_FOD$evap_lastMonth <- evap_lastMonth

# Future environment
FPA_FOD$windSpeed_MonthAfter <- windSpeed_MonthAfter
FPA_FOD$tp_MonthAfter        <- tp_MonthAfter
FPA_FOD$t2m_monthAfter       <- t2m_monthAfter

# static in time 
FPA_FOD$elevation <- elev_assigned

save(FPA_FOD, file = paste0("Data/FPA_FOD/FPA_FOD_ecmwf_",year1,"_",year2,"_new.RData"))
# The end. 