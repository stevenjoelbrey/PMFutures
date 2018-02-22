#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  year1 <- 1992
  year2 <- 2015
  grid  <- "ecmwf" # "gfed4s"
  sanityCheck <- FALSE
} else{
  year1 <- args[1]
  year2 <- args[2]
  grid  <- args[3]
  sanityCheck <- as.logical(args[4])
}

# grid_FPA_FOD_burn_area.R
# e.g. for running from command line
# Rscript --vanilla R/grid_FPA_FOD_burn_area.R 1992 2015 ecmwf FALSE

# ------------------------- Description ---------------------------------------
# This script will be used to grid monthly burn area of FPA_FOD data. The grid 
# will match GFED4s or ecmwf grid

# NOTE: Short et. al. says that fire size are stored in acres. 

print(paste("Working on year:", year1, "to year", year2))

library(stringr)
library(maps)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)
library(raster)
library(ncdf4)
library(geosphere)
library(lubridate) # For handy "month()" and "year()" functions
library(fields)

if(grid=="gfed4s"){
  
  # Get the quarter degree grid that all data lives on (GFED4s grid). All attribute
  # layers will live on the same lat lon grid as this. 
  fgrid <- "Data/GFED4s/GFED4.1s_burned_area_2004.nc"
  nc <- nc_open(fgrid)
  grid_area <- ncvar_get(nc, "grid_area")
  longitude <- ncvar_get(nc, "longitude")
  latitude  <- ncvar_get(nc, "latitude")
  nc_close(nc)
  
} else if(grid=="ecmwf"){
  
  fgrid <- "/Volumes/Brey_external/era_interim_nc_daily/t2m_1992.nc"
  nc <- nc_open(fgrid)
  longitude <- ncvar_get(nc, "longitude")
  latitude  <- ncvar_get(nc, "latitude")
  nc_close(nc)
  
  # Grid area needs to come from attribute file
  nc <- nc_open("Data/grid_attributes/grid_attributes_75x75.nc")
  grid_area <- ncvar_get(nc, "grid_area")
  nc_close(nc)
  
}

# Create the monthly time dimension, represented by a single date of that month
t1 <- as.POSIXct(paste0(year1,"-01-01"), tz="UTC")
t2 <- as.POSIXct(paste0(year2,"-12-31"), tz="UTC")
time <- seq(t1,t2, by="month")

# Get date and month arrays
timeLT <- as.POSIXlt(time)
gridMonth <- month(time) 
gridYear  <- year(time)

# Define the length of the monthly time dimension, [years X months]
nTime <- length(year1:year2) * length(1:12)

# Load the "FPA_FOD" dataframe. These data have all desired appended information
load(paste0("Data/FPA_FOD/FPA_FOD_1992_2015.RData"))

# Get fire info 
fireLon  <- FPA_FOD$LONGITUDE
fireLat  <- FPA_FOD$LATITUDE
fireTime <- FPA_FOD$DISCOVERY_DATE
fireStartMonth  <- FPA_FOD$START_MONTH
fireStartYear   <- FPA_FOD$FIRE_YEAR # came prepackaged with original FPA-FOD
containmentDate <-FPA_FOD$CONT_DATE
fireEndMonth    <- month(as.POSIXlt(FPA_FOD$CONT_DATE))
cause <- FPA_FOD$STAT_CAUSE_DESCR
nFires   <- dim(FPA_FOD)[1]

# When applying the WESTERN HEMISPHERE ONLY fire locations to ecmwf grid, add
# 360 to the fire longitude values
if(grid=="ecmwf"){
  fireLon <- fireLon + 360
}

# ------------------------- fireStats ------------------------------------------
# For fires that have containment dates, what % end in the same month they 
# start? This is an important piece of information for describing monthly burn
# area product you are about to make. 
noContainmentDate <- is.na(containmentDate)
deltaMonthContain <- (fireEndMonth - fireStartMonth)[!noContainmentDate]

quartz()
h <- hist(abs(deltaMonthContain), breaks=0:12, plot=FALSE)
plot(h$counts, log="y", type='h', lwd=10, lend=2, xlim=c(0, 11))

percentSameMonth <- sum(deltaMonthContain==0, na.rm=T) / sum(!noContainmentDate)*100
print(paste("Percent of fires with containment dates that end same as start month:", percentSameMonth))

# Of the fires that occur in one month, what % of burn area do these account for? 
totalBurnArea <- sum(FPA_FOD$FIRE_SIZE[!noContainmentDate])
sameMonthBurnArea <- sum(FPA_FOD$FIRE_SIZE[!noContainmentDate][deltaMonthContain==0])
fSame <- sameMonthBurnArea / totalBurnArea * 100
print(paste("Percent of burn area from fires that have same month containment dates:", fSame))

# ------------------------- fireCause sorting ----------------------------------
# Now handle the cause, which I will track for monthly gridded burn area
# product

causeSimple   <- rep("", nFires)
lightningMask <- "Lightning" == cause
unkownMask    <- "Missing/Undefined" == cause
humanMask     <- !lightningMask & !unkownMask

# These are the only three categories we want to track 
causeSimple[lightningMask] <- "Lightning"
causeSimple[unkownMask] <- "Unknown"
causeSimple[humanMask] <- "Human"

# Create array to save monthly burn area at each location. 
nLon <- length(longitude)
nLat <- length(latitude)
human_burn_area     <- array(data=0, dim=c(nTime, nLon, nLat))
lightning_burn_area <- array(data=0, dim=c(nTime, nLon, nLat))
unknown_burn_area   <- array(data=0, dim=c(nTime, nLon, nLat))

# ------------------------- Loop -----------------------------------------------
# Loop through each fire, placing it in the correct grid (space & time)

metersSquaredPerAcre <- 4046.86
df_burn_area_sum <- sum(FPA_FOD$FIRE_SIZE * 4046.86)

for (i in 1:nFires){# nFires
  
  # Get the burn area for this fire
  BA <- FPA_FOD$FIRE_SIZE[i] * metersSquaredPerAcre 
  
  # Find this fires place in time, i.e., correct year and month 
  tMask <- (fireStartMonth[i] == gridMonth) & (fireStartYear[i] == gridYear)
  
  # Difference in array distances
  dx <- abs(longitude - fireLon[i])
  dy <- abs(latitude - fireLat[i])
  xi <- which.min(dx)
  yi <- which.min(dy)
  
  # Make sure distance is not too large
  if(grid=="ecmwf"){
    d <- distHaversine( c(fireLon[i]-180, fireLat[i]), c(longitude[xi]-180, latitude[yi]) )/1000
    if(d > 60){
      stop(paste("60 km tolerance from grid point center for gridding exceeded at i =", i))
    }
  }
  
  # Map where I am gridding the point when I am checking in development mode
  if(sanityCheck){
    plot(longitude[xi], latitude[yi], col="black", pch=19)
    points(fireLon[i], fireLat[i], pch=3, col="red")
    map("state", add=T)
  }
  
  # Place the burn area in the correct grid box based on ignition cause 
  thisCause <- causeSimple[i]
  if(thisCause == "Lightning"){
    lightning_burn_area[tMask, xi, yi] <- lightning_burn_area[tMask, xi, yi] + BA 
  } else if(thisCause == "Unknown"){
    unknown_burn_area[tMask, xi, yi] <- unknown_burn_area[tMask, xi, yi] + BA 
  } else if(thisCause == "Human"){
    human_burn_area[tMask, xi, yi] <- human_burn_area[tMask, xi, yi] + BA
  } else{
    stop("Ignition category of this fire is not known. Find cause and rewrite.")
  }
  
  # Print progress to console
  if(i %% 100 == 0){
    print(paste("Percent complete:", i/nFires*100))
  }
  
}

# Add these three types together for fast access to a total burn_area variable.
burn_area <- lightning_burn_area + unknown_burn_area + human_burn_area

# Make a sanity plot to investigate
quartz()
f <- length(latitude):1
burn_area_total <- apply(burn_area, 2:3, sum)
# do not plot where equal zero
burn_area_total[burn_area_total==0] <- NA
image.plot(longitude, latitude[f], burn_area_total[, f])
title("Total burn area as gridded by this script")

# Now see if the total burn area was conserved or not
gridded_burn_area_sum <- sum(burn_area)
dArea <- gridded_burn_area_sum - df_burn_area_sum
print(paste("The difference in meters squared before and after gridding is:", dArea))

# The nc file needs a nice descriptive time dimension
t0 <- as.POSIXct("1900-01-01 00:00:00", tz="UTC")
t  <- as.POSIXct(paste(gridYear, gridMonth,"01", sep="-"), tz="UTC" ) 
dt_hours <- (as.numeric(t) - as.numeric(t0)) / 60^2

print("t array was created.")

# Save easy to read alternative time dimension for clarity 
monthString <- as.character(gridMonth)
l_s <- str_length(monthString)
monthString[l_s==1] <- paste0("0", monthString[l_s==1])
YYYYMM <- as.numeric(paste0(as.character(gridYear), monthString))

print("YYYYMM array was created.")


# path and file name, set dname
ncpath <- "Data/FPA_FOD/"
if(grid=="gfed"){
  ncname <- paste0("burn_area_monthly_25x25_", year1,"_",year2)
} else if(grid=="ecmwf"){
  ncname <- paste0("burn_area_monthly_75x75_", year1,"_",year2)
}
ncfname <- paste(ncpath, ncname, ".nc", sep="")

print("ncfname array was created.")


####################################################
# create and write the netCDF file -- ncdf4 version
####################################################
# define dimensions
londim <- ncdim_def("longitude","degrees_east", as.double(longitude)) 
latdim <- ncdim_def("latitude","degrees_north", as.double(latitude)) 
timedim <- ncdim_def("time", "hours since 1900-01-01 00:00:00", dt_hours, TRUE) # unlimited == TRUE

# define variables
fillvalue <- 1e32
dlname <- "grid area meters squared"
grid_area.def <- ncvar_def("grid_area","m**2",list(londim,latdim),fillvalue,dlname,prec="float")

dlname <- "YYYYMM"
YYYYMM.def <- ncvar_def("YYYYMM","month",list(timedim),fillvalue,dlname,prec="float")

dimList <- list(timedim,londim,latdim)
dlname <- "burn_area"
burn_area.def <- ncvar_def("burn_area","m**2",dimList,fillvalue,dlname,prec="float")

dlname <- "lightning_burn_area"
lightning_burn_area.def <- ncvar_def("lightning_burn_area","m**2",dimList,fillvalue,dlname,prec="float")

dlname <- "unknown_burn_area"
unknown_burn_area.def <- ncvar_def("unknown_burn_area","m**2",dimList,fillvalue,dlname, prec="float")

dlname <- "human_burn_area"
human_burn_area.def <- ncvar_def("human_burn_area","m**2",dimList,fillvalue,dlname,prec="float")

print("vars created.")


# create netCDF file and put arrays. All sperarate working. Together not. 
# Likes grid_area.df + elevation.def + ecoregion.def
ncout <- nc_create(ncfname, list(grid_area.def, 
                                 YYYYMM.def,
                                 burn_area.def, 
                                 lightning_burn_area.def, 
                                 unknown_burn_area.def,
                                 human_burn_area.def), 
                   force_v4=TRUE)

print("ncout created.")


# put variables
ncvar_put(ncout, grid_area.def, grid_area)
ncvar_put(ncout, YYYYMM.def, YYYYMM)
ncvar_put(ncout, burn_area.def, burn_area)
ncvar_put(ncout, lightning_burn_area.def, lightning_burn_area)
ncvar_put(ncout, unknown_burn_area.def, unknown_burn_area)
ncvar_put(ncout, human_burn_area.def, human_burn_area)

print("ncout variables put")

# put additional attributes into dimension and data variables
ncatt_put(ncout,"longitude","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"latitude","axis","Y")

# add global attributes
ncatt_put(ncout,0,"title",paste0("grid attributes for", grid," grid. Created by grid_FPA_FOD_burn_area.R"))
ncatt_put(ncout,0,"institution","Colorado State University")
history <- paste("Steven J. Brey", today(), sep=", ")

# Get a summary of the created file:
ncout

# close the file, writing data to disk
nc_close(ncout)