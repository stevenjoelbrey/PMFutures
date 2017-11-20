#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  year1 <- 2003
  year2 <- 2013
} else{
  year1 <- args[1]
  year2 <- args[2]
}

# grid_FPA_FOD_burn_area.R

###############################################################################
# ------------------------- Description ---------------------------------------
###############################################################################
# This script will be used to grid monthly burn area of FPA_FOD data. The grid 
# will match GFED4s.

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


# Get the quarter degree grid that all data lives on (GFED4s grid). All attribute
# layers will live on the same lat lon grid as this. 
fgrid <- "Data/GFED4s/GFED4.1s_burned_area_2004.nc"
nc <- nc_open(fgrid)
grid_area <- ncvar_get(nc, "grid_area")
longitude <- ncvar_get(nc, "longitude")
latitude  <- ncvar_get(nc, "latitude")
nc_close(nc)

# Create the time dimension
t1 <- as.POSIXct(paste0(year1,"-01-01"), tz="UTC")
t2 <- as.POSIXct(paste0(year2,"-12-31"), tz="UTC")
time <- seq(t1,t2, by="month")

# Get date and month arrays
timeLT <- as.POSIXlt(time)
gridMonth <- timeLT$mon + 1 
gridYear  <- timeLT$year + 1900

# Define the length of the monthly time dimension, [years X months]
nTime <- length(year1:year2) * length(1:12)

# Load the FPA_FOD_data
load(paste0("Data/FPA_FOD/FPA_FOD.RData"))

# Get fire info 
fireLon  <- FPA_FOD$LONGITUDE
fireLat  <- FPA_FOD$LATITUDE
fireTime <- FPA_FOD$DISCOVERY_DATE
fireStartMonth  <- as.POSIXlt(FPA_FOD$DISCOVERY_DATE)$mon + 1
fireStartYear <- FPA_FOD$FIRE_YEAR
containmentDate <-FPA_FOD$CONT_DATE
fireEndMonth    <- as.POSIXlt(FPA_FOD$CONT_DATE)$mon + 1
cause <- FPA_FOD$STAT_CAUSE_DESCR
nFires   <- dim(FPA_FOD)[1]

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


# Now handle the cause, which I will track for monthly gridded burn area
# product
causeSimple <- rep("", nFires)
lightningMask <- "Lightning" == cause
unkownMask <- "Missing/Undefined" == cause
humanMask <- !lightningMask & !unkownMask

# These are the only three categories we want to track 
causeSimple[lightningMask] <- "Lightning"
causeSimple[unkownMask] <- "Unknown"
causeSimple[humanMask] <- "Human"

# Create array to save monthly burn area at each location. 
human_burn_area <- array(data=0, dim=c(nTime, length(longitude), length(latitude)))
lightning_burn_aera <- human_burn_area
unknown_burn_area <- human_burn_area

################################################################################
# Loop through each fire, placing it in the correct grid (space & time)
################################################################################
metersSquaredPerAcre <- 4046.86
for (i in 1:nFires){# nFires
  
  # Get the burn area
  BA <- FPA_FOD$FIRE_SIZE[i] * metersSquaredPerAcre 
  
  # Find this fires place in time 
  tMask <- (fireStartMonth[i] == gridMonth) & (fireStartYear[i] == gridYear)
  
  # Difference in array distances
  dx <- abs(longitude - fireLon[i])
  dy <- abs(latitude - fireLat[i])
  xi <- which.min(dx)
  yi <- which.max(dy)
  
  # Place the burn area in the correct grid box based on ignition cause 
  cause <- causeSimple[i]
  if(cause == "Lightning"){
    lightning_burn_aera[tMask, xi, yi] <- lightning_burn_aera[tMask, xi, yi] + BA 
  } else if(cause == "Unknown"){
    unknown_burn_area[tMask, xi, yi] <- unknown_burn_area[tMask, xi, yi] + BA 
  } else if(cause == "Human"){
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
burn_area <- lightning_burn_aera + unknown_burn_area + human_burn_area

# The nc file needs a nice descriptive time dimension
t0 <- as.POSIXct("1900-01-01 00:00:00", tz="UTC")
t  <- as.POSIXct(paste(gridYear, gridMonth,"01",sep="-"), tz="UTC" ) 
dt_hours <- (as.numeric(t) - as.numeric(t0)) / 60^2

# Save easy to read alternative time dimension for clarity 
monthString <- as.character(gridMonth)
l_s <- str_length(monthString)
monthString[l_s==1] <- paste0("0", monthString[l_s==1])
YYYYMM <- as.numeric(paste0(as.character(gridYear), monthString))

# path and file name, set dname
ncpath <- "Data/FPA_FOD/"
#ncpath <- "/Volumes/Brey_external/FPA_FOD/"
ncname <- paste0("burn_area_monthly_25x25_", year1,"_",year2)
ncfname <- paste(ncpath, ncname, ".nc", sep="")

####################################################
# create and write the netCDF file -- ncdf4 version
####################################################
# define dimensions
londim <- ncdim_def("longitude","degrees_east", as.double(longitude)) 
latdim <- ncdim_def("latitude","degrees_north", as.double(latitude)) 
timedim <- ncdim_def("time", "hours since 1900-01-01 00:00:0.0", dt_hours)

# define variables
fillvalue <- 1e32
dlname <- "grid area meters squared"
grid_area.def <- ncvar_def("grid_area","m**2",list(londim,latdim),fillvalue,dlname,prec="float")

dlname <- "YYYYMM"
YYYYMM.def <- ncvar_def("YYYYMM","month",list(timedim),fillvalue,dlname,prec="float")

dlname <- "burn_area"
burn_area.def <- ncvar_def("burn_aera","m**2",list(timedim,londim,latdim),fillvalue,dlname,prec="float")

dlname <- "lightning_burn_aera"
lightning_burn_aera.def <- ncvar_def("lightning_burn_aera","m**2",list(timedim,londim,latdim),fillvalue,dlname,prec="float")

dlname <- "unknown_burn_area"
unknown_burn_area.def <- ncvar_def("unknown_burn_area","m**2",list(timedim,londim,latdim),fillvalue,dlname, prec="float")

dlname <- "human_burn_area"
human_burn_area.def <- ncvar_def("human_burn_area","m**2",list(timedim,londim,latdim),fillvalue,dlname,prec="float")


# create netCDF file and put arrays. All sperarate working. Together not. 
# Likes grid_area.df + elevation.def + ecoregion.def
ncout <- nc_create(ncfname, list(grid_area.def, 
                                 YYYYMM.def,
                                 burn_area.def, 
                                 lightning_burn_aera.def, 
                                 unknown_burn_area.def,
                                 human_burn_area.def), 
                   force_v4=TRUE)

# put variables
ncvar_put(ncout, grid_area.def, grid_area)
ncvar_put(ncout, YYYYMM.def, YYYYMM)
ncvar_put(ncout, burn_area.def, burn_area)
ncvar_put(ncout, lightning_burn_aera.def, lightning_burn_aera)
ncvar_put(ncout, unknown_burn_area.def, unknown_burn_area)
ncvar_put(ncout, human_burn_area.def, human_burn_area)

# put additional attributes into dimension and data variables
ncatt_put(ncout,"longitude","axis","X") #,verbose=FALSE) #,definemode=FALSE)
ncatt_put(ncout,"latitude","axis","Y")

# add global attributes
ncatt_put(ncout,0,"title","grid attributes for GFED4s grid. Created by grid_FPA_FOD_burn_area.R")
ncatt_put(ncout,0,"institution","Colorado State University")
history <- paste("Steven J. Brey", date(), sep=", ")

# Get a summary of the created file:
ncout

# close the file, writing data to disk
nc_close(ncout)