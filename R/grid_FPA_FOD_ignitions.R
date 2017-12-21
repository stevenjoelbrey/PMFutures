#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args=2003
}
year <- args

# TODO: Rerun for early years.

print(paste("Working on year:", year))

# grid_FPA_FOD_ignitions.R
# Rscript --vanilla R/grid_FPA_FOD_ignitions.R 2003
# This script will create daily gridded data of FPA_FOD iginition counts. 

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
fgrid <- "/Volumes/Brey_external/GFED4s/GFED4.1s_burned_area_2004.nc"
nc <- nc_open(fgrid)
grid_area <- ncvar_get(nc, "grid_area")
longitude <- ncvar_get(nc, "longitude")
latitude  <- ncvar_get(nc, "latitude")
nc_close(nc)

t1 <- as.POSIXct(paste0(year,"-01-01"), tz="UTC")
t2 <- as.POSIXct(paste0(year,"-12-31"), tz="UTC")
time <- seq(t1,t2, by="day")

nTime <- length(time)

# Load the FPA_FOD_data
load(paste0("Data/FPA_FOD/FPA_FOD_",year,".RData"))

# three categories to hold on to. Human, unknown, and lightning
humanIgnition <- array(data=0, dim=c(nTime, length(longitude), length(latitude)))
unknownIgnition <- humanIgnition
lightningIgnition <- humanIgnition

# Get fire info 
fireLon <- FPA_FOD$LONGITUDE
fireLat <- FPA_FOD$LATITUDE
fireTime <- FPA_FOD$DISCOVERY_DATE
nFires <- dim(FPA_FOD)[1]

cause <- FPA_FOD$STAT_CAUSE_DESCR

causeSimple <- rep("", nFires)
lightningMask <- "Lightning" == cause
unkownMask <- "Missing/Undefined" == cause
humanMask <- !lightningMask & !unkownMask

causeSimple[lightningMask] <- "Lightning"
causeSimple[unkownMask] <- "Unknown"
causeSimple[humanMask] <- "Human"


for (i in 1:nFires){
  
  # Find this fires place in time 
  tMask <- fireTime[i] == time
  
  # Difference in array distances
  dx <- abs(longitude - fireLon[i])
  dy <- abs(latitude - fireLat[i])

  xi <- which.min(dx)
  yi <- which.max(dy)
  
  cause <- causeSimple[i]
  if(cause == "Lightning"){
    lightningIgnition[tMask, xi, yi] <- lightningIgnition[tMask, xi, yi] + 1 
  }else if(cause == "Unknown"){
    unknownIgnition[tMask, xi, yi] <- unknownIgnition[tMask, xi, yi] + 1 
  } else{
    humanIgnition[tMask, xi, yi] <- humanIgnition[tMask, xi, yi] + 1 
  }
  
  if(i %% 100 == 0){
    print(paste("Percent complete:", i/nFires*100))
  }
    
}

# Write these as RData, then combine into a large nc file later. The nc file 
# needs to match gfed grid attributes exactly. 
t0 <- as.POSIXct("1900-01-01 00:00:00", tz="UTC")
dt_hours <- (as.numeric(time) - as.numeric(t0)) / 60^2

# path and file name, set dname
ncpath <- "Data/FPA_FOD/"
ncpath <- "/Volumes/Brey_external/FPA_FOD/"
ncname <- paste0("ignition_counts_25x25_",year)
ncfname <- paste(ncpath, ncname, ".nc", sep="")

# create and write the netCDF file -- ncdf4 version
# define dimensions
londim <- ncdim_def("longitude","degrees_east", as.double(longitude)) 
latdim <- ncdim_def("latitude","degrees_north", as.double(latitude)) 
timedim <- ncdim_def("time", "hours since 1900-01-01 00:00:0.0", dt_hours)

# define variables
fillvalue <- 1e32
dlname <- "grid area meters squared"
grid_area.def <- ncvar_def("grid_area","m**2",list(londim,latdim),fillvalue,dlname,prec="float")

dlname <- "lightningIgnition"
lightningIgnition.def <- ncvar_def("lightningIgnition","count",list(timedim,londim,latdim),fillvalue,dlname,prec="float")

dlname <- "unknownIgnition"
unknownIgnition.def <- ncvar_def("unknownIgnition","none",list(timedim,londim,latdim),fillvalue,dlname, prec="float")

dlname <- "humanIgnition"
humanIgnition.def <- ncvar_def("humanIgnition","m",list(timedim,londim,latdim),fillvalue,dlname,prec="float")


# create netCDF file and put arrays. All sperarate working. Together not. 
# Likes grid_area.df + elevation.def + ecoregion.def
ncout <- nc_create(ncfname, list(grid_area.def, 
                                 lightningIgnition.def, 
                                 unknownIgnition.def, 
                                 humanIgnition.def), 
                   force_v4=TRUE)

# put variables
ncvar_put(ncout, grid_area.def, grid_area)
ncvar_put(ncout, lightningIgnition.def, lightningIgnition)
ncvar_put(ncout, unknownIgnition.def, unknownIgnition)
ncvar_put(ncout, humanIgnition.def, humanIgnition)


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