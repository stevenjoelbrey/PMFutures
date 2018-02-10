#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  year1 <- 1992
  year2 <- 2015
}else{
  year1 <- args[1]
  year2 <- args[2]
}

# TODO: Figure out why those 2 small fires in Wisconsin do not get gridMet vars
# TODO: assigned. 

# TODO: Assign temperature, wind speed, and RH using these data. 

#------------------------------- Description -----------------------------------
# This script takes the merged FPA FOD fire data created by 
# R/merge_ecoregion_year_Data.R and assigns gridMET data.
#
# Data description: http://www.climatologylab.org/gridmet.html
#
# DataURL: https://www.northwestknowledge.net/metdata/data/
#
# NOTE: Great care was taken to make sure assignment method is optimal in 
# NOTE: "R/assign_ecmwf_reanalysis_to_FPA_FOD.R"


library(stringr)
library(maps)
library(ncdf4)
library(fields)
library(geosphere)
library(lubridate) 

# Load the fire data (TODO: could also be version that has no ecmwf?)
load("Data/FPA_FOD/FPA_FOD_ecmwf_1992_2015.RData")

# NOTE: The spatial extent of these data does not cover Alaska, HI, PR, or 
# NOTE: anything eklse outside of conus, therefor, skip them in the loop that
# NOTE: makes assignments to fires. 
ncDir <- "/Volumes/Brey_external/gridMET/"

# Function for getting gridMet data
load_gridMET_nc <- function(ncDir, year){
  
  # TODO: Give this function name arguments so it works with other variables,
  # TODO: just just 1000-hr dead fuel moisture %
  ncFile <- paste0(ncDir, "fm1000_", year, ".nc")
  nc <- nc_open(ncFile)
  day <- ncvar_get(nc, "day") # units: "days since 1900-01-01 00:00:00", local or GMT?
  lat <- ncvar_get(nc, "lat")
  lon <- ncvar_get(nc, "lon")
  
  # dead_fuel_moisture_1000hr[lat,lon,day], 2.2 Gb in memory when loaded
  values <- ncvar_get(nc, "dead_fuel_moisture_1000hr") # units: "Percent"
  
  # Close the connection 
  nc_close(nc)
  
  # Now convert day to R friendly units
  dayToSeconds <- 24 * 60 * 60 # hr min sec
  deltaSeconds <- day * dayToSeconds
  t0 <- as.POSIXct("1900-01-01 00:00:00", tz="UTC")
  t <- t0 + deltaSeconds
  
  returnList <- list()
  returnList[["lat"]] <- lat
  returnList[["lon"]] <- lon
  returnList[["t"]] <- t
  returnList[["values"]] <- values
  
  return(returnList)
  
}

# Create an array to store all fuel moisture values 
fm1000_assigned <- rep(NA, dim(FPA_FOD)[1])

# Loop through each metGRID year of data. The files are big so it is best if
# we load them one at a time, this is easy to accomplish in a loop. 
for (year in year1:year2){

  print(paste("Making assignments for fires year:", year))
  
  # Open this years netCDF file to work with 
  l <- load_gridMET_nc(ncDir, year)
  lat <- l[["lat"]]
  lon <- l[["lon"]]
  t   <- l[["t"]]
  fm1000 <- l[["values"]]
  rm(l)
  
  # Get a subset the fire data to make assignments for 
  yearMask <- lubridate::year(FPA_FOD$DISCOVERY_DATE) == year 
  df <- FPA_FOD[yearMask, ]
  
  # Get the variables required for matching in space and time
  fireLat  <- df$LATITUDE
  fireLon  <- df$LONGITUDE
  fireTime <- df$DISCOVERY_DATE
  
  # Create array to store the years assigned values
  nFire <- dim(df)[1]
  fm1000_yearly <- rep(NA, nFire)
  for (i in 1:nFire){
    
    ti <- which(fireTime[i] == t)
    xi <- which.min(abs( fireLon[i] - lon  ))
    yi <- which.min(abs( fireLat[i] - lat ))
    
    fm1000_yearly[i] <- fm1000[yi, xi, ti]
  
  }
  
  # Place yearly values into main storage array based on the yearMask 
  fm1000_assigned[yearMask] <- fm1000_yearly

}

# make sure PR, HI, AK, excluded
STATE <- FPA_FOD$STATE
noMetGRID <- STATE == "AK" | STATE == "HI" | STATE == "PR"
fm1000_assigned[noMetGRID] <- NA

# Append the information to the dataframe
FPA_FOD$fuel_moisture_1000hr <- fm1000_assigned

# Save the appended dataframe
saveName <- paste0("Data/FPA_FOD/FPA_FOD_gridMet_",year1, "-", year2, ".RData")
save(FPA_FOD, file=saveName)

# plot a slize of these assigned data, for July, to see if they make sense  
library(ggplot2)

# Make fire location color depend on fuel_moisture_1000hr
pdf(file="fire_fm1000_mapped_julys.pdf", width=15, height=10)

df_plot <- FPA_FOD[FPA_FOD$START_MONTH==7 & !noMetGRID,]
ggplot(df_plot, aes(x=LONGITUDE, y=LATITUDE, color=fuel_moisture_1000hr)) + 
  geom_point(shape=".") + 
  theme_bw()+
  ggtitle("July 1000 hr fuel moisture %")

dev.off()


# Make fire location color depend on fuel_moisture_1000hr
pdf(file="Figures/fire_fm1000_mapped_january.pdf", width=15, height=10)

df_plot <- FPA_FOD[FPA_FOD$START_MONTH==1 & !noMetGRID,]
ggplot(df_plot, aes(x=LONGITUDE, y=LATITUDE, color=fuel_moisture_1000hr)) + 
  geom_point(shape=".") + 
  theme_bw()+
  ggtitle("January 1000 hr fuel moisture %")

dev.off()





