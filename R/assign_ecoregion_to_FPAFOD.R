#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args=2003
}

# readFPAFODFireFeatures.R
# execute via command line 
# Rscript --vanilla R/assign_ecoregion_to_FPAFOD.R 2003
################################################################################
# readUSFSFireOccurance.R

# This script is used to convert spatial wildfire occurance data into a R 
# friendly easy to load dataframe and csv. This script will laod a text file of this
# data. In truth, this script hides some of the ugly work that had to be done
# to get this data into a workable format. I downloaded a .dbf file from the 
# link listed below, I then had to convert that to a text file using ArcCatalog
# and ArcMap. ESRI really does not want anyone to look at this data without 
# thier software...

# DataSource: https://www.fs.usda.gov/rds/archive/Product/RDS-2013-0009.3/

# Data Citation:
# Short, Karen C. 2015. Spatial wildfire occurrence data for the United States, 
# 1992-2013 [FPA_FOD_20150323]. 3rd Edition. Fort Collins, CO: Forest Service 
# Research Data Archive. https://doi.org/10.2737/RDS-2013-0009.3
#
#
# Excellent overview of these data available:
# http://geog.uoregon.edu/bartlein/FireStarts/fpa-fod_RODBC_01.html
################################################################################

library(stringr)
library(maps)
library(sp)
library(rgdal)
library(maptools)

year <- as.numeric(args)


print("-----------------------------------------------------------------------")
print("Reading in the really big fire .RData file")
print(paste("Working on year: ", year))
print("-----------------------------------------------------------------------")

drive <- "Data/FPA_FOD/"
dataFile <- paste0(drive,"fire_featureclass.RData")
df <- get(load(dataFile))

# Temporal subset of the FPA-FOD
timeMask <- df$FIRE_YEAR == year
df <- df[timeMask,]

# Names of the df
# [1] "OBJECTID"                   "FOD_ID"                     "FPA_ID"                     "SOURCE_SYSTEM_TYPE"        
# [5] "SOURCE_SYSTEM"              "NWCG_REPORTING_AGENCY"      "NWCG_REPORTING_UNIT_ID"     "NWCG_REPORTING_UNIT_NAME"  
# [9] "SOURCE_REPORTING_UNIT"      "SOURCE_REPORTING_UNIT_NAME" "LOCAL_FIRE_REPORT_ID"       "LOCAL_INCIDENT_ID"         
# [13] "FIRE_CODE"                  "FIRE_NAME"                  "ICS_209_INCIDENT_NUMBER"    "ICS_209_NAME"              
# [17] "MTBS_ID"                    "MTBS_FIRE_NAME"             "COMPLEX_NAME"               "FIRE_YEAR"                 
# [21] "DISCOVERY_DATE"             "DISCOVERY_DOY"              "DISCOVERY_TIME"             "STAT_CAUSE_CODE"           
# [25] "STAT_CAUSE_DESCR"           "CONT_DATE"                  "CONT_DOY"                   "CONT_TIME"                 
# [29] "FIRE_SIZE"                  "FIRE_SIZE_CLASS"            "LATITUDE"                   "LONGITUDE"                 
# [33] "OWNER_CODE"                 "OWNER_DESCR"                "STATE"                      "COUNTY"                    
# [37] "FIPS_CODE"                  "FIPS_NAME"                 


# Get start time, containment time, and out time into POSIXct format
disc_mmddyy_TIME    <- df$DISCOVERY_DATE
spaceLocation  <- str_locate(disc_mmddyy_TIME, " ")

# The time of these incidents is not needed. Just the date. These data are not 
# that precise. 
disc_mmddyy    <- str_sub(disc_mmddyy_TIME, 1, spaceLocation[,1])

# Get rid of the space that sometimes shows up at the end 
disc_mmddyy_noSpace    <- str_replace(disc_mmddyy, " ", "")
discoverd_date <- as.POSIXct(disc_mmddyy_noSpace, format="%m/%d/%y", tz="UTC")

# Add to dataframe
df$DISCOVERY_DATE <- discoverd_date

# Repeat this method for the containment date
con_mmddyy_TIME <- df$CONT_DATE
spaceLocation <- str_locate(con_mmddyy_TIME, " ")
con_mmddyy <- str_sub(con_mmddyy_TIME, 1, spaceLocation[,1])
con_mmddyy_noSpace  <- str_replace(con_mmddyy, " ", "")
con_date <- as.POSIXct(con_mmddyy_noSpace, format="%m/%d/%y", tz="UTC")

df$CONT_DATE <- con_date

# Make an array of the month of the start date 
t_LT <- as.POSIXlt(discoverd_date)
month <- t_LT$mon + 1

# Add it to the dataframe
df$START_MONTH <- month

# Now we need to assign the eco-region of each individual fire. This maintains
# highest level of precision before gridding to 0.25 degree grid. 
lat <- df$LATITUDE
lon <- df$LONGITUDE
coordinates <- cbind(lon, lat)

fireLocations <- SpatialPoints(coords=coordinates, proj4string=CRS("+proj=longlat"))

# Load the level II ecoregions of the U.S. 
# NOTE: These borders have been merged and simplified by R/simplyfyECORegions.R
layername <- "na_cec_eco_level_2.RData"
layerDir  <- "Data/GIS/na_cec_eco_l2/"
SPDF <- get(load(paste0(layerDir, layername)))

# The class of the NA_L2CODE column is wierd, making better
SPDF@data$NA_L2CODE <- as.numeric(as.character(SPDF@data$NA_L2CODE))

# We want to assign the NA_L2CODE and NA_L2NAME to each fire (row) in the fpafod
# This will require overlap analysis. 
nRow <- dim(df)[1]
NA_L2CODE <- rep(NA, nRow)

# Loop through each fire and assign the level II ecoregion
print("-----------------------------------------------------------------------")
print(paste("Beginning work on the large for loop, n =",nRow,"iterations"))
print("-----------------------------------------------------------------------")

t1 <- Sys.time()
for (i in 1:nRow){

  # Get this index fire location
  p <- fireLocations[i,]

  # Figure out which polygon this point falls inside of
  p_ecoregion <- over(p, SPDF, returnList = FALSE)
  
  # Assign level 2 identification
  NA_L2CODE[i] <- p_ecoregion$NA_L2CODE

  # Output progress to the screen
  if(i %% 1000 == 0){
    print(paste("Percent Complete: ", i/nRow*100))
  }
    
}
tf <- Sys.time()

print(paste("It took", tf - t1, "to complete loop"))

# Assign all overlap analysis information to the df and save it! 
df$NA_L2CODE <- NA_L2CODE

# Write out this data as a csv for easy viewing and sharing
# write.csv(df, file=paste0(drive, "FPA_FOD.csv"),
#           row.names = FALSE)
FPA_FOD <- df
save(FPA_FOD, file=paste0(drive, "FPA_FOD_", year, ".RData"))

print(paste("Complete:", year))


# TODO: Source border assignment script? 


# # What percent of the data have DISCOVERY_DATE
# nRecords <- length(con_date)
# print(paste(sum(is.na(discoverd_date))/nRecords * 100, "% missing end date"))
# print(paste(sum(is.na(con_date))/nRecords * 100, "% missing end date"))

# 
# 
# # Lets only look at this relationship in the western U.S. Duration relationships
# # explored by region. 
# latMask <- df$LATITUDE >= 31 & df$LATITUDE <= 49
# lonMask <- df$LONGITUDE >= -125 & df$LONGITUDE <= -104
# 
# df_subset <- df[latMask & lonMask & monthMask, ]
# 
# # Now handle lightning and human started fires differently
# lightningMask <- df_subset$STAT_CAUSE_DESCR == "Lightning"
# df_lightning <- df_subset[lightningMask,]


# sizeDuration <- function(df){
# 
#   # To make daily estimates of burn area, we need a way to spread burn area 
#   # between the start and end date of each incident. For the incidents that have
#   # start and end dates, is there a relationship between size and duration, such 
#   # that we could estimate duration from size? 
#   secondsInDay <- 60^2*24
#   fireDuration <- (as.numeric(df$CONT_DATE - df$DISCOVERY_DATE) + secondsInDay) / secondsInDay
#   fireSize <- df$FIRE_SIZE
#   
#   # Where do the fires that do not have end dates tend to occur? Pattern? 
#   missingContainment <- is.na(df$CONT_DATE)
#   
#   # TODO: consider using leaflet interactive map here. 
#   map("state", xlim=c(-125, -104))
#   points(df$LONGITUDE[missingContainment],  df$LATITUDE[missingContainment], pch=".")
#   title("Location of fires with no containment date")
#   
#   plot(fireSize, fireDuration, pch=".")
#   linearModel <- lm(fireDuration ~ fireSize)
#   title("Fire Duration as a function of fire Size")
#   abline(linearModel)
#   summary(linearModel)
#   
#   # # How big are fires that last only one day?
#   # hist(fireSize[fireDuration <=1], log="x")
#   # 
#   # plot(fireSize[fireDuration <=1], fireDuration[fireDuration <=1])
#   
#   bigFireMask <- fireSize > 1000
#   bigFireSize <- fireSize[bigFireMask]
#   bigFireDuration <- fireDuration[bigFireMask]
#   
#   plot(bigFireSize, bigFireDuration)
#   title("fire size vs. fire duration fires > 1000 acres")
#   
#   bigFireMask <- fireSize > 10000
#   bigFireSize <- fireSize[bigFireMask]
#   bigFireDuration <- fireDuration[bigFireMask]
#   
#   plot(bigFireSize, bigFireDuration)
#   title("fire size vs. fire duration fires > 10000 acres")
#   
#   
#   # 
#   # # What is the size distribution for fires that do not have containment dates? 
#   # noEndDateMask <- is.na(con_date)
#   # df_subset <- df[noEndDateMask,]
#   # 
#   # hist(log(df_subset$FIRE_SIZE), main="")
#   # title(paste("Mean size of fires with no end date:", mean(df_subset$FIRE_SIZE)))
#   # 
#   # # These are very small fires...
#   # # Subset again, if the fire size is less than 10 acres we assume that the fire 
#   # # only lasts one day.
#   # mask <- noEndDateMask & (df$FIRE_SIZE > 10)
#   # df_subset <- df[mask,]
#   # 
#   # hist(df_subset$FIRE_SIZE)
# 
# }
# 
# sizeDuration(df_lightning)


# ################################################################################
# # apply observed relationship, by region, and season, for fire size and lendth
# ################################################################################
# 
# # Subset the data by the rows of interest
# columnsToKeep <- c("FPA_ID", "FIRE_NAME","COMPLEX_NAME","DISCOVERY_DATE",
#                    "DISCOVERY_TIME",
#                    "STAT_CAUSE_DESCR", "CONT_DATE","CONT_TIME","FIRE_SIZE",
#                    "LATITUDE", "LONGITUDE","OWNER_CODE","FIPS_NAME")
# 
# # Where we do not have a containment data, leave as NA, no length calculation
# # will be possible
# # Get rid of any rows that do not have good date information
# # rowMask <- !is.na(con_date) & !is.na(discoverd_date)
# columnMask <- names(df) %in% columnsToKeep
# df_subset <- df[,columnMask]
# 
# print("-----------------------------------------------------------------------")
# print("Writing the really big (but subset for dired columns) fire csv")
# print("-----------------------------------------------------------------------")
# 
# # Add a color column, one color for human started, another for natural
# assign("fire_data", df_subset)
# save(fire_data, file=paste0("data/fireOccurrence.RData"))
# 
# # Write the data as csv as well for easy use for pandas.
# 
# write.csv(fire_data, file=paste0("data/fireOccurrence.csv"), row.names=FALSE)

