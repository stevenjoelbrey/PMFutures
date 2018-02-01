#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args=2003
}

# assign_ecoregion_to_FPAFOD.R 
# execute via command line 
# cd GoogleDrive/sharedProjects/PMFutures/
# Rscript --vanilla R/assign_ecoregion_to_FPAFOD.R 2003

#------------------------------- Desscription ---------------------------------_

# The purpose of this script is to assign ecoregions to FPA FOD fires. 
# It should be run after all FPA data have been formatted and times are stored as
# POSIXct objects. That work is performed by readUSFSFireOccurance.R


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
# #dataFile <- paste0(drive,"fire_featureclass.RData")

# UPDATED SEE read2017FPAFOD.R for details
dataFile <- paste0(drive, "FPA_FOD_1992_2015.RData") 
df <- get(load(dataFile))
rm(FPA_FOD)

# Temporal subset of the FPA-FOD so that many of these can run at once. This 
# is only done for speed, so different years can run at the same time. 
timeMask <- df$FIRE_YEAR == year
df <- df[timeMask,]

# Names of the df
# [1] "OBJECTID"                   "FOD_ID"                     "FPA_ID"                     "SOURCE_SYSTEM_TYPE"        
# [5] "SOURCE_SYSTEM"              "NWCG_REPORTING_AGENCY"      "NWCG_REPORTING_UNIT_ID"     "NWCG_REPORTING_UNIT_NAME"  
# [9] "SOURCE_REPORTING_UNIT"      "SOURCE_REPORTING_UNIT_NAME" "LOCAL_FIRE_REPORT_ID"       "LOCAL_INCIDENT_ID"         
# [13] "FIRE_CODE"                  "FIRE_NAME"                  "ICS_209_INCIDENT_NUMBER"   "ICS_209_NAME"              
# [17] "MTBS_ID"                    "MTBS_FIRE_NAME"             "COMPLEX_NAME"              "FIRE_YEAR"                 
# [21] "DISCOVERY_DATE"             "DISCOVERY_DOY"              "DISCOVERY_TIME"            "STAT_CAUSE_CODE"           
# [25] "STAT_CAUSE_DESCR"           "CONT_DATE"                  "CONT_DOY"                  "CONT_TIME"                 
# [29] "FIRE_SIZE"                  "FIRE_SIZE_CLASS"            "LATITUDE"                  "LONGITUDE"                 
# [33] "OWNER_CODE"                 "OWNER_DESCR"                "STATE"                     "COUNTY"                    
# [37] "FIPS_CODE"                  "FIPS_NAME"                               

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

# We want to assign the NA_L2CODE and NA_L2NAME to each fire (row) in the FPA FOD.
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
  assignment <- p_ecoregion$NA_L2CODE
  
  # If the fire does not overlap an ecoregion polygon the result will be NA. We
  # need to know where this occurs. Because the borders of polygons are not 
  # perfect there are gaps. We now want to assign the ecoregion with a border
  # that is closest to the fire. 
  if(is.na(assignment)){
    
    print(paste("Using distances Matrix for: ", i))
    
    # Now this calculation can be done on many fewer polygons
    dist.mat <- geosphere::dist2Line(p = p, line = SPDF)
    
    # Assign the closest polygon in the ecoregions polygons
    NA_L2CODE[i] <- SPDF@data$NA_L2CODE[dist.mat[,4]]
    
  } else{
    
    # The point fell within a polygon and we are already set! 
    NA_L2CODE[i] <- assignment
    
  }
  
  # Output progress to the screen
  if(i %% 1000 == 0){
    print(paste("Percent Complete: ", i/nRow*100))
  }
    
}
tf <- Sys.time()

print(paste("It took", tf - t1, "to complete loop"))

# Assign all overlap analysis information to the df and save it! 
df$NA_L2CODE <- NA_L2CODE

# Save the data as .RData, but only do that if you have preserved the number of
# fires in this year. 
if(dim(df)[1] == nRow & i == nRow){

  FPA_FOD <- df
  save(FPA_FOD, file=paste0(drive, "FPA_FOD_", year, ".RData"))
  
  print(paste("Complete:", year))
} else{
  
  stop("The number of fires in this year is not the same after processing.")
  
}

