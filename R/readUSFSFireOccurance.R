#!/usr/bin/env Rscript


################################################################################
# readUSFSFireOccurance.R

# This script is used to convert spatial wildfire occurance data into a R 
# friendly easy to load dataframe and csv. This script will laod a text file of this
# data. That text file was created using microft access, where the "fires"
# field was exported. The main points of this script are too
# 1) convert the time format to an R friendly format.
# 2) create a "START_MONTH" column. 
# 3) save as a csv and .RData format. 

# NOTE: All fire (1.88 million) burn area should be about 140 million acres.
# NOTE: When I sum quantity here I get 139.7639 million. Yay! 

# DataSource: https://www.fs.usda.gov/rds/archive/Product/RDS-2013-0009.4//
# .accdb file link: https://www.fs.usda.gov/rds/fedora/objects/RDS:RDS-2013-0009.4/datastreams/RDS-2013-0009.4_ACCDB/content  

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
library(lubridate) # for month()


print("-----------------------------------------------------------------------")
print("Reading in the really big fire .RData file")
print("-----------------------------------------------------------------------")

drive    <- "Data/FPA_FOD/"
dataFile <- paste0(drive,"Fires_from_accdb.txt")
df       <- read.delim(dataFile, header=FALSE, sep=",")

# Now load the file that has the headers (could only export a few rows as csv
# because Microsoft Access is garbage.)
headers_csv <- read.csv(paste0(drive, "Fires_from_accdb_headers.csv")) 

# assign these names to df. They came from the same place. 
names(df) <- names(headers_csv)

# Names of the df. For some reason I can't get stupid microsoft Access to export
# the headers. 
# [1] "FOD_ID"                     "FPA_ID"                     "SOURCE_SYSTEM_TYPE"        
# [4] "SOURCE_SYSTEM"              "NWCG_REPORTING_AGENCY"      "NWCG_REPORTING_UNIT_ID"    
# [7] "NWCG_REPORTING_UNIT_NAME"   "SOURCE_REPORTING_UNIT"      "SOURCE_REPORTING_UNIT_NAME"
# [10] "LOCAL_FIRE_REPORT_ID"       "LOCAL_INCIDENT_ID"          "FIRE_CODE"                 
# [13] "FIRE_NAME"                  "ICS_209_INCIDENT_NUMBER"    "ICS_209_NAME"              
# [16] "MTBS_ID"                    "MTBS_FIRE_NAME"             "COMPLEX_NAME"              
# [19] "FIRE_YEAR"                  "DISCOVERY_DATE"             "DISCOVERY_DOY"             
# [22] "DISCOVERY_TIME"             "STAT_CAUSE_CODE"            "STAT_CAUSE_DESCR"          
# [25] "CONT_DATE"                  "CONT_DOY"                   "CONT_TIME"                 
# [28] "FIRE_SIZE"                  "FIRE_SIZE_CLASS"            "LATITUDE"                  
# [31] "LONGITUDE"                  "OWNER_CODE"                 "OWNER_DESCR"               
# [34] "STATE"                      "COUNTY"                     "FIPS_CODE"                 
# [37] "FIPS_NAME"                  

# Format time values for use in R and place back into dataframe 

DISCOVERY_DATE <- as.character(df$DISCOVERY_DATE)
discovery_date <- as.POSIXct(DISCOVERY_DATE, format="%m/%d/%Y %H:%M:%S", tz="UTC")
df$DISCOVERY_DATE <- discovery_date

# Repeat this method for the containment date
CONT_DATE <- as.character(df$CONT_DATE)
cont_date <- as.POSIXct(CONT_DATE, format="%m/%d/%Y %H:%M:%S", tz="UTC")
df$CONT_DATE <- cont_date

# Make an array of the month of the start date 
library(lubridate) # for month()
start_month <- lubridate::month(DISCOVERY_DATE)

# Add it to the dataframe
df$START_MONTH <- start_month

print(paste("There are ", dim(df)[1]/10^6, "million in this record. ~1.88 expected"))

# For savename, we need to know the year range of these data
minYear <- min(df$FIRE_YEAR)
maxYear <- max(df$FIRE_YEAR)

# Make sure the name is descrptive of the years in these data. 
print("The data have been formatted yay!")
FPA_FOD <- df
save(FPA_FOD, file=paste0(drive, "FPA_FOD_", minYear,"_",maxYear, ".RData"))

print("script executed, output saved without error")
