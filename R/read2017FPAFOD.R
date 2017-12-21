# readUSGSFire.R

# New, probably the same data, with more metadata, to explore! 

# For all wildland fire 
# DataSource: https://www.fs.usda.gov/rds/archive/Product/RDS-2013-0009.4/
# DataURL: https://www.fs.usda.gov/rds/fedora/objects/RDS:RDS-2013-0009.4/datastreams/RDS-2013-0009.4_GDB/content

# NOTE: I did not convert these datas to a shapefile. That was done with the
# NOTE: very generous help of Jill Hubley. The data were converted from .DBF to
# NOTE: shapefiles using QGIS. Thanks Jill! 

# load the required libraries 
library(sp) 
library(rgdal)

externalDir <- "/Volumes/Brey_external/FPAFOD20170508"
layerName   <- "FPAFOD20170508"
data_projected <- readOGR(dsn=externalDir, layer=layerName) 

df <- data_projected@data
nCol <- length(names(df))

# Give columns the same names as previous, existing working version 
load("Data/FPA_FOD/FPA_FOD_2003_2013.RData")
# For loaded file "FPA_FOD" object, I added the last 2 columns. START_MONTH and
# NA_L2CODE
names(df) <- names(FPA_FOD)[1:nCol]

write.csv(df, file="TEMP.csv", row.names=FALSE)

# This will automatically handle the pesky factors that come out of the projected
# data
df <- read.csv("TEMP.csv", 
               stringsAsFactors = FALSE) 

################################################################################
# ------- Now that we have a non-factor version HANDLE THE DATES ---------------
################################################################################

# Get start time, containment time, and out time into POSIXct format
disc_mmddyy_TIME    <- df$DISCOVERY_DATE # "YYYY/MM/DD 00:00:00.000" <- format 
spaceLocation       <- str_locate(disc_mmddyy_TIME, " ")

# The time of these incidents is not needed. Just the date. These data are not 
# that precise. 
disc_mmddyy    <- str_sub(disc_mmddyy_TIME, 1, spaceLocation[,1])

# Get rid of the space that sometimes shows up at the end 
disc_mmddyy_noSpace <- str_replace(disc_mmddyy, " ", "")
discoverd_date <- as.POSIXct(disc_mmddyy_noSpace, format="%Y/%m/%d", tz="UTC")

# Add to dataframe
df$DISCOVERY_DATE <- discoverd_date

# Repeat this method for the containment date
con_mmddyy_TIME <- df$CONT_DATE
spaceLocation <- str_locate(con_mmddyy_TIME, " ")
con_mmddyy <- str_sub(con_mmddyy_TIME, 1, spaceLocation[,1])
con_mmddyy_noSpace  <- str_replace(con_mmddyy, " ", "")
con_date <- as.POSIXct(con_mmddyy_noSpace, format="%Y/%m/%d", tz="UTC")

# Store in the dataframe
df$CONT_DATE <- con_date

# Make an array of the month of the start date 
t_LT <- as.POSIXlt(discoverd_date)
month <- t_LT$mon + 1

# Add it to the dataframe
df$START_MONTH <- month


# Now write the factor free version of the dataframe as .RData file
FPA_FOD <- df
save(FPA_FOD, file="Data/FPA_FOD/FPA_FOD_1992_2015_noEcoregion.RData")
write.csv(FPA_FOD, file="Data/FPA_FOD/FPA_FOD_1992_2015_noEcoregion.csv")
