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

localDir <- "Data/GIS/FPAFOD20170508/"
layerName <- "FPAFOD20170508"
data_projected <- readOGR(dsn=localDir, layer=layerName) 

df <- data_projected@data
nCol <- length(names(df))

# Give columns the same names as previous, existing working version 
load("Data/FPA_FOD/FPA_FOD_2003_2013.RData")
# For loaded file "FPA_FOD" object, I added the last 2 columns. START_MONTH and
# NA_L2CODE
names(df) <- names(FPA_FOD)[1:nCol]

write.csv(df, file="Data/FPA_FOD/FPA_FOD_1992_2015.csv", row.names=FALSE)

# This will automatically handle the pesky factors that come out of the projected
# data
df <- read.csv("Data/FPA_FOD/FPA_FOD_1992_2015.csv", 
               stringsAsFactors = FALSE) 

# Now write the factor free version of the dataframe as .RData file
FPA_FOD <- df
save(FPA_FOD, file="Data/FPA_FOD/FPA_FOD_1992_2015.RData")
