#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  year <- 2003
}else{
  year <-   args[1]
  # year <-   args[2]
  # minLat <- args[3]
  # maxLat <- args[4]
  # minLon <- args[5]
  # maxLon <- args[6]
  # regionName <- args[7]
}

print(paste("Loading data for year:", year))

minLat <- 24
maxLat <- 49
minLon <- -125
maxLon <- -65
regionName <- "CONUS"

# subset_FINN_txt_data.R
# ------------------------- Description ---------------------------------------
# This script subsets FINN txt data to retain only wildfires in CONUS. 
# ------------------------------------------------------------------------------

dataDir <- "/Users/sbrey/GoogleDrive/Data/FINN/"
f <- paste0(dataDir, "FINN_", year, ".txt")
# Read in the data
df <- read.csv(f, stringsAsFactors = F)

latMask <- (df$LATI >= minLat) & (df$LATI <= maxLat)
lonMask <- (df$LONGI >= minLon) & (df$LONGI <= maxLon)
spatialMask <- latMask & lonMask

df_CONUS <- df[spatialMask, ]

saveName <- paste0(dataDir, "FINN_",regionName,"_", year, ".txt")
write.csv(df_CONUS, file=saveName, row.names = FALSE)

# Now make an R friendly version where time has been converted to a friendly
# analysis unit 

