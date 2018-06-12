# process_FPA_FOD_for_FINN.R

# ------------------------- Description ---------------------------------------
# This script takes FPA FOD fire data and creates .txt files as well as 
# "shape-point" files. The desired datum for the .shp files is WGS84
# Required fields for fires include: 
#   LON, LAT, AREA, LCT, TREE, HERB, BARE

################################################################################
# Arguments
################################################################################
dataInDir  <- "Data/FPA_FOD/"
dataOutDir <- "Data/FINN/"
m2perAcre  <- 4046.86 # m2/acre
blank      <- 0 # Value to fill columns that are needed but that I do not know
year1      <- 2002
year2      <- 2015
crs_string <- "+proj=longlat +datum=WGS84"

library(rgdal) # For loading required libraries. 


# Select years to create files for
years <- year1:year2

# Loop through each year
for (year in years){

  print(paste("Processing year:", year))
  
  loadFile <- paste0(dataInDir, "FPA_FOD_", year, ".RData")
  print(loadFile)
  
  # Loads a dataframe labeled "FPA_FOD""
  load(loadFile)
  
  # Mask data to be CONUS only. Do this using STATES attribute 
  statesToKeep <- !(FPA_FOD$STATE %in% c("PR", "HI", "AK"))
  df   <- FPA_FOD[statesToKeep,]
  
  # Calculate area in SI units (m2)
  AREA <- df$FIRE_SIZE * m2perAcre # acres  * (m2/acres) = m2
  
  # Create the dataframe to save out. This takes attributes from FPA FOD for 
  # this year and creates columns for FINN to add information. 
  df_out <- data.frame(LON=df$LONGITUDE, LAT=df$LATITUDE, AREA=AREA, 
                       LCT=blank, TREE=blank, HERB=blank, BARE=blank, 
                       FOD_ID=df$FOD_ID)
  
  saveFile <- paste0(dataOutDir, "FPA_FOD_", year, ".csv")
  write.csv(df_out, file=saveFile, row.names=FALSE)
  
  #######################################################
  # Now create and write these data as a shape point file
  #######################################################
  coords <- cbind(df$LONGITUDE, df$LATITUDE)
  
  spdf <- SpatialPointsDataFrame(coords=coords, 
                                 data=df_out, 
                                 proj4string=CRS(crs_string))
  
  writeOGR(spdf, 
           dsn = paste0(dataOutDir, "/shape/"), 
           layer=paste0("FPA_FOD_", year), 
           driver="ESRI Shapefile")
  
  # Clear out the variables before we come back through here. 
  rm(spdf, df_out, FPA_FOD, df, coords, AREA)

} # End of looping through years 