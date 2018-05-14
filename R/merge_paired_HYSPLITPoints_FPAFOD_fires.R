# merge_paired_HYSPLITPoints_FPAFOD_fires.R

# ------------------------- Description ---------------------------------------

year1 <- 2007 
year2 <- 2015

# Load and append the yearly wildfire data 
dataDir <- "Data/HMS/" 

# If the merged datasets do not exist, make them. 

# Merge the yearly files 
df <- get(load( paste0(dataDir, "FPA_FOD_with_HP_", year1, ".RData") ))
df_hms <- get(load( paste0(dataDir, "hysplitPoints_with_FPAFOD_", year1, ".RData") ))
for(y in (year1+1):year2){
  
  print(paste("Appending:", y, "to", year1))
  
  # Load "FPA_FOD" dataframe 
  load( paste0(dataDir, "FPA_FOD_with_HP_", y, ".RData") )
  df <- rbind(df, FPA_FOD)
  rm(FPA_FOD)
  
  # Load "hysplitPoints" data frame
  load( paste0(dataDir, "hysplitPoints_with_FPAFOD_", y, ".RData") )
  df_hms <- rbind(df_hms, hysplitPoints)
  rm(hysplitPoints)
  
}

# Give them the original load nice name now that they are merged 
FPA_FOD <- df
hysplitPoints <- df_hms
rm(df, df_hms) # So there is no confusion about the data floating around

# Make a "paired" variable for making factors needed for plotting easy later
# for each dataset. This will help with plotting later when re-doing this 
# analysis with AQ context. 
FPA_FOD$pairedWithHMS <- FPA_FOD$n_HP > 0
hysplitPoints$pairedWithFPAFOD <- hysplitPoints$nFPAFODPaired > 0 

################################################################################
# Now, assign large regions to the dataframes. I want a factor that tells me if
# the FPA FOD wildfire is in the SE or West. 
################################################################################
region <- rep("", dim(hysplitPoints)[1])

load("Data/GIS/southeast_bounds.RData") # loads "minLat"... 
EMask <- (hysplitPoints$Lon > minLon) & (hysplitPoints$Lon < maxLon) & 
         (hysplitPoints$Lat>minLat) & (hysplitPoints$Lat < maxLat)

load("Data/GIS/west_bounds.RData") # loads "minLat"... 
WMask <- (hysplitPoints$Lon > minLon) & (hysplitPoints$Lon < maxLon) & 
  (hysplitPoints$Lat>minLat) & (hysplitPoints$Lat < maxLat)

region[WMask] <- "West"
region[EMask] <- "Southeast"
hysplitPoints$region <- factor(region, levels = c("West", "Southeast", ""))

################################################################################
# Now FPA FOD
region <- rep("", dim(FPA_FOD)[1])

load("Data/GIS/southeast_bounds.RData") # loads "minLat"... 
EMask <- (FPA_FOD$LONGITUDE > minLon) & (FPA_FOD$LONGITUDE < maxLon) & 
         (FPA_FOD$LATITUDE > minLat) & (FPA_FOD$LATITUDE < maxLat)

load("Data/GIS/west_bounds.RData") # loads "minLat"... 
WMask <- (FPA_FOD$LONGITUDE > minLon) & (FPA_FOD$LONGITUDE < maxLon) & 
         (FPA_FOD$LATITUDE > minLat) & (FPA_FOD$LATITUDE < maxLat)

region[WMask] <- "West"
region[EMask] <- "Southeast"
FPA_FOD$region <- factor(region, levels = c("West", "Southeast", ""))


save(FPA_FOD, file=paste0(dataDir, "FPA_FOD_with_HP_", year1, "_", year2, ".RData"))
save(hysplitPoints, file=paste0(dataDir, "hysplitPoints_with_FPA_FOD_", 
                                year1,"_", year2, ".RData"))
