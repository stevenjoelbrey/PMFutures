# sanity_check_FINN_processed_data.R

# ------------------------- Description ---------------------------------------
# This script is used to perform basic checks on the FPA FOD wildfires that 
# were put through the FINN framework to have PM25 emissions estimated. It will
# save out the processed data with FPA FOD attributes appended. 

# ------------------------- Arguments ------------------------------------------

year <- 2002


library(rgdal)
library(maps)
library(ggplot2)

# Read in the data that was given to be processed
data_projected <- readOGR(dsn="Data/FINN/shape/", layer=paste0("FPA_FOD_", year))
df_in <- data_projected@data

load(paste0("Data/FINN/FPA_FOD_forFINN_", year, ".RData"))

if(!sum(df_in$FOD_ID == df$FOD_ID)/length(df$FOD_ID)==1){
  stop("The data do not have the same number of rows, attributes lost")
} else{
  # If the code is still running it means all the rows matched. Now see if they
  # match the FINN processed data. 
  rm(df_in)
}

# # Read in the data processed by Christine
df_processed <- read.csv("Data/FINN/FINNv1.5_BREY_2006_061220182.txt")

# headers <- names(df_processed)
# # Kg of PM2.5
# df_processed <- read.delim("Data/FINN/FINNv1.5_BREY_2002_06122018_test2.txt",
#                            sep=",", header = FALSE, skip=1)
# d <- dim(df_processed)
# df_processed <- df_processed[,1:length(headers)]
# names(df_processed) <- headers
#   
# print(names(df_processed))

# Do the dimensions (rows) match the file you gave? 
print(paste("Original dimensions", dim(df)))
print(paste("Processed dimensions", dim(df_processed)))

# Do the FOD_ID match for every row? If yes cbind them
test <- unique(df_processed$FOD_ID == df$FOD_ID)

IDinOriginal <- match(df_processed$FOD_ID, df$FOD_ID)

# Make sure the match aligns
unique(df_processed$FOD_ID == df$FOD_ID[IDinOriginal])

if(test){
  df_combined <- cbind(df[IDinOriginal,], df_processed)
}
print(names(df_combined))

# Map the locations of the data you gave and what you got back. 
# Are they exactly the same? 
unique(df_combined$lat - df_combined$LATITUDE)
unique(df_combined$long - df_combined$LONGITUDE)

quartz()
plot(df_combined$long, df_combined$lat, pch=".", col=df_combined$genLC)
map("state", add=T)

# Map and color code the points by land cover
hist(df_combined$genLC)

# Plot emissions [PM2.5] vs. fire size for each land cover type
plot(df_combined$PM25, df_combined$AREA_m2, col=df_combined$genLC, pch=19)

# Look at human and lightning emission totals
lightningMask <- df_combined$STAT_CAUSE_DESCR == "Lightning"

# Total emissions of each. 
lightningPM25 <- sum(df_combined$PM25[lightningMask])
humanPM25     <- sum(df_combined$PM25[!lightningMask])




