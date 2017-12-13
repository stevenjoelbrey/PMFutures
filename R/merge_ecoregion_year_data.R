# merge_ecoregion_year_data.R

# assign_ecoregion_to_FPAFOD.R assigns each fire event the ecoregion it falls
# into. To make this way faster, this is done one year at a time. This script
# will be used to sanity check the assignments and also to merge these yearly 
# files together. Really this whole script and seperation of the former is
# only needed because I do not know how to parallel program in R. 

# TODO: Eventually this will  not need to draw from the "_borderfix_" version of
# TODO: the yearly data. 

library(maps)
library(sp)
library(rgdal)
library(maptools)
library(rgeos)

dataDir   <- "Data/FPA_FOD/"
startYear <- 1992
endYear   <- 2013
years <- startYear:endYear

# Loop through and load desired files. 
for (year in years){
  
  f <- paste0(dataDir, "FPA_FOD_borderfix_",year,".RData")
  
  if(year == years[1]){
    FPA_FOD_merged <- get(load(f))
  } else{
    FPA_FOD        <- get(load(f))
    FPA_FOD_merged <- rbind(FPA_FOD_merged, FPA_FOD)
  }
  
}

print("Unique years in the merged dataset include..")
print(unique(FPA_FOD_merged$FIRE_YEAR))


FPA_FOD_merged$NA_L2CODE <- as.numeric(as.character(FPA_FOD_merged$NA_L2CODE))

# Need to handle the points that fall inbetween polygons. Assignments must be
# made for these ~ 33,000 fires. 

# TODO: Remove Puerto Rico fires, these are outside our scope and do not have
# TODO: ecoregions. AND HAWAII
unique(FPA_FOD_merged$STATE)

print("Number of fires with no ecoregion assignment made:")
print(sum(is.na(FPA_FOD_merged$NA_L2CODE)))


coordinates <- cbind(FPA_FOD_merged$LONGITUDE, FPA_FOD_merged$LATITUDE)
fireLocations <- SpatialPoints(coordinates)

map("world", ylim=c(10,80), xlim=c(-170, -50))
# Any NAs? Where are they
plot(fireLocations[is.na(FPA_FOD_merged$NA_L2CODE),], add=T, col="red", pch=".")
FPA_FOD_merged$NA_L2CODE[is.na(FPA_FOD_merged$NA_L2CODE)] <- 0
plot(fireLocations[FPA_FOD_merged$NA_L2CODE==6.2,], add=T, col="forestgreen",pch=".")

# Alaska and Northwest coast
plot(fireLocations[FPA_FOD_merged$NA_L2CODE==7.1,], add=T, col="lightblue",pch=".")
# Great Basin and Plains of Western Wyoming and four corners. 
plot(fireLocations[FPA_FOD_merged$NA_L2CODE==10.1,], add=T, col="tan",pch=".")
# South tip of Florida
plot(fireLocations[FPA_FOD_merged$NA_L2CODE==15.4,], add=T, col="purple",pch=".")

FPA_FOD <- FPA_FOD_merged
save(FPA_FOD, file=paste0("Data/FPA_FOD/FPA_FOD_",startYear,"_",endYear,".RData"))
write.csv(FPA_FOD, file=paste0("Data/FPA_FOD/FPA_FOD_",startYear,"_",endYear,".csv"))

