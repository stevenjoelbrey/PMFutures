# plot_ecoregion_burn_totals.R

# ------------------------- Description ---------------------------------------
# This script is meant to make a plot that shows that states (or eco-regions) 
# with large burn area are dominated by lightning ignitions. The script creates
# a bar plot and a map. Together they make the point. We are also going to show
# GFED4s emission estimate totals for each ecoregion. This will demonstrate that
# not all burn_area is equal. 
# NOTE: When comparing GFED4s emissions and FPA-Burn area can only go back to
# NOTE: 1997 for GFED. 

# TODO: rerun GFED4s ecoregion totals when createECMWFGridAttributes.R has been
# TODO: re-created such that assignments are made for all interior grid cells
# TODO: and coastal grid cells. 

library(rgdal)
library(ggplot2)
library(rgeos)
library(maptools)
library(grid)
library(gridExtra)
library(sfsmisc)
library(ncdf4)

startYear <- 1997
endYear   <- 2015

# What region are you investigating? 
minLat <- 31
maxLat <- 49
minLon <- -125
maxLon <- -100

keepRegions <- c(6.2, 10.1, 7.1, 9.3, 9.4, 13.1, 12.1, 11.1, 10.2)
regionColors <- c("#65B657", "#E7ED90", "#60BAAF", "#F2DBA1", "#EDCB9B",
                  "#B5D67C",  "#D6D292", "#D2E5B1", "#F3DB70")
nRegions <- length(keepRegions)
regionBA_L <- rep(NA, nRegions) # lightning
regionBA_H <- rep(NA, nRegions) # human 
  
# Load FPA-FOD data and subset by latitude and longitude 
load("Data/FPA_FOD/FPA_FOD_1992_2015.RData")

# Spatial subset 
latMask <- FPA_FOD$LATITUDE > minLat & FPA_FOD$LATITUDE < maxLat
lonMask <- FPA_FOD$LONGITUDE > minLon & FPA_FOD$LONGITUDE < maxLon
regionMask <- FPA_FOD$NA_L2CODE %in% keepRegions

# spatially subset the data and get a mask for the fires that are started by
# lightning 
FPA_FOD <- FPA_FOD[latMask & lonMask & regionMask, ]
lightningMask <- FPA_FOD$STAT_CAUSE_DESCR == "Lightning"

# Count ecoregion total burn area from FPA-FOD. 
for (i in 1:nRegions){
  
  # Create region Mask
  rMask <- keepRegions[i] == FPA_FOD$NA_L2CODE
  
  regionBA_L[i] <- sum(FPA_FOD$FIRE_SIZE[rMask & lightningMask])
  regionBA_H[i] <- sum(FPA_FOD$FIRE_SIZE[rMask & !lightningMask])
  
}

################################################################################
# Now sum the total GFED4s emissions by ecoregion for the same spatial and 
# temporal domain. 
################################################################################
if(startYear == 1997 & endYear == 2015){
  
  gfedDir <- "/Volumes/Brey_external/GFED4s/"
  ncFile <- paste0(gfedDir, "GFED4.1s_ecmwf_monthly_DM_1997.nc")
  
  nc <- nc_open(ncFile)
  grid_lat <- ncvar_get(nc, "latitude") 
  grid_lon <- ncvar_get(nc, "longitude")
  time <- ncvar_get(nc, "time")
  monthly_DM <- ncvar_get(nc, "monthly_DM")
  nc_close(nc)
  
  # Sum the first years 12 months of data
  yearSum <- apply(monthly_DM, 1:2, sum)
 
  
  # Add each year sum to yearSum
  for(y in (startYear + 1):endYear){
    
    ncFile <- paste0(gfedDir, "GFED4.1s_ecmwf_monthly_DM_1997.nc")
    nc <- nc_open(ncFile)
    monthly_DM <- ncvar_get(nc, "monthly_DM")
    nc_close(nc)
    
    # append the total burn area on this grid. Take all months. 
    yearSum <- yearSum + apply(monthly_DM, 1:2, sum)
    
  }
  
  # Plot the totals to make sure it looks normal 
  quartz() 
  f <- length(grid_lat):1
  yearSumPlot <- yearSum
  yearSumPlot[yearSum==0] <- NA
  image.plot(grid_lon, grid_lat[f], yearSumPlot[,f])
  title("Sum of DM")
  
  # Now we need to subset by the western US coords and then load grid attributes
  # to subset. 
  ncFile <- "Data/grid_attributes/grid_attributes_75x75.nc"
  nc <- nc_open(ncFile)
  ecoregion <- round(ncvar_get(nc, "ecoregion"),1)
  elevation <- ncvar_get(nc, "elevation")
  nc_close(nc)
  
  # For finding North America CONUS
  psuedo_lon <- grid_lon - 360
  
  quartz()
  image.plot(psuedo_lon, grid_lat[f], elevation[,f])
  title("elevation")
    
  quartz()
  image.plot(psuedo_lon, grid_lat[f], ecoregion[,f])
  title("ecoregions")
  
  # Create the spatial mask subset
  lonMask <- psuedo_lon >= minLon & psuedo_lon <= maxLon
  latMask <- grid_lat >= minLat & grid_lat <= maxLat
  
  # subset all the variables that these of these coords
  yearSum_subset <- yearSum[lonMask, latMask]
  ecoregion_subset <- ecoregion[lonMask, latMask]
  elevation_subset <- elevation[lonMask, latMask]
  grid_lat_subset <- grid_lat[latMask]
  grid_lon_subset <- grid_lon[lonMask]
  psuedo_lon_subset <- psuedo_lon[lonMask]
  
  quartz()
  f <- length(grid_lat_subset):1
  image.plot(psuedo_lon_subset, grid_lat_subset[f], yearSum_subset[,f])
  title("Total DM subset")
  map("state", add=T)
  
  quartz()
  f <- length(grid_lat_subset):1
  image.plot(psuedo_lon_subset, grid_lat_subset[f], elevation_subset[,f])
  title("And now elevation it has been subset")
  map("state", add=T)
  
  quartz()
  f <- length(grid_lat_subset):1
  image.plot(psuedo_lon_subset, grid_lat_subset[f], ecoregion_subset[,f])
  title("And now ecoregion has been subset. Make sure no white space in middle of country")
  map("state", add=T)
  
  # Count keeper ecoregion total burn area from GFED4s. 
  regionDM <- rep(NA, nRegions) # lightning & human. GFED does not know difference 
  for (i in 1:nRegions){
    
    # Create region Mask
    rMask_2D <- keepRegions[i] == ecoregion_subset
    rMask_2D[is.na(rMask_2D)] <- FALSE
    
    regionDM[i] <- sum(yearSum_subset[rMask_2D])

  }
  
  
  
}
################################################################################
# Bar plot of the area burned totals! 
# For nice labels 
# https://jangorecki.gitlab.io/data.table/library/sfsmisc/html/axTexpr.html
################################################################################
BA <- regionBA_L + regionBA_H

# sort 
ORDER <- order(BA, decreasing = TRUE)
BA <- BA[ORDER]
keepRegions <- keepRegions[ORDER]
regionBA_L  <- regionBA_L[ORDER]
regionBA_H  <- regionBA_H[ORDER]
regionColors <- regionColors[ORDER]

png(filename="Figures/summary/ecoregio_burn_area_bar.png", 
    res=250, height=1600, width=3300)
par(las=1, mar=c(4,14,4,8))

bp <- barplot(height=BA, 
              names.arg=keepRegions, yaxt="n", col="gray", 
              cex.names = 1.8)
#axis(side=1, at = bp, labels=keepRegions, col=)

barplot(height=regionBA_H, add=T, yaxt="n", col="orange")

# Nice axis label 
aY <- axTicks(2); axis(2, at=aY, label= axTexpr(2, aY), cex.axis=1.5)
mtext("Acres \n Burned", side=2, line=7, cex=2)

x = par("usr")

oldPar <- par()


# Add emissions from GFED4s if we are looking at the right years 
if(startYear == 1997 & endYear == 2015){
  par(new = TRUE)

  bp_new <- barplot(height=regionDM, pch="*", border="transparent", col="transparent",
                    axes = FALSE, bty = "n", xlab = "", ylab = "")
  text(bp_new, regionDM[ORDER], "*", col="red", cex=4)
  eaxis(side=4, at = pretty(range(regionDM)), col="red")
  axis(side=4, at = pretty(range(regionDM)), col="red", labels=rep("",7) )
  
  mtext("", side=4, line=3)
}

dev.off()


# Load spatial data 
ecoRegionFile <- "Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData"
map.det <- get(load(ecoRegionFile))

# Make the annoying factor a numeric array 
map.det$NA_L2CODE <- as.numeric(as.character(map.det$NA_L2CODE))

png(filename="Figures/summary/ecoregio_burn_area_map.png", 
    res=250, height=2350, width=3000)
par(las=1, mar=c(0,0,0,0))

# Make a simple map showing the ecoregions 
map("state", xlim=c(minLon-7, maxLon), ylim=c(minLat, maxLat), lwd=0.1)

for (i in 1:nRegions){
  m <- as.numeric(as.character(map.det$NA_L2CODE)) == keepRegions[i]
    plot(map.det[m,], col=regionColors[i], add=T)
}
map("state", add=T, lty=2)
legend("topleft",
       legend=keepRegions,
       fill=regionColors,
       bty="n",
       cex=3)

dev.off()
# # get centroids
# map.test.centroids <- gCentroid(map.test, byid=T)
# map.test.centroids <- as.data.frame(map.test.centroids)
# map.test.centroids$OBJECTID <- row.names(map.test.centroids)


# TODO: optional cool barplot on map version of figure
# map the data! 
# https://stackoverflow.com/questions/36063043/how-to-plot-barchart-onto-ggplot2-map



