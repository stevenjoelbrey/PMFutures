# plot_ecoregion_burn_totals.R

# ------------------------- Description ---------------------------------------
# This script is meant to make a plot that shows that states (or eco-regions) 
# with large burn area are dominated by lightning ignitions. The script creates
# a bar plot and a map. Together they make the point. 

library(rgdal)
library(ggplot2)
library(rgeos)
library(maptools)
library(grid)
library(gridExtra)
library(sfsmisc)

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

# Count ecoregion totals
for (i in 1:nRegions){
  
  # Create region Mask
  rMask <- keepRegions[i] == FPA_FOD$NA_L2CODE
  
  regionBA_L[i] <- sum(FPA_FOD$FIRE_SIZE[rMask & lightningMask])
  regionBA_H[i] <- sum(FPA_FOD$FIRE_SIZE[rMask & !lightningMask])
  
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
    res=250, height=1600, width=3000)
par(las=1, mar=c(4,14,4,2))

bp <- barplot(height=BA, names.arg=keepRegions, yaxt="n", col="gray", 
        cex.names = 1.8)
#axis(side=1, at = bp, labels=keepRegions, col=)

barplot(height=regionBA_H, add=T, yaxt="n", col="orange")

# Nice axis label 
aY <- axTicks(2); axis(2, at=aY, label= axTexpr(2, aY), cex.axis=1.5)
mtext("Acres \n Burned", side=2, line=7, cex=2)

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



