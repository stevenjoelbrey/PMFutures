# ------------------------- Description ---------------------------------------
# This script will be used to show the regional differences in inteannual 
# variability in total burn area across inventories. 
# sister function for spatial component is called Python/compare_burn_area_maps.py

# What region are you investigating? 
minLat <- 31
maxLat <- 49
minLon <- -125
maxLon <- -100

includedMonths <- 5:10

# for labelling files
monthSpan <- paste0("months_",min(includedMonths), "_", max(includedMonths))

library(ncdf4)
library(ggplot2)
library(fields)
library(RColorBrewer)

# Load ecoregions spatial shapes
load("Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData")


# Load the needed burn area products and grid attributes
grid_nc_file <- "Data/grid_attributes/grid_attributes_25x25.nc"
nc_grid   <- nc_open(grid_nc_file)
ecoregion <- round(ncvar_get(nc_grid, "ecoregion"),3) # annoying precision issues
elevation <- ncvar_get(nc_grid, "elevation")
state     <- ncvar_get(nc_grid, "state")
latitude  <- nc_grid$dim[["latitude"]][["vals"]]
longitude <- nc_grid$dim[["latitude"]][["vals"]]
nc_close(nc_grid)

# Get GFED4s monthly burn area from DM file
GFED_nc_file <- "Data/GFED4s/GFED4.1s_monthly_DM_2003_2016.nc"
nc_GFED <- nc_open(GFED_nc_file)
GFED_YYYYMM <- ncvar_get(nc_GFED, "time")
GFED_latitude <- ncvar_get(nc_GFED, "latitude")
GFED_longitude <- ncvar_get(nc_GFED, "longitude")
GFED_BAF <- ncvar_get(nc_GFED, "burn_area_fraction")
grid_area <- ncvar_get(nc_GFED, "grid_area")
nc_close(nc_GFED)

# Need total burn area, not burn area fraction as is stored in nc files. 
# NOTE: We want the time dimension to come first, to match the FPA-FOD data. so
# NOTE: we will rewrite this array one time dimension at a time. 

GFED_BA_t <- array(NA, dim=c(168, 1440, 720)) # No NA in array to be transformed so check after
for (i in 1:dim(GFED_BAF)[3]){
  GFED_BA_t[i,,] <- GFED_BAF[,,i] * grid_area  # [1440  720] dot [1440  720]
}

# Sanity check, need to make sure every position of NA was overwritten by data
# write
if(sum(is.na(GFED_BA_t))){
  stop("Data not placed correctly into transposed matrix")
} else{
  GFED_BA <- GFED_BA_t 
  rm(GFED_BA_t)
}

print("Dim of GFED_BA")
print(dim(GFED_BA))

# Get FPA_FOD gridded burn area
FPA_nc_file <- "Data/FPA_FOD/burn_area_monthly_25x25_2003_2013.nc"
nc_FPA <- nc_open(FPA_nc_file)
FPA_YYYYMM <- ncvar_get(nc_FPA, "YYYYMM")
FPA_latitude <- ncvar_get(nc_FPA,"latitude")
FPA_longitude <- ncvar_get(nc_FPA,"longitude")
FPA_BA <- ncvar_get(nc_FPA, "burn_aera")
FPA_BA_human <- ncvar_get(nc_FPA, "human_burn_area")
FPA_BA_lightning <- ncvar_get(nc_FPA, "lightning_burn_aera")
nc_close(nc_FPA)

print("Dim of FPA_BA")
print(dim(FPA_BA))

# Sanity check, make sure the spatial dimensions of the data agree
if (unique(GFED_latitude - FPA_latitude)==0 & 
    unique(GFED_longitude - FPA_longitude)==0){
  print("dimensions match, good to go.")
} else{
  stop("Spatial dimensions do not match. Something went wrong.")
}

# ------------------------- visualSanityCheck ----------------------------------
# visual sanity check, turn zeros into nans so it is easier to see gradients.
# NOTE: The time dimension does not actually match here because GFED has not 
# NOTE: been subset yet. 

FPA_BA_NA <- FPA_BA
FPA_BA_NA[FPA_BA==0] <- NA

GFED_BA_NA <- GFED_BA
GFED_BA_NA[GFED_BA==0] <- NA

quartz( width=10, height=5)
par(mfrow=c(1,2))
image(FPA_BA_NA[7,,])
image(GFED_BA_NA[7,,])

# Save memory 
rm(FPA_BA_NA, GFED_BA_NA)


# ------------------------- spatialSubsets -------------------------------------
# All grids need to be subset to fit the CONUS or western US, depending on the
# region of interest. Lets start with US states. 

# Create large scale spatial domain masks! 
latMask <- GFED_latitude >= minLat & GFED_latitude <= maxLat 
lonMask <- GFED_longitude >= minLon & GFED_longitude <= maxLon

state_dict <- read.csv("Data/grid_attributes/state_dict.csv")
geographyMask <- (state != 0) & (state != 2) # no alaska or water/none
geographyMask[, GFED_latitude > 49] <- FALSE

# TODO: Make this variable!
ecoregionMask <- 0 != ecoregion 

# ------------------------- temporalSubsets ------------------------------------
# First make sure that the datasets time dimensions match. This will be limited
# by FPA-FOD when it is used, since it only goes up to 2013. 
timeMatch <- GFED_YYYYMM %in% FPA_YYYYMM
GFED_BA <- GFED_BA[timeMatch,,]

# Make the official time array for this scripts plotting and analysis
YYYYMM <- FPA_YYYYMM

# Get the year and month pieces of the time for subsetting 
YYYYMM_str <- as.character(YYYYMM)
mm         <- as.numeric(str_sub(YYYYMM_str, 5,6))
yyyy       <- as.numeric(str_sub(YYYYMM_str, 1,4))

# This will be a season or month mask application
timeMask <- mm %in% includedMonths

# ------------------------- subsets data check ---------------------------------
# make sure the dimensions allign across datasets 

GFED_BA <- GFED_BA[timeMask,,]
FPA_BA  <- FPA_BA[timeMask,,]
FPA_BA_human <- FPA_BA_human[timeMask,,]
FPA_BA_lightning <- FPA_BA_lightning[timeMask,,]

YYYYMM <- YYYYMM[timeMask]
mm     <- mm[timeMask]
yyyy   <- yyyy[timeMask]

nTime <- length(YYYYMM)

# ------------------------- conversion to acres ---------------------------------
# m**2 per acre 4046.86

GFED_BA <- GFED_BA / 4046.86
FPA_BA  <- FPA_BA / 4046.86
FPA_BA_human <- FPA_BA_human / 4046.86
FPA_BA_lightning <- FPA_BA_lightning /   4046.86

# ------------------------- create3DAttributeGrids -----------------------------
# Now, make a month and ecoregion mask to apply to each time dim. This will be
# nice labelling tools when it comes to making scatterplots. 

ecoregion_grid <- GFED_BA; ecoregion_grid[] <- NA
state_grid     <- ecoregion_grid

# ------------------------- monthlyBurnTotalsMaps ------------------------------
# Plot each grid, side by side, total burn area, for each month. 
# Store and compare. Do this for western CONUS. This means subsetting lat and 
# lon arrays

figureDir <- "Figures/burn_area_comparison/monthly/maps/"

latSubset <- GFED_latitude[latMask]
lonSubset <- GFED_longitude[lonMask]

# remove zero values, they are not interesting and clog the visuals. 
GFED_BA_NA <- GFED_BA[,lonMask, latMask]
FPA_BA_NA <- FPA_BA[,lonMask, latMask]
ecoregion_NA <- ecoregion[lonMask, latMask]

GFED_BA_NA[GFED_BA_NA==0] <- NA 
FPA_BA_NA[FPA_BA_NA==0]   <- NA

# Now flip the rows 
flipper <- length(latSubset):1

latSubset_flip <- latSubset[flipper]
GFED_BA_flip   <- GFED_BA_NA[,,flipper]
FPA_BA_NA_flip <- FPA_BA_NA[,,flipper]

# Create save-worthy image
for (i in 1:length(YYYYMM)){
  
  saveName <- paste0(figureDir, "GFED4s_FPA_FOD_",YYYYMM[i], ".png")
  png(filename=saveName, width=3500, height=1500, res=200)
  #quartz(width=5, height=5)
  par(mfrow=c(1,2), las=1, mar=c(4,4,4,8))
  
  # SETUP
  # shade burn area by these colors 
  logGFED <- log10(GFED_BA_flip[i,,])
  logFPA  <- log10(FPA_BA_NA_flip[i,,])
  
  # Take the log10 of the burn areas to setup a common colorbar for the two. 
  maxVal <- max(c(logGFED,logFPA), na.rm=T)
  minVal <- min(c(logGFED,logFPA), na.rm=T)
  
  zMax <- round(maxVal)
  zMin <- round(minVal)  
  x <- zMin:zMax
  xLabels <- 10^x
  
  fields::image.plot(lonSubset, latSubset_flip, logGFED, 
                     zlim=c(zMin, zMax),
                     axis.args=list( at=x, labels=xLabels, cex.axis = 1.5),
                     xlab="",
                     ylab=""
                     )
  map("state", add=T)
  map("world", add=T)
  title(paste("GFED4s:", YYYYMM[i]))
  
  fields::image.plot(lonSubset, latSubset_flip, logFPA, 
                     zlim=c(zMin, zMax),
                     axis.args=list( at=x, labels=xLabels, cex.axis = 1.5),
                     xlab="",
                     ylab=""
                     )
  
  map("state", add=T)
  map("world", add=T)
  title(paste("FPA-FOD:", YYYYMM[i]))
  
  dev.off()
  
}


# ----------------------------- mapTimespanSum ---------------------------------
# Plot total burn area for each! Over the whole time period! Show ecoregion
# perimeters? No, totals in each ecoregion! 


GFED_BA_flip_total <- apply(GFED_BA_flip, 2:3, sum, na.rm=T)
FPA_BA_flip_total  <- apply(FPA_BA_NA_flip, 2:3, sum, na.rm=T)

# remove zeros from plotting again
GFED_BA_flip_total[GFED_BA_flip_total==0] <- NA
FPA_BA_flip_total[FPA_BA_flip_total==0] <- NA

eco_mask <- as.character(SPDF$NA_L2CODE) == "6.2"
eco_6.2 <- SPDF[eco_mask,]

quartz(width=10, height=5)

png(filename="Figures/burn_area_comparison/summary/maps/western_US_2003_2013_burn_area.png", 
    res=200,
    width=3500, 
    height=1500)

par(mfrow=c(1,2), las=1, mar=c(4,4,4,9))

# TODO: Add ecoregions? 

# Label the legend. Span the orders of magnitude needed. 
x <- -2:6
acreLabel <- 10^x

image.plot(lonSubset, latSubset_flip, log10(GFED_BA_flip_total), 
           zlim=c(-2, 6),
           axis.args=list( at=x, labels=acreLabel, cex.axis = 1.5, bty="n"),
           xlab="", 
           ylab="",
           cex.axis=1.5,
           legend.width=3,
           cex.axis=1.5,
           xaxt="n",
           yaxt="n",
           bty="n"
           )
map("state", add=T, lty=1)
#plot(eco_6.2, add=T, border="black", lwd=3)
title(paste("GFED4s burn area 2003-2013"))

image.plot(lonSubset, latSubset_flip, log10(FPA_BA_flip_total), 
           zlim=c(-2, 6),
           axis.args=list( at=x, labels=acreLabel, cex.axis = 1.5),
           xlab="", 
           ylab="",
           cex.axis=1.5,
           legend.width=3,
           cex.axis=1.5,
           xaxt="n",
           yaxt="n",
           bty="n"
           )
map("state", add=T, lty=1)
#plot(eco_6.2, add=T, border="black", lwd=3)
title(paste("FPA-FOD burn area 2003-2013"))


dev.off()

################################################################################
# Show the totals again only this time only show region 6.2
################################################################################
nTime <- dim(GFED_BA_NA)[1]
GFED_BA_6.2 <- GFED_BA_NA
FPA_BA_6.2  <- FPA_BA_NA

spaceMask <- ecoregion_NA!=6.2

# Should plot blue over the area of the ecoregion polygon 
quartz()
image.plot(lonSubset, latSubset_flip, ecoregion_NA[,flipper])
plot(eco_6.2, add=T)

for (i in 1:nTime){
  # Plot this spatial mask to make sure that it incliudes southern colorado
  GFED_BA_6.2[i,,][spaceMask] <- NA 
  FPA_BA_6.2[i,,][spaceMask]  <- NA
}

# Flip
GFED_6.2_flip   <- GFED_BA_6.2[,,flipper]
FPA_BA_6.2_flip <- FPA_BA_6.2[,,flipper]

# Sum
GFED_6.2_total <- apply(GFED_6.2_flip, 2:3, sum, na.rm=T)
FPA_6.2_total  <- apply(FPA_BA_6.2_flip, 2:3, sum, na.rm=T)


png(filename="Figures/burn_area_comparison/summary/maps/western_US_2003_2013_6.2_burn_area.png", 
    res=200,
    width=3500, 
    height=1500)

par(mfrow=c(1,2), las=1, mar=c(4,4,4,9))

# Label the legend. Span the orders of magnitude needed. 
x <- -2:6
acreLabel <- 10^x

image.plot(lonSubset, latSubset_flip, log10(GFED_6.2_total), 
           zlim=c(-2, 6),
           axis.args=list( at=x, labels=acreLabel, cex.axis = 1.5, bty="n"),
           xlab="", 
           ylab="",
           cex.axis=1.5,
           legend.width=3,
           cex.axis=1.5,
           xaxt="n",
           yaxt="n",
           bty="n"
)
map("state", add=T, lty=3)
plot(eco_6.2, add=T, border="black", lwd=3)
title(paste("GFED4s burn area 2003-2013"))

image.plot(lonSubset, latSubset_flip, log10(FPA_6.2_total), 
           zlim=c(-2, 6),
           axis.args=list( at=x, labels=acreLabel, cex.axis = 1.5),
           xlab="", 
           ylab="",
           cex.axis=1.5,
           legend.width=3,
           cex.axis=1.5,
           xaxt="n",
           yaxt="n",
           bty="n"
)
map("state", add=T, lty=3)
plot(eco_6.2, add=T, border="black", lwd=3)
title(paste("FPA-FOD burn area 2003-2013"))


dev.off()

################################################################################
# Now show human ignition totals vs. non - human totals 
################################################################################
human     <- FPA_BA_human[,lonMask, latMask][,,flipper]
lightning <- FPA_BA_lightning[,lonMask, latMask][,,flipper]

human_sum <- apply(human, 2:3, sum, na.rm=T)
lightning_sum <- apply(lightning, 2:3, sum, na.rm=T)

png(filename="Figures/burn_area_comparison/summary/maps/western_US_2003_2013_6.2_burn_area_ignition.png", 
    res=200,
    width=3500, 
    height=1500)

par(mfrow=c(1,2), las=1, mar=c(4,4,4,9))

# Label the legend. Span the orders of magnitude needed. 
x <- -2:6
acreLabel <- 10^x

image.plot(lonSubset, latSubset_flip, log10(lightning_sum), 
           zlim=c(-2, 6),
           axis.args=list( at=x, labels=acreLabel, cex.axis = 1.5, bty="n"),
           xlab="", 
           ylab="",
           cex.axis=1.5,
           legend.width=3,
           cex.axis=1.5,
           xaxt="n",
           yaxt="n",
           bty="n"
)
map("state", add=T, lty=3)
plot(eco_6.2, add=T, border="black", lwd=3)
title(paste("Lightning burn area 2003-2013"))

image.plot(lonSubset, latSubset_flip, log10(human_sum), 
           zlim=c(-2, 6),
           axis.args=list( at=x, labels=acreLabel, cex.axis = 1.5),
           xlab="", 
           ylab="",
           cex.axis=1.5,
           legend.width=3,
           cex.axis=1.5,
           xaxt="n",
           yaxt="n",
           bty="n"
)
map("state", add=T, lty=3)
plot(eco_6.2, add=T, border="black", lwd=3)
title(paste("Human burn area 2003-2013"))


dev.off()


# Log scale makes for a tough comparison. Now try doing % 

png(filename="Figures/burn_area_comparison/summary/maps/percent_human_burn_area.png", 
    res=200,
    width=3500/2, 
    height=1500)

human_plus_lightning <- (human_sum + lightning_sum)
percentHuman <- human_sum / human_plus_lightning * 100

# Mask out Nans, this is where both numbers are zero
percentHuman[is.nan(percentHuman)] <- NA

# TODO: compare this addition sum to the overall value stored in the nc file 


image.plot(lonSubset, latSubset_flip, percentHuman, 
           #zlim=c(-2, 6),
           col=tim.colors(100),
           #axis.args=list( at=x, labels=acreLabel, cex.axis = 1.5),
           xlab="", 
           ylab="",
           cex.axis=1.5,
           legend.width=3,
           cex.axis=1.5,
           xaxt="n",
           yaxt="n",
           bty="n"
)
map("state", add=T, lwd=3)
plot(eco_6.2, add=T, border="white", lwd=3)

dev.off()

# ----------------------------- monthlyScatterPlots -----------------------------------

# Create domain mask
nLon <- length(GFED_longitude)
nLat <- length(GFED_latitude)
domainMask <- matrix(FALSE, nrow = nLon, ncol = nLat)
domainMask[lonMask, latMask] <- TRUE

# Create ecoregion mask, default is no unlabelled areas
ecoregionMask <- ecoregion != 0 # ecoregion == 6.2 

# Use both to create a spatial mask
spatialMask <- ecoregionMask & domainMask

# interannual variability, total acres burned western U.S. summer. 
years  <- sort(unique(yyyy))
months <- sort(unique(mm))
nYear  <- length(unique(yyyy))
nMonth <- length(unique(mm))

# storate arrays for the subsets of ignition burn area data we are interested in 
a <- rep(NA, nYear*nMonth)
GFED_month_total   <- a
FPAFOD_month_total <- a
human_month_total  <- a
lightning_month_total <- a

a_annual <- rep(NA, nYear)
GFED_summer_total      <- a_annual
FPAFOD_summer_total    <- a_annual
human_summer_total     <- a_annual
lightning_summer_total <- a_annual

# Sum monthly area burnedfor specified region
for (y in 1:nYear){
  
  # Count up the annual total (for given month range).
  # These lines are a bit tricky. First I sum over the time (year) direction
  # then apply a spatial mask to that 2D array of month sums for a year, then
  # sum the spatailly subset data.
  yearMask <- years[y] == yyyy
  
  GFED_summer_total[y]      <- sum(apply(GFED_BA[yearMask,,], 2:3, sum)[spatialMask])
  FPAFOD_summer_total[y]    <- sum(apply(FPA_BA[yearMask,,], 2:3, sum)[spatialMask])
  human_summer_total[y]     <- sum(apply(FPA_BA_human[yearMask,,], 2:3, sum)[spatialMask])
  lightning_summer_total[y] <- sum(apply(FPA_BA_lightning[yearMask,,], 2:3, sum)[spatialMask])
  
  # Save out the monthly totals too
  for (m in 1:nMonth){
    
    monthsMask <- years[y] == yyyy & months[m] == mm
    if(sum(monthsMask) > 1){
      stop("Time mask mistake.")
    }
    # Count up the monthlty total in specified spatial mask 
    GFED_month_total[monthsMask]   <- sum(GFED_BA[monthsMask,,][spatialMask])
    FPAFOD_month_total[monthsMask] <- sum(FPA_BA[monthsMask,,][spatialMask])
    human_month_total[monthsMask]  <- sum(FPA_BA_human[monthsMask,,][spatialMask])
    lightning_month_total[monthsMask] <- sum(FPA_BA_lightning[monthsMask,,][spatialMask])

  }
  
}


################################################################################
# plot monthly correlation between the monthly region totals
################################################################################

# Set up month color pallete 
# http://colorbrewer2.org/?type=qualitative&scheme=Set3&n=12

colorPallete <- c('#8dd3c7', '#33a02c','#bebada', '#fb8072', '#80b1d3', '#fdb462', 
                 '#b3de69','#fccde5','#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f')
monthColors <- colorPallete[1:12 %in% months]

# Month colors needs to be repped the number of years that are being used
monthColors <- rep(monthColors, nYear)

# Now set up year colors, no need to repeat 
yearColors  <- colorPallete[1:12 %in% 1:nYear]


figureDir <- "Figures/burn_area_comparison/monthly/scatterplots/"


plot_burn_area_ignition <- function(GFED_month_total=GFED_month_total, 
                                    FPAFOD_month_total=FPAFOD_month_total,
                                    human_month_total=human_month_total,
                                    lightning_month_total=lightning_month_total,
                                    dotLabel = dotLabel, # months | years
                                    dotColors=dotColors, # monthColors | yearColors
                                    dotTitle = "Month", # "Year"
                                    timeIntegral="Monthly", # "Monthly" or "Summer"
                                    fileName=fileName){


  png(filename=fileName, height=1000, width=2000, res=200)
  
  par(las=1, mar=c(5,5,5,5), cex=0.5, mfrow=c(1,2))
  
  # Create the plot frame that covers data maxes
  maxVal <- max(c(max(GFED_month_total), max(FPAFOD_month_total)))
  minVal <- 0 # Note, we may want to use data min, depending on how this looks 
  a <- c(minVal, maxVal)
  
  plot(a,a, 
       lwd=0.5, type="l", 
       xlab="GFED4s (acres)", 
       ylab="") 
  
  mtext("FPA-FOD (acres)", side=2, las=0, line=4)
  
  # All ecoregions all igntion types
  points(GFED_month_total, 
         FPAFOD_month_total, 
         col=dotColors, 
         pch=19) # Solid circle
  #points(GFED_month_total, FPAFOD_month_total)
  
  # Now show the same but for human ignitions 
  points(GFED_month_total, 
         human_month_total, 
         col=dotColors,
         pch=17,
         cex=1.5)
  
  # Now show the same but for Lightning ignitions 
  points(GFED_month_total, 
         lightning_month_total, 
         col=dotColors,
         pch=15)
  
  legend("topright",
         inset=c(-0.2,0),
         legend=dotLabel,
         col=unique(dotColors),
         pch=19,
         bty="n", 
         title=dotTitle,
         cex=0.7,
         xpd=TRUE
         )
  
  legend("top",
         inset=c(0,-0.1),
         legend=c("All FPA-FOD", "Human Ignitions", "Lightning Ignitions"),
         col="black",
         pch=c(19,17,15),
         bty="n",
         horiz = T,
         cex=0.5,
         xpd=T
  )
  
  title(paste("Western U.S.", timeIntegral,"Burn Area"))
  
  lim <- 0.25 * maxVal
  # Place lines showing where the subset plot will be
  lines(c(-1e36, lim), c(lim,lim), col="blue", xpd=F, lwd=2)
  lines(c(lim, lim), c(-1e36,lim), col="blue", xpd=F, lwd=2)
  
  
  # Inmap! Create inset map! 
  # https://gis.stackexchange.com/questions/222799/create-an-inset-map-in-r
  par(usr=c(-236, -83, 22, 130))
  rect(xleft = minLon, ybottom = minLat, xright = maxLon, ytop = maxLat, col = "white")
  map("state", xlim=c(minLon, maxLon), ylim=c(minLat, maxLat), add=T, boundary = T, lty=1)
  
  # Inset 
  aa <- c(0,lim)
  plot(aa, aa, 
       lwd=0.5, type="l", 
       xlab="GFED4s (acres)", 
       ylab="",
       axes=F) 
  box(col = 'blue')
  axis(1, col="blue")
  axis(2, col="blue")
  
  mtext("FPA-FOD (acres)", side=2, las=0, line=4)
  
  # All ecoregions all igntion types
  points(GFED_month_total, 
         FPAFOD_month_total, 
         col=monthColors, 
         pch=19) # Solid circle
  title("Inset", col="blue")
  
  # All ecoregions all igntion types
  points(GFED_month_total, 
         FPAFOD_month_total, 
         col=monthColors, 
         pch=19) # Solid circle
  
  # Now show the same but for human ignitions 
  points(GFED_month_total, 
         human_month_total, 
         col=monthColors,
         pch=17,
         cex=1.5)
  
  # Now show the same but for Lightning ignitions 
  points(GFED_month_total, 
         lightning_month_total, 
         col=monthColors,
         pch=15)
  
  legend("topright",
         inset=c(-0.2,0),
         legend=dotLabel,
         col=unique(monthColors),
         pch=19,
         bty="n", 
         title=dotTitle,
         cex=0.7, 
         xpd=TRUE
  )
  
  legend("top",
         inset=c(0,-0.1),
         legend=c("All FPA-FOD", "Human Ignitions", "Lightning Ignitions"),
         col="black",
         pch=c(19,17,15),
         bty="n",
         horiz = T,
         cex=0.5,
         xpd=T
  )
  
  
  dev.off()

}

# Make annual plot
fileName <- paste0(figureDir, "western_US_burn_area_by_ignition_summer_",monthSpan,".png")
timeIntegral <- "Summer"
plot_burn_area_ignition(GFED_summer_total,
                        FPAFOD_summer_total,
                        human_summer_total,
                        lightning_summer_total,
                        dotLabel = years,
                        dotColors = yearColors,
                        dotTitle = "Year",
                        timeIntegral = "Summer",
                        fileName=fileName)

# Make monthly plot
fileName <- paste0(figureDir, "western_US_burn_area_by_ignition_",monthSpan,".png")
timeIntegral <- "Monthly"
plot_burn_area_ignition(GFED_month_total,
                        FPAFOD_month_total,
                        human_month_total,
                        lightning_month_total,
                        dotLabel = months,
                        dotColors = monthColors,
                        dotTitle = "Month",
                        timeIntegral = "Monthly",
                        fileName=fileName)


#------------------------------- ecoRegionScatter ------------------------------

# TODO: Definitly don't need to go through ecoregions that do not show up 
# outside wesdtern US 

# make a scatter plot of the relation for EVERY ecoregion. 
figureDir <- "Figures/burn_area_comparison/monthly/scatterplots/"


ecoregionLabels <- sort(unique(as.numeric(ecoregion[spatialMask])))

# Get the burn area for both products in the domain, and make a scatterplot 
for (eco in ecoregionLabels){
  
  # Sum monthly area burnedfor specified region
  eco_GFED <- rep(NA, nYear*nMonth)
  eco_FPA  <- rep(NA, nYear*nMonth)
  
  for (y in 1:nYear){
    for (m in 1:nMonth){
      
      tMask <- years[y] == yyyy & months[m] == mm
      eco_GFED[tMask] <- sum(GFED_BA[tMask,,][spatialMask & ecoregion == eco])    
      eco_FPA[tMask] <- sum(FPA_BA[tMask,,][spatialMask & ecoregion == eco])
      
    }
  }
  
  # You have data for selected ecoregion, get to plotting 
  fileName <- paste0(figureDir, "western_US_ecoregion_",eco,".png")
  png(filename=fileName, height=1000, width=2000, res=200)
  par(las=1, mar=c(5,5,5,5), cex=0.5, fig = c(0,1,0,1), mfrow=c(1,2) )
  
  # Create the plot frame that covers data maxes
  maxVal <- max(c(max(eco_GFED), max(eco_FPA)))
  minVal <- 0 # Note, we may want to use data min, depending on how this looks 
  a <- c(minVal, maxVal)
  
  plot(a,a, 
       lwd=0.5, type="l", 
       xlab="GFED4s (acres)", 
       ylab="") 
  
  mtext("FPA-FOD (acres)", side=2, las=0, line=4)
  
  # All ecoregions all igntion types
  points(eco_GFED, 
         eco_FPA, 
         col=monthColors, 
         pch=19,
         cex=1.5) # Solid circle
  points(eco_GFED,eco_FPA, lwd=0.3, cex=1.5)
  
  title(paste("Monthly Burn Area Ecoregion", eco))

  # I want a month legend, outside right!
  legend("topright", 
         inset=c(-0.2,0),
         legend=months,
         pch=19,
         col=monthColors,
         title="Month",
         xpd=TRUE,
         bty="n"
         )
  
  # Map of the plotted and shaded region
  ecoregion_highlight <- (ecoregion==eco & spatialMask)[,nLat:1]
  ecoregion_highlight[ecoregion_highlight==0] <- NA
  
  # Now make the same figure, only an inset of the lower left hand corner. This
  # will be the boxed portion of the previous figure.
  lim <- 0.25 * maxVal
  # Place lines showing where the subset plot will be
  lines(c(-1e36, lim), c(lim,lim), col="blue", xpd=F)
  lines(c(lim, lim), c(-1e36,lim), col="blue", xpd=F)
  

  
  # Set up map inset figure space 
  par(usr=c(-236, -83, 22, 130))
  rect(xleft = minLon, ybottom = minLat, xright = maxLon, ytop = maxLat, col = "white")
  map("state", xlim=c(minLon, maxLon), ylim=c(minLat, maxLat), add=T, boundary = T, lty=1)
  image(GFED_longitude, GFED_latitude[nLat:1], ecoregion_highlight, add=T, col="forestgreen")
  
  # TODO: inset of lower left part of plot with arrows indicating that is what 
  # TODO: is being shown
  
  # Inset panal right
  #par(new=T)
  aa <- c(0,lim)
  plot(aa, aa, 
       lwd=0.5, type="l", 
       xlab="GFED4s (acres)", 
       ylab="",
       axes=F) 
  box(col = 'blue')
  axis(1, col="blue")
  axis(2, col="blue")
  
  mtext("FPA-FOD (acres)", side=2, las=0, line=4)
  
  # All ecoregions all igntion types
  points(eco_GFED, 
         eco_FPA, 
         col=monthColors, 
         pch=19,
         cex=1.5) # Solid circle
  points(eco_GFED,eco_FPA, lwd=0.3, cex=1.5)
  title("Inset", col="blue")
  
  dev.off()
  
 
  
  
}



#------------------------------- gridPointScatter ------------------------------

# Point for point comparison of inventories. 
# Now do the spatial subset across the time dimension. We will lose all spatial
# information, we will only know that we are looking at the region and ecoregion
# selected.

nSpatialPoints <- sum(spatialMask)
GFED_BA_ <- numeric(0)
FPA_BA_ <- numeric(0)
FPA_BA_human_ <- numeric(0)
FPA_BA_lightning_ <- numeric(0)

# Subset the spatial dimension one time index at a time
# NOTE: THERE HAS TO BE A BETTER WAY TO DO THIS WITH apply() 
# NOTE: These data will be flattened to 1D arrays
for (i in 1:dim(GFED_BA)[1]){
  
  # These are flat 
  GFED_BA_ <- append(GFED_BA_, GFED_BA[i,,][spatialMask])
  FPA_BA_ <- append(FPA_BA_, FPA_BA[i,,][spatialMask])
  FPA_BA_human_ <- append(FPA_BA_human_, FPA_BA_human[i,,][spatialMask])
  FPA_BA_lightning_ <- append(FPA_BA_lightning_, FPA_BA_lightning[i,,][spatialMask])
  
}

png(filename=paste0(figureDir, "gridPoints_US.png"), height=800, width=800)
par(las=1, mar=c(5,5,5,5), cex=2)

plot(GFED_BA_,
     GFED_BA_, 
     type="l", xlab="GFED4s (acres)", ylab="")
points(GFED_BA_, FPA_BA_, pch=19)

mtext("FPA-FOD (acres)", side=2, las=0, line=4, cex=2)
title("Western U.S. Monthly grid point burn area comparison")

# Inmap! Create inset map! 
# https://gis.stackexchange.com/questions/222799/create-an-inset-map-in-r
par(usr=c(-216-20, -63-20, 22, 144))
rect(xleft =-126.2, ybottom = 30, xright = -100, ytop = 50.1, col = "white")
map("state", xlim=c(-126.2,-100), ylim=c(30,50.1), add=T, boundary = T, lty=1,
    main="plotted area")

dev.off()




