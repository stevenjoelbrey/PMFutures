# plot_interannual_variability_vs_met.R

# ------------------------- Description ---------------------------------------
# This script will be used to plot regional interannual variability in burn area
# segregated by ecoregion and other grid attributes. 

# The main goal is to show if weather variables explain variance in interannual 
# variability and if relationships differ between ignition types. This is done 
# using ecmwf reanalysis fields and FPA-FOD fire occurance data.
# Timescales of interest include season (e.g. 5-10), monthly, and lagged monthly. 

# TODO: Download ecmwf reanalysis on monthly scale, or average what I have already

library(stringr)
library(maps)
library(ncdf4)
library(fields)
library(sfsmisc)
library(lubridate) # for month()

# What region are you investigating? 
minLat <- 31
maxLat <- 49
minLon <- -125
maxLon <- -100

year1 <- 1992
year2 <- 2015 # extent of current FPA-FOD, will be 2015 soon. 
ecoregion_select <- 11.1
state_select <- 7 #c("Washington", "Oregon", "Montana") # coming soon
month_select  <- 5:12

years  <- year1:year2
nYears <- length(years)

################################################################################
# Save figures in directory specific to how the data are subset
################################################################################
experimentDir <- paste0("eco=", ecoregion_select, "_", 
                        "months=", min(month_select),"_", max(month_select), 
                        "/" )
figureDir     <- paste0("Figures/regional_met_relations/", experimentDir)

if(!dir.exists(figureDir)){
  dir.create(figureDir)
} else{
  print("Saving figures into an already existing directory")
}


# Load FPA-FOD data, the one with all attributes (ecoregions) assgined.
# TODO: load the one with weather and make sure they match! 
load(paste0("Data/FPA_FOD/FPA_FOD_ecmwf_", year1,"_", year2,".RData"))

# Get fire parameters too be used in subsetting the data 
fireDate      <- FPA_FOD$DISCOVERY_DATE
fireLat       <- FPA_FOD$LATITUDE
fireLon       <- FPA_FOD$LONGITUDE
fireYear      <- FPA_FOD$FIRE_YEAR
fireEcoregion <- FPA_FOD$NA_L2CODE
fireCause     <- FPA_FOD$STAT_CAUSE_DESCR
fireMonth     <- FPA_FOD$START_MONTH      
fireSize      <- FPA_FOD$FIRE_SIZE

# Handy arrays
humanStart <- fireCause != "Lightning"
nFires     <- length(fireSize)

# Load ecoregion borders, for plotting, analysis etc. 
load("Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData") # <- SPDF
ecoregion_polygon <- SPDF[SPDF@data$NA_L2CODE==ecoregion_select,]

#################################################################
# create a fire size attribute bin that will be used for plotting 
#################################################################
sizeBins <- c(0, 10^c(1:6)) # a clean span of 6 orders of magnitude 
cexBins  <- c(0, 0.1, 0.4, 1, 2, 3.8)
nBins <- length(sizeBins)
sizeClass <- rep(NA, nFires)

# Assign the fires to size bins
for (b in 1:(nBins-1)){
  
  # Mask them 
  m <- fireSize >= sizeBins[b] & fireSize < sizeBins[b+1]
  sizeClass[m] <- cexBins[b]
  
}

# TODO: Make this less terrible... Use scientific notation
legendText <- c("(10 - 99]", "(100 - 1,000]", 
                "(1,000 - 10,000]", "(10,000 - 100,000]", 
                "(100,000 - 1,100,000]")

# Create annual totals for specified ecoregion and months 
FPA_BA_lightning <- rep(NA, nYears)
FPA_BA_human     <- rep(NA, nYears)

# Keep track of monthly totals too. ( "_m" indicates monthly variable )
nMonths <- nYears * length(month_select)
monthlyTimeArray   <- rep(as.POSIXct("1750-04-26", tz="UTC"), nMonths)
FPA_BA_lightning_m <- rep(NA, nMonths)
FPA_BA_human_m     <- rep(NA, nMonths)

monthCount <- 0 # Advance the month counter 
for (i in 1:nYears){
  
  yearMask <- years[i] == fireYear
  monthMask <- fireMonth %in% month_select
  ecoRegionMask <-  fireEcoregion %in% ecoregion_select
  #stateMask <- fire_state
  m <- yearMask & monthMask & ecoRegionMask
  
  # Sum to burn area for the mask made above! 
  FPA_BA_human[i]     <- sum(fireSize[m & humanStart])
  FPA_BA_lightning[i] <- sum(fireSize[m & !humanStart])
  
  # Keep track of monthly totals too. Go within the year
  for (m in month_select){

    # Advance the month index 
    monthCount <- monthCount + 1

    # Make a nice time array for plotting monthly data
    if(m < 10){
      mm <- paste0("0", m)
    } else{
      mm <- as.character(m)
    }
    dateString <- paste0(years[i], "-", mm, "-15")
    monthlyTimeArray[monthCount] <- as.POSIXct(dateString, tz="UTC")
    
    # Mask isloated month in year timeframe    
    oneMonthMask <- fireMonth == m
    m2 <- oneMonthMask & yearMask & ecoRegionMask
    
    FPA_BA_human_m[monthCount]     <- sum(fireSize[m2 & humanStart])
    FPA_BA_lightning_m[monthCount] <- sum(fireSize[m2 & !humanStart])
    
  }
  
}

################################################################################
# Plot Monthly time series of burn area
################################################################################
if(FALSE){
png(filename=paste0(figureDir,"FPA_FOD_monthly_timeSeries_",
                ecoregion_select,".png"),
    height=2000, width=3000, res=250)

par(mar=c(4,8,4,4), lty=1, cex=2)

# Set the limits on the y-axis for the plot so we can create a correct sized
# blank. 
yMin <- min( c(FPA_BA_human_m, FPA_BA_lightning_m) )
yMax <- max( c(FPA_BA_human_m, FPA_BA_lightning_m) )

# Create the blank space (baby, I'll right your name)
plot(monthlyTimeArray, FPA_BA_lightning_m, col="gray", 
     bty="n", yaxt="n", xaxt="n",
     ylab="", xlab="", cex=1,
     ylim = c(yMin, yMax),
     xlim = c(min(monthlyTimeArray), max(monthlyTimeArray)),
     pch="")
lines(monthlyTimeArray, FPA_BA_lightning_m, col="gray")

# Label the x-axis
BY <- seq(1, length(monthlyTimeArray), by = 12)
axis(1, at = monthlyTimeArray[BY], labels = years)


# vertical axis
eaxis(2)
mtext("Acres Burned", side=2, line=5, cex=2)

# Human
#points(monthlyTimeArray, FPA_BA_human_m, col="orange", pch=1)
lines(monthlyTimeArray, FPA_BA_human_m, col="orange")

title(paste("Monthly Burn area in ecoregion", ecoregion_select),
      cex.main=2)


legend("topleft",
       bty="n",
       legend=c("Human", "Lightning"),
       pch= 19,
       col=c("Orange", "gray"),
       cex=1.3
       )

dev.off()
}

################################################################################
# Plot Interannual variability over time for selected months 
################################################################################
png(filename=paste0(figureDir,"FPA_FOD_interannual_variability_mapped_",
                    ecoregion_select,
                    ".png"),
    height=2000, width=4800, res=250)

par(mfrow=c(1,2), xpd=T, mar=c(4,10,4,0))

maxValue <- max(c(FPA_BA_lightning/10^6, FPA_BA_human/10^6))


plot(years, FPA_BA_lightning/10^6, col="gray", 
     pch=19, bty="n", yaxt="n", xaxt="n",
     ylab="", xlab="", cex=2,
     ylim=c(0, maxValue)) 
lines(years, FPA_BA_lightning/10^6, col="gray", lty=2)

# Add axis to the plot 
eaxis(1, cex.axis=2)
eaxis(2, sub10=5, cex.axis=2)
mtext("Acres Burned (millions)", side=2, line=4, cex=2)

points(years, FPA_BA_human/10^6,pch=19, col="orange", cex=2)
lines(years, FPA_BA_human/10^6, col="orange", lty=2)


legend("topleft",
       legend=c("Lightning", "Human"),
       pch=19,
       pt.cex=2,
       col=c("gray", "orange"),
       cex=2,
       bty = "n"
       )

# Map the fires in the time series
# TODO: make sure to plot the fires in ascending size! That way we can see small
# TODO: fires next to big fires. 

par(xpd=F, mar=c(4,4,4,4))

monthMask      <- fireMonth %in% month_select
ecoRegionMask  <- fireEcoregion == ecoregion_select

humanLocations <- monthMask & ecoRegionMask & humanStart

plot(FPA_FOD$LONGITUDE[humanLocations], FPA_FOD$LATITUDE[humanLocations], 
     cex=sizeClass[humanLocations],
     col=adjustcolor("orange", 0.5), pch=1,
     bty="n", xaxt="n", yaxt="n",
     ylab="",
     xlab="")

map("state", add=T)

lightningLocations <- monthMask & ecoRegionMask & !humanStart

points(FPA_FOD$LONGITUDE[lightningLocations], FPA_FOD$LATITUDE[lightningLocations], 
     cex=sizeClass[lightningLocations],
     col=adjustcolor("gray",0.5), pch=1,
     bty="n", xaxt="n", yaxt="n",
     ylab="",
     xlab="")

# Show the ecoregion border 
plot(ecoregion_polygon, add=T, lty=1, border="lightgray")

legend("bottomleft", bty="n",
       #title="Acres",
       inset=c(0,0),
       xpd=TRUE,
       legendText,
       pt.cex=cexBins[2:(nBins-1)],
       pch=1,
       col="black",
       horiz=F)

dev.off()

# NOTE: This only works since all fires are treated as points, and only exist
# NOTE: in the western hemisphere.
fireLonAdjusted <- fireLon + 360

################################################################################
# Get corrosponding environmental data to make seasonal comparisons for
# specific variables.
# ERA-interim ecmwf averaging definitions can be found at the link below:
# https://software.ecmwf.int/wiki/display/CKB/ERA-Interim%3A+monthly+means
################################################################################

# These are nice 0:360 lon bounds that include what we need from North America 
# when loading ecmwf data 
minLon     <- 230.0
maxLon     <- 263

################################################################################
# Apply ecmwf gridded ecoregion mask. Want to only look at meteorology in 
# locations where we are looking at fire! 
################################################################################
gridAttributesFile <- "Data/grid_attributes/grid_attributes_75x75.nc"
grid_attributes    <- nc_open(gridAttributesFile)

# Spatial first
grid_latitude <- grid_attributes$dim$latitude$vals
grid_longitude <- grid_attributes$dim$longitude$vals

# Actual attributes next 
grid_ecoregion     <- round(ncvar_get(grid_attributes, "ecoregion"),2)
grid_state         <- ncvar_get(grid_attributes, "state")
grid_ecoregion[is.na(grid_ecoregion)] <- -1
nc_close(grid_attributes)

# Get rid of pesky NA values, they make it hard to make masks 

grid_latMask       <- grid_latitude > minLat & grid_latitude < maxLat
grid_lonMask       <- grid_longitude > minLon & grid_longitude < maxLon

lati <- min(which(grid_latMask)); 
lat_count <- max(which(grid_latMask)) - lati

loni <- min(which(grid_lonMask));
lon_count <- max(which(grid_lonMask)) - loni

# Location of local nc data batch
ncDir <- "/Volumes/Brey_external/era_interim_nc_daily_merged/"

# The gridded array are very large. So we only want to extact exactly what
# we need. This needs to be done for time, lat, and lon. 
nc_file <- paste0(ncDir, "t2m_", year1, "_", 2016, ".nc")
nc      <- nc_open(nc_file)
nc_time <- ncvar_get(nc, "time") # hours since 1900-01-01 00:00:0.0
nc_close(nc)

# create time extraction array
t <- as.POSIXct(nc_time*60^2, origin="1900-01-01 00:00:0.0", tz="utc")
tMask <- year(t) <= 2015
ti <- 1 # first extact index
tf <- max(which(tMask==1))

# Get month and year arrays from ecmwf grid. Also subset t such that it only 
# goes through the desired year 
t <- t[tMask]
tMon <- month(t)
tYear<- year(t)

# Now actually load the data
nc_file <- paste0(ncDir, "t2m_", year1, "_", 2016, ".nc")
nc      <- nc_open(nc_file)
t2m <- ncvar_get(nc, "t2m", start=c(loni, lati, ti), count=c(lon_count, lat_count, tf))
ecmwf_latitude  <- ncvar_get(nc, "latitude", start=lati, count=lat_count)
ecmwf_longitude <- ncvar_get(nc, "longitude", start=loni, count=lon_count)

# t2m <- ncvar_get(nc, "t2m", start=c(1, 1, 180), count=c(length(grid_longitude), length(grid_latitude), 1))
# ecmwf_latitude  <- ncvar_get(nc, "latitude")
# ecmwf_longitude <- ncvar_get(nc, "longitude")
nc_close(nc)

# Make sure these coordinates match the fire coordinates by ensuring that
# fires are plotted over this western US snapshot
flipper <- length(ecmwf_latitude):1
quartz(width=8, height=5)
image.plot(ecmwf_longitude, ecmwf_latitude[flipper], t2m[,flipper, 180])
# Add fire
points(fireLonAdjusted, fireLat, pch=".")
map("state", add=T)

#Function for loading files
getMetVar <- function(VAR,loni, lati, ti, lon_count, lat_count, tf){
  nc      <- nc_open(paste0(ncDir, VAR, "_", year1, "_", 2016, ".nc"))
  VAR_ <- ncvar_get(nc, VAR, start=c(loni, lati, ti), count=c(lon_count, lat_count, tf))
  nc_close(nc)
  return(VAR_)
}

print("Loading really big met files into workspace")
tp <- getMetVar("tp", loni, lati, ti, lon_count, lat_count, tf)
rh2m <- getMetVar("rh2m", loni, lati, ti, lon_count, lat_count, tf)
e <- getMetVar("e", loni, lati, ti, lon_count, lat_count, tf)
d2m <- getMetVar("d2m", loni, lati, ti, lon_count, lat_count, tf)
  
# wind speed
u10 <- getMetVar("u10", loni, lati, ti, lon_count, lat_count, tf)
v10 <- getMetVar("v10", loni, lati, ti, lon_count, lat_count, tf)
ws <- sqrt(u10*u10 + v10*v10)
rm(v10, u10)

# wind gust 
wg <- getMetVar("fg10", loni, lati, ti, lon_count, lat_count, tf)


# Now we need the attributes to be subset the same way spatially
grid_lat_keeps <- grid_latitude %in% ecmwf_latitude
grid_lon_keeps <- grid_longitude %in% ecmwf_longitude

grid_state_subset     <- grid_state[grid_lon_keeps, grid_lat_keeps]
grid_ecoregion_subset <- grid_ecoregion[grid_lon_keeps, grid_lat_keeps]

grid_ecoregion_mask <- grid_ecoregion_subset == ecoregion_select

# Make summer totals / means depending on the paramter
v <- rep(NA, nYears)
t2m_mean <- v
tp_total <- v
e_total  <- v
rh2m_mean<- v 
d2m_mean <- v
ws_mean  <- v
wg_mean  <- v

# This fuction handles the multi-step hasle of subsetting in both time and
# a non regular grid space 
spaceTimeStat <- function(x, tMask, FUN="mean"){
  
  # Time subset first 
  spatialSubset <- x[,, tMask]
  # stats on time (called mean but could be "sum")
  gridMeans <- apply(spatialSubset, 1:2, FUN)
  # Spatial subset
  # TODO: Think about weather mean of precip grid box totals in more valueble
  # TODO: than the total precip of all grid boxes. For now going with mean. 
  spatialMean <- mean(gridMeans[grid_ecoregion_mask])
  
  return(spatialMean)
  
}

# TODO: Make a spatial correlation version to be plotted! Make correlations
# TODO: on the domain grid, then plot white over non-ecoregion or states 
# TODO: desired later. 

# TODO: Save out monthly relationships 
for (i in 1:nYears){

  # Time masks 
  yearMask  <- years[i] == tYear
  monthMask <- tMon %in% month_select
  tMask     <- yearMask & monthMask 
  
  # Take the spatial temporal subset average or sum and store 
  t2m_mean[i] <- spaceTimeStat(t2m, tMask, FUN="mean")
  tp_total[i] <- spaceTimeStat(tp, tMask, FUN="sum")
  e_total[i]  <- spaceTimeStat(e, tMask, FUN="sum")
  rh2m_mean[i]<- spaceTimeStat(rh2m, tMask, FUN="mean")  
  d2m_mean[i] <- spaceTimeStat(d2m, tMask, FUN="mean")  
  ws_mean[i]  <- spaceTimeStat(ws, tMask, FUN="mean") 
  wg_mean[i]  <- spaceTimeStat(wg, tMask, FUN="mean")
  
}

################################################################################
# Plot the relationships 
################################################################################
# TODO: Save "experiment" directories that document time and spatial extent 
# TODO: correlations being plotted. 
scatter_plot <- function(x=t2m_mean, y1=FPA_BA_lightning, y2=FPA_BA_human,
                         xLab="", yLab=""){

  f <- paste0(figureDir, xLab, "_vs_", yLab, ".png")

  yLab <- str_replace(yLab, " ", "\n")
  
  png(filename=f, width=1700, height=1400, res=300)
  par(mar=c(4,8,4,4))
  
  plot(x, y1, pch=19, cex=1.8,
       col="gray", bty="n", yaxt="n",
       ylab="", xlab="")
  mtext(xLab, side=1, line = 2)
  mtext(yLab, side=2, line = 4.5, las=1)
  
  minVal <- min(c(y1,y2)); maxVal <- max(c(y1,y2))
  ylabValues <- seq(minVal, maxVal, length.out = 5)
  Labels <- pretty10exp(ylabValues, digits=2)
  eaxis(side=2, at=ylabValues, labels=Labels)
  
  
  points(x, y2, col="orange", pch=19, cex=1.4)
  
  dev.off()

}


scatter_plot(x=t2m_mean, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="Temperature", yLab="Acres Burned")
scatter_plot(x=tp_total, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="Precip", yLab="Acres Burned")
scatter_plot(x=rh2m_mean, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="Relative Humidity", yLab="Acres Burned")
scatter_plot(x=d2m_mean, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="Dew Point Temperature", yLab="Acres Burned")
scatter_plot(x=ws_mean, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="Wind Speed",  yLab="Acres Burned")
scatter_plot(x=wg_mean, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="Wind Gust",  yLab="Acres Burned")


################################################################################
# These relationships should happen on a grid by grid basis. Show correlations
# between these grid boxes burn area. Call upon monthly gridded burn area
# created by grid_FPA_FOD_burn_area.R
################################################################################

# Here is where you throw out all of the spatially subsetting you previously
# did! Just plot your desired spatial domain. 

ncFile <- "Data/FPA_FOD/burn_area_monthly_75x75_1992_2015.nc"
nc <- nc_open(ncFile)
grid_area <- ncvar_get(nc, "grid_area")
longitude <- ncvar_get(nc, "longitude")
latitude  <- ncvar_get(nc, "latitude")
burn_area <- ncvar_get(nc, "burn_area")
t_hours   <- ncvar_get(nc, "time")
nc_close(nc)

# Time in hours from origin to POSIXct
t_mon <- as.POSIXct(t_hours*60^2, origin="1900-01-01 00:00:0.0", tz="utc")

# Spatially subset the data to match the loaded ecmwf data
burn_area <- burn_area[, longitude %in% ecmwf_longitude, latitude %in% ecmwf_latitude]

# Plot the spatial total as a sanity check
burn_area_total <- apply(burn_area, 2:3, sum)

# NOTE: Note the shifted coordinates so I can plot the U.S. map in this sanity
# NOTE: check plot 
quartz()
f <- length(ecmwf_latitude):1
image.plot(ecmwf_longitude-360, ecmwf_latitude[f], burn_area_total[,f])
map("state", add=T)

################################################################################
# Create Monthly ecmwf means and totals. Remember, pearson correlations in case
# the relations are not linear. 
################################################################################

# t2m subset by t made up of tMon and tYear
# TODO: 
nMonths <- length(unique(tYear)) * length(unique(tMon))
monthly <- array(NA, dim = c(length(ecmwf_longitude), 
                             length(ecmwf_latitude), 
                             nMonths))

# Give monthly array dimensions to all variables 
t2m_monthly <- monthly
tp_monthly <- monthly
e_monthly  <- monthly
rh2m_monthly<- monthly 
d2m_monthly <- monthly
ws_monthly  <- monthly
wg_monthly  <- monthly

i <- 0
for (y in sort(unique(tYear))) {
  
  for (m in sort(unique(tMon))){
    i <- i + 1
    tMask <- y == tYear & m == tMon
    
    #print(paste( min(t[tMask]) ,"-", max(t[tMask])))
    
    t2m_monthly[,,i] <- apply(t2m[,,tMask], 1:2, mean)
    tp_monthly[,,i] <- apply(tp[,,tMask], 1:2, sum)
    e_monthly[,,i] <- apply(e[,,tMask], 1:2, sum)
    rh2m_monthly[,,i] <- apply(rh2m[,,tMask], 1:2, mean)
    d2m_monthly[,,i] <- apply(d2m[,,tMask], 1:2, mean)
    ws_monthly[,,i] <- apply(ws[,,tMask], 1:2, mean)
    wg_monthly[,,i] <- apply(wg[,,tMask], 1:2, mean)
    
  }
  
}

# Make spatial correlations with gridded burn area and display as sketched in 
# notebook. 





