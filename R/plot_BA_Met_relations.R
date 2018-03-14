# plot_BA_Met_relations.R

# ------------------------- Description ---------------------------------------
# This script will be used to compare the interannual variability of 
# environmental variables and burn area, segregated by ignition type. 

library(stringr)
library(maps)
library(ncdf4)
library(fields)
library(sfsmisc)
library(lubridate) # for month()

# What region are you investigating?
regionName <- "southeast" # west

if (regionName == "west"){
  
  minLat <- 31.2
  maxLat <- 49
  minLon <- -125
  maxLon <- -104
  
} else if(regionName == "southeast"){
  
  minLat <- 25 
  maxLat <- 39 
  minLon <- -91   
  maxLon <- -75.5  
  
}

# TODO: Assign ecoregion nice names to all FPA_FOD df 
year1 <- 1992
year2 <- 2015 
ecoregion_select <- 11.1
ecoregion_name   <- "Mediterranean California" #"Mississippi Alluvial"#"Ozark" #"Southern Plains" # "Mediterranean California" "Forested mountains" "High deserts"
month_select     <- 1:12 # THIS MAY BE VERY WRONG FOR HUMAN and it does vary by region. 
minSize <- 0     # acres
maxSize <- Inf   # acres
years  <- year1:year2
nYears <- length(years)

print(paste("Working on:", ecoregion_select, ecoregion_name,"minSize:",minSize))

################################################################################
# Save figures in directory specific to how the data are subset
################################################################################
experimentDir <- paste0("eco=", ecoregion_select, "_", 
                        "months=", min(month_select),"_", max(month_select), "_",
                        "minSize=", minSize, "_maxSize=", maxSize, 
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

# Get rid of fires that have "Missing/Undefined" start types, cleaner analysis
# with fewer assumptions. 
HasStartInfo <- FPA_FOD$STAT_CAUSE_DESCR != "Missing/Undefined"
bigEnough    <- FPA_FOD$FIRE_SIZE >= minSize
smallEnough  <- FPA_FOD$FIRE_SIZE <= maxSize
FPA_month_mask <- lubridate::month(FPA_FOD$DISCOVERY_DATE) %in% month_select
FPA_ecoregion_mask <- FPA_FOD$NA_L2CODE %in% ecoregion_select

m <- HasStartInfo & bigEnough & smallEnough & FPA_month_mask & FPA_ecoregion_mask

# Subset the data
FPA_FOD <- FPA_FOD[m,]

print("The dim of subset FPA_FOD is:")
print(dim(FPA_FOD))

# Get fire parameters too be used in subsetting the data 
fireDate      <- FPA_FOD$DISCOVERY_DATE
fireLat       <- FPA_FOD$LATITUDE
fireLon       <- FPA_FOD$LONGITUDE
fireYear      <- lubridate::year(fireDate)
fireEcoregion <- FPA_FOD$NA_L2CODE
fireCause     <- FPA_FOD$STAT_CAUSE_DESCR
fireMonth     <- lubridate::month(fireDate)      
fireSize      <- FPA_FOD$FIRE_SIZE

# Handy arrays
humanStart <- fireCause != "Lightning" # because "Missing" have been removed
nFires     <- length(fireSize)

# Create annual sums for lightning and human-started fire burn area. Exclude
# fires less than "minSize"
FPA_BA_lightning <- rep(NA, nYears)
FPA_BA_human     <- rep(NA, nYears)

for (i in 1:nYears){
  yMask <- years[i] == fireYear
  FPA_BA_lightning[i] <- sum(FPA_FOD$FIRE_SIZE[(fireCause == "Lightning") & yMask])
  FPA_BA_human[i] <- sum(FPA_FOD$FIRE_SIZE[(fireCause != "Lightning") & yMask])
}

print(paste("Total FPA_BA_lightning:", sum(FPA_BA_lightning)))
print(paste("Total FPA_BA_human:", sum(FPA_BA_human)))

# Load ecoregion borders, for plotting, analysis etc. 
load("Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData") # <- SPDF
ecoregion_polygon <- SPDF[SPDF@data$NA_L2CODE==ecoregion_select,]

# NOTE: This only works since all fires are treated as points, and only exist
# NOTE: in the western hemisphere.
fireLonAdjusted <- fireLon + 360

################################################################################
# Get corrosponding environmental data to make seasonal comparisons for
# specific variables.
# ERA-interim ecmwf averaging definitions can be found at the link below:
# https://software.ecmwf.int/wiki/display/CKB/ERA-Interim%3A+monthly+means
################################################################################

# These are nice 0:360 lon bounds that include what we need from CONUS 
# when loading ecmwf data 
minLon     <- 230.0
maxLon     <- 290

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

# Where in grid are the indicies that cover my limits?
lati <- min(which(grid_latMask)); 
lat_count <- max(which(grid_latMask)) - lati

loni <- min(which(grid_lonMask));
lon_count <- max(which(grid_lonMask)) - loni

# Location of local nc data batch
ncDir <- "/Volumes/Brey_external/era_interim_nc_daily_merged/"

# The gridded array are very large. So we only want to extract exactly what
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

nc_close(nc)

# Make sure these coordinates match the fire coordinates by ensuring that
# fires are plotted over this western US snapshot
flipper <- length(ecmwf_latitude):1
quartz(width=8, height=5)
image.plot(ecmwf_longitude, ecmwf_latitude[flipper], t2m[,flipper, 180])
# Add fires
points(fireLonAdjusted, fireLat, pch=".")
map("state", add=T)
title("This is where fireLon is adjusted + 360")


quartz(width=8, height=5)
map("state")
image.plot(ecmwf_longitude-360, ecmwf_latitude[flipper], t2m[,flipper, 180])
# Add fires
points(fireLon, fireLat, pch=".")
map("state", add=T)
title("This is where ecmwf is adjusted - 360")


# Function for loading large .nc files disired latitude longitude extent as 
# demonstrated by t2m above
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

# Now we need the grid_attributes to be subset by the same latitude longitude 
# we extracted from nc. 
grid_lat_keeps <- grid_latitude %in% ecmwf_latitude
grid_lon_keeps <- grid_longitude %in% ecmwf_longitude

# Do the subsetting with the masks created above
grid_state_subset     <- grid_state[grid_lon_keeps, grid_lat_keeps]
grid_ecoregion_subset <- grid_ecoregion[grid_lon_keeps, grid_lat_keeps]
grid_ecoregion_mask   <- grid_ecoregion_subset == ecoregion_select

# Plot the masked ecoregion with a state and ecoregion map 
quartz()
f <- length(ecmwf_latitude):1
image.plot(x=(ecmwf_longitude-360), y=ecmwf_latitude[f], grid_ecoregion_mask[,f])
map("state", add=T, lwd=0.3)
plot(ecoregion_polygon, add=T)
title("Red should be over plotted ecoregion")

# Make summer totals / means depending on the paramter. This is only done for
# selected eco_region and months. 
v <- rep(NA, nYears)
t2m_mean <- v
tp_total <- v
e_total  <- v
rh2m_mean<- v 
d2m_mean <- v
ws_mean  <- v
wg_mean  <- v

# This fuction handles the multi-step hasle of subsetting in both time and
# a non regular grid space (ecoregion) 
spaceTimeStat <- function(x, tMask, FUN="mean", grid_ecoregion_mask){
  
  # Time subset first, preserve only desired months 
  spatialSubset <- x[,, tMask]
  # stats on time (called mean but could be "sum")
  gridMeans <- apply(spatialSubset, 1:2, FUN)
  
  # Take the final stat of the spatial subset
#  if(FUN=="mean"){
    spatialMean <- mean(gridMeans[grid_ecoregion_mask])
  # } else{
  #   # FUN == "sum"
  #   spatialMean <- sum(gridMeans[grid_ecoregion_mask])
  # }
  return(spatialMean)
  
}

# TODO: Save out monthly relationships 
for (i in 1:nYears){
  
  # Time masks 
  yearMask  <- years[i] == tYear
  monthMask <- tMon %in% month_select
  tMask     <- yearMask & monthMask 
  
  # Take the spatial temporal subset average or sum and store 
  t2m_mean[i] <- spaceTimeStat(t2m, tMask, FUN="mean", grid_ecoregion_mask)
  tp_total[i] <- spaceTimeStat(tp, tMask, FUN="sum", grid_ecoregion_mask)
  e_total[i]  <- spaceTimeStat(e, tMask, FUN="sum", grid_ecoregion_mask)
  rh2m_mean[i]<- spaceTimeStat(rh2m, tMask, FUN="mean", grid_ecoregion_mask)  
  d2m_mean[i] <- spaceTimeStat(d2m, tMask, FUN="mean", grid_ecoregion_mask)  
  ws_mean[i]  <- spaceTimeStat(ws, tMask, FUN="mean", grid_ecoregion_mask) 
  wg_mean[i]  <- spaceTimeStat(wg, tMask, FUN="mean", grid_ecoregion_mask)
  
}

################################################################################
# Plot the relationships 
################################################################################
# TODO: Save "experiment" directories that document time and spatial extent 
# TODO: correlations being plotted. 
# TODO: ggplopt, shade by year? 
scatter_plot <- function(x=t2m_mean, y1=FPA_BA_lightning, y2=FPA_BA_human,
                         xLab="", yLab="", regionName=""){
  
  # Calculate the correlations, these are being calculated with respect to 
  # predictive ability so now use pearson correlation coef. 
  r1 <- cor(x, y1, method="pearson")
  r2 <- cor(x, y2, method="pearson")
  
  r1_show <- round(r1, 2)
  r2_show <- round(r2, 2)
  
  
  print(paste("r2=",summary(lm(y1 ~ x))$r.squared))
  
  # Set up the file name based on plotted variable
  f <- paste0(figureDir, xLab, "_vs_", yLab, ".png")
  yLab <- str_replace(yLab, " ", "\n")
  
  # Set up the figure
  png(filename=f, width=3000, height=1400, res=300)
  par(mar=c(4,12,4,12), font=2)
  
  plot(x, y1,
       pch=19, cex=1.5,
       col="gray", bty="n", 
       yaxt="n", xaxt="n",
       ylab="", xlab="",
       axes=F
  )
  mtext(xLab, side=1, line = 2)
  mtext(paste("Lightning-ignited\n acres burned"), 
        side=2, line = 4.5, las=1, col="darkgray")
  
  # Handle axis labels 
  ylabValues <- axis(2, labels = FALSE, lwd=0) # Get defual values without labelling
  Labels <- pretty10exp(ylabValues, digits=1)
  eaxis(side=2, f.smalltcl=0, drop.1=T, at=ylabValues, labels=Labels)
  eaxis(side=1, f.smalltcl=0)
  
  # Plot and label the lightning area burned
  par(new=TRUE)
  plot(x, y2, col="orange", pch=19, cex=1.5,
       axes=FALSE, ann=FALSE )
  eaxis(side=4, f.smalltcl=0)#, max.at=5)
  #axis(4, col="orange", labels = FALSE)
  mtext("Human-ignited\n acres burned", 
        side=4, line = 4.5, las=1, col="orange", font=2)
  
  # Put the linear-correlations on the y-axis
  legend("topleft",
         legend=c(paste0("r=",r1_show), paste0("r=",r2_show) ),
         title=" Pearson correlation: ",
         pch=19,
         box.lwd=0.7,
         box.lty=3,
         cex=0.8,
         pt.cex=1.3,
         col=c("gray", "orange")
         )
  
  title(regionName, line=0.3)
  
  
  dev.off()
  
  return(c(r1, r2))
  
}

inPerM <- 39.3701

# TODO: make total precip mean of all boxes totals.
r_T2m <- scatter_plot(x=t2m_mean-273.15, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="Temperature [C]", yLab="Acres Burned", regionName = ecoregion_name)
r_TP <- scatter_plot(x=tp_total*inPerM, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="Total precipitation [in]", yLab="Acres Burned", regionName = ecoregion_name)
r_RH <- scatter_plot(x=rh2m_mean, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="Relative humidity", yLab="Acres Burned", regionName = ecoregion_name)
scatter_plot(x=d2m_mean-273.15, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="Dew point temperature [C]", yLab="Acres Burned", regionName = ecoregion_name)
scatter_plot(x=ws_mean, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="10-meter wind speed [meters per second]",  yLab="Acres Burned", regionName = ecoregion_name)
scatter_plot(x=wg_mean, y1=FPA_BA_lightning, y2=FPA_BA_human, xLab="10-meter wind gust [meters per second]",  yLab="Acres Burned", regionName = ecoregion_name)

x <- rep(NA, 6)
r_df <- data.frame(r=x, variable=x, ignitionType=x, ecoregion=x)
r_df[1,] <-c(r_T2m[1], "T2M", "Lightning", ecoregion_select)
r_df[2, ] <- c(r_T2m[2], "T2M", "Human", ecoregion_select)
r_df[3,] <-c(r_TP[1], "TP", "Lightning", ecoregion_select)
r_df[4, ] <- c(r_TP[2], "TP", "Human", ecoregion_select)
r_df[5,] <-c(r_RH[1], "RH", "Lightning", ecoregion_select)
r_df[6, ] <- c(r_RH[2], "RH", "Human", ecoregion_select)

saveName <- paste0("Data/correlations/ecoregion_", ecoregion_select, ".RData")
save(r_df, file=saveName)
 

# ################################################################################
# # These relationships should happen on a grid by grid basis. Show correlations
# # between these grid boxes burn area. Call upon monthly gridded burn area
# # created by grid_FPA_FOD_burn_area.R
# ################################################################################
# 
# # Here is where you throw out all of the spatially subsetting you previously
# # did! Just plot your desired spatial domain. 
# 
# ncFile <- "Data/FPA_FOD/burn_area_monthly_75x75_1992_2015.nc"
# nc <- nc_open(ncFile)
# grid_area <- ncvar_get(nc, "grid_area")
# BA_longitude <- ncvar_get(nc, "longitude") # "BA" for burn_area. to make sure coordinates align with ecmwf
# BA_latitude  <- ncvar_get(nc, "latitude")
# burn_area <- ncvar_get(nc, "burn_area")
# human_burn_area <- ncvar_get(nc, "human_burn_area")
# lightning_burn_area <- ncvar_get(nc, "lightning_burn_area")
# unknown_burn_area <- ncvar_get(nc, "unknown_burn_area")
# t_hours   <- ncvar_get(nc, "time")
# nc_close(nc)
# 
# 
# # Time in hours from origin to POSIXct
# t_mon_FPA <- as.POSIXct(t_hours*60^2, origin="1900-01-01 00:00:0.0", tz="utc")
# 
# # Spatially subset the data to match the loaded ecmwf data
# lonMask <- BA_longitude %in% ecmwf_longitude
# latMask <- BA_latitude %in% ecmwf_latitude
# 
# if(unique(BA_longitude[lonMask] - ecmwf_longitude) != 0 |
#    unique(BA_longitude[lonMask] - ecmwf_longitude) != 0){
#   stop("Weather dimensions and burn area dimensions do not match.")
# } else{
#   BA_longitude <- BA_longitude[lonMask]
#   BA_latitude  <- BA_latitude[latMask]
# }
# 
# # If the code is still running then the following subset of the data works! 
# burn_area           <- burn_area[, lonMask, latMask]
# human_burn_area     <- human_burn_area[, lonMask, latMask]
# lightning_burn_area <- lightning_burn_area[, lonMask, latMask]
# unknown_burn_area   <- unknown_burn_area[, lonMask, latMask]
# 
# # Make sure the three types of burn area add up to the "burn_area" quantity, 
# # which is supposed to be the sum of the three.
# total <- human_burn_area + lightning_burn_area + unknown_burn_area
# dArea <- burn_area - total
# 
# # Plot the spatial total as a sanity check
# burn_area_total <- apply(burn_area, 2:3, sum)
# dArea_total <- apply(dArea, 2:3, sum)
# 
# # Makes plotting clearer by removing locations == 0 
# burn_area_total[burn_area_total==0] <- NA 
# dArea_total[dArea_total==0] <- NA
# 
# # NOTE: Note the shifted coordinates so I can plot the U.S. map in this sanity
# # NOTE: check plot 
# pdf(file="Figures/gridded_fires_test_grid.pdf", width=24, height=20)
# par(mfrow=c(1,1))
# f <- length(ecmwf_latitude):1
# image.plot(ecmwf_longitude-360, ecmwf_latitude[f], burn_area_total[,f])
# map("state", add=T)
# title("Total burn area in dataset. No time subsetting.")
# 
# # Now add the fires that were gridded to the plot 
# points(FPA_FOD$LONGITUDE, FPA_FOD$LATITUDE, pch=".", col="black")
# 
# dev.off()
# 
# # image.plot(ecmwf_longitude-360, ecmwf_latitude[f], dArea_total[,f])
# # map("state", add=T)
# # title("difference from variable saved as burn area and sum of its parts")
# 
# # NOTE: R/grid_FPA_FOD_burn_area.R grids burn area of three types of ignitions
# # NOTE: human, lightning, and unknown. At the end these three are added to create
# # NOTE: the "burn_area" i.e. total burn area variable. Any differenes between
# # NOTE: this variable and the sum of its parts is due to floating point inprecision
# # NOTE: in saving. 
# 
# # Thus far I have made the chioce that "unknown" fire start cuases are human, so
# # I am going to merge these two totals. 
# human_burn_area <- human_burn_area + unknown_burn_area
# rm(unknown_burn_area)
# 
# ################################################################################
# # Create Monthly ecmwf means and totals. Remember, spearman correlations in case
# # the relations are not linear. 
# ################################################################################
# 
# # t2m subset by t made up of tMon and tYear
# # TODO: 
# nMonths <- length(unique(tYear)) * length(unique(tMon))
# monthly <- array(NA, dim = c(length(ecmwf_longitude), 
#                              length(ecmwf_latitude), 
#                              nMonths))
# 
# # Give monthly array dimensions to all variables 
# t2m_monthly <- monthly
# tp_monthly  <- monthly
# e_monthly   <- monthly
# rh2m_monthly<- monthly 
# d2m_monthly <- monthly
# ws_monthly  <- monthly
# wg_monthly  <- monthly
# 
# # CHECK ME
# i <- 0
# for (y in sort(unique(tYear))) {
#   
#   for (m in sort(unique(tMon))){
#     i <- i + 1
#     tMask <- y == tYear & m == tMon
#     
#     #print(paste( min(t[tMask]) ,"-", max(t[tMask])))
#     
#     t2m_monthly[,,i] <- apply(t2m[,,tMask], 1:2, mean)
#     tp_monthly[,,i] <- apply(tp[,,tMask], 1:2, sum)
#     e_monthly[,,i] <- apply(e[,,tMask], 1:2, sum)
#     rh2m_monthly[,,i] <- apply(rh2m[,,tMask], 1:2, mean)
#     d2m_monthly[,,i] <- apply(d2m[,,tMask], 1:2, mean)
#     ws_monthly[,,i] <- apply(ws[,,tMask], 1:2, mean)
#     wg_monthly[,,i] <- apply(wg[,,tMask], 1:2, mean)
#     
#   }
#   
# }
# 
# print(paste("The dimension of e.g. t2m_montly is:", 
#             dim(t2m_monthly)[1],
#             dim(t2m_monthly)[2], 
#             dim(t2m_monthly)[3]))
# 
# # Make spatial correlations with gridded burn area and display as sketched in 
# # notebook. Plot the indivisual grid correlations to see what these correlation
# # values really mean. Make sure to label with r value, lon, and lat. 
# # TODO: Interannual variability of this relation. This means summing the desired
# # TODO: months yearly totals. 
# # TODO: Lagged correlations? 
# # TODO: Consider consistent 0:1 -1:1 colorbars so different variables easier to
# # TODO: compare. 
# 
# spatialCorrelationMap <- function(varName="2-meter temperature", 
#                                   varValues=t2m_monthly,
#                                   includedMonths=1:12,
#                                   annual=FALSE, 
#                                   colorScheme="heat",
#                                   makePlot=FALSE){
#   
#   # If we are interested in burn area correlation for select months, year totals
#   if(annual){
#     
#     nYear <- length(years)
#     l_BA  <- array(data=NA, dim=c(nYear, nLon, nLat))
#     h_BA  <- array(data=NA, dim=c(nYear, nLon, nLat))
#     varValuesAnnual <- array(data=NA, dim=c(nLon, nLat, nYear))
#     
#     # Month mask does not change within year loop 
#     month_mask <- lubridate::month(t_mon_FPA) %in% includedMonths
#     for (i in 1:nYears){
#       
#       # Mask out desired months for this year only 
#       seasonMask <- years[i] == lubridate::year(t_mon_FPA) & month_mask
#       
#       # Sum the burn area in these months for this year, per spatial dim  
#       l_BA[i,,] <- apply(lightning_burn_area[seasonMask,,], 2:3, sum)
#       h_BA[i,,] <- apply(human_burn_area[seasonMask,,], 2:3, sum)
#       
#       # For the met variable now too. 
#       # TODO: Think about if mean vs. sum matters for variables like precip. For
#       # TODO: a spearnman correlation I do not think it does, as mean is simple
#       # TODO: sum scaled by N. 
#       varValuesAnnual[,,i] <- apply(varValues[,,seasonMask], 1:2, sum)
#       
#     }
#     # Now change the name to be consistent with the rest of this function. 
#     varValues <- varValuesAnnual
#     
#   } else{
#     
#     # Perform correlations for the subset monthly time array, not summing the 
#     # yearly values. 
#     
#     # Mask out months we do not want to include in temporal correlation 
#     month_mask <- lubridate::month(t_mon_FPA) %in% includedMonths
#     print(paste("The number of months used per box is:", sum(month_mask)))
#     
#     # Subset burn area only by the months
#     l_BA <- lightning_burn_area[month_mask,,]
#     h_BA <- human_burn_area[month_mask,,]
#     
#     # Subset varValues
#     varValues <- varValues[,,month_mask]
#     
#   }
#   
#   # Time the correlation loop, which can be slop when making scatterplots for
#   # each grid point. 
#   t0 <- Sys.time()
#   
#   nLat <- length(BA_latitude)
#   nLon <- length(BA_longitude)
#   lightning_corMat <- array(NA, dim=c(nLon,nLat))
#   human_corMat     <- array(NA, dim=c(nLon,nLat))
#   for(i in 1:nLat){
#     
#     print(paste("Correlation calculation & plotting % complete:",i/nLat*100))
#     
#     for(j in 1:nLon){
#       
#       x <- varValues[ j, i,] # time last for met var
#       # NOTE: the time dimension in burn area and weather arrays are in different
#       # NOTE: places, and this is annoying...
#       y_lightning <- l_BA[, j, i] 
#       y_human     <- h_BA[, j, i]
#       
#       # sometimes there is zero burn area all the time. That means the sd = 0 
#       # and that means you cannot take a spearman correlation. So check the data
#       # for some variation. If no variation, well NA says it all. 
#       if(sd(y_lightning) != 0 & sd(x) != 0){
#         r_lightning <- cor(x, y_lightning, method="spearman")
#         lightning_corMat[j,i] <- r_lightning
#       }
#       if(sd(y_human) != 0 & sd(x) != 0){
#         r_human <- cor(x, y_human, method="spearman")
#         human_corMat[j,i] <- r_human
#       }
#       
#       # Only plot the raw data that make the corrrelation when there is 
#       # variability for both and we have said, yes plase plot. 
#       if(sd(y_lightning) != 0 & sd(y_human) != 0 & makePlot){  
#         
#         # Get lat lon for labeling
#         lat <- BA_latitude[i] 
#         lon <- BA_longitude[j]
#         
#         # Round the correlation coef for labeling purposes
#         r_lightning_pretty <- round(r_lightning, 3)
#         r_human_pretty     <- round(r_human, 3)
#         
#         fileName <- paste0("Figures/var_correlation_by_grid_box/",
#                            varName, "_", lat, "_", lon, "_annual=",annual,"_",
#                            "months=", min(includedMonths), "-", max(includedMonths),
#                            "_rh=",r_human_pretty,"_rl=", r_lightning_pretty, 
#                            ".png")
#         
#         # Create the figure 
#         png(filename=fileName, res=100, width=1200, height=500)
#         
#         par(mar=c(4,7,4,4), mfrow=c(1,2))
#         
#         yMax <- max(c(y_lightning, y_human))
#         
#         plot(x, y_lightning, col="gray", ylim=c(0, yMax),
#              xlab="", ylab="", yaxt="n")
#         eaxis(2)
#         points(x, y_human, col="orange")
#         mtext(varName, side=1, line=3)
#         mtext("Burn Area", side=2, line=4.5)
#         
#         legend("topleft",
#                legend=c(r_lightning_pretty, r_human_pretty),
#                fill=c("gray", "orange"), 
#                bty="n")
#         
#         title(paste(varName, "coords:(", lat, ",", lon,")"))
#         
#         map("state", xlim=c(-125, -98))
#         points(lon-360, lat, pch=3, col="red")
#         
#         dev.off()
#         
#       } # End of plotting if statement
#       
#     }
#   }
#   
#   dt <- Sys.time() - t0
#   print(paste("It took ", dt/60, "minutes to run the loop"))
#   
#   ################################################################################
#   # Plot the correlation maps side by side 
#   ################################################################################
#   
#   # Get max and min values present in these correlation matrixs 
#   minCor <- min( c(min(lightning_corMat, na.rm=T), min(human_corMat, na.rm=T)) )
#   maxCor <- max( c(max(lightning_corMat, na.rm=T), max(human_corMat, na.rm=T)) )
#   
#   # Set the breaks and colorbar. 49 unique colors should be plenty for these data
#   # and correlations. 
#   if(colorScheme == "heat"){
#     
#     colorpallete <- rev(heat.colors(99))
#     breaks       <- seq(minCor, maxCor, length.out=100)
#     
#   }else if(colorScheme == "cool"){
#     
#     rwb <- colorRampPalette(colors = c("white", "blue"))
#     colorpallete <-rwb(99)
#     breaks       <- seq(minCor, maxCor, length.out=100)    
#     
#   }else if(colorScheme == "diverging"){
#     
#     # Set up diverging colorbar with center at zero (white)
#     rwb <- colorRampPalette(colors = c("blue", "white", "red"))
#     colorpallete <- rwb(99)
#     
#     # This now needs to be symetric for fair assessment 
#     maxAbs <- max(abs( c( (minCor), abs(maxCor) ) ))
#     #maxAbs <- 1 # would make colorbar consistent across all diverging. 
#     breaks       <- seq(maxAbs*-1, maxAbs, length.out=100)
#     
#   }
#   
#   
#   # Plot the correlation maps side by side 
#   plotLon <- ecmwf_longitude - 360
#   f <- nLat:1
#   plotLat <- ecmwf_latitude[nLat:1]
#   
#   # Name of the saved correlation map needs to indicate what temporal information
#   # was used to create it
#   if(annual){
#     saveName <- paste0("Figures/spatial_correlation_maps/",
#                        varName, "_burned_area_cor_",
#                        "months=", min(includedMonths), "-", max(includedMonths),
#                        "_annual_sums.png")
#     
#   }else{
#     # Show included months in monthly time series 
#     saveName <- paste0("Figures/spatial_correlation_maps/",
#                        varName,
#                        "_burned_area_cor_",
#                        "months=", min(includedMonths), "-", max(includedMonths),
#                        ".png")
#   }
#   
#   
#   png(filename=saveName, res=250, width=2350, height=1000)
#   
#   par(mfrow=c(1,2))
#   
#   par(mar=c(2,5.2,2,1))
#   image(plotLon, plotLat, lightning_corMat[, f], 
#         col=colorpallete, breaks=breaks,
#         xlab="", ylab="", bty="n", xaxt="n", yaxt="n",
#         xlim=c(-125, -100))
#   map("state", col="black", add=T)
#   plot(SPDF, add=T, border="black", lty=3)
#   mtext("Lightning ignited fires", side=1, cex=1.5, font=2, line=0.5)
#   
#   par(mar=c(2,1,2,5.2))
#   image.plot(plotLon, plotLat, human_corMat[, f], 
#              col=colorpallete, breaks=breaks,
#              xlab="", ylab="", bty="n", xaxt="n", yaxt="n",
#              xlim=c(-125, -100),
#              legend.width=2)
#   map("state", col="black", add=T)
#   plot(SPDF, add=T, border="black", lty=3)
#   mtext("Human ignited fires", side=1, cex=1.5, font=2, line=0.5)
#   
#   # plot over the center of both 
#   if(annual){
#     timeText <- "yearly correlation"
#   }else{
#     timeText <- "monthly correlation"
#   }
#   
#   mtext( paste(varName, "& area burned |", timeText), 
#          outer = TRUE , line=-1.8, cex=1.6, font=2)
#   
#   dev.off()
#   
# }
# includedMonths <- 6:9
# annual <- F
# spatialCorrelationMap(varName="2-meter temperature",  varValues=t2m_monthly, includedMonths, annual, colorScheme="diverging", makePlot=F)
# spatialCorrelationMap(varName="total precip", varValues=tp_monthly, includedMonths, annual, colorScheme = "diverging", makePlot=F)
# spatialCorrelationMap(varName="2-meter RH",  varValues=rh2m_monthly, includedMonths, annual, colorScheme = "diverging", makePlot=F)
# spatialCorrelationMap(varName="2-meter dew point",  varValues=d2m_monthly, includedMonths, annual, colorScheme = "diverging", makePlot=F)
# spatialCorrelationMap(varName="total evaporation",  varValues=e_monthly, includedMonths, annual, colorScheme = "diverging", makePlot=F)
# spatialCorrelationMap(varName="2-meter wind gusts",  varValues=wg_monthly, includedMonths, annual, colorScheme = "diverging", makePlot=F)
# spatialCorrelationMap(varName="2-meter wind speed",  varValues=ws_monthly, includedMonths, annual, colorScheme = "diverging", makePlot=F)