# plot_ignition_weather.R 

# ------------------------- Description ---------------------------------------
# This script will be used to determine if the weather on ignition days is 
# different for different ignition types. 

# Load the required libraries
library(maps)
library(sfsmisc)

# ----------------------- Subset arguments -------------------------------------

# #Lat lon extent for FPA FOD data in the western US
# minLat <- 30
# maxLat <- 50
# minLon <- -125
# maxLon <- -100

# ecoregion
# ecoregion <- 11.1

# Lat lon extent for FPA FOD data in contiguous US (CONUS)
minLat <- 25
maxLat <- 50
minLon <- -125
maxLon <- -60
# ecoregion
ecoregion <- c(6.2,  9.2,  9.3, 10.1, 9.4, 13.1, 12.1, 10.2, 11.1,  7.1,  8.4,  8.3,  8.5,  9.5,  5.2,  8.1, 5.3,  3.1,  6.1, 15.4,  0.0 , 2.2,  3.2,  2.3,  8.2,  9.6)



# Select which environmental variable assigned to fires will be subset by
# elevation and compared between ignition types
weatherVAR <- "fuel_moisture_1000hr" # "fuel_moisture_1000hr" "t2m"

# Season
includeMonths <- c(1:12)

# Min fire size. Recal, even though fire managers view every small fire as a 
# potentially large fire, in terms on when real emissions occur, there is a 
# good reason to ignore lots of really small fires. 
minFireSize <- 0

# ----------------------- Setup Figure Dirs ------------------------------------
# THESE ARGUMENTS SET UP A UNIQUE EXPERIMENT. MAKE A DIRECTORUY TO STORE EACH
# EXISTING EXPERIMENT. 
if(length(ecoregion)>1){
  ecoName <- "all"
}else{
  ecoName <- ecoregion
}
experimentVARS <- paste0("ecoregion=", ecoName,"_", #"west_",
                         "months=", min(includeMonths), "-", max(includeMonths), "_",
                         "sizeMin=", minFireSize
                         )

experimentDir <- paste0("Figures/ignitionCompare_", experimentVARS, "/")
dailyDir <- paste0(experimentDir, "JDayIgnitionTvsElevation/")
binnedDir <- paste0(experimentDir, "binnedElevationDistributions/")
if(!dir.exists(experimentDir)){

  dir.create(experimentDir)
  dir.create(dailyDir)
  dir.create(binnedDir)
  
} else{
  
  print("The Direcory exists. Overwriting files within")
  
}

# Load the fire occurence data, "FPA_FOD" dataframe
#dataFile <- paste0("Data/FPA_FOD/FPA_FOD_ecmwf_1992_2015.RData")
dataFile <- paste0("Data/FPA_FOD/FPA_FOD_gridMet_1992-2015.RData")
load(dataFile)

# C is much nicer than K. What am I a chemist?
FPA_FOD$t2m <- FPA_FOD$t2m - 273.15  

# Get the spatial coordinates of the fires 
lat <- FPA_FOD$LATITUDE
lon <- FPA_FOD$LONGITUDE
startMonth <- FPA_FOD$START_MONTH

# spatial subset 
latMask <- lat > minLat & lat < maxLat
lonMask <- lon > minLon & lon < maxLon

# time subset 
timeMask <- FPA_FOD$START_MONTH %in% includeMonths

# Ecoregion ecoregion
ecoregionMask <- FPA_FOD$NA_L2CODE %in% ecoregion

# Fire size mask 
sizeMask <- FPA_FOD$FIRE_SIZE >= minFireSize

# Mask out the fires that do not have a known cause (NEW)
causeMask <- FPA_FOD$STAT_CAUSE_DESCR != "Missing/Undefined"

# Elevation Mask
# NOTE: Because the data are not always precise for location of fire, I have at
# NOTE: at least one fire that is just off coast of San Bern CA in the ocean and 
# NOTE: was assigned an elevation of -557m (FOD_ID: 1276313). For fires less than
# NOTE: 100m I will assume something similar happenned, the data are suspect,
# NOTE: and throw them out in this analysis where elevation are used. I know
# NOTE: that Death Valley is lowest point, at about -85m.
ELEV_MASK <- FPA_FOD$elevation <= -100 # meters

# Subset the data, one mask to rule them all. 
m <- timeMask & latMask & lonMask & ecoregionMask & sizeMask & causeMask & !ELEV_MASK
FPA_FOD <- FPA_FOD[m,]

# Get the cuase so we can create comparisons 
fireCause     <- FPA_FOD$STAT_CAUSE_DESCR
lightningMask <- fireCause == "Lightning"
humanMask     <- fireCause != "Lightning" 

# Function for plotting two distributions
plot_ignition_weather_distribution <- function(VAR="t2m", varLabel, humanMask){
  
  # Create the density estimates for the two types. 
  dHuman     <- density(FPA_FOD[[VAR]][humanMask], na.rm=T)
  dLightning <- density(FPA_FOD[[VAR]][!humanMask])
  
  par(las=1)
  
  # Get the needed max value for the y-axis
  maxVal <- max(c(dHuman$y,  dLightning$y))
  
  plot(dHuman$x, dHuman$y, 
       type="l", col='orange',
       xlab=varLabel,
       ylim=c(0, maxVal),
       lwd=4,
       bty="n", 
       cex.axis=2,
       cex.lab=2,
       yaxt="n",
       ylab="")
  eaxis(2, cex.axis=2)
  
  mtext("Probability\nDensity   ", 2,  line=5, cex=2)
  
  lines(dLightning$x, dLightning$y, col="darkgray", lwd=4)
  
  title(paste("Mean", varLabel, "at Ignition Locations"), cex.main=1.4)
  
  legend_human <- paste("Human, n=", sum(humanMask))
  legend_lightning <- paste0("Lightning, n=", sum(!humanMask))
  
  legend("topright",
         legend=c(legend_human, legend_lightning),
         col=c("orange", "darkgray"),
         lty=1,
         cex=1,
         lwd=4,
         bty="n")
  
  
}

# Show the distribution of environmental variables for two ignition types
# for all elevations and days to motivated subset comparison
png(filename=paste0(experimentDir, "t2m_distribution_all.png"),
    width=2500, height=2200, res=250)
par(mar=c(5,8,5,5))
plot_ignition_weather_distribution("t2m", "Temperature [C]", humanMask)
dev.off()

png(filename=paste0(experimentDir, "fm1000_distribution_all.png"),
    width=2900, height=2200, res=250)
par(mar=c(5,15,5,5))
plot_ignition_weather_distribution("fuel_moisture_1000hr", "1000-hr fuel moisture %", humanMask)
dev.off()


################################################################################
# Sample mean differences vs. Julain day
################################################################################


# Time difference tolerance in JDays for fires to be comparable
# Sets the size for a batch of days where I consider season to be accounted for 
JDayTol       <- 20 

# How many fires of each ignition type are required for us to compare them at a 
# given JDay window and elevation bin.
minSampleSize <- 10 

# The elevation bins. We make the assumption that locations within these bins
# will not have different temperatures that can be entirely explained by 
# differences in elevation
elevationBins <- seq(-100, 4000, by=100)
nBins         <- length(elevationBins) - 1 # bins are values between these

# Get the fire occurence day of year for ignitions information 
JDay       <- FPA_FOD$DISCOVERY_DOY
JDaysArray <- sort(unique(JDay)) # limited by selected months
nDays      <- length(JDaysArray)

# Do the same for burn area of fires started by Julian day. 
# TODO: Do this with dplr or aggreate 
# TODO: https://stackoverflow.com/questions/7560671/aggregate-data-in-one-column-based-on-values-in-another-column
JDayBurnArea_H <- rep(NA, nDays)
JDayBurnArea_L <- rep(NA, nDays)
for (j in 1:nDays){
  # Mask the data by each JDay and count the burn area that occurs from fires
  # that start on that JDay. 
  jDayMask        <- FPA_FOD$DISCOVERY_DOY == JDaysArray[j] 
  JDayBurnArea_H[j] <- sum(FPA_FOD$FIRE_SIZE[jDayMask & humanMask])
  JDayBurnArea_L[j] <- sum(FPA_FOD$FIRE_SIZE[jDayMask & lightningMask])
  
}

cumulative_BA_H <- cumsum(JDayBurnArea_H)
cumulative_BA_L <- cumsum(JDayBurnArea_L)
cumulative_BA <- cumulative_BA_H + cumulative_BA_L

# Estimate the 10th and 90th percent accumulation of burn area JDays
get_span <- function(BAC){
  percentOfBurn <- BAC/max(BAC)
  x10 <- which.min(abs(percentOfBurn-.10))
  x90 <- which.min(abs(percentOfBurn-.90))
  SPAN <- c(JDaysArray[x10], JDaysArray[x90])
  return(SPAN)
}



# Find the edges of the julain day experiment using a histogram of ignition counts
# by JDay 
png(filename=paste0(experimentDir, "JDayIgnitionCounts.png"),
    res=250, height=2000, width=3700)

par(mar=c(4,5.5,4,8))

# Histrogram of all ignitions by Julian day
h <- hist(JDay, breaks = JDaysArray, 
          main="", las=1, xaxt="n", 
          xlab="", ylab="Ignition Count",
          cex.lab=1.8,
          col="gray")
title(paste("Ecoregion", ecoName,"ignition counts: 1992-2015"))

# Non-lightning only 
hist(JDay[!lightningMask], breaks = JDaysArray, 
     col="orange", add=T, xlim=c(0, 370), xaxt="n")
mtext("Day of year", 1, line=2.4, cex=1.5)
# TODO: Make this by burn area? That way we are choosing the area based on 
# TODO: significance? 
span <- as.numeric(quantile(JDay, c(0.10, 0.90) ) ) 
#abline(v=span, col="black", lwd=3, lty=3)

# Add a date axis
t0 <- as.POSIXct("2012-01-01", tz="UTC")
DOY <- t0 + JDaysArray*24*60^2

axis(1, at=JDaysArray, labels=format(DOY, "%m/%d"))
# TODO: Also consider axis.Date()
# TODO: burn area in background on second axis would make compelling argument..
par(new=TRUE)

# Set up the ylim, because sometimes there is more burn area from humanm ignited
yMax <- max(c(cumulative_BA_L, cumulative_BA_H))

plot(JDaysArray, cumulative_BA_L, 
     ylim=c(0, yMax),
     pch=19, col="gray", bty="n", yaxt="n", xaxt="n",
     ylab="", xlab="", type="l", lwd=5)
lines(JDaysArray, cumulative_BA_H, col="orange", pch=19, lwd=5)

eaxis(4)
mtext("Cumulative Burn Area (Acres) -", side=4, col="black", line=-2, cex=1.5)

# Show date of 10% and 90th percent of burn area met for both ignitions
abline(v=get_span(cumulative_BA_L), col="gray",lty=3, lwd=3)
abline(v=get_span(cumulative_BA_H), col="orange",lty=3, lwd=3)

# Give this a white background so vertical lines do not go through this text
legend("topleft", 
       legend=c("Lightning Ignitions", "Human Ignitions"),
       fill=c("gray", "orange"), 
       #bty="n", 
       box.col="white",
       cex=2,
       bg="white"
)

dev.off()

# # print(paste("80% of ignitions are between JDays", span[1], "-", span[2]))
# # print(paste("80% of BA is between JDays", span[1], "-", BA_span[2]))
# 
# # Define the range explicitly
# # NOTE: I am not sure it makes sense to do this, as compared to the whole year 
# # NOTE: story, since we account for JDay as it is. 
# minJDay <- min(JDaysArray) #span[1]
# maxJDay <- max(JDaysArray) #span[2]
# 
# # Create expiment data storage arrays
# # These are what we loop. These datys and tol gets inbetween 
# loopJDays <- seq(minJDay, maxJDay, by = JDayTol) 
# nLoop <- length(loopJDays)
# JDayFireCount <- rep(NA, nLoop)
# 
# # Empty storage dataframes creation. Lots of calculations per Julian day chunk
# # needs to be saved for plotting. 
# nCol <- nBins
# a <- array(NA, dim= c(nLoop,  nCol) )
# df_dummy <- data.frame(a)
# row.names(df_dummy) <- loopJDays
# 
# # For storing samples
# df_difference_of_means <- df_dummy
# df_sd <- df_dummy
# 
# # Shapiro-wilks tests a distribution for normality
# df_H_wilkes_p <- df_dummy
# df_H_wilkes_W <- df_dummy
# 
# df_L_wilkes_p <- df_dummy
# df_L_wilkes_W <- df_dummy
# 
# # Difference of the means test. 
# df_t_test  <- df_dummy
# df_t_upper <- df_dummy
# df_t_lower <- df_dummy
# 
# 
# # Loop through each day and perform experiment. Loop by JDayTol so no data is 
# # used more than once.  
# for (i in 1:nLoop){ 
#   
#   # Subset fires to a given julain day (or maybe range of +-10 or something?)
#   JDayMask <- abs(JDay - loopJDays[i]) <= (JDayTol/2)
#   JDayFireCount[i] <- sum(JDayMask)
#   
#   day_left_bound  <- min(JDay[JDayMask])
#   day_right_bound <- max(JDay[JDayMask])
#   
#   print(paste("These data are coming from JDays:", 
#               day_left_bound, "-", day_right_bound))
#   
#   # Subset the dataframe by those that are within the JDay tolerance, l==lightning
#   df    <- FPA_FOD[JDayMask,]
#   lMask <- df$STAT_CAUSE_DESCR == "Lightning"
#   
#   # Default color is orange until seen that it is lightning 
#   COL <- rep("orange", sum(JDayMask))
#   COL[lMask] <- "gray"
#   
#   # Plot the fire on an elevation vs. met variable set of axis. 
#   # TODO: boxplot of these data by elevation bin! 
#   # TODO: Also show the 
#   if(TRUE){
#     jDayFilename <- paste0(dailyDir, "JDay_", loopJDays[i],
#                            "_elevation_vs_",weatherVAR,".png")
#     png(filename=jDayFilename, res=250, height=2000, width=4000)
#     
#     plot(df$elevation, df[[weatherVAR]], 
#          xlab="Elevation [m]", ylab=weatherVAR,
#          col=COL, pch=19, cex=0.5, 
#          ylim=c(0, 30))
#     abline(v=elevationBins, lty=2)
#     title(paste("Julain Days in these data:", 
#                 day_left_bound, "-", day_right_bound))
#   
#     # All elevations difference of means 
#     meanLT <- mean(df[[weatherVAR]][lMask])
#     meanHT <- mean(df[[weatherVAR]][!lMask])
#     
#     abline(h=meanLT, col="gray", lwd=3)
#     abline(h=meanHT, col="orange", lwd=3)
#     
#     allDT <- meanLT - meanHT
#     
#     dev.off()
#   }
#   
#   # Now account for elevation by separating data by the chosen elevation bins. 
#   # IDEA: Adjust temperature based on elevation using standard lapse rate?
#   sampleElavations <- df$elevation
#   
#   # png(filename = paste0(dailyDir, "JDay_", loopJDays[i],"_sample_elevations.png"),
#   #     width=1000, height=1000, res=200)
#   # 
#   # # TODO: change this to a density curve and make it easier to read
#   # hist(sampleElavations[lMask], col=adjustcolor("gray", 0.4) )
#   # hist(sampleElavations[!lMask], col=adjustcolor("orange", 0.4), add=T)
#   # 
#   # dev.off()
#   
#   # Count total in each elevation bin for each ignition type. We will exclude
#   # the bins where there are not enough of each type
#   lightningBinCount <- hist(sampleElavations[lMask], breaks=elevationBins, plot=F)$counts
#   
#   humanBinCount <- hist(sampleElavations[!lMask], breaks=elevationBins, plot=F)$counts
#   
#   # Mask out (later) the elevation bins that do not have enough data
#   fullBins <- humanBinCount > minSampleSize & 
#               lightningBinCount > minSampleSize
#   
#   #################################################
#   # Stats for each elev bin calculated in this loop 
#   #################################################
#   for ( m in 1:nBins ){
#     
#     # Only bother with these calculations when there is enough data in each bin. 
#     if(fullBins[m]){
#     
#       elevMask <- sampleElavations >= elevationBins[m] & 
#                   sampleElavations < elevationBins[m+1]
#       
#       # # We need to know if the eleation distributions within this bin look 
#       # # alike. Elevation differences could still explain small differences
#       # # in temperature or other variable. 
#       # png(filename=paste0(binnedDir, "JDay_",loopJDays[i], "_",
#       #                     elevationBins[m],"-", elevationBins[m],".png"),
#       #     res=250, height=1000, width=1000)
#       # 
#       # hist(sampleElavations[elevMask & lMask], col=adjustcolor("gray", 0.4),
#       #      xlab="Ignition Elevations")
#       # hist(sampleElavations[elevMask & !lMask], add=T, col=adjustcolor("orange", 0.4))
#       # 
#       # dev.off()
#       
#       # Subset the data by elevation. Make these arrays ready to use. From here
#       # forward in JDayLoop:
#       # H = Human
#       # L = Lightning
#       H <- df[[weatherVAR]][elevMask & !lMask]
#       L <- df[[weatherVAR]][elevMask & lMask]
#       
#       # mean human ignited location weatherVAR sample  
#       meanHTs <- mean(H)
#       sdHTs   <- sd(H)
#       
#       # mean lightning weatherVAR sample 
#       meanLTs <- mean(L)
#       sdLTs   <- sd(L)
#       
#       # Difference in this elevation/time bin means 
#       df_difference_of_means[i, m] <- meanLTs - meanHTs
#       
#       # Are the means different? 
#       # 1) Wilkes test to see if the distributions are normal. 
#       # 2) If they are normal, t-test for difference of sample means. 
#       
#       # NOTE for wilkes test:
#       # Ho: The distribution is normal.
#       # h1: reject thre null, (evidense distribution is not normal)
#       
#       # Shapiro-Wilkes test requires sample size between 3-5000. If in
#       # this loop there are more than 3. If more than 5000, randomly draw 5000
#       # from the array. 
#       # TODO: Make into a function so that you do not have two ugly if statements
#       # TODO: back to back. 
#       if(length(L) > 5000){
#         L <- sample(L, size=5000)
#       }
#       if(length(H) > 5000){
#         H <- sample(H, size=5000)
#       }
#       
#       df_H_wilkes_p[i, m] <- shapiro.test(H)$p.value
#       df_H_wilkes_W[i, m] <- shapiro.test(H)$statistic
#         
#       df_L_wilkes_p[i, m] <- shapiro.test(L)$p.value
#       df_L_wilkes_W[i, m] <- shapiro.test(L)$statistic
#       
#       # Test the two distributions to see if the means are different. Two sided. 
#       # no strong a-priori for what value should be larger. 
#       t_test_results <- t.test(L, H, alternative = "two.sided")
#       df_t_test[i,m] <- t_test_results$p.value
#       
#       # store the 95% confidense bounds 
#       df_t_lower[i,m] <- t_test_results$conf.int[1]
#       df_t_upper[i,m] <- t_test_results$conf.int[2]
#       
#     }
#     
#   } # end of elevation bin looping 
#   
#   
# }
# 
# ################################################################################
# # Plot up these differences over Julain dates, one elevation curve at a time
# # TODO: ggplot2, size of dot related to sample size!!!
# ################################################################################
# library(fields)
# 
# # Label colums with bin mid-point elevation
# elev_labs <- elevationBins[1:nBins] + mean(diff(elevationBins))/2
# 
# df <- df_difference_of_means
# colnames(df) <- as.character(elev_labs)
# 
# 
# png(filename = paste0(experimentDir, "JDay_diff_Series_",weatherVAR,".png"), 
#     width=3600, height=1500, res=250)
# par(mar=c(4,6.5,4,6))
# 
# # Absolute biggest value in data
# maxValue <- max(abs(range(df, na.rm = T)))
# minValue <- min(range(df, na.rm = T))
# ylim <- c(minValue, maxValue)
#   
# 
# if(weatherVAR == "fuel_moisture_1000hr"){
#   niceYLab <- "1000-hr dead fuel moisture %"
# } else if (weatherVAR=="t2m"){
#   niceYLab <- "2-meter T [C]"
# }
# 
# plot(loopJDays, df[,8], 
#      pch="", bty="n", 
#      ylim=ylim, 
#      ylab=paste(niceYLab,"\n(Lightning - Human)"),
#      xlab="Day of Year", 
#      xaxt="n",
#      cex.lab=2, 
#      las=1)
# 
# axis(1, at=JDaysArray, labels=format(DOY, "%m/%d"), cex=2)
# 
# 
# # title(paste("(Mean Lightning Ignition", weatherVAR, 
# #             ") - (Mean Human Ignition", weatherVAR,")"),
# #       cex.main=2)
# 
# # Place lines for different elevation bands on blank plot, elevation bands
# # are represented by the number of columns. Time is rows. 
# nLines <- dim(df)[2]
# lineColors <- terrain.colors(nLines)
# 
# for (r in 1:nLines){ # 
#   
#   # Place points where significantly different from zero 
#   sigMask <- df_t_test[,r] < 0.05
#   points(loopJDays[sigMask], df[,r][sigMask], col=lineColors[r], pch=19, cex=2)
#   
#   # Now plot the points where the sigMask has p-value above 0.05
#   points(loopJDays[!sigMask], df[,r][!sigMask], col=lineColors[r], pch=1, cex=2)
#   
#   # # Add lines to make connexting the dots easier, also messing looking. 
#   # lines(loopJDays, df[,r], col=lineColors[r], lwd=2, lty=2)
#   
# }
# 
# abline(h=0, lwd=1, lty=2)
# 
# # Legend to understand the colors relation to elevation
# fields::image.plot( legend.only=TRUE, 
#                     zlim= range(elevationBins), 
#                     col=lineColors,
#                     #legend.lab="Elevation [m]",
#                     line=3,
#                     bg="transparent") 
# 
# mtext("Elevation [m]", side=4, line=-1, xpd=T, cex=2)
# 
# dev.off()

# ################################################################################
# # Histogram of differences chunked by Julain day and elevation
# ################################################################################
# 
# # # Skip columns (height bins) that are all NA, i.e. sample sizes were never 50
# # # for each type of ignition on each day. 
# nCol   <- dim(df_difference_of_means)[2]
# # useCol <- rep(TRUE, nCol)
# # for (x in 1:nCol){
# #   is.na(unique(df_difference_of_means[, 1]))
# # }
# 
# # TODO: Show the total number of fires that went into the calculation for
# # TODO: each elevation bin. Segregated by ignition type. 
# # TODO: Need objective way to determine which columns (height bins) are to be
# # TODO: looped and plotted. 
# 
# # What elevation bins have enough values to plot? 
# for (i in 1:nCol){
#   
#   DT <- df_difference_of_means[, i]
#   
#   if (sum(is.na(DT)) < nLoop){
#     
#     # Figure out the correct elevation bin associated with each row. 
#     binLabel <- paste(elevationBins[i], "-", elevationBins[i+1], "meters")
#     
#     # There is non NA data to make the histogram 
#     png(filename=paste0(experimentDir,"diffence_of_means_elevation_bin=",
#                         binLabel,".png"), 
#         width=1000, height=1000, res=250)
#     par(mfrow=c(1,1), las=1)
#     
#     # Show consistent size of hist 
#     hist(DT, xlim=c(-5,5), ylim = c(0, 40),
#          main="", ylab="", xlab="(T_lightning - T_human) [C]", 
#          col=lineColors[i])
#     abline(v=0, lty=2)
#     title(binLabel)
#     
#     dev.off()
#     
#   }
#   
# }


# ################################################################################
# # Regular Distributions 
# ################################################################################
# 
# 
# quartz(width=15, height=5)
# 
# png(filename=paste0("Figures/summary/ignitionTemperature_6_2_Fires>1000.png"),
#     width=3000, height=1000, res=200)
# par(mfrow=c(1,3), mar=c(4,5,4,4))
# plot_ignition_weather_distribution("t2m", "Temperature [K]", humanMask)
# 
# # Show the locations of the fires
# # TODO: Dot size is size of fire!
# map("state", xlim=c(-125,-102), ylim=c(33,50))
# points(FPA_FOD$LONGITUDE[humanMask], FPA_FOD$LATITUDE[humanMask], 
#        col="orange", pch=".")
# points(FPA_FOD$LONGITUDE[!humanMask], FPA_FOD$LATITUDE[!humanMask], 
#        col="darkgray", pch=".")
# title(paste("Fires locations"), cex.main=2)
# 
# 
# # TODO: histogram of the JDay that these fires are occuring on. 
# minDay <- min(FPA_FOD_subset$DISCOVERY_DOY)
# maxDay <- max(FPA_FOD_subset$DISCOVERY_DOY)
# 
# hist(FPA_FOD$DISCOVERY_DOY[humanMask], 
#      col=adjustcolor("orange",0.5),
#      breaks=c(minDay:maxDay), las=1, xlab="Day of Year",
#      ylim=c(0,30),
#      border="transparent",
#      main="Day of year ignition counts",
#      cex.main=2)
# hist(FPA_FOD$DISCOVERY_DOY[!humanMask], 
#      col=adjustcolor("darkgray", 0.5), breaks=c(minDay:maxDay), add=T,
#      border="transparent")
# dev.off()
# 
# plot_ignition_weather_distribution("tp_lastMonth", "Total Precip", humanMask)
# 
# 
# 
# quartz()
# plot(FPA_FOD_subset$windSpeed_MonthAfter, FPA_FOD_subset$FIRE_SIZE)

#windSpeed_MonthAfter



