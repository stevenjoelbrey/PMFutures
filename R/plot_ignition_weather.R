# assign_ignition_weather.R 

# ------------------------- Description ---------------------------------------
# This script will be used to determine if the weather on ignition days is 
# different for different ignition types. 

# Load the required libraries
library(maps)
library(sfsmisc)

# ----------------------- Subset arguments -------------------------------------

# Lat lon extent
minLat <- 30
maxLat <- 50
minLon <- -125
maxLon <- -100

# Season
includeMonths <- c(5:10)

# ecoregion 
ecoregion <- 6.2

# Min fire size. Recal, even though fire managers view every small fire as a 
# potentially large fire, in terms on when real emissions occur, there is a 
# good reason to ignore lots of really small fires. 
minFireSize <- 0

# ----------------------- setup figure dirs ------------------------------------
# THESE ARGUMENTS SET UP A UNIQUE EXPERIMENT. MAKE A DIRECTORUY TO STORE EACH
# EXISTING EXPERIMENT. 
experimentVARS <- paste0("ecoregion=", ecoregion,"_",
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

# Load the fire occurence data 
dataFile <- paste0("Data/FPA_FOD/FPA_FOD_ecmwf_1992_2015.RData")
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

# Subset the data, one mask to rule them all. 
m <- timeMask & latMask & lonMask & ecoregionMask & sizeMask
FPA_FOD <- FPA_FOD[m,]

# Get the cuase so we can create comparisons 
fireCause     <- FPA_FOD$STAT_CAUSE_DESCR
lightningMask <- fireCause == "Lightning"

# Function for plotting two distributions
plot_ignition_weather_distribution <- function(VAR="t2m", varLabel, humanMask){
  
  # Create the density estimates for the two types. 
  dHuman     <- density(FPA_FOD[[VAR]][humanMask])
  dLightning <- density(FPA_FOD[[VAR]][!humanMask])
  
  par(las=1)
  
  # Get the needed max value for the y-axis
  maxVal <- max(c(dHuman$y,  dLightning$y))
  
  plot(dHuman$x, dHuman$y, 
       type="l", col='orange',
       xlab=varLabel,
       ylab="", 
       ylim=c(0, maxVal),
       lwd=4,
       bty="n", 
       cex.axis=2,
       cex.lab=2,
       yaxt="n")
  eaxis(2, cex.axis=2)
  
  mtext("Probability Density", 2,  line=5, cex=2)
  
  lines(dLightning$x, dLightning$y, col="darkgray", lwd=4)
  
  title(paste("Mean", varLabel, "at Ignition Locations"), cex.main=1.4)
  
  legend_human <- paste("Human, n=", sum(humanMask))
  legend_lightning <- paste0("Lightning, n=", sum(!humanMask))
  
  legend("topleft",
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
par(mar=c(5,5,5,5))
plot_ignition_weather_distribution("t2m", "Temperature [C]", !lightningMask)
dev.off()

################################################################################
# Sample mean differences vs. Julain day
################################################################################

# Select which environmental variable assigned to fires will be subset by
# elevation and compared between ignition types
weatherVAR <- "t2m"

# Time difference tolerance in JDays for fires to be comparable
# Sets the size for a batch of days where I consider season to be accounted for 
JDayTol       <- 10 

# How many fires of each ignition type are required for us to compare them at a 
# given JDay window and elevation bin.
minSampleSize <- 50 

# The elevation bins. We make the assumption that locations within these bins
# will not have different temperatures that can be entirely explained by 
# differences in elevation
elevationBins <- seq(-300, 3900, by=100)
nBins         <- length(elevationBins) - 1 # bins are values between these

# Get the fire occurence day of year for ignitions information 
JDay       <- FPA_FOD$DISCOVERY_DOY
JDaysArray <- sort(unique(JDay)) # limited by selected months
nDays      <- length(JDaysArray)

# Find the edges of the julain day experiment using a histogram of ignition counts
# by JDay 
png(filename=paste0(experimentDir, "JDayIgnitionCounts.png"),
    res=250, height=2000, width=3700)
# All
hist(JDay, breaks = JDaysArray, main="Julian Day Ignition Counts: 1992-2015", las=1)
# Non-lightning only 
hist(JDay[!lightningMask], breaks = JDaysArray, col="orange", add=T, xlim=c(0, 370))
span <- as.numeric(quantile(JDay, c(0.10, 0.90) ) ) 
abline(v=span, col="red", lwd=3)
legend("topleft", legend=c("Total Ignitions", "Human Ignitions"),
       fill=c("white", "orange"), bty="n", cex=2)

dev.off()

print(paste("80% of data is between JDays", span[1], "-", span[2]))


# Define the range explicitly 
minJDay <- span[1]
maxJDay <- span[2]

# Create expiment data storage arrays
loopJDays <- seq(minJDay, maxJDay, by = JDayTol) # This is what we loop. These datys and tol gets inbetween 
nLoop <- length(loopJDays)
JDayFireCount <- rep(NA, nLoop)

# Empty storage dataframes creation. Lots of calculations per Julain day chunk
# needs to be saved for plotting. 
nCol <- nBins
a <- array(NA, dim= c(nLoop,  nCol) )
df_dummy <- data.frame(a)
row.names(df_dummy) <- loopJDays

# For samples
df_difference_of_means <- df_dummy
df_sd <- df_dummy

# Shapiro-wilks tests a distribution for normality
df_H_wilkes_p <- df_dummy
df_H_wilkes_W <- df_dummy

df_L_wilkes_p <- df_dummy
df_L_wilkes_W <- df_dummy

# Difference of the means test. 
df_t_test  <- df_dummy
df_t_upper <- df_dummy
df_t_lower <- df_dummy
# TODO: Create array to store 95% confidense level values. 

# Loop through each day and perform experiment. Loop by JDayTol so no data is 
# used more than once.  
for (i in 1:nLoop){ 
  
  # Subset fires to a given julain day (or maybe range of +-10 or something?)
  JDayMask <- abs(JDay - loopJDays[i]) <= JDayTol
  JDayFireCount[i] <- sum(JDayMask)
  
  # Subset the dataframe by those that are within the JDay tolerance
  df <- FPA_FOD[JDayMask,]
  lMask <- df$STAT_CAUSE_DESCR == "Lightning"
  
  # Default color is orange until seen that it is lightning 
  COL <- rep("orange", sum(JDayMask))
  COL[lMask] <- "gray"
  
  # TODO: boxplot of these data by elevation bin! 
  # TODO: Also show the 
  jDayFilename <- paste0(dailyDir, "JDay_", loopJDays[i],
                         "_elevation_vs_",weatherVAR,".png")
  png(filename=jDayFilename, res=250, height=2000, width=4000)
  
  plot(df$elevation, df[[weatherVAR]], 
       xlab="Elevation [m]", ylab=weatherVAR,
       col=COL, pch=19, cex=0.5, 
       ylim=c(0, 30))
  abline(v=elevationBins, lty=2)
  title(paste("Julain Day Center:", loopJDays[i], "+-", JDayTol))

  # All elevations difference of means 
  meanLT <- mean(df[[weatherVAR]][lMask])
  meanHT <- mean(df[[weatherVAR]][!lMask])
  
  abline(h=meanLT, col="gray", lwd=3)
  abline(h=meanHT, col="orange", lwd=3)
  
  allDT <- meanLT - meanHT
  
  dev.off()
  
  # Now account for elevation by separating data by the chosen elevation bins. 
  # IDEA: Adjust temperature based on elevation using standard lapse rate?
  sampleElavations <- df$elevation
  
  png(filename = paste0(dailyDir, "JDay_", loopJDays[i],"_sample_elevations.png"),
      width=1000, height=1000, res=200)
  # TODO: change this to a density curve and make it easier to read
  hist(sampleElavations[lMask], col=adjustcolor("gray", 0.4) )
  hist(sampleElavations[!lMask], col=adjustcolor("orange", 0.4), add=T)
  
  dev.off()
  
  # Count total in each elevation bin for each ignition type. We will exclude
  # the bins where there are not enough of each type
  lightningBinCount <- hist(sampleElavations[lMask], breaks=elevationBins, plot=F)$counts
  
  humanBinCount <- hist(sampleElavations[!lMask], breaks=elevationBins, plot=F)$counts
  
  # Mask out (later) the elevation bins that do not have enough data
  fullBins <- humanBinCount > minSampleSize & 
              lightningBinCount > minSampleSize
  
  #################################################
  # Stats for each elev bin calculated in this loop 
  #################################################
  for ( m in 1:nBins ){
    
    # Only bother with these calculations when there is enough data in each bin. 
    if(fullBins[m]){
    
      elevMask <- sampleElavations >= elevationBins[m] & 
                  sampleElavations < elevationBins[m+1]
      
      # We need to know if the eleation distributions within this bin look 
      # alike. Elevation differences could still explain small differences
      # in temperature or other variable. 
      png(filename=paste0(binnedDir, "JDay_",loopJDays[i], "_",
                          elevationBins[m],"-", elevationBins[m],".png"),
          res=250, height=1000, width=1000)
      
      hist(sampleElavations[elevMask & lMask], col=adjustcolor("gray", 0.4),
           xlab="Ignition Elevations")
      hist(sampleElavations[elevMask & !lMask], add=T, col=adjustcolor("orange", 0.4))
      
      dev.off()
      
      # Subset the data by elevation. Make these arrays ready to use. From here
      # forward in JDayLoop:
      # H = Human
      # L = Lightning
      H <- df[[weatherVAR]][elevMask & !lMask]
      L <- df[[weatherVAR]][elevMask & lMask]
      
      # mean human Temperature sample  
      meanHTs <- mean(H)
      sdHTs   <- sd(H)
      
      # mean lightning Temperature sample 
      meanLTs <- mean(L)
      sdLTs   <- sd(L)
      
      # Difference in this elevation/time bin means 
      df_difference_of_means[i, m] <- meanLTs - meanHTs
      
      # Are the means different? 
      # 1) Wilkes test to see if the distributions are normal. 
      # 2) If they are normal, t-test for difference of sample means. 
      
      # NOTE for wilkes test:
      # Ho: The distribution is normal.
      # h1: reject thre null, (evidense distribution is not normal)
      
      # Shapiro-Wilkes test requires sample size between 3-5000. If in
      # this loop there are more than 3. If more than 5000, randomly draw 5000
      # from the array. 
      # TODO: Make into a function so that you do not have two ugly if statements
      # TODO: back to back. 
      if(length(L) > 5000){
        L <- sample(L, size=5000)
      }
      if(length(H) > 5000){
        H <- sample(H, size=5000)
      }
      
      df_H_wilkes_p[i, m] <- shapiro.test(H)$p.value
      df_H_wilkes_W[i, m] <- shapiro.test(H)$statistic
        
      df_L_wilkes_p[i, m] <- shapiro.test(L)$p.value
      df_L_wilkes_W[i, m] <- shapiro.test(L)$statistic
      
      # Test the two distributions to see if the means are different. Two sided. 
      # no strong a-priori for what value should be larger. 
      t_test_results <- t.test(L, H, alternative = "two.sided")
      df_t_test[i,m] <- t_test_results$p.value
      
      # store the 95% confidense bounds 
      df_t_lower[i,m] <- t_test_results$conf.int[1]
      df_t_upper[i,m] <- t_test_results$conf.int[2]
      
    }
    
  } # end of elvation bin looping 
  
  
}

################################################################################
# Plot up these differences over Julain dates, one elevation curve at a time
################################################################################
library(fields)
df <- df_difference_of_means

png(filename = paste0(experimentDir, "JDay_diff_Series.png"), 
    width=3600, height=1500, res=250)
par(mar=c(4,6,4,6))

# Absolute biggest value in data
maxValue <- max(abs(range(df, na.rm = T)))
ylim <- c(maxValue * -1, maxValue)
  
plot(loopJDays, df[,8], pch="", bty="n", ylim=ylim, 
     ylab="Detla Temperature [C]",
     xlab="Julain Day", 
     cex.lab=2)
title(paste("(Mean Lightning Ignition", weatherVAR, ") - (Mean Human Ignition Temperature)"),
      cex.main=2)

# # Add a nice date axis
# DATES <- as.Date( (loopJDays-1), origin = "2014-01-01")
# axis(1, at=loopJDays, labels=DATES, lines=4)

# Place lines for different elevation bands on blank plot 
nLines <- dim(df)[2]
lineColors <- terrain.colors(nLines)

for (r in 1:nLines){ # 
  
  lines(loopJDays, df[,r], col=lineColors[r], lwd=3)
  
  # Place points where significantly different from zero 
  sigMask <- df_t_test[,r] < 0.05
  points(loopJDays[sigMask], df[,r][sigMask], col=lineColors[r], pch=1)
  
  # # Add points with errorbars 
  # arrows(x0=loopJDays, y0= df_t_lower[,r], 
  #        x1=loopJDays, y1= df_t_upper[,r],
  #        length=0.05, angle=90, code=3,
  #        col=lineColors[r])
  
}

abline(h=0, lwd=1, lty=2)

# Legend to understand the colors relation to elevation
image.plot( legend.only=TRUE, zlim= range(elevationBins), col=lineColors) 

dev.off()

################################################################################
# Histogram of differences chunked by Julain day and elevation
################################################################################

# # Skip columns (height bins) that are all NA, i.e. sample sizes were never 50
# # for each type of ignition on each day. 
nCol   <- dim(df_difference_of_means)[2]
# useCol <- rep(TRUE, nCol)
# for (x in 1:nCol){
#   is.na(unique(df_difference_of_means[, 1]))
# }

# TODO: Show the total number of fires that went into the calculation for
# TODO: each elevation bin. Segregated by ignition type. 
# TODO: Need objective way to determine which columns (height bins) are to be
# TODO: looped and plotted. 

# What elevation bins have enough values to plot? 
for (i in 1:nCol){
  
  DT <- df_difference_of_means[, i]
  
  if (sum(is.na(DT)) < nLoop){
    
    # Figure out the correct elevation bin associated with each row. 
    binLabel <- paste(elevationBins[i], "-", elevationBins[i+1], "meters")
    
    # There is non NA data to make the histogram 
    png(filename=paste0(experimentDir,"diffence_of_means_elevation_bin=",
                        binLabel,".png"), 
        width=1000, height=1000, res=250)
    par(mfrow=c(1,1), las=1)
    
    # Show consistent size of hist 
    hist(DT, xlim=c(-5,5), ylim = c(0, 40),
         main="", ylab="", xlab="(T_lightning - T_human) [C]", 
         col=lineColors[i])
    abline(v=0, lty=2)
    title(binLabel)
    
    dev.off()
    
  }
  
}


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



