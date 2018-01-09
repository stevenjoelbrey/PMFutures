# assign_ignition_weather.R 

# ------------------------- Description ---------------------------------------
# This script will be used to determine if the weather on ignition days is 
# different for different ignition types. 

# TODO: eliminate small fires. (e.g. less than 1 acre). What do we learn?

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

# THESE ARGUMENTS SET UP A UNIQUE EXPERIMENT. MAKE A DIRECTORUY TO STORE EACH
# EXISTING EXPERIMENT. 
experimentVARS <- paste0("ecoregion=", ecoregion,"_",
                         "months=", min(includeMonths), "-", max(includeMonths), "_",
                         "sizeMin=", minFireSize
                         )

experimentDir <- paste0("Figures/ignitionCompare_", experimentVARS, "/")
if(!dir.exists(experimentDir)){
  dir.create(experimentDir)
}


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

# Time difference tolerance in JDays for fires to be comparable
JDayTol       <- 10
minSampleSize <- 50
elevationBins <- seq(-300, 3900, by=300)
nBins         <- length(elevationBins)

# Get the fire occurence day of year for ignitions information 
JDay       <- FPA_FOD$DISCOVERY_DOY
JDaysArray <- sort(unique(JDay))
nDays      <- length(JDaysArray)

# Find the edges of the julain day experiment using a histogram of JDay counts
# TODO: Do this by ignition type 
png(filename="Figures/summary/JDayIgnitionCounts.png",
    res=250, height=2000, width=3700)
hist(JDay, breaks = JDaysArray, main="Julian Day Ignition Counts: 1992-2015", las=1)
hist(JDay[!lightningMask], breaks = JDaysArray, col="orange", add=T, xlim=c(0, 370))
span <- as.numeric(quantile(JDay, c(0.10, 0.90) ) ) 
abline(v=span, col="red", lwd=3)
legend("topleft", legend=c("Total Ignitions", "Human Ignitions"),
       fill=c("white", "orange"), bty="n", cex=2)

dev.off()

# Define the range explicitly 
minJDay <- span[1]
maxJDay <- span[2]

# Create expiment data storage arrays
loopJDays <- minJDay:maxJDay
nLoop <- length(loopJDays)
JDayFireCount <- rep(NA, nLoop)

# Empty storage dataframes creation
nCol <- (nBins-1)
a <- array(NA, dim= c(nLoop,  nCol) )
df_dummy <- data.frame(a)
row.names(df_dummy) <- loopJDays

df_difference_of_means <- df_dummy

# Space for other will be stored here
df_H_wilkes_p <- df_dummy
df_H_wilkes_W <- df_dummy

df_L_wilkes_p <- df_dummy
df_L_wilkes_W <- df_dummy

# Difference of the means test. 
df_t_test <- df_dummy

# Loop through each day and perform experiment
# TODO: Should I loop over every JDay or by the tolerance? Probably Tol.
# TODO: otherwise we are reusing a lot of data. 
for (i in 1:nLoop){ 
  
  # Subset fires to a given julain day (or maybe range of +-10 or something?)
  JDayMask <- abs(JDay - loopJDays[i]) < JDayTol
  JDayFireCount[i] <- sum(JDayMask)
  
  # Subset the dataframe by those that are within the JDay tolerance
  df <- FPA_FOD[JDayMask,]
  lMask <- df$STAT_CAUSE_DESCR == "Lightning"
  
  # Default cover is orange until seen that it is lightning 
  COL <- rep("orange", sum(JDayMask))
  COL[lMask] <- "gray"
  
  # TODO: boxplot of these data by elevation bin! 
  # TODO: Also show the 
  jDayFilename <- paste0("Figures/JDayIgnitionTvsElevation/JDay_",
                         loopJDays[i],".png")
  png(filename=jDayFilename, res=250, height=2000, width=4000)
  
  plot(df$elevation, df$t2m, 
       xlab="Elevation [m]", ylab="Temperature [C]",
       col=COL, pch=19, cex=0.5, 
       ylim=c(0, 30))
  abline(v=elevationBins, lty=2)
  title(paste("Julain Day Center:", loopJDays[i], "+-", JDayTol))
  
  # All elevations difference of means 
  meanLT <- mean(df$t2m[lMask])
  meanHT <- mean(df$t2m[!lMask])
  
  abline(h=meanLT, col="gray", lwd=3)
  abline(h=meanHT, col="orange", lwd=3)
  
  allDT <- meanLT - meanHT
  
  dev.off()
  
  # Now account for elevation 
  # IDEA: Adjust temperature based on elevation using standard lapse rate?
  sampleElavations <- df$elevation
  
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
  for ( m in 1:(nBins-1) ){
    
    # Only bother with these calculations when there is enough data in each bin. 
    if(fullBins[m]){
    
      elevMask <- sampleElavations >= elevationBins[m] & 
                  sampleElavations < elevationBins[m+1]
      
      # Subset the data by elevation. Make these arrays ready to use
      H <- df$t2m[elevMask & !lMask]
      L <- df$t2m[elevMask & lMask]
      
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
      
      # Test the two distributions to see if the means are different. 
      df_t_test[i,m] <- t.test(H, L, alternative = "two.sided")$p.value
    }
    
  } # end of elvation bin looping 
  
  
}


################################################################################
# Histogram of differences chunked by Julain day and elevation
################################################################################
pdf(file="Figures/summary/diffence_of_means_by_elevation_bin.pdf", 
    width=12, height=20)
par(mfrow=c(5,2))
# TODO: Show the total number of fires that went into the calculation for
# TODO: each elevation bin. Segregated by ignition type. 
# columns 3:12
for (i in 3:12){
  
  # Figure out the correct elevation bin associated with each row. 
  binLabel <- paste(elevationBins[i], "-", elevationBins[i+1], "meters")
  
  # Show consistent size of hist 
  hist(df_difference_of_means[, i], xlim=c(-5,5), ylim = c(0, 40),
       main="", ylab="", xlab="Difference of sample means")
  abline(v=0, lty=2)
  title(binLabel)
  
}

dev.off()

################################################################################
# Plot up these differences over Julain dates, one elevation curve at a time
################################################################################
png(filename = "JDay_diff_Series.png", width=3600, height=1500, res=250)
par(mar=c(4,6,4,6))

library(fields)
maxValue <- max(abs(df_mean), na.rm=T)
ylim <- c(maxValue * -1, maxValue)
  
plot(loopJDays, df_mean[,8], pch="", bty="n", ylim=ylim, 
     ylab="Detla Temperature [C]",
     xlab="Julain Day", 
     cex.lab=2)
title("(Mean Lightning Ignition Temperature) - (Mean Human Ignition Temperature)",
      cex.main=2)

# Place lines for different elevation bands on blank plot 
nLines <- dim(df_mean)[2]
lineColors <- terrain.colors(nLines)

for (r in 1:nLines){
  lines(loopJDays, df_mean[,r], col=lineColors[r], lwd=3)
}

abline(h=0, lwd=1, lty=2)

# Legend to understand the colors relation to elevation
image.plot( legend.only=TRUE, zlim= range(elevationBins), col=lineColors) 

dev.off()

################################################################################
# Regular Distributions 
################################################################################


quartz(width=15, height=5)

png(filename=paste0("Figures/summary/ignitionTemperature_6_2_Fires>1000.png"),
    width=3000, height=1000, res=200)
par(mfrow=c(1,3), mar=c(4,5,4,4))
plot_ignition_weather_distribution("t2m", "Temperature [K]", humanMask)

# Show the locations of the fires
# TODO: Dot size is size of fire!
map("state", xlim=c(-125,-102), ylim=c(33,50))
points(FPA_FOD$LONGITUDE[humanMask], FPA_FOD$LATITUDE[humanMask], 
       col="orange", pch=".")
points(FPA_FOD$LONGITUDE[!humanMask], FPA_FOD$LATITUDE[!humanMask], 
       col="darkgray", pch=".")
title(paste("Fires locations"), cex.main=2)


# TODO: histogram of the JDay that these fires are occuring on. 
minDay <- min(FPA_FOD_subset$DISCOVERY_DOY)
maxDay <- max(FPA_FOD_subset$DISCOVERY_DOY)

hist(FPA_FOD$DISCOVERY_DOY[humanMask], 
     col=adjustcolor("orange",0.5),
     breaks=c(minDay:maxDay), las=1, xlab="Day of Year",
     ylim=c(0,30),
     border="transparent",
     main="Day of year ignition counts",
     cex.main=2)
hist(FPA_FOD$DISCOVERY_DOY[!humanMask], 
     col=adjustcolor("darkgray", 0.5), breaks=c(minDay:maxDay), add=T,
     border="transparent")
dev.off()

plot_ignition_weather_distribution("tp_lastMonth", "Total Precip", humanMask)



quartz()
plot(FPA_FOD_subset$windSpeed_MonthAfter, FPA_FOD_subset$FIRE_SIZE)

#windSpeed_MonthAfter



