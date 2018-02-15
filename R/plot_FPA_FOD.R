# plot_FPA_FOD.R

# ------------------------- Description ---------------------------------------
# This script produces a summary plot of FPA-FOD. 

load("Data/FPA_FOD/FPA_FOD_1992_2015.RData")

library(raster)
library(randomcoloR)

# # What region are you investigating? 
# minLat <- 31
# maxLat <- 49
# minLon <- -125
# maxLon <- -100

# Mask out non-CONUS
noThanks <- c("PR", "HI", "AK") 
stateMask <- FPA_FOD$STATE %in% noThanks
FPA_FOD <- FPA_FOD[!stateMask,]


lat <- FPA_FOD$LATITUDE
lon <- FPA_FOD$LONGITUDE
cause <- as.character(FPA_FOD$STAT_CAUSE_DESCR)
dot_color <- rep("blue", length(lat))
dot_color[cause == "Lightning"] <- "gray"
dot_color[cause != "Missing/Undefined" & cause != "Lightning"] <- "orange"


png(file="Figures/CONUS_FPA_FOD_fires.png", width=6100, height=4000, res=250)
map("state", xlim=c(-125, -65))
points(lon, lat, pch=".", col=dot_color)
map("state", add=T)
title("US Wildfires | 1992 - 2015", cex.main=4)
abline(v=-100, lty=3, lwd=3)

# # Cover part of line we do not want with Canada and Mexico
# canada <- getData("GADM",country="CAN", level=1)
# mexico <- getData("GADM",country="MEX", level=1)
# plot(canada, add=T, col="white")
# plot(mexico, add=T, col="white")
legend("bottomleft", 
       legend =c("Lightning", "Human", "Unknown"), 
       fill=c("gray", "orange", "blue"),
       bty="n",
       cex=3,
       title="Ignition Type"
       )
dev.off()

################################################################################
# Now I am interested in only showing the locations of class C fires. Probably
# only west of 100 deg W.
################################################################################

# https://www.nwcg.gov/term/glossary/size-class-of-fire
# Class A - one-fourth acre or less;
# Class B - more than one-fourth acre, but less than 10 acres;
# Class C - 10 acres or more, but less than 100 acres;
# Class D - 100 acres or more, but less than 300 acres;
# Class E - 300 acres or more, but less than 1,000 acres;
# Class F - 1,000 acres or more, but less than 5,000 acres;
# Class G - 5,000 acres or more.

# Make charactor not factor
FPA_FOD$FIRE_SIZE_CLASS <- as.character(FPA_FOD$FIRE_SIZE_CLASS)

# Store the array 
size_class <- FPA_FOD$FIRE_SIZE_CLASS

# Save the unique, in alphabetical order
unique_classes <- sort(unique(FPA_FOD$FIRE_SIZE_CLASS))
nClasses <- length(unique_classes)

# Assign a color to each month of the year
monthColorPalette <- distinctColorPalette(12)
#colorRampPalette(c("blue", "red") )(12)

# Loop through months of year and assigned color 
monthColor <- rep("", length(lat))
for(i in 1:12){
  monthColor[i == FPA_FOD$START_MONTH] <- monthColorPalette[i]
}

classCex <- rev(c(0.4, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1))

# Show each panel of fire locations, by class, by ignition. 
png(file="Figures/CONUS_FPA_FOD_fires_w_class.png", width=1800, height=3700, res=250)
par(mfrow=c(7,2))

# Plot the large fires first
for(i in nClasses:1){
  
  m <- unique_classes[i] == size_class
  
  # Lightning 
  map("state", xlim=c(-125, -65), mar=c(0,0,3,0))
  lightningMask <- (cause == "Lightning") & m 
  points(lon[lightningMask], lat[lightningMask], 
         col=monthColor[lightningMask], 
         pch=19, cex=classCex[i])
  title(paste("Lightning-ignited class", unique_classes[i], "wildfires. n = ",
              sum(lightningMask)), 
        line=0)
  
  legend("bottomright",
         title="month",
         legend=c(1:12),
         col=monthColorPalette,
         pch=19,
         bty="n",
         cex=0.7)
  
  # Human 
  map("state", xlim=c(-125, -65), mar=c(0,0,3,0))
  humanMask <- (cause != "Missing/Undefined") & (cause != "Lightning") & m
  points(lon[humanMask], lat[humanMask], 
         col=monthColor[humanMask], 
         pch=19, cex=classCex[i])
  title(paste("Human-ignited class", unique_classes[i], "wildfires. n=", 
              sum(humanMask)), 
        line=0)

}

dev.off()


################################################################################
# Plot the same for western US only and 3 biggest size class only
################################################################################
# TODO: THIS NEEDS TO BE CHECKED BUT DON'T BOTHER IF YOU DONT USE THE FIGURE

# We want to investigate the largest fires only. Use this mask to also get rid
# of fires in the east. 
lonMask <- lon <= -100
classMask <- size_class == "G" | size_class == "F" 

# Mask lightning ignitions of G and F Class in the west
lightningMask <- (cause == "Lightning") 
m1            <-  classMask & lightningMask & lonMask

# Mask human ignitions of G and F Class in the west
humanMask <- (cause != "Missing/Undefined") & (cause != "Lightning") 
m2        <-  classMask & humanMask & lonMask

# Now I want to make a histogram showing the months these fires occur. First
# count the number of fires in each month for each ignition
m1_month_count <- rep(NA, 12)
m2_month_count <- rep(NA, 12)
for (i in 1:12){
  m1_month_count[i] <- sum(FPA_FOD$START_MONTH[m1]==i)
  m2_month_count[i] <- sum(FPA_FOD$START_MONTH[m2]==i)
}

YLIM <- c(0, max(c(m1_month_count, m2_month_count)))

# Plot the locations
png(file="Figures/large_western_US_FPA_FOD_fires.png", 
    width=2300, height=2000, res=250)
par(mfrow=c(2,2))

map("state", xlim=c(-125, -100), mar=c(0,0,4,0))
points(lon[m1], lat[m1], pch=19, col=monthColor[m1], cex=0.5)
title(paste0("Lightning-ignited fires > 1000 acres"), line=0, cex.main=1.5)

par(mar=c(4,6,4,0))
barplot(height = m1_month_count, names.arg=1:12, ylim=YLIM,
        xlab="Month", ylab="Number of fires",
        space=0, col=monthColorPalette,
        bty="n", las=1, cex.lab=2)
title( paste0("n=", sum(m1)))#, line=0)

# Plot the large western US fires started by people
map("state", xlim=c(-125, -100), mar=c(0,0,4,0))
points(lon[m2], lat[m2], pch=19, col=monthColor[m2], cex=0.5)
title(paste0("Human-ignited fires > 1000 acres"), line=0, cex.main=1.5)

par(mar=c(4,6,4,0))
barplot(height = m2_month_count, names.arg=1:12, ylim=YLIM, 
        xlab="Month", ylab="Number of fires",
        space=0, col=monthColorPalette,
        bty="n", las=1, cex.lab=2)
title( paste0("n=", sum(m2)))#, line=0)

dev.off()

