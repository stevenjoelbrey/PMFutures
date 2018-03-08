# ------------------------- Description ---------------------------------------
# This script is used to compare wildfires by ignition in a DOY, elevation, size
# fuel moisture parameter space. 
# NOTE: cowplot sometimes does not save figure when source() in RStudio is used. 
# NOTE: run all code by selecting and executing till you figure this out. 

library(ggplot2)
library(RColorBrewer)
library(scales) 
library(ggthemes)

regionName <- "east"

if(regionName == "east_and_west"){
  ecoregionSelect <- c(6.2, 10.1, 11.1, 8.5, 8.4, 8.3, 15.4) 
}else if(regionName == "east"){
  ecoregionSelect <- c(8.5, 8.4, 8.3, 15.4) 
} else if(regionName == "west"){
  ecoregionSelect <- c(6.2, 10.1, 11.1) 
}

minSize         <- 0 # minimum fire size to include

# Load the data that has griodMET appended. 
load("Data/FPA_FOD/FPA_FOD_gridMet_1992-2015.RData")

# Mask by ecoregion and lat and lon. 
# TODO: Do we even need to subset by lat and lon since we are doing ecoregions? 
latMask <- (FPA_FOD$LATITUDE <= 50) & (FPA_FOD$LATITUDE > 25)
lonMask <- (FPA_FOD$LONGITUDE <= -70) & (FPA_FOD$LONGITUDE > -127)
ecoregionMask <- FPA_FOD$NA_L2CODE %in% ecoregionSelect
hasStartMask <- FPA_FOD$STAT_CAUSE_DESCR != "Missing/Undefined"
sizeMask <- FPA_FOD$FIRE_SIZE >= minSize
elevationMask <- FPA_FOD$elevation > -80 # this is ~lowest elevation in death valley
m <- latMask & lonMask & ecoregionMask & hasStartMask & sizeMask & elevationMask

# Subset the fire data 
df <- FPA_FOD[m, ]
nRow <- dim(df)[1]

# Give the ecoregions nice names
ecoRegionNames <- rep("", nRow)

ecoRegionNames[df$NA_L2CODE==6.2]  <- "Forested mountains"
ecoRegionNames[df$NA_L2CODE==10.1] <- "High deserts"
ecoRegionNames[df$NA_L2CODE==11.1] <- "Mediterranean California"

ecoRegionNames[df$NA_L2CODE==8.3] <- "Southern Plains"
ecoRegionNames[df$NA_L2CODE==8.4] <- "Ozark"
ecoRegionNames[df$NA_L2CODE==8.5] <- "Mississippi Alluvial" 
ecoRegionNames[df$NA_L2CODE==15.4] <- "Everglades"


ecoRegionNames <-  base::factor(ecoRegionNames, 
                                levels = c("Forested mountains",
                                           "High deserts",
                                           "Mediterranean California",
                                           "Southern Plains",
                                           "Ozark",
                                           "Mississippi Alluvial",
                                           "Everglades"))
df$ecoRegionNames <- ecoRegionNames

# Add ignitionType column as a factor
ignition <- rep("Human-ignition", length(df$STAT_CAUSE_DESCR) )
ignition[df$STAT_CAUSE_DESCR=="Lightning"] <- "Lightning-ignition"
df$ignitionType <- base::factor(ignition, levels = c("Lightning-ignition", 
                                                     "Human-ignition"))


# Consistent breaks, fire size class
sizeClassBreaks <- c(0 , 10, 100, 300, 1000, 5000, 10000, 100000)
sizeTheme <- scale_size_area(labels = comma, 
                             breaks=sizeClassBreaks, 
                             max_size = 14)


# Create a color palette (NOT USED)
myPalette <- colorRampPalette(colors = c("blue", "orange"))

# Create a colorscale that is shared
sc <- scale_colour_gradient( labels=comma, limits=c(2, 38),
                            low = "orange", 
                            high = "blue")
if(minSize > 100){
  sc <- scale_colour_gradient( labels=comma, limits=c(min(df$fuel_moisture_1000hr), 
                                                      max(df$fuel_moisture_1000hr)),
                               low = "orange", 
                               high = "blue")
}

# Plot the space for the desired ecoregions 
p <- ggplot(df, aes(x=DISCOVERY_DOY, y=elevation, 
                    color=fuel_moisture_1000hr,
                    size=FIRE_SIZE)) +
  geom_point(aes(x=DISCOVERY_DOY, y=elevation, 
                 color=fuel_moisture_1000hr,
                 size=FIRE_SIZE), 
             alpha=0.75)+
  facet_grid(ignitionType~ecoRegionNames)+ # Note was facet_wrap for west 
  sizeTheme+
  sc +  
  theme_bw()+
  xlab("Day of year")+
  ylab("Wildfire elevation [m]")+
  labs(color='1000 hr \nfuel moisture %') +
  labs(size='Fire size \n(acres)')+
  guides(size = guide_legend(order=2))+
  theme_tufte(ticks=T, base_size = 25)
  
print("Getting to code to save figure.")
saveName <- paste0("Figures/elevation_DOY_space/minSize=",
                   minSize,"_",regionName,".png")
png(filename=saveName, width=4000, height=2700, res=250)
p
print("Ran save figiure code")
dev.off()

