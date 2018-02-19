# ------------------------- Description ---------------------------------------
# This script is used to compare wildfires by ignition in a DOY, elevation, size
# fuel moisture parameter space. 
# NOTE: cowplot sometimes does not save figure when source() in RStudio is used. 
# NOTE: run all code by selecting and executing till you figure this out. 

library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(scales) 


ecoregionSelect <- 6.2
ecoRegionName   <- "Forested mountains" # "Mediterranean CA" # mediterranean
minSize         <- 1000 # minimum fire size to include
# Figure out maximum and minimum values for colorbar, make it consistent across
# ecoregions. 
x1 <- 0
x2 <- 25 # Ok to saturate at 40 when including the very small fires. 25 good for 1000

# Load the data that has griodMET appended. 
load("Data/FPA_FOD/FPA_FOD_gridMet_1992-2015.RData")

# Mask by ecoregion and lat and lon. 
latMask <- (FPA_FOD$LATITUDE <= 50) & (FPA_FOD$LATITUDE > 31)
lonMask <- (FPA_FOD$LONGITUDE <= -100) & (FPA_FOD$LONGITUDE > -127)
ecoregionMask <- FPA_FOD$NA_L2CODE == ecoregionSelect
hasStartMask <- FPA_FOD$STAT_CAUSE_DESCR != "Missing/Undefined"
sizeMask <- FPA_FOD$FIRE_SIZE >= minSize
elevationMask <- FPA_FOD$elevation > -80 # this is ~lowest elevation in death valley
m <- latMask & lonMask & ecoregionMask & hasStartMask & sizeMask & elevationMask

# In order for elevation scale to be consistent we need to set the ylim before
# creating the sibset dataframe
yMin <- -10
yMax <- max(FPA_FOD$elevation)

df <- FPA_FOD[m, ]
# Try making FIRE_SIZE_CLASS numeric via factor
df$numeric_SIZE_CLASS <- as.numeric(as.factor(df$FIRE_SIZE_CLASS ))



# maxSize <- max(df$FIRE_SIZE)
# if(minSize==0){
#   sizeClassBreaks <- c(0 , 10, 100, 300, 1000, 5000, 10000, 100000)
#   sizeTheme <- scale_size_area(labels = comma, 
#                                breaks=sizeClassBreaks, max_size = 10)
# }else{
#   sizeTheme <- scale_size_area(labels = comma, 
#                                limits=c(minSize, maxSize), max_size = 10)
# }

# Consistent breaks
sizeClassBreaks <- c(0 , 10, 100, 300, 1000, 5000, 10000, 100000)
sizeTheme <- scale_size_area(labels = comma, 
                             breaks=sizeClassBreaks, max_size = 10)

# Make the seperate dataframes based on ignition
df_L <- df[df$STAT_CAUSE_DESCR == "Lightning",]
df_H <- df[df$STAT_CAUSE_DESCR != "Lightning",]

# Create a color palette (NOT USED)
myPalette <- colorRampPalette(colors = c("blue", "orange"))

# Create a colorscale that is shared
sc <- scale_colour_gradient(limits=c(x1, x2), labels=comma, 
                            low = alpha("orange", 0.75), high = alpha("blue", 0.75) )

# Plot the lightning-ignited 
lightning <- ggplot(df_L, aes(x=DISCOVERY_DOY, y=elevation, 
                              color=fuel_moisture_1000hr,
                              size=FIRE_SIZE)) +
  geom_point(aes(x=DISCOVERY_DOY, y=elevation, 
                 color=fuel_moisture_1000hr,
                 size=FIRE_SIZE), 
             alpha=0.75)+
  sc +  
  sizeTheme+
  theme_bw()+
  ggtitle(paste(ecoRegionName, "lightning-ignited wildfires"))+
  xlab("Day of year")+
  ylab("Elevation [m]")+
  ylim(yMin,yMax)+
  labs(color='1000 hr \nfuel moisture %') +
  labs(size='Fire size \n(acres)')+
  theme(axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=16, face="bold"),
        legend.title = element_text(size=16, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.key.width=unit(2,"line"),
        title = element_text(size=15, face="bold")
  )+
  guides(size = guide_legend(order=2))

print("Ran lightning figure")

# Plot the human-ignited 
human <- ggplot(df_H, aes(x=DISCOVERY_DOY, y=elevation, 
                          color=fuel_moisture_1000hr,
                          size=FIRE_SIZE)) + 
  geom_point(aes(x=DISCOVERY_DOY, y=elevation, 
                 color=fuel_moisture_1000hr,
                 size=FIRE_SIZE),
             alpha=0.75)+
  sc +  
  sizeTheme+
  theme_bw()+
  ggtitle(paste(ecoRegionName, "human-ignited wildfires"))+
  xlab("Day of year")+
  ylab("Elevation [m]")+
  ylim(yMin,yMax)+
  labs(color='1000 hr \nfuel moisture %') +
  labs(size='Fire size \n(acres)')+
  theme(axis.text=element_text(size=12, face="bold"),
        axis.title=element_text(size=16, face="bold"),
        legend.title = element_text(size=16, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.key.width=unit(2,"line"),
        title = element_text(size=15, face="bold")
  )+
  guides(size = guide_legend(order=2))

print("Getting to code to save figure.")
saveName <- paste0("Figures/elevation_DOY_space/minSize_",
                   minSize,"_ecoregion_", ecoregionSelect,"_acres.png")
png(filename=saveName, width=4000, height=1500, res=250)
plot_grid(lightning, human, labels = "AUTO")
print("Ran save figiure code")
dev.off()