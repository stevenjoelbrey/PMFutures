# plot_air_quality_wildfires.R

# ------------------------- Description ---------------------------------------
# This script is used to take a look at the wildfires that qualify as air quality
# relevent via assignment of HYSPLIT points with assign_HYSPLITPoints_to_FPAFOD.R
# This script is also responsible for creating dataframes merging the yearly 
# analysis files output by R/assign_HYSPLITPoints_to_FPAFOD.R

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(ggthemes)

dataDir <- "Data/HMS/" 
figureDir <- "Figures/HMS_assignments/"

# Load the merged yearly paired fire data created by R/merge_paired_fires.R
print("Loading data file that exists")
load(paste0(dataDir, "FPA_FOD_with_HP_2007_2015.RData"))
load(paste0(dataDir, "hysplitPoints_with_FPA_FOD_2007_2015.RData"))

# I do not want Alaska, HI, or PR, for plotting purposes. 
fire_states <- FPA_FOD$STATE
stateMask <- fire_states != "HI" & fire_states != "AK" & fire_states != "PR"
FPA_FOD <- FPA_FOD[stateMask, ]

fire_cause <- rep("", dim(FPA_FOD)[1])
fire_cause[FPA_FOD$STAT_CAUSE_DESCR != "Lightning" & FPA_FOD$STAT_CAUSE_DESCR != "Missing/Undefined"] <- "Human"
fire_cause[FPA_FOD$STAT_CAUSE_DESCR == "Missing/Undefined"] <- "Unknown"
fire_cause[FPA_FOD$STAT_CAUSE_DESCR == "Lightning"] <- "Lightning"

FPA_FOD$fire_cause <- fire_cause

# Lets get some basics out of the way. Wildfires per year, proportion associated
# with HMS analysis
# TODO: Make monthly
years <- unique(FPA_FOD$FIRE_YEAR)
nYears <- length(years)
n_fires     <- rep(NA, nYears)
n_HMS_fires <- rep(NA, nYears)
for (i in 1:nYears){
  yearMask <- FPA_FOD$FIRE_YEAR == years[i]
  n_fires[i] <- sum(yearMask)
  n_HMS_fires[i] <- sum(yearMask & (FPA_FOD$n_HP > 0) ) 
}

# Save the figure that shows if there is a trend in fires associated with AQ
# forecasts. Based on the way the data are generated it is not clear what this
# tells us. 
png(filename=paste0(figureDir, "percent_AQ_wildfires.png"), 
    res=250, height=1000, 
    width=2000)

plot(years, n_HMS_fires/n_fires*100, 
     xlab="year", ylab="% wildfires paired with AQ forecast", las=1)
lines(years, n_HMS_fires/n_fires*100)
title("Annual percent of FPA FOD wildfires paired with HYSPLIT Points")

dev.off()


################################################################################
# Show the number of FPA wildfires that are associated with air quality 
# forecasts in each region. 
################################################################################
png(filename=paste0(figureDir, "FPA_FOD_AQfires_count_by_region.png"), 
    res=250, height=1000, 
    width=1500)

# Mask out non west or southeast region, also wildfires with no hysplits
m <- (FPA_FOD$region != "") & (FPA_FOD$n_HP > 0)
df <- FPA_FOD[m,]

ggplot(df, aes(x=region, fill=fire_cause))+
  geom_bar(stat="count", width=0.5, position = "dodge")+
  scale_fill_manual(values = c("Lightning"="gray", "Human"="orange", "Unknown"="blue"))+
  guides(fill=guide_legend(title="Fire Cause"))+
  theme_tufte(ticks=T, base_size = 20)+
  ylab("Wildfire count")+
  ggtitle("FPA FOD wildires associated\nwith air quality forecasts")


dev.off()

################################################################################
# Show the number of Hysplit points that are associated with FPA FOD wildfires
# in each region. 
################################################################################
png(filename=paste0(figureDir, "HysplitPoints_wFPA_FOD_count_by_region.png"), 
    res=250, height=1000, 
    width=1500)

# Mask out non west or southeast region, also wildfires with no hysplits
m <- (hysplitPoints$region != "") 
df <- hysplitPoints[m,]

# Make an association factor 
wild <- rep("", sum(m))
wild[df$pairedWithFPAFOD] <- "Yes"
wild[!df$pairedWithFPAFOD] <- "No"
df$wild <- wild

ggplot(df, aes(x=region, fill=wild))+
  geom_bar(stat="count", width=0.5, position = "dodge")+
  #scale_fill_manual(values = c("Lightning"="gray", "Human"="orange", "Unknown"="blue"))+
  guides(fill=guide_legend(title="Accounted for in\nFPA FOD"))+
  theme_tufte(ticks=T, base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5))+ # Center the title
  ylab("Hysplit point count")+
  ggtitle("HYSPLIT points analysed")

dev.off()

################################################################################
# Map when and where the paired fire data occur in the US. 
# http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
################################################################################
AQ_fires <- FPA_FOD[FPA_FOD$n_HP > 0, ] 

png(filename=paste0(figureDir, "AQ_fires_mapped.png"), res=100, 
    width=2000, height=1200)

states <- map_data("state")
ggplot() + geom_polygon(data = states, aes(x=long, y = lat, group = group),
                        fill = NA, color = "gray") + 
  coord_fixed(1.3)+ # fixes the relationship between one unit in the y direction and one unit in the x direction.
  theme_tufte(ticks=T, base_size = 20)+
  geom_point(data = AQ_fires, aes(x = LONGITUDE, y = LATITUDE, 
                                  color = fire_cause, size = FIRE_SIZE)
                                  #shape=fire_cause)
             )+
  scale_colour_manual(values = c("Lightning"="gray", "Human"="orange", "Unknown"="blue"))
  
dev.off()  

################################################################################
# Plot the locations of HYSPLIT points that have been paired with FPA FOD and 
# those that have not. 
################################################################################

world <- map_data("world")


hysplitPoints$paried  <- as.factor(hysplitPoints$nFPAFODPaired > 0)

png(filename=paste0(figureDir, "paired_HYSPLITPoints_mapped.png"), res=100, 
    width=2000, height=1200)

ggplot() + geom_polygon(data = states, aes(x=long, y = lat, group = group),
                        fill = NA, color = "gray") + 
  coord_fixed(1.3)+ 
  theme_tufte(ticks=T, base_size = 20)+
  geom_point(data = hysplitPoints, aes(x = Lon, y = Lat, 
                                  color = paried), shape=19, size=0.5)+
  xlim(c(-125,-60))+
  ylim(c(27,49))+
  ggtitle("Hysplit points shaded by pairing with FPA FOD")
  # geom_polygon(data = world, aes(x=long, y = lat, group = group),
  #              fill = "white", color = NA) 

dev.off()  

################################################################################
# Plot the paired datasets on the same map. Color by dataset. Dots off alone
# in the wild represent errors as there should be no more than 10 km of seperation.
################################################################################

HP_paired <- hysplitPoints[hysplitPoints$nFPAFODPaired > 0, ]

png(filename=paste0(figureDir, "paired_HYSPLITPoints_mapped.png"), 
    width=2000, height=1200)

ggplot() + geom_polygon(data = states, aes(x=long, y = lat, group = group),
                        fill = NA, color = "gray") + 
  coord_fixed(1.3)+
  geom_point(data = HP_paired, aes(x = Lon, y = Lat), shape=19, size=0.5, col="blue")+
  geom_point(data = AQ_fires, aes(x = LONGITUDE, y = LATITUDE), shape=19, size=0.5, col="red")+
  theme_tufte(ticks=T, base_size = 20)+
  ylim(c(25,49.5))+
  xlim(c(-125,-60))


dev.off()

################################################################################
# HMS smoke duration hours as a function of 
################################################################################
png(filename=paste0(figureDir, "fire_size_vs.png"), width=1700, height=1400,
    res=200)

p <- ggplot(AQ_fires, aes(x=FIRE_SIZE, y=n_HP, col=fire_cause, size=totalDurationHours))
p + geom_point() +  
  scale_colour_manual(values = c("Lightning"="gray", "Human"="orange", "Unknown"="blue"))+
  theme_tufte(ticks=T, base_size = 20)+
  xlab("Fire Size (acres)") + ylab("Hysplit points associated with fire")

dev.off()


################################################################################
# Plot bar chart of total AQ wildfires and the proportion by ignition type in 
# the different ecoregions. 
################################################################################
df <- FPA_FOD[FPA_FOD$n_HP==0, ] # FPA_FOD[FPA_FOD$n_HP==0, ] AQ_fires
df$NA_L2CODE <- as.factor(df$NA_L2CODE)

# Make associatedion with air quality forecast a factor in the data 
fireType <- rep("AQ forecast", dim(df)[1])
fireType[df$n_HP == 0] <- "No AQ forecast"
fireType <- factor(fireType, levels = c("No AQ forecast", "AQ forecast"))
df$fireType <- fireType

png(filename=paste0(figureDir, "regional_counts.png"), width=2000, height=800,
    res=200)

g <- ggplot(df, aes(NA_L2CODE,fill=fire_cause))
g + geom_bar()+
  scale_fill_manual(values = c("Lightning"="gray", "Human"="orange", "Unknown"="blue"))+
  theme_tufte(ticks=T, base_size = 15)+
  #facet_wrap(~fireType)+
  #ylab("Wildfires associated with AQ forecasts")+
  xlab("Level II ecoregion")

dev.off()


  