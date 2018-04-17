# plot_ignition_types.R

# ------------------------- Description ---------------------------------------
# This script shows the number of wildfires by large region (west and southeast)
# and color codes them by the ignition count assigned in the FPA FOD. 

complexColors <- FALSE

library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(gridExtra)
library(grid)

dataFile <- paste0("Data/FPA_FOD/FPA_FOD_gridMet_1992-2015.RData")
load(dataFile)

# Create a column in the dataframe that will be able to store what large region
# a wildfire occurs in. 
FPA_FOD$REGION <- "none" # mark as none until a smart assignment is made

# Mask the wildfires in the west
load("Data/GIS/west_bounds.RData")
latMask <- (FPA_FOD$LATITUDE <= maxLat) & (FPA_FOD$LATITUDE >= minLat)
lonMask <- (FPA_FOD$LONGITUDE <= maxLon) & (FPA_FOD$LONGITUDE >= minLon)
mWest       <- latMask & lonMask
FPA_FOD$REGION[mWest] <- "West"
rm(minLat, maxLat, minLon, maxLon) # we are loading these again, clear for safety

# Mask the wildfires in the southeast 
load("Data/GIS/southeast_bounds.RData")
latMask <- (FPA_FOD$LATITUDE <= maxLat) & (FPA_FOD$LATITUDE >= minLat)
lonMask <- (FPA_FOD$LONGITUDE <= maxLon) & (FPA_FOD$LONGITUDE >= minLon)
mSoutheast <- latMask & lonMask
FPA_FOD$REGION[mSoutheast] <- "Southeast"

# Now get rid of any wildfire that is not in the SE or the west 
df <- FPA_FOD[FPA_FOD$REGION != "none", ]
rm(FPA_FOD) # remove the main array, df contains all the wildfires of interest

df$REGION <- factor(df$REGION, levels = c("West", "Southeast") ) 
df$month  <- factor(df$START_MONTH, levels=c(1,2,3,4,5,6,7,8,9,10,11,12)) 


# Make a bar chart of the monthly occurence of wildfires by region and show
# the types of igntiions responsible. 
q <- ggplot(df, aes(x=STAT_CAUSE_DESCR, fill = month))+
  geom_bar() +
  facet_wrap(~REGION)+
  theme_tufte(ticks=F, base_size = 18)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("FPA FOD wildfire causes 1992-2015")+
  xlab("Wildfire cause")

png(filename="Figures/regional_ignition_counts.png",
    width=2200, height=1100, res=180)
print(q)
dev.off()

# here is an idea for how to communicate burn area at the same time
# https://stackoverflow.com/questions/8612920/pie-charts-in-ggplot2-with-variable-pie-sizes


# Make colors for these cause factors custom. I need lightning to stand out. 
if(complexColors){
  
  cause <- factor(df$STAT_CAUSE_DESCR, levels=c("Campfire", "Equipment Use", 
                                                "Debris Burning", "Smoking", "Arson", 
                                                "Miscellaneous", "Railroad", "Children",
                                                "Powerline", "Fireworks", "Structure",
                                                "Missing/Undefined", "Lightning"))
  
  df$cause <- cause
  # A unique color for each cause but make lightning gray like always
  causeColors <- c(brewer.pal(n=12, name="Paired"), "gray")
  
} else{
  
  cause <- df$STAT_CAUSE_DESCR
  cause[cause != "Lightning"] <- "Human"
  
  # Make cause a factor
  cause <- factor(cause, levels=c("Human", "Lightning"))
  df$cause <- cause
  
  # Set the colors for the causes, in this case there are only two.
  causeColors <- c("orange", "gray")
  
  # Save this dataframe, this will be handy for other plotting scripts that 
  # require comparing these large regions. 
  save(df, file="Data/FPA_FOD/west_and_southeast_FPA_FOD.RData")
  
}

p <- ggplot(df, aes(x=DISCOVERY_DOY, fill = cause))+
  geom_bar(width=1) +
  scale_fill_manual(values=causeColors)+
  facet_grid(~REGION)+
  theme(panel.spacing = unit(1, "lines"))+
  theme_tufte(ticks=T, base_size = 18)+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill="Ignition Type")+
  #ggtitle("FPA FOD day of year count and causes | 1992-2015")+
  xlab("Day of year")+
  ylab("Wildfire Count")+
  theme(legend.position="top")
  #scale_x_discrete(expand = c(0,0))

png(filename="Figures/DOY_ignition_counts_by_cause.png",
    width=2200, height=1100, res=250)
print(p)
dev.off()

# Make a bar plot that shows the burn area by Month but is nicely aligned with 
# DOY
g <- ggplot(df, aes(x=month, y=FIRE_SIZE, fill = cause))+
  geom_col(width=1) +
  scale_fill_manual(values=causeColors)+
  facet_grid(~REGION)+
  theme(panel.spacing = unit(1, "lines"))+
  theme_tufte(ticks=T, base_size = 18)+
  theme(plot.title = element_text(hjust = 0.5))+
  labs(fill="Ignition Type")+
  #ggtitle("FPA FOD day of year count and causes | 1992-2015")+
  xlab("Month")+
  ylab("Burn Area (acres)")+
  guides(fill=FALSE)
  #scale_x_discrete(expand = c(0,0))

# Use grid.arrange function to combine the day of year igntion counts with the
# monthly burn area bar plot. 
png(filename="Figures/DOY_ignition_counts_and_burn_area.png",
    width=2000, height=2000, res=250)
grid.arrange(p, g, nrow = 2)
dev.off()
