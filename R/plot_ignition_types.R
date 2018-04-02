# plot_ignition_types.R

library(ggplot2)
library(ggthemes)

dataFile <- paste0("Data/FPA_FOD/FPA_FOD_gridMet_1992-2015.RData")
load(dataFile)

FPA_FOD$REGION <- "none"

# Mask the wildfires in the west
load("Data/GIS/west_bounds.RData")
latMask <- (FPA_FOD$LATITUDE) <= maxLat & (FPA_FOD$LATITUDE >= minLat)
lonMask <- (FPA_FOD$LONGITUDE) <= maxLon & (FPA_FOD$LONGITUDE >= minLon)
m       <- latMask & lonMask
FPA_FOD$REGION[m] <- "West"
rm(minLat, maxLat, minLon, maxLon)

# Mask the wildfires in the southeast 
load("Data/GIS/southeast_bounds.RData")
latMask <- (FPA_FOD$LATITUDE) <= maxLat & (FPA_FOD$LATITUDE >= minLat)
lonMask <- (FPA_FOD$LONGITUDE) <= maxLon & (FPA_FOD$LONGITUDE >= minLon)
m       <- latMask & lonMask
FPA_FOD$REGION[m] <- "Southeast"

# Now get rid of any wildfire that is not in the SE or the west 
df <- FPA_FOD[FPA_FOD$REGION != "none", ]
rm(FPA_FOD)

df$REGION <- factor(df$REGION, levels = c("West", "Southeast") ) 
df$month <- factor(df$START_MONTH) 

g <- ggplot(df, aes(x=STAT_CAUSE_DESCR, fill = month))+
  geom_bar() +
  facet_wrap(~REGION)+
  theme_tufte(ticks=F, base_size = 18)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5))+
  ggtitle("FPA FOD wildfire causes 1992-2015")+
  xlab("Wildfire cause")

png(filename="Figures/regional_ignition_counts.png",
    width=2200, height=1100, res=180)
print(g)
dev.off()

# here is an idea for how to communicate burn area at the same time
# https://stackoverflow.com/questions/8612920/pie-charts-in-ggplot2-with-variable-pie-sizes






