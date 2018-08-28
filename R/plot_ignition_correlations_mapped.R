# plot_ignition_correlations_mapped.R

# plot the correlation in interannual burn area by ignition type for all US 
# ecoregions.

library(sp)
library(maps)
require(rgdal)
require(ggplot2)
library(ggthemes)
library(rgeos)
library(maptools)

# Load the ecoregion spatial polygons into the work space
load("Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData")

# Assign human and lightning ignition burn area correlations
# TODO: Do this dynamically for all ecoregions with FPA FOD wildfires occuring.
# TODO: Currently these numbers have been hard coded by looking at the figures
# TODO: that show the spearman r values for these burn area quantities by ecoregion,
# TODO: Those correlations are calculated and plotted by R/plot_interannual_variability.
# TODO: The numbers typed in the arrays below were taken from the table of correlations
# TODO: in the SI and typed in here on 4/18/2018. 
ecoregion             <- c(6.2,  10.1, 11.1, 8.3,  8.4,  8.5,  15.4)
ignition_correlations <- c(0.79, 0.83, 0.19, 0.51, 0.67, 0.66, 0.3)

# Figure out where the ecoregions listed above live in the array of ecoregions 
# stored in the polygon object.
SPDF@data$NA_L2CODE[match(ecoregion, SPDF@data$NA_L2CODE)] == ecoregion

# URL for notes on how to convert to ggplot workability
# http://mazamascience.com/WorkingWithData/?p=1494

# Load map_data needed to complete a nice map. 
states <- map_data("state") # for showing USA states
world <- map_data("world") # needs to have USA removed or data over CONUS will be hidden
notUSA <- world[world$region != "USA",] # for covering ecoregions that go outside of CONUS

# I am going to want to place the west and southeast region borders on the map
# for full context of the areas data being used in the correlation plot. Load
# and save the corners here for easy segment plotting later. 
load("Data/GIS/southeast_bounds.RData")
SEBorder <- data.frame(minLat=minLat, maxLat=maxLat,
                       minLon=minLon, maxLon=maxLon)
load("Data/GIS/west_bounds.RData")
WBorder <- data.frame(minLat=minLat, maxLat=maxLat,
                      minLon=minLon, maxLon=maxLon)  


# Convert the spatialpolygondataframe to an ordinary dataframe that retains 
# spatial info. Regular dataframe required for plotting in ggplot2. 
spdf.map <- ggplot2::fortify(SPDF, region = "NA_L2CODE")

# Initialize correltions to NA, that way ecoregions were they are not assigned will be left blank
spdf.map$cor <- NA 

# Assign ecoregion correlations by matching the ID of the spatialpolygon object
# that has been converted to a long lat dataframe for ggplot2 mapping. 
for (eco in ecoregion){
  m <- eco == spdf.map$id
  spdf.map$cor[m] <- ignition_correlations[eco==ecoregion]
}

# Create a color pallete Jeff likes..
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(100), na.value="white")#, limits=c(1, 8))

cor_map <- ggplot(spdf.map, aes(x = long, y = lat, group = group, fill = cor)) +
  coord_fixed(1.3) +
  geom_polygon(colour = "white", size = 0.5, aes(group = group)) +
  #sc+
  scale_fill_continuous(na.value="white")+
  geom_polygon(data=states, aes(x = long, y = lat, group = group),  
               fill="transparent", color = "black", linetype = 1, size = 0.2)+
  # Cover the parts of the ecoregions that go outside the US
  geom_polygon(data=notUSA, aes(x = long, y = lat, group = group),  
               fill="white", color = "transparent")+
  # Add line segmengs that show the SE
  geom_segment(aes(x = SEBorder$minLon, y = SEBorder$minLat, xend = SEBorder$minLon, yend = SEBorder$maxLat), linetype=3, size=0.2)+
  geom_segment(aes(x = SEBorder$minLon, y = SEBorder$maxLat, xend = SEBorder$maxLon, yend = SEBorder$maxLat), linetype=3, size=0.2)+
  geom_segment(aes(x = SEBorder$maxLon, y = SEBorder$maxLat, xend = SEBorder$maxLon, yend = SEBorder$minLat), linetype=3, size=0.2)+
  geom_segment(aes(x = SEBorder$maxLon, y = SEBorder$minLat, xend = SEBorder$minLon, yend = SEBorder$minLat), linetype=3, size=0.2)+
  # Add line segments showing western border
  geom_segment(aes(x = WBorder$minLon, y = WBorder$minLat, xend = WBorder$minLon, yend = WBorder$maxLat), linetype=3, size=0.2)+
  geom_segment(aes(x = WBorder$minLon, y = WBorder$maxLat, xend = WBorder$maxLon, yend = WBorder$maxLat), linetype=3, size=0.2)+
  geom_segment(aes(x = WBorder$maxLon, y = WBorder$maxLat, xend = WBorder$maxLon, yend = WBorder$minLat), linetype=3, size=0.2)+
  geom_segment(aes(x = WBorder$maxLon, y = WBorder$minLat, xend = WBorder$minLon, yend = WBorder$minLat), linetype=3, size=0.2)+
  
  labs(fill="r")+
  
  theme_tufte(ticks=T, base_size = 18)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  theme(legend.position=c(0.92, .3))+
  coord_cartesian(xlim=c(-123, -68), ylim=c(26, 50))

# Save the map that shows the correlations spatially
png(filename="Figures/Figure_7.png",
    width=1200*1.3, height=980, res=250)
print(cor_map)
dev.off()
