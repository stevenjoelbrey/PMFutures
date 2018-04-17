# plot_figure_4.R

# This plots figure 4. 

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(ggthemes)

wildfireSize <- 0.4 # size of plotting dot

# Load the west and SE data only
load("Data/FPA_FOD/west_and_southeast_FPA_FOD.RData")

# Get rid of wildfires with unknown start
hasStart <- df$STAT_CAUSE_DESCR != "Missing/Undefined"

# We only want to see the large wildfires
sizeMask <- df$FIRE_SIZE >= 1000


df <- df[hasStart & sizeMask,]

# Make a map of the wildfires and regions
# http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
states <- map_data("state")

# Create subsets of states for the different regions
west <- subset(states, region %in% c("california", "oregon", "washington",
                                     "montana", "wyoming", "idaho", "nevada",
                                     "utah", "colorado", "new mexico", "arizona"))


load("Data/GIS/west_bounds.RData")
west_map <- ggplot(data = west) + 
  geom_polygon(aes(x = long, y = lat, group = group),  fill="white", color = "black") +
  coord_fixed(1.3) +
  guides(fill=FALSE)+
  theme_tufte(ticks=T, base_size = 18)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  geom_point(data = df, 
             mapping = aes(x = LONGITUDE, y = LATITUDE, color = cause),
             size=wildfireSize)+
  scale_color_manual(values=c("orange", "gray"))+
  guides(color=FALSE)+
  coord_cartesian(xlim=c(minLon+1, maxLon+2), 
                  ylim=c(minLat, maxLat+0.1))

load("Data/GIS/southeast_bounds.RData")
SE <- subset(states, region %in% c("florida", "georgia", "south carolina",
                                   "north carolina",
                                   "virginia", "west virginia", "kentucky",
                                   "tennessee", "alabama", "mississippi",
                                   "louisiana", "arkansas", "indiana",
                                   "illinois", "maryland"))
SE_map <- ggplot(data = SE) + 
  geom_polygon(aes(x = long, y = lat, group = group),  fill="white", color = "black") +
  coord_fixed(1.3) +
  guides(fill=FALSE)+
  theme_tufte(ticks=T, base_size = 18)+
  theme(axis.title=element_blank(),
        axis.text=element_blank(),
        axis.ticks=element_blank())+
  geom_point(data = df, 
             mapping = aes(x = LONGITUDE, y = LATITUDE, color = cause),
             size=wildfireSize)+
  scale_color_manual(values=c("orange", "gray"))+
  guides(color=FALSE)+
  coord_cartesian(xlim=c(minLon-1, maxLon+2), ylim=c(minLat, maxLat))


# Show the number of large wildfires by region starting with the west
count_plot <- ggplot(df, aes(x=month, fill = cause))+
  geom_bar() +
  scale_fill_manual(values=c("orange","gray"))+
  facet_wrap(~REGION)+
  theme_tufte(ticks=T, base_size = 18)+
  #ggtitle("FPA FOD wildfire causes 1992-2015")+
  xlab("Wildfire start month")+
  ylab("Number of wildfires")+
  theme(legend.position="top")

lay <- rbind(c(1,2),
             c(3,3))

png(filename="Figures/Figure_4.png",
    width=2000, height=2000, res=250)
grid.arrange(west_map, SE_map, count_plot, layout_matrix = lay)
dev.off()
