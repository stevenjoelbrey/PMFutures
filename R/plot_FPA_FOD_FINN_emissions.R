# plot_FPA_FOD_FINN_emissions.R
#
# ------------------------- Description ---------------------------------------
# This script plots FPA FOD PM25 emissions as estimated by the FINN modeling 
# framework. FPA FOD fires were used instead of MODIS dections and PM25 
# emissions for individual fires were estimated by Christine Weidinmyer. 
#
# Loads data created by:
#   R/process_FPA_FOD_for_FINN.R
#   R/sanity_check_FINN_procssed_data.R

# ------------------------- Arguments ------------------------------------------
dataInFile   <- "Data/FINN/processed/FPA_FOD_wFINN-PM25_2002-2015.RData" 
figureDir <- paste0("Figures/FINN_emissions/")

a <- 1 # alpha for color in plots 
lightningCol <- adjustcolor("cyan4", alpha.f = a)
humanCol     <- adjustcolor("orange", alpha.f = a)
unknownCol   <- adjustcolor("gray", alpha.f = a)

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(ggthemes)
library(dplyr)

# If the figure dirctory does not exist, create it
if(!dir.exists(figureDir)){
  dir.create(figureDir)
}

# Load the CONUS FPA FOD wildfires where FINN PM2.5 emission estimates have
# been made. 
load(dataInFile)

# Handle wildfire cause
fire_cause <- rep("-", dim(FPA_FOD)[1])
fire_cause[FPA_FOD$STAT_CAUSE_DESCR != "Lightning" & FPA_FOD$STAT_CAUSE_DESCR != "Missing/Undefined"] <- "Human"
fire_cause[FPA_FOD$STAT_CAUSE_DESCR == "Missing/Undefined"] <- "Missing/Undefined"
fire_cause[FPA_FOD$STAT_CAUSE_DESCR == "Lightning"] <- "Lightning"
# Make this a factor and assign it
fire_cause <- factor(fire_cause, levels=c("Lightning", "Human", "Missing/Undefined"))
FPA_FOD$fire_cause <- fire_cause


# ----------------------------Make regions mask---------------------------------
region <- rep("Other", dim(FPA_FOD)[1])
load("Data/GIS/west_bounds.RData")
latMask <- minLat <= FPA_FOD$LATITUDE  & maxLat >= FPA_FOD$LATITUDE
lonMask <- minLon <= FPA_FOD$LONGITUDE & maxLon >= FPA_FOD$LONGITUDE
region[latMask & lonMask] <- "West"

load("Data/GIS/southeast_bounds.RData")
latMask <- minLat <= FPA_FOD$LATITUDE  & maxLat >= FPA_FOD$LATITUDE
lonMask <- minLon <= FPA_FOD$LONGITUDE & maxLon >= FPA_FOD$LONGITUDE
region[latMask & lonMask] <- "Southeast"

# Make it a factor for plotting 
region <- factor(region, levels=c("West", "Southeast", "Other"))
FPA_FOD$region <- region

# Mask out non west or southeast region areas for plotting "df"
m <- (FPA_FOD$region != "Other") 
df <- FPA_FOD[m,]

# Sanity checks
sum(df$PM25[df$region=="West" & df$fire_cause=="Lightning"])
sum(df$PM25[ (df$region=="Southeast") & (df$fire_cause=="Human")])

# Create a new category, sum PM25 for a given region ignitionType combo
g_sum <- ggplot(df %>%
                  group_by(region, fire_cause) %>% 
                  summarise(PM25_SUMMED = sum(PM25)),
                aes(x = region, y = PM25_SUMMED, fill = fire_cause))

g <- g_sum + geom_bar(stat = "identity", position = "dodge")+ 
  scale_fill_manual(values = c("Lightning"=lightningCol, 
                               "Human"=humanCol, 
                               "Missing/Undefined"=unknownCol))+
  guides(fill=guide_legend(title="Ignition type"))+
  theme_tufte(ticks=T, base_size = 20)+
  #theme(legend.position="top")+
  theme(plot.title = element_text(hjust = 0.5))+ # Center the title
  ylab(expression("FINN PM"[2.5]*" [kg]"))+
  xlab("")
  #ggtitle(expression("Wildfire PM"[2.5]*" emissions"))

png(filename=paste0(figureDir, "FPA_FOD_PM25_emissions_by_ignition.png"), 
    res=250, height=1300, width=1700)
print(g)
dev.off()


library(rgdal)
library(rasterVis)
library(RColorBrewer)

# Create a grid to store PM25 emissions
lon  <- seq(-125, -65, by=0.25)
lat  <- seq(24, 50, by=0.25)
nLon <- length(lon)
nLat <- length(lat)

nFires <- dim(FPA_FOD)[1]

# Create a grid to store wildfire PM25 emissions
PM25_h <- matrix(data=0, nrow=nLon, ncol=nLat)
PM25_L <- PM25_h
PM25_M <- PM25_L

for (i in 1:nFires){

  # Get fire attributes
  PM25    <- FPA_FOD$PM25[i]
  fireLat <- FPA_FOD$LATITUDE[i]
  fireLon <- FPA_FOD$LONGITUDE[i]
  cause   <- FPA_FOD$fire_cause[i]

  xi <- which.min(abs(fireLon - lon))
  yi <- which.min(abs(fireLat - lat))

  if(cause=="Human"){
    PM25_h[xi, yi] <- PM25_h[xi, yi] + PM25
  }else if(cause=="Lightning"){
    PM25_L[xi, yi] <- PM25_L[xi, yi] + PM25
  } else{
    PM25_M[xi, yi] <- PM25_M[xi, yi] + PM25
  }

  if(i %%1000 == 0){
    print(paste("Percent complete:", i/nFires*100))
  }

}

PM25[PM25==0] <- NA
r_h <- raster(PM25_h)
extent(r_h) <- c(min(lon), max(lon), min(lat), max(lat))

plot(r_h)


# https://stackoverflow.com/questions/33227182/how-to-set-use-ggplot2-to-map-a-raster






