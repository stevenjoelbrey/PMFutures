# plot_FINN_FPAFOD_wildfires.R

# ------------------------- Description ---------------------------------------
# This script plots FINN wildfires and thier pairing with FPA FOD as well as
# FPA FOD and how they are paired to FINN.
#
# Loads data created by:
#   R/assign_FINNFires_to_FPAFOD.R
#   R/merge_paired_FINN_FPAFOD_fires.R

# TODO: arguments for loading data with difference time and distance tolerances
# TODO: for pairing 

lightningCol <- adjustcolor("cyan4", alpha.f = 1)
humanCol     <- adjustcolor("orange", alpha.f = 1)
unknownCol   <- adjustcolor("gray", alpha.f = 1)

library(ggplot2)
library(ggmap)
library(maps)
library(mapdata)
library(ggthemes)
library(dplyr)

dataDir   <- "Data/FINN/" 
figureDir <- "Figures/FINN_assignments/"

# Load the merged yearly paired fire data created by R/merge_paired_fires.R
print("Loading data file that exists")
load(paste0(dataDir, "FINN_with_FPAFOD_dxdy=10_DT=1_7_2002_2015.RData"))
load(paste0(dataDir, "FPA_FOD_with_FINN_dxdy=10_DT=1_7_2002_2015.RData"))

# # I do not want Alaska, HI, or PR, for plotting purposes. 
# fire_states <- FPA_FOD$STATE
# stateMask <- fire_states != "HI" & fire_states != "AK" & fire_states != "PR"
# FPA_FOD <- FPA_FOD[stateMask, ]

fire_cause <- rep("", dim(FPA_FOD)[1])
fire_cause[FPA_FOD$STAT_CAUSE_DESCR != "Lightning" & FPA_FOD$STAT_CAUSE_DESCR != "Missing/Undefined"] <- "Human"
fire_cause[FPA_FOD$STAT_CAUSE_DESCR == "Missing/Undefined"] <- "Unknown"
fire_cause[FPA_FOD$STAT_CAUSE_DESCR == "Lightning"] <- "Lightning"
# Make this a factor 
fire_cause <- factor(fire_cause, levels=c("Lightning", "Human", "Unknown"))

FPA_FOD$fire_cause <- fire_cause

# Lets get some basics out of the way. Wildfires per year, proportion associated
# with FINN analysis
# TODO: Make monthly
years <- unique(FPA_FOD$FIRE_YEAR)
nYears <- length(years)
n_fires     <- rep(NA, nYears)
n_FINN_fires <- rep(NA, nYears)
for (i in 1:nYears){
  yearMask <- FPA_FOD$FIRE_YEAR == years[i]
  n_fires[i] <- sum(yearMask)
  n_FINN_fires[i] <- sum(yearMask & (FPA_FOD$n_FINN > 0) ) 
}

percentPaired <- round(sum(FPA_FOD$n_FINN>0)/length(FPA_FOD$n_FINN)*100,2)

n_Finn_west <- FPA_FOD$n_FINN[FPA_FOD$region == "West"]
percentPaired_west <- round(sum(n_Finn_west>0)/length(n_Finn_west)*100,2)

n_Finn_se <- FPA_FOD$n_FINN[FPA_FOD$region == "Southeast"]
percentPaired_se <- round(sum(n_Finn_se>0)/length(n_Finn_se)*100,2)


# Save the figure that shows if there is a trend in fires associated with AQ
# forecasts. Based on the way the data are generated it is not clear what this
# tells us. 
png(filename=paste0(figureDir, "percent_AQ_wildfires.png"), 
    res=250, height=1000, 
    width=2000)

plot(years, n_FINN_fires/n_fires*100, 
     xlab="year", ylab="% wildfires paired with FINN", las=1)
lines(years, n_FINN_fires/n_fires*100)
title("Annual percent of FPA FOD wildfires paired with FINN",
      sub=paste(percentPaired, "% of wildfires paired with FINN"))

dev.off()


################################################################################
# Show the number of FPA wildfires that are associated with FINN fires in
# in each region. 
################################################################################

# Mask out non west or southeast region
m <- (FPA_FOD$region != "") & (FPA_FOD$n_FINN > 0)
df <- FPA_FOD[m,]

sum((df$region=="Southeast") & (df$fire_cause=="Human"))
sum((df$region=="West") & (df$fire_cause=="Human"))


g <- ggplot(df, aes(x=region, fill=fire_cause))+
  geom_bar(stat="count", width=0.5, position = "dodge")+
  scale_fill_manual(values = c("Lightning"=lightningCol, "Human"=humanCol, "Unknown"=unknownCol))+
  guides(fill=guide_legend(title="Fire Cause"))+
  theme_tufte(ticks=T, base_size = 20)+
  ylab("FPA FOD fires linked\nto FINN fires")+
  ggtitle("FPA FOD wildires co-located\nwith FINN fires")

png(filename=paste0(figureDir, "FPA_FOD_with_FINN_count_by_region.png"), 
    res=250, height=1000, 
    width=1500)
print(g)
dev.off()

################################################################################
# Show the total PM2.5 emissions estimated by FINN for these paired wildfires
# NOTE: Be careful when using dodge
# https://stackoverflow.com/questions/46279720/using-dodge-position-in-ggplot-changing-column-values
################################################################################
west_lightning_PM25 <- sum(df$TOTAL_PM25[ (df$region=="West") & (df$fire_cause=="Lightning")])
east_human_PM25 <- sum(df$TOTAL_PM25[ (df$region=="Southeast") & (df$fire_cause=="Human")])

# Create a new category, sum PM25 for a given region fire_cause combo
p <- ggplot(df %>%
            group_by(region, fire_cause) %>%
            summarise(TOTAL_PM25_SUMMED = sum(TOTAL_PM25)),
            aes(x = region, y = TOTAL_PM25_SUMMED, fill = fire_cause))

# Now set up the dodged bars by color and other styling needed for plot 
p2 <- p + geom_bar(stat = "identity", position = "dodge")+ 
  scale_fill_manual(values = c("Lightning"=lightningCol, 
                               "Human"=humanCol, 
                               "Unknown"=unknownCol))+
  guides(fill=guide_legend(title="Fire Cause"))+
  theme_tufte(ticks=T, base_size = 20)+
  ylab(expression("PM"[2.5]*" [kg]"))+
  ggtitle(expression("PM"[2.5]*" emissions from FPA FOD wildires"))

# Save the image 
png(filename=paste0(figureDir, "Total_PM25_FPA_FOD_paired_to_FINN.png"), 
    res=250, height=1300, width=2000)
print(p2)
dev.off()

################################################################################
# Show the number of FINN fires that are associated with FPA FOD wildfires
# in each region. 
# TODO: This needs an ignition attribution and then the y-axis needs to be 
# TODO: changed to total PM25.
# TODO: Use the summary teqnique shown above! 
################################################################################

# Mask out non west or southeast region, also wildfires with no hysplits
m <- (FINN$region != "") 
df <- FINN[m,]

# Make an association factor 
wild <- rep("", sum(m))
wild[df$pairedWithFPAFOD] <- "Yes"
wild[!df$pairedWithFPAFOD] <- "No"
df$wild <- wild

sum(df$region=="West" & df$wild=="No")

p <- ggplot(df, aes(x=region, fill=wild))+
  geom_bar(stat="count", width=0.5, position = "dodge")+
  #scale_fill_manual(values = c("Lightning"=lightningCol, "Human"=humanCol, "Unknown"=unknownCol))+
  guides(fill=guide_legend(title="Accounted for in\nFPA FOD"))+
  theme_tufte(ticks=T, base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5))+ # Center the title
  ylab("FINN detection count")+
  ggtitle("FINN detections")

png(filename=paste0(figureDir, "FINN_wFPA_FOD_count_by_region.png"), 
    res=250, height=1300, width=2000)
print(p)
dev.off()

################################################################################
# Map when and where the paired fire data occur in the US. 
# http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
################################################################################
FINN_fires <- FPA_FOD[FPA_FOD$n_FINN > 0, ] 

states <- map_data("state")
P <- ggplot() + geom_polygon(data = states, aes(x=long, y = lat, group = group),
                        fill = NA, color = lightningCol) + 
  coord_fixed(1.3)+ # fixes the relationship between one unit in the y direction and one unit in the x direction.
  theme_tufte(ticks=T, base_size = 20)+
  geom_point(data = FINN_fires, aes(x = LONGITUDE, y = LATITUDE, 
                                  color = fire_cause, size = FIRE_SIZE)
             #shape=fire_cause)
  )+
  scale_colour_manual(values = c("Lightning"=lightningCol, "Human"=humanCol, "Unknown"=unknownCol))


png(filename=paste0(figureDir, "FPA_FOD_paired_wFINN_fires_mapped.png"), res=100, 
    width=2000, height=1200)
print(P)
dev.off()  

################################################################################
# Plot the locations of FINN fires that have been paired with FPA FOD and 
# those that have not. 
################################################################################

world <- map_data("world")


h <- ggplot() + geom_polygon(data = states, aes(x=long, y = lat, group = group),
                        fill = NA, color = lightningCol) + 
  coord_fixed(1.3)+ 
  theme_tufte(ticks=T, base_size = 20)+
  geom_point(data = FINN, aes(x = LONGI, y = LATI, 
                                       color = pairedWithFPAFOD), 
             shape=19, size=0.5)+
  xlim(c(-125,-60))+
  ylim(c(27,49))+
  ggtitle("FINN fires shaded by pairing with FPA FOD")


png(filename=paste0(figureDir, "paired_FINNfires_mapped.png"), res=100, 
    width=2000, height=1200)
print(h)
dev.off()  

# ################################################################################
# # Plot the paired datasets on the same map. Color by dataset. Dots off alone
# # in the wild represent errors as there should be no more than 10 km of seperation.
# ################################################################################
# 
# HP_paired <- hysplitPoints[hysplitPoints$nFPAFODPaired > 0, ]
# 
# png(filename=paste0(figureDir, "paired_HYSPLITPoints_mapped.png"), 
#     width=2000, height=1200)
# 
# ggplot() + geom_polygon(data = states, aes(x=long, y = lat, group = group),
#                         fill = NA, color = lightningCol) + 
#   coord_fixed(1.3)+
#   geom_point(data = HP_paired, aes(x = Lon, y = Lat), shape=19, size=0.5, col=unknownCol)+
#   geom_point(data = AQ_fires, aes(x = LONGITUDE, y = LATITUDE), shape=19, size=0.5, col="red")+
#   theme_tufte(ticks=T, base_size = 20)+
#   ylim(c(25,49.5))+
#   xlim(c(-125,-60))
# 
# 
# dev.off()
# 
# ################################################################################
# # HMS smoke duration hours as a function of 
# ################################################################################
# png(filename=paste0(figureDir, "fire_size_vs.png"), width=1700, height=1400,
#     res=200)
# 
# p <- ggplot(AQ_fires, aes(x=FIRE_SIZE, y=n_HP, col=fire_cause, size=totalDurationHours))
# p + geom_point() +  
#   scale_colour_manual(values = c("Lightning"=lightningCol, "Human"=humanCol, "Unknown"=unknownCol))+
#   theme_tufte(ticks=T, base_size = 20)+
#   xlab("Fire Size (acres)") + ylab("Hysplit points associated with fire")
# 
# dev.off()
# 
# 
# ################################################################################
# # Plot bar chart of total AQ wildfires and the proportion by ignition type in 
# # the different ecoregions. 
# ################################################################################
# df <- FPA_FOD[FPA_FOD$n_HP==0, ] # FPA_FOD[FPA_FOD$n_HP==0, ] AQ_fires
# df$NA_L2CODE <- as.factor(df$NA_L2CODE)
# 
# # Make associatedion with air quality forecast a factor in the data 
# fireType <- rep("AQ forecast", dim(df)[1])
# fireType[df$n_HP == 0] <- "No AQ forecast"
# fireType <- factor(fireType, levels = c("No AQ forecast", "AQ forecast"))
# df$fireType <- fireType
# 
# png(filename=paste0(figureDir, "regional_counts.png"), width=2000, height=800,
#     res=200)
# 
# g <- ggplot(df, aes(NA_L2CODE,fill=fire_cause))
# g + geom_bar()+
#   scale_fill_manual(values = c("Lightning"=lightningCol, "Human"=humanCol, "Unknown"=unknownCol))+
#   theme_tufte(ticks=T, base_size = 15)+
#   #facet_wrap(~fireType)+
#   #ylab("Wildfires associated with AQ forecasts")+
#   xlab("Level II ecoregion")
# 
# dev.off()


