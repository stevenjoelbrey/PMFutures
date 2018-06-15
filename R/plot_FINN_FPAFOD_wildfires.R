# plot_FINN_FPAFOD_wildfires.R
#
# ------------------------- Description ---------------------------------------
# This script plots FPA FOD PM25 emissions as estimated by the FINN modeling 
# franework. FPA FOD fires were used instead of MODIS dections and PM25 
# emissions were estimated by Christine. 
#
# Loads data created by:
#   R/process_FPA_FOD_for_FINN.R
#   R/sanity_check_FINN_procssed_data.R

# ------------------------- Arguments ------------------------------------------
year1 <- 2002 
year2 <- 2015
dataInFile   <- "Data/FINN/processed/FPA_FOD_wFINN-PM25_2002-2015.RData" 
figureDir <- paste0("Figures/FINN_emissions/")

a <- 1 # alpha for color in plots 
lightningCol <- adjustcolor("cyan4", alpha.f = a)
humanCol     <- adjustcolor("orange", alpha.f = a)
unknownCol   <- adjustcolor("gray", alpha.f = a)
noPairCol    <- adjustcolor("black", alpha.f = a)

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
fire_cause[FPA_FOD$STAT_CAUSE_DESCR == "Missing/Undefined"] <- "Unknown"
fire_cause[FPA_FOD$STAT_CAUSE_DESCR == "Lightning"] <- "Lightning"
# Make this a factor and assign it
fire_cause <- factor(fire_cause, levels=c("Lightning", "Human", "Unknown"))
FPA_FOD$fire_cause <- fire_cause


# Mask out non west or southeast region
m  <- (FPA_FOD$region != "") & (FPA_FOD$n_FINN > 0)
df <- FPA_FOD[m,]


################################################################################
# Show the number of FPA wildfires that are associated with FINN fires in
# in each region. 
################################################################################

# These sums are created and shown so that I can be sure I know what ggplot
# code is doing/showing
print(paste("Southeast human total:", sum( (df$region=="Southeast") & (df$fire_cause=="Human") ) ))
print(paste("Southeast lightning total:", sum( (df$region=="Southeast") & (df$fire_cause=="Lightning") ) ))

print(paste( "West human total:" ,sum((df$region=="West") & (df$fire_cause=="Human")) ))
print(paste( "West lightning total:" ,sum((df$region=="West") & (df$fire_cause=="Lightning")) ))

g <- ggplot(df, aes(x=region, fill=fire_cause))+
  geom_bar(stat="count", width=0.5, position = "dodge")+
  scale_fill_manual(values = c("Lightning"=lightningCol, "Human"=humanCol, "Unknown"=unknownCol))+
  guides(fill=guide_legend(title="Fire Cause"))+
  theme_tufte(ticks=T, base_size = 16)+
  ylab("FPA FOD fires linked\nto FINN fires")+
  ggtitle("FPA FOD wildires co-located\nwith FINN fires")

png(filename=paste0(figureDir, "FPA_FOD_with_FINN_count_by_region.png"), 
    res=250, height=1000, 
    width=1300)
print(g)
dev.off()

################################################################################
# Show the total PM2.5 emissions estimated by FINN for these paired wildfires
# NOTE: Be careful when using dodge. It can do strange things when summing a 
# NOTE: numeric variable for a specific subset of the data 
# https://stackoverflow.com/questions/46279720/using-dodge-position-in-ggplot-changing-column-values
################################################################################
# west_lightning_PM25 <- sum(df$TOTAL_PM25[ (df$region=="West") & (df$fire_cause=="Lightning")])
# east_human_PM25    <- sum(df$TOTAL_PM25[ (df$region=="Southeast") & (df$fire_cause=="Human")])
# 
# # Create a new category, sum PM25 for a given region fire_cause combo
# p <- ggplot(df %>%
#             group_by(region, fire_cause) %>%
#             summarise(TOTAL_PM25_SUMMED = sum(TOTAL_PM25)), # i.e. sums TOTAL_PM25 by region and fire cause
#             aes(x = region, y = TOTAL_PM25_SUMMED, fill = fire_cause))
# 
# # Now set up the dodged bars by color and other styling needed for plot
# p2 <- p + geom_bar(stat = "identity", position = "dodge")+
#   scale_fill_manual(values = c("Lightning"=lightningCol,
#                                "Human"=humanCol,
#                                "Unknown"=unknownCol))+
#   guides(fill=guide_legend(title="Fire Cause"))+
#   theme_tufte(ticks=T, base_size = 20)+
#   ylab(expression("PM"[2.5]*" [kg]"))+
#   ggtitle(expression("PM"[2.5]*" emissions from FPA FOD wildires"))
# 
# # Save the image
# png(filename=paste0(figureDir, "Total_PM25_FPA_FOD_paired_to_FINN.png"),
#     res=250, height=1300, width=2000)
# print(p2)
# dev.off()
# 
# rm(df)
################################################################################
# Show the number of FINN fires that are associated with FPA FOD wildfires
# in each region. 
################################################################################

# Mask out non west or southeast region areas
m <- (FINN$region != "") 
df <- FINN[m,]

# Create a datafrane if the number of wildfires for each ignition type associated
# with a given FINN (row) detection
noFPAPaired <- rep(0, length(df$nFPAFODPaired))
noFPAPaired[df$nFPAFODPaired==0] <- Inf # Want this to be the category when no FPA FOD

ignition <- data.frame(nFPALightningPaired=df$nFPALightningPaired, 
                       nFPAHumanPaired=df$nFPAHumanPaired, 
                       nMissingPaired=df$nMissingPaired,
                       noFPAPaired=noFPAPaired)
# Which column has the biggest number?
# If the row is all zeros than which.max returns the first index. This is why
# Inf assigned to noFPAPaired rows that are all zero.
mode_column <- apply(ignition, 1, which.max)

# Assign the wildfire responsible and make this a factor
ignitionType <- names(ignition)[mode_column]
ignitionType <- factor(ignitionType, 
                       levels=c("nFPALightningPaired", "nFPAHumanPaired", "nMissingPaired", "noFPAPaired"), 
                       labels=c("Lightning", "Human", "Unknown/Missing", "No FPA FOD paired"))

# Put this factor on the dataframe to be plotted
df$ignitionType <- ignitionType

# Sanity checks
sum(df$PM25[df$region=="West" & df$ignitionType=="Lightning"])
sum(df$PM25[ (df$region=="Southeast") & (df$ignitionType=="Human")])
sum(df$PM25[df$region=="West" & (df$nFPALightningPaired>0 | df$nFPAHumanPaired>0)])

# Create a new category, sum PM25 for a given region ignitionType combo
g_sum <- ggplot(df %>%
                group_by(region, ignitionType) %>% 
                summarise(PM25_SUMMED = sum(PM25)),
                aes(x = region, y = PM25_SUMMED, fill = ignitionType))

g <- g_sum + geom_bar(stat = "identity", position = "dodge")+ 
  scale_fill_manual(values = c("Lightning"=lightningCol, "Human"=humanCol, "Unknown/Missing"=unknownCol, "No FPA FOD paired"=noPairCol))+
  guides(fill=guide_legend(title="Ignition type"))+
  theme_tufte(ticks=T, base_size = 20)+
  theme(plot.title = element_text(hjust = 0.5))+ # Center the title
  ylab(expression("PM"[2.5]*" [kg]"))+
  ggtitle(expression("FINN PM"[2.5]*" emissions"))

png(filename=paste0(figureDir, "FINN_PM25_emissions_by_ignition.png"), 
    res=250, height=1300, width=1700)
print(g)
dev.off()

# ################################################################################
# # Map when and where the paired fire data occur in the US. 
# # http://eriqande.github.io/rep-res-web/lectures/making-maps-with-R.html
# ################################################################################
# FINN_fires <- FPA_FOD[FPA_FOD$n_FINN > 0, ] 
# 
# states <- map_data("state")
# P <- ggplot() + geom_polygon(data = states, aes(x=long, y = lat, group = group),
#                         fill = NA, color = lightningCol) + 
#   coord_fixed(1.3)+ # fixes the relationship between one unit in the y direction and one unit in the x direction.
#   theme_tufte(ticks=T, base_size = 20)+
#   geom_point(data = FINN_fires, aes(x = LONGITUDE, y = LATITUDE, 
#                                   color = fire_cause, size = FIRE_SIZE)
#              #shape=fire_cause)
#   )+
#   scale_colour_manual(values = c("Lightning"=lightningCol, "Human"=humanCol, "Unknown"=unknownCol))
# 
# 
# png(filename=paste0(figureDir, "FPA_FOD_paired_wFINN_fires_mapped.png"), res=100, 
#     width=2000, height=1200)
# print(P)
# dev.off()  

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


