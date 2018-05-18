# plot_annual_burn_area.R

# ------------------------- Description ---------------------------------------
# This script is used to plot the burn area total, color coded by wildfire
# ignition type, for 1992 - 2015. This is a very simple plot that demonstrates
# interannual variability. 
# ------------------------------------------------------------------------------


library(sfsmisc)

load("Data/FPA_FOD/FPA_FOD_1992_2015.RData")

df <- FPA_FOD[!(FPA_FOD$STATE %in% c("AK", "HI", "PR") ),]

years <- 1992:2015
nYears <- length(years)

h_BA <- rep(NA, nYears)
l_BA <- rep(NA, nYears)

# NOTE: This is simplified for this schematic, since we are not color coding 
# NOTE: wildfires in unknown causes. 
humanMask <- df$STAT_CAUSE_DESCR != "Lightning"

for (i in 1:nYears){
  
  yearMask <- years[i] == df$FIRE_YEAR
  
  h_BA[i] <- sum(df$FIRE_SIZE[yearMask & humanMask])
  l_BA[i] <- sum(df$FIRE_SIZE[yearMask & !humanMask])
  
}

png(filename="Figures/annual_BA_bar.png", height=1300, width=2000, res=200)
par(mar=c(4,8,4,4))
barplot(height=(h_BA+l_BA), names.arg = years, yaxt='n', col="cyan4")
eaxis(2)
dev.off()

png(filename="Figures/annual_BA_bar_ignitions.png", height=1300, width=2000, res=200)
par(mar=c(4,8,4,4))
barplot(height=(h_BA+l_BA), names.arg = years, yaxt='n', col="cyan4")
barplot(height=h_BA, names.arg = "", yaxt='n', xaxt='n', col="orange", add=T)
eaxis(2)
dev.off()

# Make plot for the SE
load("Data/GIS/southeast_bounds.RData")
m <- (df$LONGITUDE >= minLon) & (df$LONGITUDE<=maxLon) & (df$LATITUDE >= minLat) & (df$LATITUDE <= maxLat)
df_SE <- df[m,]

humanMask <- df_SE$STAT_CAUSE_DESCR != "Lightning"

h_BA_SE <- rep(NA, nYears)
l_BA_SE <- rep(NA, nYears)
for (i in 1:nYears){
  
  yearMask <- years[i] == df_SE$FIRE_YEAR
  
  h_BA_SE[i] <- sum(df_SE$FIRE_SIZE[yearMask & humanMask])
  l_BA_SE[i] <- sum(df_SE$FIRE_SIZE[yearMask & !humanMask])
  
}

png(filename="Figures/annual_BA_bar_ignitions_southeast.png", height=1300, width=2000, res=200)
par(mar=c(4,8,4,4))
barplot(height=(h_BA_SE+l_BA_SE), names.arg = years, yaxt='n', col="cyan4")
barplot(height=h_BA_SE, names.arg = "", yaxt='n', xaxt='n', col="orange", add=T)
eaxis(2)
dev.off()

# Make plot for the West
load("Data/GIS/west_bounds.RData")
m <- (df$LONGITUDE >= minLon) & (df$LONGITUDE<=maxLon) & (df$LATITUDE >= minLat) & (df$LATITUDE <= maxLat)
df_W <- df[m,]

humanMask <- df_W$STAT_CAUSE_DESCR != "Lightning"

h_BA_W <- rep(NA, nYears)
l_BA_W <- rep(NA, nYears)
for (i in 1:nYears){
  
  yearMask <- years[i] == df_W$FIRE_YEAR
  
  h_BA_W[i] <- sum(df_W$FIRE_SIZE[yearMask & humanMask])
  l_BA_W[i] <- sum(df_W$FIRE_SIZE[yearMask & !humanMask])
  
}

png(filename="Figures/annual_BA_bar_ignitions_west.png", height=1300, width=2000, res=200)
par(mar=c(4,8,4,4))
barplot(height=(h_BA_W+l_BA_W), names.arg = years, yaxt='n', col="cyan4")
barplot(height=h_BA_W, names.arg = "", yaxt='n', xaxt='n', col="orange", add=T)
eaxis(2)
dev.off()



