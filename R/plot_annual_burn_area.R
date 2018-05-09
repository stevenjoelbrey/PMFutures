# plot_annual_burn_area.R

library(sfsmisc)

load("Data/FPA_FOD/FPA_FOD_1992_2015.RData")

df <- FPA_FOD[!(FPA_FOD$STATE %in% c("AK", "HI", "PR") ),]

years <- 1992:2015
nYears <- length(years)

h_BA <- rep(NA, nYears)
l_BA <- rep(NA, nYears)

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
