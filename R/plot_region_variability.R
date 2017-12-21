# plot_region_varaibility.R

# This script will be used for comparing the interannual variability, and 
# relations to mean synoptic meteorology for human and lightning started fires. 

# Are human-started fires independent of inter-annual variability? 
################################################################################


library(sp)
library(maps)
library(ncdf4)
library(stringr)
library(sfsmisc)

# Make sketched figure 1. 

# What region are you investigating? 
minLat <- 31
maxLat <- 49
minLon <- -125
maxLon <- -100

includedMonths <- c(5:10)

# Load required datasets
# 1) Burn area from FPA and GFED4s (FINN)
# 2) GIS borders of ecoregions. 

load("Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData")

# Make ecoregion color map, I got these colors using colorzilla
includedRegions <- c(6.2, 10.1, 7.1, 9.3, 9.4, 13.1, 12.1, 11.1, 10.2)
regionColors <- c("#65B657", "#E7ED90", "#60BAAF", "#F2DBA1", "#EDCB9B",
                  "#B5D67C",  "#D6D292", "#D2E5B1", "#F3DB70")

SPDF$NA_L2CODE <- as.numeric(as.character(SPDF$NA_L2CODE))

regionMask <- SPDF$NA_L2CODE %in% includedRegions
SPDF_subset <- SPDF[regionMask, ]


ORDER <- match(includedRegions, SPDF_subset$NA_L2CODE)
SPDF_subset <- SPDF_subset[ORDER, ]


# load burn area grids, subset spatial domain right away. 
ncFile <- "Data/FPA_FOD/burn_area_monthly_25x25_2003_2013.nc"
nc <- nc_open(ncFile)
fpaLat <- ncvar_get(nc, "latitude")
fpaLon <- ncvar_get(nc, "longitude")
ba_fpa <- ncvar_get(nc, "burn_aera")
fpaTime<- ncvar_get(nc, "YYYYMM")

fpaYear <- as.numeric(str_sub(fpaTime, 1,4))
fpaMonth <- as.numeric(str_sub(fpaTime, 5,6))

nc_close(nc)

ncFile <- "Data/GFED4s/GFED4.1s_monthly_DM_2003_2016.nc"
nc <- nc_open(ncFile)
gfedLat <- ncvar_get(nc, "latitude")
gfedLon <- ncvar_get(nc, "longitude")
gfedTime <- ncvar_get(nc, "time")
grid_area <- ncvar_get(nc, "grid_area")
burn_area_fraction <- ncvar_get(nc, "burn_area_fraction")
nc_close(nc)

gfedYear <- as.numeric(str_sub(gfedTime, 1,4))
gfedMonth <- as.numeric(str_sub(gfedTime, 5,6))


# Need total burn area, not burn area fraction as is stored in nc files. 
# NOTE: We want the time dimension to come first, to match the FPA-FOD data. so
# NOTE: we will rewrite this array one time dimension at a time. 

ba_gfed <- array(NA, dim=c(168, 1440, 720)) # No NA in array to be transformed so check after
for (i in 1:dim(burn_area_fraction)[3]){
  ba_gfed[i,,] <- burn_area_fraction[,,i] * grid_area  # [1440  720] dot [1440  720]
}

# Load the needed grid attributes
grid_nc_file <- "Data/grid_attributes/grid_attributes_25x25.nc"
nc_grid   <- nc_open(grid_nc_file)
ecoregion <- round(ncvar_get(nc_grid, "ecoregion"), 2) # annoying precision issues
#elevation <- ncvar_get(nc_grid, "elevation")
#state     <- ncvar_get(nc_grid, "state")
nc_close(nc_grid)

# Spatialy subset the burn area data 
acrePerMsquared <- 4046.86

latMask <- fpaLat >= minLat & fpaLat <= maxLat
lonMask <- fpaLon >= minLon & fpaLon <= maxLon

# Subset and conert to acres
ba_gfed <- ba_gfed[, lonMask, latMask] / acrePerMsquared
ba_fpa <- ba_fpa[, lonMask, latMask] / acrePerMsquared
ecoregion<- ecoregion[lonMask, latMask]

lat <- fpaLat[latMask]
lon <- fpaLon[lonMask]


# Now build the summary dataframes! 
nRow <- length(2003:2016) * length(1:12)
nCol <- length(includedRegions)
m <- matrix(NA, nrow=nRow, ncol=nCol)
df <- data.frame(m)
names(df) <- as.character(includedRegions)
row.names(df) <- gfedTime

fpa_monthly_df <- df
gfed_monthly_df <- df

# Assigm burn area by month to the repsective ecoregion! 
for (i in 1:nRow){
  
  for (eco in includedRegions){
    
    # spatial mask
    ecoMask <- eco == ecoregion 
    
    # sum of this months spatial mask value
    # if statement for fpa, which does not pspan the same time 
    if (i <= length(fpaTime)){
      fpa_monthly_df[i, eco == includedRegions] <- sum(ba_fpa[i,,][ecoMask]) 
    }
    
    gfed_monthly_df[i, eco == includedRegions] <- sum(ba_gfed[i,,][ecoMask]) 
    
  }
  
}


# Quickly great the annual versions 
# Now build the summary dataframes! 
nRow <- length(2003:2016)
nCol <- length(includedRegions)
m <- matrix(NA, nrow=nRow, ncol=nCol)
df <- data.frame(m)
names(df) <- as.character(includedRegions)
row.names(df) <- as.character(2003:2016)

fpa_yearly_df  <- df
gfed_yearly_df <- df

years <- 2003:2016
for (y in 1:length(years) ){
  
  yearMask <- years[y] == gfedYear
  
  for (eco in includedRegions){
   
    eco <- as.character(eco)
    fpa_yearly_df[[eco]][y] <- sum(fpa_monthly_df[[eco]][yearMask], na.rm=T)
    gfed_yearly_df[[eco]][y] <- sum(gfed_monthly_df[[eco]][yearMask])
     
    
  }
  
}

################################################################################
# Make the figure, two panal with the map on the left.
################################################################################
png(filename="Figures/summary/figure_1.png",
    width=2200, 
    height=1000,
    res=200)

par(mfrow=c(1,2), mar=c(5,8,4,11))

map("state",  lty=3, col="darkgray", 
    xlim=c(minLon, maxLon), ylim=c(minLat, maxLat))
plot(SPDF_subset, col=regionColors, add=T)
map("state",  lty=3, col="black",  add=T)

legend("topleft",inset=c(-0.2, 0), xpd=T,
       legend=c(includedRegions),
       fill=regionColors,
       bty="n"
       )

# Weave the two annual datasets together
fpa <- t(as.matrix(fpa_yearly_df))
gfed <- t(as.matrix(gfed_yearly_df))

both <- matrix(NA, nrow=9, ncol=28)
j <- 0
for (i in 1:14){
  
  j <- j + 1
  both[,j] <- fpa[, i]
  j <- j + 1
  both[,j] <- gfed[, i]
  
}

# Set the spacing between bars such that the years are separated
space     <- rep(c(1,0), 14)

labels  <- as.character()
for (i in 1:length(years)){
  labels <- append(labels, c(years[i]," "))
}


bp <- barplot(both, col= regionColors, border=NA, yaxt='n', 
              space=space, width=1.5)

# Get the midpoints of bp
atValue <- bp[-length(bp)] + diff(bp)/2
# I want the 1st and 3rd, and so on
at <- atValue[seq(1,27, by=2)]

eaxis(2, cex.axis=1.0)
mtext("Burn Area\n(acres)  ", side=2, cex=1.4, las=2, line=3.5)
axis(1, at=at, labels=years, las=2, cex.axis=1.0)

dev.off()


################################################################################
# Now, a histogram of the monthly totals
################################################################################

# TODO: Make sure these draw from the same years!

# make the monthly totals by ecoregion
df  <- data.frame(matrix(NA, nrow=12, ncol=length(includedRegions)))
names(df) <- as.character(includedRegions)
gfed_seasonality <- df
fpa_seasonality  <- df


gfed_monthly_df <- gfed_monthly_df[gfedYear <= 2013,]
fpa_monthly_df <- fpa_monthly_df[gfedYear <= 2013,]
gfedMonth <- gfedMonth[gfedYear <= 2013]

for (m in 1:12){
  
  monthMask <- m == gfedMonth
  
  for (eco in includedRegions){
    
    ecoMask <- eco == includedRegions
    gfed_seasonality[m, ecoMask] <- sum(gfed_monthly_df[monthMask, ecoMask])
    fpa_seasonality[m, ecoMask] <- sum(fpa_monthly_df[monthMask, ecoMask], na.rm=T)
    
  }
}

png(filename="Figures/summary/western_US_seasonality.png", res=200,
    width=2000, height=1000)
par(mfrow=c(1,2), mar=c(4,5,4,4))

# Make a common y-axis
ylim <- c(0,1.4*10^7)

bp <- barplot(t(as.matrix(gfed_seasonality)), 
              col=regionColors, border=NA, yaxt='n',
              ylim=ylim)
axis(1, at=bp, labels=1:12)
mtext("Month", side=1, line=2)
title("GFED4s")
eaxis(2)
#mtext("Burn Area\n(acres)  ", side=2, cex=1.4, las=2, line=3.5)


bp <- barplot(t(as.matrix(fpa_seasonality)), 
              col=regionColors, border=NA, yaxt='n',
              ylim=ylim)
axis(1, at=bp, labels=1:12)
mtext("Month", side=1, line=2)
title("FPA-FOD")
#eaxis(2)

dev.off()

