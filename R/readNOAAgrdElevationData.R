# readNOAAElevationData.R

#------------------------------- Description ---------------------------------_
# This script is used to convert some 1-arc minute NOAA elevation data to an 
# nc format that will be easy to use in the existing data pipeline. These elevation
# data were donwloaded to replace a much courase 0.25 degree dataset provided
# by JISAO. 

# Data Description: https://www.ngdc.noaa.gov/mgg/global/global.html
# ETOP1 Bedrock Data URL: https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/ice_surface/grid_registered/netcdf/

library(raster)
f <- "Data/GIS/ETOPO1_Bed_g_gmt4.grd"
r <- raster(f)

# NOTE:
# 60 arc-minute (1 degree) and these data are about 1 arc-minute

print(r)
quartz()
plot(r)
data <- r@data

# NOTE: The code at this URL shows how to save a raster object as a . nc file
# http://geog.uoregon.edu/GeogR/topics/netcdf-to-raster.html

# write the raster layer as a netCDF file, so that it is used consistently as
# other data, as a regular grid. 
outfile <- "Data/GIS/ETOPO1.nc"
crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
writeRaster(r, outfile, overwrite=TRUE, format="CDF", 
            varname="elevation", varunit="meters", 
            longname="ETOPO1 Global Relief Model https://www.ngdc.noaa.gov/mgg/global/global.html", 
            xname="longitude", 
            yname="latitude")




