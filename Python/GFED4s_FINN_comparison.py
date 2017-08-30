#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to show the correlations between GFED and FINN. It
# serves as both a sanity check and a nice little analysis

# NOTE: 1 will reference GFED and 2 will reference FINN, thoughout


import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import datetime
import matplotlib.ticker as tkr
import cesm_nc_manager as cnm
import scipy.stats as stats
import statsmodels
import matplotlib.colors as colors


figureDir = '../Figures/GFED4s_vs_FINN/'
region    = "_"

# Load both emission inventories, get time dimensions of each, make sure all 
# dimensions agree for comparisons moving forward.
GFEDFile = '/barnes-scratch/sbrey/GFED4s/GFED4.1s_ecmwf_C_2003_2013.nc'
GFED_nc = Dataset(GFEDFile, 'r')

FINNFile = '/barnes-scratch/sbrey/FINN/FINN_ecmwf_CO2_2003_2013.nc'
FINN_nc  = Dataset(FINNFile, 'r')

gramsCarbonPerCO2 = 12./(12.+32.)

###############################################################################
# Make sure the time and space dimensions match exactly 
###############################################################################
t1 = GFED_nc.variables['time']
t2 = FINN_nc.variables['time']

# Time
if np.unique(t1[:]-t2[:])[0] == 0:
	time, month, year = cnm.get_era_interim_time(t1)
	nTime = len(time)
	yday = np.zeros(nTime)
	for i in range(nTime):
		yday[i] = time[i].timetuple().tm_yday
		

# y and x 
lat1 = GFED_nc.variables['latitude'][:]
lat2 = FINN_nc.variables['latitude'][:]

if np.unique(lat1 - lat2)[0] == 0:
	lat = lat1
	
lon1 = GFED_nc.variables['longitude'][:]
lon2 = FINN_nc.variables['longitude'][:]

if np.unique(lon1 - lon2)[0] == 0:
	lon = lon1


###############################################################################
# Handle making the spatial subset
###############################################################################
x1 = GFED_nc.variables['C'][:]
x2 = FINN_nc.variables['CO2'][:] * gramsCarbonPerCO2

if region  != '_':
	
	# Get the chosen region bounds
	minLat, maxLat, minLon, maxLon, resolution = cnm.getRegionBounds(region)
	
	# Subset the data based on the bounds of this region. 
	x1, latNew, lonNew = cnm.mask2dims(x1, lon, lat, 0,\
							     xmin=minLon, xmax=maxLon,\
							     ymin=minLat, ymax=maxLat)
							     
	x2, latNew, lonNew = cnm.mask2dims(x2, lon, lat, 0,\
							           xmin=minLon, xmax=maxLon,\
							           ymin=minLat, ymax=maxLat)
	
	# Need to use the subset dimensions moving forward						           
	lon = lonNew
	lat = latNew						           

###############################################################################
# Begin comparisons
###############################################################################

# x1_ = x1.ravel()
# x2_ = x2.ravel()

# Global total daily correlation plot 
x1_daily = np.sum(x1, axis=(1,2))
x2_daily = np.sum(x2, axis=(1,2))

# Show a basic time series of these data on the same axis
fig = plt.figure(figsize=(12, 5))
plt.plot(time, x1_daily, label="GFED4s")
plt.plot(time, x2_daily, label="FINN")
plt.legend()
plt.title('Daily domain sum time series 2003-2013')
plt.savefig(figureDir + "daily_timeSeries" + region + ".png")

fig = plt.figure(figsize=(12, 9))

plt.plot(x1_daily, x1_daily, color='k') # identity line 
c = plt.scatter(x1_daily, x2_daily, s=4, c=yday)
plt.xlabel('GFED4s [grams carbon]', fontsize=18)
plt.ylabel('FINN [grams carbon]', fontsize=18)
cbar = plt.colorbar(c, label='Day of year')
plt.title('GFED4s vs. FINN daily'+ region[1:-1] + 'emissions 2003-2013', fontsize=18)

plt.savefig(figureDir + "daily" + region + "summed_emissions.png")
plt.close()



if region != "_":

	fig = plt.figure(figsize=(12, 9))

	plt.scatter(x1.ravel(), x2.ravel(), s=4)
	plt.xlabel('GFED4s [grams carbon]', fontsize=18)
	plt.ylabel('FINN [grams carbon]', fontsize=18)
	plt.title('GFED4s vs. FINN daily'+ region[1:-1] + 'emissions 2003-2013', fontsize=18)

	plt.savefig(figureDir + "grid_emission_correlation" + region + ".png")
	plt.close()

###############################################################################
# Show maps of totals (and normalized), Global and US each product
###############################################################################

x1_total = np.sum(x1, axis=0)
x2_total = np.sum(x2, axis=0)

# Mask where zero
x1_total = ma.masked_where(x1_total == 0, x1_total)
x2_total = ma.masked_where(x2_total == 0, x2_total)

vMin = np.array( (x1_total.min(), x2_total.min()) ).min()
vMax = np.array( (x1_total.max(), x2_total.max()) ).max()

# Set up the map
m = Basemap(projection='merc',llcrnrlat=minLat, urcrnrlat=maxLat,\
            llcrnrlon=minLon, urcrnrlon=maxLon, resolution=resolution,\
	    	lon_0=0, lat_0=-90)

# grid coords for mesh plotting of values. 
lons, lats = np.meshgrid(lon, lat)
x, y = m(lons, lats)

###############################################################################
# Compare total emissions, log scale colorbar
# Shows us which one has higher emissions, how this compares by location
###############################################################################
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

plt.subplot(121)	
#axes[0].set_title("first")

m.drawcoastlines()
m.drawstates()
m.drawcountries()
#c1 = m.pcolor(x, y,  x1_total, norm=colors.LogNorm(vmin=vMin, vmax=vMax))
c1 = m.pcolor(x, y,  x1_total, vmin=vMin, vmax=vMax)
cbar = plt.colorbar(c1, pad=0.01, orientation='horizontal', extend='both',\
			        norm=colors.LogNorm())
cbar.set_label('grams carbon', fontsize=25)


plt.subplot(122)	

m.drawcoastlines()
m.drawstates()
m.drawcountries()
#c2 = m.pcolor(x, y,  x2_total, norm=colors.LogNorm(vmin=vMin, vmax=vMax))
c2 = m.pcolor(x, y,  x2_total, vmin=vMin, vmax=vMax)
cbar = plt.colorbar(c2, pad=0.01, orientation='horizontal', extend='both')
cbar.set_label('grams carbon', fontsize=25)

fig.tight_layout()
plt.savefig(figureDir + "total_emissions_map" +region[0:-1] +'.png')
plt.close()

###############################################################################
# Plot normalized comparison
# This shows us where each inventory has the most emissions, highlights the
# inventory heavy locations 
###############################################################################
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12,6))

plt.subplot(121)	
#axes[0].set_title("first")

m.drawcoastlines()
m.drawstates()
m.drawcountries()
c1 = m.pcolor(x, y,  x1_total / x1_total.max() )
cbar = plt.colorbar(c1, pad=0.01, orientation='horizontal', extend='both',\
			        norm=colors.LogNorm())
cbar.set_label('normalized carbon emissions', fontsize=18)
plt.title('GFED4s', fontsize=22)

plt.subplot(122)	

m.drawcoastlines()
m.drawstates()
m.drawcountries()
c2 = m.pcolor(x, y,  x2_total / x2_total.max())
cbar = plt.colorbar(c2, pad=0.01, orientation='horizontal', extend='both')
cbar.set_label('normalized carbon emissions', fontsize=18)
plt.title('FINN', fontsize=22)

fig.tight_layout()
plt.savefig(figureDir + "total_emissions_map_norm" + region[0:-1] +'.png')
plt.close()


# TODO: compare the seasonality 
x1_month_total = np.zeros(12)
x2_month_total = np.zeros(12)

for m in range(12):
	month_mask = (m+1) == month
	x1_month_total[m] = np.sum( x1_daily[month_mask] )
	x2_month_total[m] = np.sum( x2_daily[month_mask] )
	
months = np.arange(1,13)
width = 0.4	


fig = plt.figure()
ax = fig.add_subplot(111)
rects1 = ax.bar(months, x1_month_total, width, label="GFED4s")
rects2 = ax.bar(months + width, x2_month_total, width,label="FINN")
plt.legend()

plt.savefig(figureDir + "seasonal_totals" + region[0:-1] +'.png')



