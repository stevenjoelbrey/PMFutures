#!/usr/bin/env python2


###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to make sure that regridding of the GFED4s emissions
# worked. 


import os
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import matplotlib.ticker as tkr

def mask2dims(data, x, y, tDim,xmin, xmax, ymin, ymax):
	"""Takes a data array of shape [time, lat, lon] and subsets
           the data to the domain specified by x and y min and max arguments
           tDim specifies the location of the time dimension. Use tDim=-1
	   when the data array has no time dimension.

	   return: subset, ynew, xnew
	"""

	xMask = np.where((x >= xmin) & (x <= xmax))
	xnew = x[xMask]
	xmin_i = np.min(xMask)
	xmax_i = np.max(xMask) + 1

	yMask = np.where((y >= ymin) & (y <= ymax))
	ynew = y[yMask]
	ymin_i = np.min(yMask)
	ymax_i = np.max(yMask) + 1

	if (tDim==0):
		subset = data[:, ymin_i:ymax_i, xmin_i:xmax_i]
	elif(tDim==2):
		subset = data[ymin_i:ymax_i, xmin_i:xmax_i, :]
	elif(tDim==-1):
		# No time dimension
		subset = data[ymin_i:ymax_i, xmin_i:xmax_i]

	return subset, ynew, xnew

################################################################################
# Select and load particulate species emissions and AirQualityData Masks
################################################################################
dataDirBase  = '/barnes-scratch/sbrey/'
emissionsDir = 'GFED4s/'
species      = 'C'
year         = '2016'

emissionFileNew = dataDirBase + emissionsDir + 'GFED4.1s_METGrid_' + species + '_' + year + '.nc'

nc = Dataset(emissionFileNew, 'r')
E_new   = nc.variables[species]
lon_new = nc.variables['longitude']
lon_new = lon_new[:]-180.
lat_new = nc.variables['latitude'][:]
deltaMass=nc.variables['deltaMass']


# load the original data for a comparison of total emissions
emissionFileOld = dataDirBase + emissionsDir + 'GFED4.1s_' + species + '_' + year + '.nc'
nc      = Dataset(emissionFileOld, 'r')
E_old   = nc.variables[species]
lon_old = nc.variables['longitude'][:]
lat_old = nc.variables['latitude'][:]

# Calculate deltaMass as a percent of the mass of the original grid
deltaMassPercent = deltaMass[:] / np.sum(E_old, axis=(1,2))*100.


################################################################################
#------------------ Subset emissions in space (and time?) ------------------
################################################################################
# For report use western U.S. only

ymin   = 10    #30.    # 10.
ymax   = 90.   #50.    # 90.
xmin   = -130.  #234.   # 190.
xmax   = -60.  #259.   # 320


E_new, lat_new, lon_new = mask2dims(E_new, lon_new, lat_new, 0, xmin, xmax, ymin, ymax)
E_old, lat_old, lon_old = mask2dims(E_old, lon_old, lat_old, 0, xmin, xmax, ymin, ymax)

# Sum over the daily axis to make emissions annual
E_new_all = np.sum(E_new[:], axis=0)
E_old_all = np.sum(E_old[:], axis=0)

E_new_norm = E_new_all / E_new_all.max()
E_old_norm = E_old_all / E_old_all.max()


# Set up a relevant map to use later
#m = Basemap(projection='cyl')

m = Basemap(width=12000000,height=8000000,
            resolution='l',projection='stere',\
            lat_ts=50,lat_0=50,lon_0=-107.)


# Set up projection on map for regridded (new) data 
lons, lats = np.meshgrid(lon_new, lat_new)
x_new, y_new = m(lons, lats)

lons, lats = np.meshgrid(lon_old, lat_old)
x_old, y_old = m(lons, lats)




################################################################################
#----------------------- Load desired met event masks --------------------------
################################################################################

#highWindMask = cnm.getSelf(AScenario, 'highWindMask_9')
#highTMask    = cnm.getSelf(AScenario, 'highTMask_298')
#lowPrecMask  = cnm.getSelf(AScenario, 'lowPRECTMask_0.01')
#stagMask     = cnm.getSelf('2000Base', 'stagnation_mask') # needs to be rerun for 2000Base

#ncFile = cnm.makeAQNCFile('stagnationMask',AScenario, 'daily')
#nc     = Dataset(ncFile, 'r')
#stagMask = nc.variables['stagnationMask'][:]
#AQLat  = nc.variables['lat'][:]
#AQLon  = nc.variables['lon'][:]
#nc.close()


# cyclone & blocking days coming soon (I hope)


################################################################################
# Map the occurance of each type of met event 
################################################################################

fig = plt.figure(figsize=(12,8))

m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c    = m.pcolor(x_old, y_old, E_old_all, cmap="jet")
cbar = m.colorbar(c, location='bottom', pad="1%", extend='both')
cbar.set_label('Grams C day$^{-1}$')
plt.title(year +  ' Total emissions (original grid)')
plt.draw()
plt.show(block=False)
#plt.savefig('../Figures/GFED_original_' + year + '.png')
#plt.close()


fig = plt.figure(figsize=(12,8))

m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x_new, y_new, E_new_all, cmap="jet")
cbar = m.colorbar(c, location='bottom',pad="1%")
cbar.set_label('Grams C day$^{-1}$')
plt.title(year + ' Total emissions (re-gridded)')
plt.draw()
plt.show(block=False)

#plt.savefig('../Figures/GFED_regridded_' + year + '.png')
#plt.close()

fig = plt.figure()
plt.plot(deltaMassPercent)
plt.title('Daily change in mass from regriddiding %')
plt.draw()
plt.show(block=False)

#plt.savefig('../Figures/deltaMass_' + year + '.png')
#plt.close()











