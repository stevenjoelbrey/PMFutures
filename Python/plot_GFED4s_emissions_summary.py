#!/usr/bin/env python2

# plotEmissionSummary.py

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to summarize and plot emissions from GFED4s. 
# Summaries will be met event specific so we can learn about the relationships
# between fire emissions and synoptic meteorology.


# TODO: Make the distribution + boxplot that you drafted with Emily. 

# TODO: Make this into a python notebook for better integrated analysis.

# TODO: Different lags for testing and different averaging timescales. 

# NOTE: In this scrip E will generally stand for "emissions", m for mask,
# NOTE: and _ will represent a flattened multidimensional array. 

# Load resources
import os
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import numpy.ma as ma
import datetime
from datetime import date
from datetime import timedelta
import matplotlib.ticker as tkr
import cesm_nc_manager as cnm

################################################################################
# Select and load emissions species and AirQualityData Masks
################################################################################

dataDirBase = "/barnes-scratch/sbrey/era_interim_nc_daily_merged/"
figureDir = "../Figures/GFED_era_interm_analysis/"


# Get emissions, use this to get dimensions
ncFile  = "/barnes-scratch/sbrey/GFED4s/GFED4.1s_METGrid_C_NA_2003_2016.nc"
nc = Dataset(ncFile, 'r')
latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
time = nc.variables['time'][:]
C = nc.variables['C'][:]
nc.close()

# Make time into datetime arrays
time, month, year = cnm.get_era_interim_time(time)

################################################################################
#------------------ Subset model emissions in time -----------------------------
# -------- For spatial subsetting use cesm_nc_manager.mask2dims() --------------
################################################################################
# For report use western U.S. only. That is what you get when you load a file
# that has "_NA_" in the name. 

startMonth = 6
endMonth   = 9
# Get bounds for creating a map region 
minLat     = latitude.min()  # NOTE: This may be different from the bounds on the "_NA_" data
maxLat     = latitude.max()    
minLon     = longitude.min()   
maxLon     = longitude.max()   

# Make and apply month mask 
month_mask = (month >= startMonth) & (month <= endMonth)
C = C[month_mask,:,:]
time = time[month_mask]
month = month[month_mask]
year = year[month_mask]

################################################################################
# Set up a relevant map to use later
################################################################################
m = Basemap(projection='merc',llcrnrlat=minLat, urcrnrlat=maxLat,\
            llcrnrlon=minLon, urcrnrlon=maxLon,resolution='c',\
	    lon_0=0, lat_0=-90)

# grid coords for mesh plotting of values. 
lons, lats = np.meshgrid(longitude, latitude)
x, y = m(lons, lats)

################################################################################
# Plot total emissions to sanity check the grid, also to show total emissions
# and reveal where emissions matter. 
################################################################################
C_total = np.sum(C,axis=0)
C_total_ma = np.ma.masked_where(C_total==0, C_total)

fig=plt.figure(figsize=(8,8))
m.drawcoastlines()
m.drawstates()
m.drawcountries()
#m.drawmapboundary(fill_color='aqua')
c = m.pcolor(x, y,  C_total_ma)
cbar = plt.colorbar(c, pad=0.01, orientation='horizontal', extend='both',\
			norm=matplotlib.colors.LogNorm())
cbar.set_label('Emissions (grams carbon)', fontsize=25)
cbar.ax.tick_params(labelsize=15) 
#m.drawmapboundary(fill_color='aqua')
plt.title("Total June-Sept Emissions 2003-2016", fontsize=25)
fig.tight_layout()
plt.savefig(figureDir + 'total_C_emissions_6_9_2003_2016.png')
plt.close()

################################################################################
#----------------------- Load desired met event masks --------------------------
################################################################################

maskFile = '/barnes-scratch/sbrey/era_interim_nc_daily_merged/met_event_masks_NA_2003_2016.nc'
nc = Dataset(maskFile, 'r')

# TODO: change this to calling on the function that makes the mask to make a dynamic
# TODO: exploration of the space of cutoffs.


latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
stagnation_mask = nc.variables['stagnation_mask'][month_mask,:,:]
high_T_mask     = nc.variables['high_T_mask'][month_mask,:,:]
low_precip_mask = nc.variables['low_precip_mask'][month_mask,:,:]
high_wind_mask  = nc.variables['high_wind_mask'][month_mask,:,:]
low_RH_mask     = nc.variables['low_RH_mask'][month_mask,:,:]
blocking_mask   = nc.variables['blocking_mask'][month_mask,:,:]
# TODO: cyclone days coming soon
nc.close()

################################################################################
# Count the number of events at each grid cell for each type and make units
# total per summer day. 
################################################################################

uniqueYears = np.unique(year)
nYears = len(uniqueYears)
daysPerSummer = np.sum(time < datetime.datetime(uniqueYears[1],1,1), dtype=float)
n = daysPerSummer * nYears * 1.

# Sum across time dimension TODO: and make % June-September
nHighWindGrid = np.sum(high_wind_mask, axis=0) / n
nHighTGrid    = np.sum(high_T_mask,    axis=0) / n
nLowPrecGrid  = np.sum(low_precip_mask,  axis=0) / n
nStagnationGrid = np.sum(stagnation_mask, axis=0) / n
nLowRHGrid = np.sum(low_RH_mask, axis=0) / n
nBlockingGrid = np.sum(blocking_mask, axis=0) / n

# Count total occurance
nHighWind = np.sum(high_wind_mask)
nHighT    = np.sum(high_T_mask)
nLowPrec  = np.sum(low_precip_mask)
nStagnation = np.sum(stagnation_mask)
nLowRH = np.sum(low_RH_mask)
nBlocking = np.sum(blocking_mask)

################################################################################
# Map the occurance of each type of met event 
################################################################################

fig = plt.figure(figsize=(10,6))

ax = fig.add_subplot(2,3,1)

m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nHighWindGrid, vmin=0., vmax=0.1 )
cbar = m.colorbar(c, location='bottom', pad="1%", extend='max', ticks=[0, 0.05, .1])
cbar.set_label('proportion of days')
plt.title('Days mean wind > 8 m/s')

ax = fig.add_subplot(2,3,2)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nHighTGrid, vmin=0, vmax=1 ) # cmap="Purples"
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days T > 24 C')

ax = fig.add_subplot(2,3,3)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nLowPrecGrid, vmin=0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days with < 0.01 inches precip')

ax = fig.add_subplot(2,3,4)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nStagnationGrid, vmin = 0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Stagnation Days')

ax = fig.add_subplot(2,3,5)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nLowRHGrid, vmin = 0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('RH < 25%')

ax = fig.add_subplot(2,3,6)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nBlockingGrid, vmin = 0, vmax=0.1, cmap='Reds')
# nBlockingGrid.max() = 0.0398 for sd = 1
cbar = m.colorbar(c, location='bottom', pad="1%", ticks=[0, 0.05, 0.1])
cbar.set_label('proportion of days')
plt.title('500 mb 5 day blocking event')

fig.tight_layout()

plt.savefig(figureDir + 'era_interim_MetMaskCounts.png')
plt.close()

################################################################################
# Now we want to make the exact same plot, only we want to mask out the locations
# where emissions are low. We want to look at met events where emissions for 
# fires are high. Once those locations are known. Show exact same plot again
# only where emissions are high 
################################################################################

# Show the distribution of emissions. With and without zeros. 


#bins = np.logspace(1, C.max(), 100)
fig = plt.figure(figsize=(8,7))
plt.hist(C.flatten())#, bins=bins)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Emissions (g Carbon) day$^{-1}$ grid$^{-1}$ ', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.tick_params(labelsize=15) 
plt.title('distribution of daily emissions', fontsize=24)
plt.savefig(figureDir + 'dailyEmissionsDistribution.png')
plt.show()
plt.close()

# Mask out zero values, extra care for histogram. 
C_no_zero = np.ma.masked_where(C==0., C)
C_flat = C.flatten()
C_flat_noZero = C_flat[C_flat > 0.]

# make a cutoff value, we are going to choose to ignore small emission
# events for this analysis. 
cutoff = np.percentile(C_flat_noZero, 50)

allEmissions = np.sum(C_flat_noZero)
topEmissions = np.sum(C_flat_noZero[C_flat_noZero >= cutoff])


# Make some log scale bins for counting 
bins = np.logspace(1, 12, 30)

fig = plt.figure(figsize=(8,7))
plt.hist(C_flat_noZero, bins=bins)
plt.axvline(x=cutoff, color='red', linewidth=3)
plt.xscale('log')
#plt.yscale('log')
plt.xlabel('Emissions (g Carbon) day$^{-1}$ grid$^{-1}$ ', fontsize=20)
plt.ylabel('Count', fontsize=20)
plt.tick_params(labelsize=20) 
plt.title('distribution of daily emissions', fontsize=24)
fig.tight_layout()
plt.savefig(figureDir + 'dailyEmissionsDistribution_noZeros.png')
plt.show()
plt.close()


###############################################################################
# Figure out spatial location to retain for analysis. Also switch from C to
# 'E' which will generically refer to emissions. 
###############################################################################

ETotal = np.sum(C, axis=0) # g/grid integrated over all time
E = C

# Spatially, I tink it makes the most sense to ignore boxes that have no and 
# very low emissions. 
spatial_cutoff    = np.percentile(ETotal, 75)
highEMask = ETotal > spatial_cutoff
ETotalSum = np.sum(ETotal)

# Estimate proportion of total emissions these account for after masking out
# the low emission amounts 
highETotal = ma.masked_where(highEMask==False, ETotal)
highEmittersSum = np.sum(highETotal)
percentTotalERetained = highEmittersSum / ETotalSum * 100. 
print 'Percent of domain emissions retained with chosen mask: ' + str(percentTotalERetained)


###############################################################################
# Show the locations that account for the majority of emissions during this
# time period. 
###############################################################################
fig = plt.figure(figsize=(8,8))
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, highETotal)
cbar = m.colorbar(c, location='bottom',pad="1%")
cbar.set_label('carbon emitted [g]', fontsize=20)
cbar.ax.tick_params(labelsize=18) 
plt.title('highest emission locations, June-Sept 2003-2016', fontsize=20)
plt.savefig(figureDir + 'GFED4s_HighEmitterTotal_6_9_2003_2016.png')

###############################################################################
# Mask the met event totals by this high emission mask 
###############################################################################
nHighWindGrid = ma.masked_where(highEMask==False, nHighWindGrid)
nHighTGrid    = ma.masked_where(highEMask==False, nHighTGrid)
nLowPrecGrid  = ma.masked_where(highEMask==False, nLowPrecGrid)
nStagnationGrid = ma.masked_where(highEMask==False, nStagnationGrid)
nLowRHGrid = ma.masked_where(highEMask==False, nLowRHGrid)
nBlockingGrid = ma.masked_where(highEMask==False, nBlockingGrid)


# Now do it with the emissions in time, we do not ever want to look at the 
# locations where the values are cumulatively this low, but we do want the
# continuity of days with high emissions followed by days with low emissions. 
nTime = len(time)
E_masked_spatial = np.ma.empty(shape=E.shape, dtype=float) 
for i in range(nTime):
	E_masked_spatial[i,:,:] = np.ma.masked_where(highEMask==False, E[i,:,:])

# Plot the met event counter again now that they are masked to only show
# locations with significant cumulative emissions
fig = plt.figure(figsize=(10,6))

ax = fig.add_subplot(2,3,1)

m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nHighWindGrid, vmin=0., vmax=0.1 )
cbar = m.colorbar(c, location='bottom', pad="1%", extend='max', ticks=[0, 0.05, .1])
cbar.set_label('proportion of days')
plt.title('Days mean wind > 8 m/s')

ax = fig.add_subplot(2,3,2)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nHighTGrid, vmin=0, vmax=1 ) # cmap="Purples"
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days T > 24 C')

ax = fig.add_subplot(2,3,3)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nLowPrecGrid, vmin=0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days with < 0.01 inches precip')

ax = fig.add_subplot(2,3,4)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nStagnationGrid, vmin = 0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Stagnation Days')

ax = fig.add_subplot(2,3,5)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nLowRHGrid, vmin = 0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('RH < 25%')

ax = fig.add_subplot(2,3,6)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nBlockingGrid, vmin = 0, vmax=0.1 )
cbar = m.colorbar(c, location='bottom',pad="1%", extend='max', ticks=[0, 0.05, 0.1])
cbar.set_label('proportion of days')
plt.title('500 mb 5 day blocking event')

fig.tight_layout()

plt.savefig(figureDir + 'era_interim_MetMaskCounts_highEmitters.png')

################################################################################
# Make masked arrays for each emission day type. DO NOT show where mask == 0
################################################################################
#zeroMask = E_masked == 0

# highWindE = ma.masked_where(high_wind_mask == 0, E_masked_spatial)
# HighTE    = ma.masked_where(high_T_mask == 0, E_masked_spatial)
# LowPrecE  = ma.masked_where(low_precip_mask == 0, E_masked_spatial)
# stagE     = ma.masked_where(stagnation_mask == 0, E_masked_spatial)
# blockE    = ma.masked_where(blocking_mask == 0, E_masked_spatial)
# lowRHE    = ma.masked_where(low_RH_mask == 0, E_masked_spatial)
# NOTE: This implies that E is all summer E

# We also need to mask where emissions are zero, because zero is non-interesting 
# highWindE = ma.masked_where(E_masked == 0, highWindE)
# HighTE    = ma.masked_where(E_masked == 0, HighTE)
# LowPrecE  = ma.masked_where(E_masked == 0, LowPrecE)
# stagE     = ma.masked_where(E_masked == 0, stagE)
# blockE    = ma.masked_where(E_masked == 0, blockE)
# lowRHE    = ma.masked_where(E_masked == 0, lowRHE)
# E_noZero  = ma.masked_where(E_masked == 0, E_masked)

################################################################################
# TODO: Map total emissions for each met event type 
################################################################################

# NOTE: E_masked_spatial is E masked in the spatial locations that did not meet
# NOTE: the cutoff value 

def arrayOfNonZeroMetEmission(metMask, E_masked_spatial):
	metE = E_masked_spatial[metMask == 1]
	metE = metE[metE > 0.]
	metE = ma.compressed(metE)
	return metE

# Flatten all of these values, trailing _ indicates flattened to 1D, removes 
# masked values 
highWindE_ = arrayOfNonZeroMetEmission(high_wind_mask, E_masked_spatial)
HighTE_    = arrayOfNonZeroMetEmission(high_T_mask, E_masked_spatial)
LowPrecE_  = arrayOfNonZeroMetEmission(low_precip_mask, E_masked_spatial)
stagE_     = arrayOfNonZeroMetEmission(stagnation_mask, E_masked_spatial)
blockE_    = arrayOfNonZeroMetEmission(blocking_mask, E_masked_spatial)
lowRHE_    = arrayOfNonZeroMetEmission(low_RH_mask, E_masked_spatial)
E_         = ma.compressed(E_masked_spatial)
E_         = E_[E_ > 0.]

# First make a multihist of these met events. 
bins = np.logspace(5, 12, 40)
x = np.cumsum(np.diff(bins))

# Count the amount in each emission bin for each type of met event
E_heights         = plt.hist(E_, bins, label='E | No meteorology event', alpha=0.2)[0]
HighTE_heights    = plt.hist(HighTE_, bins, label='E | high T', alpha=0.2)[0]
highWindE_heights = plt.hist(highWindE_, bins, label = 'E | high wind', alpha=0.2)[0]
LowPrecE_heights  = plt.hist(LowPrecE_, bins, label = 'E | low precip', alpha=0.2)[0]
stagE_heights     = plt.hist(stagE_, bins, label = 'E | stagnation', alpha=0.2)[0]
blockE_heights    = plt.hist(blockE_, bins, label = 'E | 500mb blocking', alpha=0.2)[0]
lowRHE_heights    = plt.hist(lowRHE_, bins, label = 'E | low RH%', alpha=0.2)[0]
plt.xscale('log')
plt.show(block=False)
plt.close()

# Create data structure to plot nice side by side by histograms
lw=2
fig = plt.figure(figsize=(10,10))
ax= fig.add_subplot(1,1,1)
#plt.subfigure(1,1)
plt.plot(x, E_heights, label='No meteorology event', linewidth=lw)
plt.plot(x, HighTE_heights, label='E | high T', linewidth=lw)
plt.plot(x, highWindE_heights, label='E | high wind', linewidth=lw)
plt.plot(x, LowPrecE_heights, label='E | low precip', linewidth=lw)
plt.plot(x, stagE_heights, label='E | stagnation', linewidth=lw)
plt.plot(x, blockE_heights , label = 'E | 500mb blocking', linewidth=lw)
plt.plot(x, lowRHE_heights, label='E | low RH%', linewidth=lw)

# Set the look of the plot
plt.xscale('log')
ax.tick_params(labelsize=20) 
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
#plt.xlim([0, 0.00004])
plt.xlabel('g day$^{-1}$ grid$^{-1}$', fontsize=26)
plt.ylabel('Total grid box days with given emission', fontsize=26)
plt.legend(loc='best', frameon=False, fontsize=18)

fig.tight_layout()

plt.savefig(figureDir + 'lineHistogram_GDFED_era_interim.png')

################################################################################
# Now normalize the occurance of these events and plot the lines
# of the heights of the histogram. _n stands for normalized.
################################################################################

def normalizeE(heights):
	heights_n = heights / np.sum(heights)
	return heights_n

E_heights_n         = normalizeE(E_heights)
HighTE_heights_n    = normalizeE(HighTE_heights)
highWindE_heights_n = normalizeE(highWindE_heights)
LowPrecE_heights_n  = normalizeE(LowPrecE_heights)
stagE_heights_n     = normalizeE(stagE_heights) 
blockE_heights_n    = normalizeE(stagE_heights) 
lowRHE_heights_n    = normalizeE(stagE_heights) 


fig = plt.figure(figsize=(10,10))
ax= fig.add_subplot(1,1,1)
#plt.subfigure(1,1)
plt.plot(x, E_heights_n, label='No meteorology event', linewidth=lw)
plt.plot(x, HighTE_heights_n, label='E | high T', linewidth=lw)
plt.plot(x, highWindE_heights_n, label='E | high wind', linewidth=lw)
plt.plot(x, LowPrecE_heights_n, label='E | low precip', linewidth=lw)
plt.plot(x, stagE_heights_n, label='E | stagnation', linewidth=lw)
plt.plot(x, blockE_heights_n, label = 'E | 500mb blocking', linewidth=lw)
plt.plot(x, lowRHE_heights_n, label='E | low RH%', linewidth=lw)

# Set the look of the plot
plt.xscale('log')
ax.tick_params(labelsize=20) 
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
#plt.xlim([0, 0.00004])
plt.xlabel('g day$^{-1}$ grid$^{-1}$', fontsize=26)
plt.ylabel('Proportion of identified days', fontsize=26)
plt.legend(loc='best', frameon=False, fontsize=18)

fig.tight_layout()

plt.savefig(figureDir + 'lineHistogram_normalized_GDFED_era_interim.png')

plt.close()















