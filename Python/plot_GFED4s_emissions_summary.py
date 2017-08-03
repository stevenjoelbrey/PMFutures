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

# TODO: Make paths dynamic 

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

drive = "/barnes-scratch/sbrey/"
dataDirBase = drive + "era_interim_nc_daily_merged/"
figureDir = "../Figures/GFED_era_interm_analysis/"

################################################################################
#------------------ Subset model emissions in time -----------------------------
# -------- For spatial subsetting use cesm_nc_manager.mask2dims() --------------
################################################################################
# For report use western U.S. only. That is what you get when you load a file
# that has "_NA_" in the name. 

startMonth = 6
endMonth   = 9
region     = "_west_" # "_west_"| "_PNW_" | "_CAL_" | "_NorthRockies_" 

# Get region lat lon range	
minLat, maxLat, minLon, maxLon  = cnm.getRegionBounds(region)

# Get emissions, use this to get dimensions
# TODO: load global emissions, let region sorting take care of the rest. 
ncFile  = drive + "GFED4s/GFED4.1s_METGrid_C_NA_2003_2016.nc"
nc = Dataset(ncFile, 'r')
latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
time = nc.variables['time'][:]
C = nc.variables['C'][:]
nc.close()

# Make time into datetime arrays
time, month, year = cnm.get_era_interim_time(time)

# Spatially subset the data 
C, ynew, xnew = cnm.mask2dims(C, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)

################################################################################
# Show emissions time series for the domain
################################################################################
C_daily_total = np.sum(C,axis=(1,2))
C_cumulative = np.cumsum(C_daily_total)

fig = plt.figure(figsize=(12,8))
ax = plt.subplot(111)
plt.plot(time, C_daily_total)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', labelsize=20)
plt.xlabel("date", fontsize=26)
plt.ylabel("grams carbon emitted", fontsize=26)
plt.title("GDED4.1s Daily Emissions", fontsize=29)
plt.savefig(figureDir + "daily_timeSeries" +region+".png")
plt.close()

# TODO: Also plot different region contributions lines! That would be dope!

################################################################################
# show emissions monthly histogram
################################################################################
uniqueMonths = np.unique(month)
monthTotal = np.zeros(12)
for i in range(12):
	monthMask = month == uniqueMonths[i]
	monthTotal[i] = np.sum(C_daily_total[monthMask])

fig = plt.figure(figsize=(12,8))
ax = plt.subplot(111)
plt.bar(uniqueMonths, monthTotal)
plt.xticks(uniqueMonths, uniqueMonths)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', labelsize=20)
plt.xlabel("Month", fontsize=26)
plt.ylabel("grams carbon emitted", fontsize=26)
plt.title("GDED4.1s Seasonality", fontsize=29)
plt.savefig(figureDir + 'emissions_seasonality' + region + ".png")
plt.close()


################################################################################
# Show total summer emissions
################################################################################
uniqueYears = np.unique(year)
nYears = len(uniqueYears)
C_summer = np.zeros(nYears)
C_june = np.zeros(nYears)
C_july = np.zeros(nYears)
C_aug = np.zeros(nYears)
C_sept = np.zeros(nYears)

for i in range(nYears):
	summerMask = (month >= 6.) & (month <= 9.) & (year == uniqueYears[i])
	juneMask = (month == 6.) & (year == uniqueYears[i])
	julyMask = (month == 7.) & (year == uniqueYears[i])
	augMask  = (month == 8.) & (year == uniqueYears[i])
	septMask  = (month == 9.) & (year == uniqueYears[i])
	
	C_summer[i] = np.sum(C_daily_total[summerMask])
	C_june[i] = np.sum(C_daily_total[juneMask])
	C_july[i] = np.sum(C_daily_total[julyMask])
	C_aug[i] = np.sum(C_daily_total[augMask])
	C_sept[i] = np.sum(C_daily_total[septMask])

# First plot with all one color 
fig = plt.figure(figsize=(12,8))
ax = plt.subplot(111)
p1 = plt.bar(uniqueYears, C_summer)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', labelsize=20)
plt.xlabel("date", fontsize=26)
plt.ylabel("grams carbon emitted", fontsize=26)
plt.title("GDED4.1s June-Sept Emissions", fontsize=29)
plt.savefig(figureDir + "summer_interannual_variability"+region+".png")
plt.close()

# Now stack the monthly contributions 
fig = plt.figure(figsize=(12,8))
ax = plt.subplot(111)
p1 = plt.bar(uniqueYears, C_june, color="blue")
p2 = plt.bar(uniqueYears, C_july, bottom=C_june, color="green")
p3 = plt.bar(uniqueYears, C_aug, bottom=(C_june+C_july), color="grey")
p4 = plt.bar(uniqueYears, C_sept, bottom=(C_june+C_july+C_aug), color="lightpink")
plt.legend( (p1[0], p2[0], p3[0], p4[0]), ("June", "July", "Aug", "Sept"), 
			frameon=False, loc="best", fontsize=27)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', labelsize=20)
plt.xlabel("date", fontsize=26)
plt.ylabel("grams carbon emitted", fontsize=26)
plt.title("GDED4.1s Summer Emissions", fontsize=29)
plt.savefig(figureDir + "summer_interannual_variability_months"+region+".png")
plt.close()

################################################################################
# Make and apply month mask to start all summer mask analysis 
################################################################################

month_mask = (month >= startMonth) & (month <= endMonth)
C = C[month_mask,:,:]
time = time[month_mask]
month = month[month_mask]
year = year[month_mask]

################################################################################
# Set up a relevant map to use later
################################################################################
m = Basemap(projection='merc',llcrnrlat=minLat, urcrnrlat=maxLat,\
            llcrnrlon=minLon, urcrnrlon=maxLon, resolution=resolution,\
	    	lon_0=0, lat_0=-90)

# grid coords for mesh plotting of values. 
lons, lats = np.meshgrid(xnew, ynew)
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
c = m.pcolor(x, y,  C_total_ma, cmap='viridis_r')
cbar = plt.colorbar(c, pad=0.01, orientation='horizontal', extend='both',\
			norm=matplotlib.colors.LogNorm())
cbar.set_label('Emissions (grams carbon)', fontsize=25)
cbar.ax.tick_params(labelsize=15) 
#m.drawmapboundary(fill_color='aqua')
plt.title("Total June-Sept Emissions 2003-2016", fontsize=25)
fig.tight_layout()
plt.savefig(figureDir + 'total_C_emissions_6_9_2003_2016' + region +'.png')
plt.close()

################################################################################
#----------------------- Load desired met event masks --------------------------
################################################################################

maskFile = drive + 'era_interim_nc_daily_merged/met_event_masks_NA_2003_2016.nc'
nc = Dataset(maskFile, 'r')

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

# Spatially subset these masks 
high_T_mask, ynew, xnew = cnm.mask2dims(high_T_mask, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
low_precip_mask, ynew, xnew = cnm.mask2dims(low_precip_mask, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
stagnation_mask, ynew, xnew = cnm.mask2dims(stagnation_mask, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
high_wind_mask, ynew, xnew = cnm.mask2dims(high_wind_mask, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
low_RH_mask, ynew, xnew = cnm.mask2dims(low_RH_mask, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
blocking_mask, ynew, xnew = cnm.mask2dims(blocking_mask, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)

# Now that we are done with using the old bounds altogether. 
longitude = xnew
latitude = ynew 

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
c = m.pcolor(x, y, nHighWindGrid, vmin=0., vmax=0.1, cmap='viridis_r', )
cbar = m.colorbar(c, location='bottom', pad="1%", extend='max', ticks=[0, 0.05, .1])
cbar.set_label('proportion of days')
plt.title('Days mean wind > 8 m/s')

ax = fig.add_subplot(2,3,2)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nHighTGrid, vmin=0, vmax=1, cmap='viridis_r') 
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days T > 24 C')

ax = fig.add_subplot(2,3,3)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nLowPrecGrid, vmin=0, vmax=1, cmap='viridis_r')
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days with < 0.01 inches precip')

ax = fig.add_subplot(2,3,4)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nStagnationGrid, vmin = 0, vmax=1, cmap='viridis_r' )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Stagnation Days')

ax = fig.add_subplot(2,3,5)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nLowRHGrid, vmin = 0, vmax=1, cmap='viridis_r' )
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

plt.savefig(figureDir + 'era_interim_MetMaskCounts'+region+'.png')
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
plt.savefig(figureDir + 'dailyEmissionsDistribution'+region+'.png')
plt.show(block=False)
plt.close()

# Mask out zero values, extra care for histogram. 
C_no_zero = np.ma.masked_where(C==0., C)
C_flat = C.flatten()
C_flat_noZero = C_flat[C_flat > 0.]

# make a cutoff value, we are going to choose to ignore small emission
# events for this analysis. 
cutoff = np.percentile(C_flat_noZero, 50) # TODO: make this cutoff an argument
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
plt.savefig(figureDir + 'dailyEmissionsDistribution_noZeros'+region+'.png')
plt.close()


###############################################################################
# Figure out spatial location to retain for analysis. Also switch from C to
# 'E' which will generically refer to emissions. 
###############################################################################
ETotal = np.sum(C, axis=0) # g/grid integrated over all time
E = C

# Spatially, I tink it makes the most sense to ignore boxes that have no and 
# very low emissions. 
spatial_cutoff    = np.percentile(ETotal, 75) # TODO: make percentile an argument
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
c = m.pcolor(x, y, highETotal, cmap='viridis_r')
cbar = m.colorbar(c, location='bottom',pad="1%")
cbar.set_label('carbon emitted [g]', fontsize=20)
cbar.ax.tick_params(labelsize=18) 
plt.title('highest emission locations, June-Sept 2003-2016', fontsize=20)
plt.savefig(figureDir + 'GFED4s_HighEmitterTotal_6_9_2003_2016'+region+'.png')
plt.close()

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
c = m.pcolor(x, y, nHighWindGrid, vmin=0., vmax=0.1, cmap='viridis_r')
cbar = m.colorbar(c, location='bottom', pad="1%", extend='max', ticks=[0, 0.05, .1])
cbar.set_label('proportion of days')
plt.title('Days mean wind > 8 m/s')

ax = fig.add_subplot(2,3,2)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nHighTGrid, vmin=0, vmax=1, cmap='viridis_r' ) 
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days T > 24 C')

ax = fig.add_subplot(2,3,3)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nLowPrecGrid, vmin=0, vmax=1, cmap='viridis_r' )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days with < 0.01 inches precip')

ax = fig.add_subplot(2,3,4)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nStagnationGrid, vmin = 0, vmax=1, cmap='viridis_r' )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Stagnation Days')

ax = fig.add_subplot(2,3,5)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nLowRHGrid, vmin = 0, vmax=1, cmap='viridis_r' )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('RH < 25%')

ax = fig.add_subplot(2,3,6)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nBlockingGrid, vmin = 0, vmax=0.1, cmap="Reds" )
cbar = m.colorbar(c, location='bottom',pad="1%", extend='max', ticks=[0, 0.05, 0.1])
cbar.set_label('proportion of days')
plt.title('500 mb 5 day blocking event')

fig.tight_layout()

plt.savefig(figureDir + 'era_interim_MetMaskCounts_highEmitters'+region+'.png')
plt.close()

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
#plt.show(block=False)
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

plt.savefig(figureDir + 'lineHistogram_GDFED_era_interim'+region+'.png')
plt.close()

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
plt.xlabel('g day$^{-1}$ grid$^{-1}$', fontsize=26)
plt.ylabel('Proportion of identified days', fontsize=26)
plt.legend(loc='best', frameon=False, fontsize=18)

fig.tight_layout()

plt.savefig(figureDir + 'lineHistogram_normalized_GDFED_era_interim'+region+'.png')
plt.close()


