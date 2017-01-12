#!/usr/bin/env python2

# plotEmissionSummary.py

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to summarize and plot emissions from CESM. 
# A Multipanel plot will sumarize the present day, 2050, 2100 emissions for a 
# selected variable. Later on, this script will be able to make the standard 
# plots for desired time and season subsets. Much later, this code will be
# incorperated into a shiny app that will allow easy exploration of the model 
# output. 

# Load resources
import cesm_nc_manager as cnm
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


###############################################################################
# ------------------------- USER INPUT ---------------------------------------- 
###############################################################################
# TODO: Make sure these can work as agruments passed by the command line. 

variable = 'BC'    # Anything listed in outputFromYellowstone/FireEmissions
scenario = 'RCP85' # Will be a loop
year     = '2100'  # Will be a loop

def getVariableData(variable, scenario, year):
	"""This function loads the emission 'variable' for the given secenario
       and year (ending decade of 11 year period of data). The time variable
	   is transformed to a python datetime.date object and assumes that there
       are no leap years in the model output (tested). 
	   
	   return value: bb, bbUnits, t, lat, lon

	"""

	# Load the data
	ncFile   = cnm.makeEmissionNCFile(variable, scenario, year)
	nc       = Dataset(ncFile, 'r')
	bb       = nc.variables['bb'][:]
	bbUnits  = nc.variables['bb'].units

	# Handle the massive fill_value for bb so math does not go crazy later
	bbMask = bb.mask
	bb.data[bbMask==1] = 0 # changes locations of data fill value from huge to 0

	# Handle time dimension. NOTE: Leap years do not exist in the model. 
	daysSincet0 = nc.variables['time'][:]   # only used for length of time dim
	 
	months      = [1,2,3,4,5,6,7,8,9,10,11,12]
	daysInMonth = [31,28,31,30,31,30,31,31,30,31,30,31] # Never a 29 day February
	nMonth      = len(daysInMonth)
	daysInYear  = np.sum(daysInMonth)
	nYears      = len(daysSincet0) / daysInYear 

	# start year is always 10 years before year selection string. 
	# start month and day are always Jan 1 for emissions
	startYear = int(year) - 10
	t0        = date(startYear, 1, 1) 

	t = [] # where time will be appended
	for deltaYear in range(nYears):
		# Advance the year
		YEAR = t0.year + deltaYear
		for m in range(nMonth):
			dim = daysInMonth[m]
			for d in np.linspace(1, dim, dim, dtype='int'):
				newDate = date(YEAR, months[m], d)
				t.append(newDate)
			
	# Make useful array
	t = np.array(t)

	# Get spacial dimensions
	lon = nc.variables['lon'][:]
	lat = nc.variables['lat'][:]

	nlon = len(lon)
	nlat = len(lat)

	nc.close()

	return bb, bbUnits, t, lat, lon

def loadEmissionGridAttributes(lat, lon):
	"""Load model grid cell atributes. These do not depend on time. Must be subset 
	   to match bb grid that is already loaded into workspace.
	
   	   Args:
   	   lon --- longitude that defines the grid of emission data.
       lat --- latitude that defines the grid of emission data. 
	
	   return value: area of grid cells in units of m^2
	"""

	dataDir = '/fischer-scratch/sbrey/outputFromYellowstone/FireEmissions/'
	moduleLayers = ['area','landfrac','landmask']
	layer = {}
	layerUnits = {}
	for f in moduleLayers:
		fileName = dataDir + 'cesm130_clm5_firemodule_'+ f + '_f09x125.nc'
		nc = Dataset(fileName, 'r')
		layer[f] = nc.variables[f][:]
		allLons = nc.variables['lon'][:]
		allLats = nc.variables['lat'][:]
		nc.close()

	# Get area from km**2 to m*2
	m2Perkm2 = 1e6 * 1.
	area = layer['area']
	area = area * m2Perkm2
	landMask = layer['landmask']

	# Get rid of the insanely large fill value that messes up multiplication. 
	area.data[landMask==0] = 0

	# Subset based on the dimensions where lat and lon match
	lonI = np.where(np.in1d(allLons, lon))[0]
	latI = np.where(np.in1d(allLats, lat))[0]

	# TODO: Handle this masking in a non horrible way. The plus one in the max
	#       end of indexing is required to include that maximum index. 
	areaSubset = area[latI.min():(latI.max()+1) , lonI.min():(lonI.max()+1)]
	#fill_value = area.fill_value
	area = areaSubset

	return area # later others may be needed as well. 

def makeTotalEmissions(bb, area, t): 
	"""This function takes change the units of emissions from kg/sec/m^2 to 
	   an equal size array that returns kg/day/gridcell and another array
	   that returns total emissions (kg) in time period, shape=(lat, lon). 	
	
	   Args:
		bb --- data in kg/sec/m^2 to be changed. 
		area --- associated grid showing the area of each cell in m^2. 
		t --- the time array. Datetime object that describes time dimension
			  of bb. 
		
	   return vaue: kgPerDay, kgTotal
	"""

	# Get emissions from kg/m2/sec to kg in time period of interest
	secondsPerDay = 24. * 60**2
	kgPerDay = np.zeros(bb.shape)
	for i in range(len(t)):
		kgPerDay[i, :, :] = bb[i,:,:] * (area.data * secondsPerDay)

	# Now sum over time to come up with maximum values
	kgTotal = np.sum(kgPerDay,0) * len(t)

	return kgPerDay, kgTotal



# TODO: add extent arg so that all have same limit
def makePcolorFig(m,x,y,z,titleText): 

	# Create the figure 
	fig = plt.figure(figsize=(4, 4))

	m.drawcoastlines()
	m.drawstates()
	m.drawcountries()

	parallels = np.arange(0.,90,15.)
	meridians = np.arange(180.,360.,15.)

	m.drawparallels(parallels,labels=[1,0,0,0], fontsize=10)
	m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

	z = ma.masked_where(z == 0, z)
	cs = m.pcolor(x, y, z, cmap='Reds')

	# add colorbar.
	cbar = m.colorbar(cs, location='bottom',pad='3%')
	cbar.set_label(variable, fontsize=18)

	# add title
	plt.title(titleText, fontsize=20)
	#plt.savefig(figureName, format='png')
		    
	return(fig) 


###############################################################################
# TODO: Build a function that subsets time and space based on a startDate, 
# TODO: endDate argument, minLon, maxLon, minLat, maxLat set of arguments. 
###############################################################################
def timeSeries(s1, s2, s1Label, s2Label, t, titleText):
	"""s1 and s2 are data time series for emission summed over a specified time
       and spatial extent. That work is done by another function. """

	# TODO: consider keeping on the same yaxis so it is easy to tell bigger values
	# TODO: trend lines maybe? Slope for decade in legend? 
	
	fig, ax1 = plt.subplots(figsize=(4,3))

	ax2 = ax1.twinx()
	ax1.plot(t, s1, 'b-', linewidth=1)
	ax1.tick_params(axis='x', labelsize=11)
	ax1.tick_params(axis='y', labelsize=11)

	ax2.plot(t, s2, 'r-', linewidth=1)
	ax2.tick_params(axis='y', labelsize=11)
	ax2.tick_params(axis='x', labelsize=11)

	# shared x-axis 
	ax1.set_xlabel('Time', fontsize=11)
	# Label y axes	
	ax1.set_ylabel(s1Label, color='b', fontsize=11)
	ax2.set_ylabel(s2Label, color='r', fontsize=11)

	plt.title(titleText, fontsize=18)

	return(fig)


def monthlyTotals(s, t):
	"""Counts monthly totals of s for all months of the year. Returns
	   an array that contains s totals for each month. 
    """

	# Create a place to store monthly total
	sMonthTotal = np.zeros(12)
	for i in range(len(t)):
		day = t[i]
		dayMonth                = day.month  # what month is this date in?
		monIndex                = dayMonth-1 # the index is month# - 1
		sMonthTotal[monIndex]  += s[i]       # add to ongoing month total 

	return sMonthTotal
		
		
def makeHist(s1MonthTotal, s2MonthTotal, s1Label, s2Label, titleText):
	"""Function makes lovely histogram (bar plot) for two series of data
	   TODO: Make histogram actually lovely
	   TODO: startMonth endMonth dynamic argument accept. 
	"""

	n_groups = len(s1MonthTotal)

	fig, ax = plt.subplots()

	index = np.arange(n_groups)
	bar_width = 0.35

	opacity = 1
	error_config = {'ecolor': '0.3'}

	rects1 = plt.bar(index, s1MonthTotal, bar_width,
                 	 alpha=opacity,
                	 color='b',
                 	 label=s1Label)

	rects2 = plt.bar(index + bar_width, s2MonthTotal, bar_width,
                 	 alpha=opacity,
                 	 color='r',
                 	 label=s2Label)

	plt.xlabel('Month')
	plt.ylabel('Emissions [passed units argument]')
	plt.title(titleText)
	plt.xticks(index + bar_width, ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'))
	plt.legend()

	plt.tight_layout()

	return(fig)



###############################################################################
# ------------------------- exe envire ---------------------------------------- 
###############################################################################
bb, bbUnits, t, lat, lon = getVariableData(variable, scenario, year)
area                     = loadEmissionGridAttributes(lat, lon)
kgPerDay, kgTotal        = makeTotalEmissions(bb, area, t)

# Start developing the actual plot. One vertical panel at a time. 


###############################################################################
# north america plot setup area 
# Set up the projection etc. Things that can be outside time loop
###############################################################################
m = Basemap(width=8000000,height=6500000,
            rsphere=(6378137.00, 6356752.3142),\
            resolution='l',area_thresh=10000.,projection='lcc',\
            lat_1=45.,lat_2=55,lat_0=50,lon_0=-105.)



lons, lats = np.meshgrid(lon, lat) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.


kgPerDayDomain = np.sum(kgPerDay, (1,2))




