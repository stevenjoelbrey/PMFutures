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
# TODO: Make this into a class. Copy the stucture of Sam Atwoods code. 
# TODO: Update function documentation to be the accepted format. 

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

	nc.close()

	return bb, bbUnits, t, lat, lon

def loadEmissionGridAttributes(lat, lon):
	"""Load model grid cell atributes. These do not depend on time. Must be subset 
	   to match bb grid that is already loaded into workspace.
   	   Parameters:
           lon: longitude that defines the grid of emission data.
           lat: latitude that defines the grid of emission data. 
	
           return value: area of grid cells in units of m^2
	"""

	dataDir = '/fischer-scratch/sbrey/outputFromYellowstone/FireEmissions/'
	moduleLayers = ['area','landfrac','landmask']
	layer = {}
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
	
	   Parameters:
            bb:      data in kg/sec/m^2 to be changed. 
            area:    associated grid showing the area of each cell in m^2. 
            t:       the time array. Datetime object that describes time 
                     dimension
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
def makePcolorFig(ax, m, x, y, z, titleText, maxVal): 
    '''Creates a simple pcolor image of the NxM array passed in z. 
           Parameters:
               m: The map object to be plotted on.
               x: The x coords of m.
               y: The y coords of m.
               z: NxM array to be plotted using pColor
               titleText: The title of the figure returned by this func. 
    '''

    # Create the figure on the passes axis
    ax
    
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    
    #parallels = np.arange(0.,90,15.)
    #meridians = np.arange(180.,360.,15.)
    
    #m.drawparallels(parallels,labels=[1,0,0,0], fontsize=10)
    #m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    
    z = ma.masked_where(z == 0, z)
    cs = m.pcolor(x, y, z, cmap='Reds', vmin=0, vmax = maxVal)
    
    # add colorbar.
    cbar = m.colorbar(cs, location='bottom', pad='5%')
    cbar.set_label(variable, fontsize=11)
    
    # add title
    ax.set_title(titleText, fontsize=16)
    #plt.savefig(figureName, format='png')
    		    
    return(ax) 



def subsetModelEmissions(data, t, lat, lon, startMonth=1, endMonth=12, 
                         minLat=11., maxLat=90., minLon=190., maxLon=320.):
    '''This function subsets a three dimensional grid (time, lat, lon) based
    on the passed spatial and temporal arguments.
    Parameters:
        data:       The array(time, lat, lon) to be subset. 
        t:          data time dimension as datetime.date object
        lat:        The latitude of data
        lon:        The longitude of data
        startMonth: The first month to be considered for all years.
        endMonth:   The last month to be consider for all years. 
        minLat:     These four parameters subset the passed data array
        maxLat:
        minLon:
        maxLon: 
        
        return Value: dataSubset, tSubset, latSubset, lonSubset

    '''
    # NOTE:masks do not work on multidimensional arrays
    latMask = (lat >= minLat) & (lat <= maxLat)
    lonMask = (lon >= minLon) & (lon <= maxLon)
    lati = np.where(latMask == True)[0]
    loni = np.where(lonMask == True)[0]
    
    # TODO: Figure out qwhy the TODO is needed for last index,
    #       this is a strange python thing I do not understand. 
    lonSubset = lon[loni[0]:(loni[-1]+1)]
    latSubset = lat[lati[0]:(lati[-1]+1)]
    
    # Make the time subset 
    nt = len(t)
    mon = np.zeros(nt, dtype='int')
    timeMask = np.zeros(nt, dtype='bool')
    for i in range(nt):
        mon[i] = t[i].month
        
    tMask = (mon >= startMonth) & (mon <= endMonth)
    ti = np.where(tMask == True)[0] # you can have whatever you like
    tSubset = t[ti[0]:(ti[-1]+1)]
    
    # Super ugly statement for getting the subset of lat and lon
    dataSubset = data[ti[0]:(ti[-1]+1), lati[0]:(lati[-1]+1), loni[0]:(loni[-1]+1)]

    return dataSubset, tSubset, latSubset, lonSubset    
    
    
    
def timeSeries(ax1, s1, s2, s1Label, s2Label, t, titleText):    
    '''kgPerDay1 and kgPerDay2 are emissions data (time, lat, lon). This 
    function will sum over time and lon to plot total emissions vs. time in
    the entire spatial domain of the passed arrays. '''
  
    ax1
        
    ax1.plot(t, s1, 'b-', linewidth=1, label=s1Label)
    ax1.tick_params(axis='x', labelsize=11)
    ax1.tick_params(axis='y', labelsize=11)
        
    ax1.plot(t, s2, 'r-', linewidth=1, label=s2Label)
    ax1.tick_params(axis='y', labelsize=11)
    ax1.tick_params(axis='x', labelsize=11)
        
    # shared x-axis 
    ax1.set_xlabel('Time', fontsize=11)
    # Label y axes	
    ax1.set_ylabel('[kg]', fontsize=11)
    
    plt.title(titleText, fontsize=18)
    
    ax1.legend(frameon=False, loc='upper left')
    
    return(ax1)


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
		
		
def makeHist(ax, s1MonthTotal, s2MonthTotal, s1Label, s2Label, titleText):
    """Function makes lovely histogram (bar plot) for two series of data
    TODO: Make histogram actually lovely
    TODO: startMonth endMonth dynamic argument accept. 
    """
    n_groups = len(s1MonthTotal)
 
    ax

    index = np.arange(n_groups)
    bar_width = 0.35

    opacity = 1

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
    plt.xticks(index + bar_width, ('1', '2', '3', '4', '5', '6', '7', '8', \
                                    '9', '10', '11', '12'))
    plt.legend()

    plt.tight_layout()

    return(ax)


###############################################################################
# ------------------------- USER INPUT ---------------------------------------- 
###############################################################################
# TODO: Make sure these can work as agruments passed by the command line. 

variable   = 'BC'    # Anything listed in outputFromYellowstone/FireEmissions
startMonth = 1
endMonth   = 12
minLat     = 37.   # 10.
maxLat     = 41.   # 90.
minLon     = -109 + 360.  # 190.
maxLon     = -102. + 360. # 320

###############################################################################
# ------------------------- exe envire ---------------------------------------- 
###############################################################################

# Create loop[ that gives us output we want for plotting and saves
# everything into python dictionaries.
kgTotals     = {}; kgTotals_max     = []
tSeries      = {}; tSeries_max      = []
sMonthTotals = {}; sMonthTotals_max = []
ts           = {}
for scenario in ['RCP45', 'RCP85']:
    for year in ['2050', '2100']:

        # Key names
        NAME = scenario+year        
        
        # Load specified variable
        bb, bbUnits, t, lat, lon = getVariableData(variable, scenario, year)

        # Subset variable
        bb, t, lat, lon = subsetModelEmissions(bb, t, lat, lon, 
                                               startMonth, endMonth, 
                                               minLat, maxLat,
                                               minLon, maxLon)
        ts[year] = t
                                              
        # Get Associated area mask
        area = loadEmissionGridAttributes(lat, lon)

        # Summarize the data for plotting 
        kgPerDay, kgTotal = makeTotalEmissions(bb, area, t)
        kgTotals[NAME]    = kgTotal
        kgTotals_max.append(kgTotal.max())
        
        s                  = np.sum(kgPerDay, (1,2)) 
        tSeries[NAME]      = s 
        tSeries_max.append(s.max())
        
        sMonthTotals[NAME] = monthlyTotals(s, t)
        sMonthTotals_max.append(sMonthTotals[NAME].max())
        
# Reduce max lists to single values
kgTotals_max     = max(kgTotals_max)
tSeries_max      = max(tSeries_max)
sMonthTotals_max = max(sMonthTotals_max)

###############################################################################
# north america plot setup area 
# Set up the projection etc. Things that can be outside time loop
# TODO: Make a map projection that adjust to tghe limts passed above.
# TODO: Will probably have to go with a square projection. 
###############################################################################
m = Basemap(width=9000000, height=6000000,
            rsphere=(6378137.00, 6356752.3142),\
            resolution='l',area_thresh=10000.,projection='lcc',\
            lat_1=45.,lat_2=55,lat_0=50,lon_0=-105.)

lons, lats = np.meshgrid(lon, lat) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

fig = plt.figure(figsize=(12,12))

ax1 = fig.add_subplot(421)
makePcolorFig(ax1, m, x, y, z=kgTotals['RCP452050'], 
              titleText='RCP45 2050', maxVal=kgTotals_max)

ax2 = fig.add_subplot(422)
makePcolorFig(ax2, m, x, y, z=kgTotals['RCP452100'], 
              titleText='RCP45 2100', maxVal=kgTotals_max)

ax3 = fig.add_subplot(423)
makePcolorFig(ax3, m, x, y, z=kgTotals['RCP852050'], 
              titleText='RCP85 2050', maxVal=kgTotals_max)

ax4 = fig.add_subplot(424)
makePcolorFig(ax4, m, x, y, z=kgTotals['RCP852100'], 
              titleText='RCP85 2100', maxVal=kgTotals_max)

ax5 = fig.add_subplot(425)
timeSeries(ax5, tSeries['RCP452050'], tSeries['RCP852050'], 
           s1Label='RCP45', s2Label='RCP85', t=ts['2050'], titleText='2050')

ax6 = fig.add_subplot(426)
timeSeries(ax6, tSeries['RCP452100'], tSeries['RCP852100'], 
           s1Label='RCP45', s2Label='RCP85', t=ts['2100'], titleText='2100')

ax7 = fig.add_subplot(427)
makeHist(ax7, sMonthTotals['RCP452050'], sMonthTotals['RCP852050'],
         s1Label='RCP45', s2Label='RCP85', titleText='2050')

ax8 = fig.add_subplot(428)
makeHist(ax8, sMonthTotals['RCP452100'], sMonthTotals['RCP852100'],
         s1Label='RCP45', s2Label='RCP85', titleText='2100')


# Handle the amount of space between plots 
plt.subplots_adjust(wspace=0.5, hspace=0.5)

plt.show()




