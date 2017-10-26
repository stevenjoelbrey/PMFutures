#!/usr/bin/env python2

###############################################################################
# ------------------------- Description ---------------------------------------
###############################################################################
# The purpose of this script is to regrid FINN or GFFED emissions data to ecmwf
# era-interim grid (or any lat lon passed) and save as a single daily nc file.
# There is a good deal of special handlign of HYSPLIT Points fire data.

# TODO: Calculate and save the grid sizes m**2 to written nc files

import os
import re # regular expression library
import sys # for reading command line arguments
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
import numpy as np
import scipy.interpolate
import time as timer
import cesm_nc_manager as cnm


print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	year      =  str(sys.argv[1])
	inventory =  str(sys.argv[2])
	METGrid   =  str(sys.argv[3])
	species   =  str(sys.argv[4])
	getBurnArea = False

else:
	# Development environment. Set variables manually here.
	year      = str(2003)
	inventory = 'GFED4s'  # FINN | GFED4s | HMS
	METGrid   = 'ecmwf'   #
	species   = 'monthly_DM'   # Either CO2 or a specific vegetation type or C for GFED, SPDH
	getBurnArea = True    # Only works for monthly data...Will get burn area from C or DM monthly files

# When this is true extra plots are made to ensure that the grids are being
# handled correctly. These were examined during development
sanityCheck = False

# Figure out what machine this code is running on. Set file paths.
dataDir = cnm.getDrive()


# Get the era-interim (or other MET grid). They ALL live on the same grid, except
# time
# TODO: Handle time matching in this script as well!
metFile = dataDir + 'era_interim_nc_daily/sp_'+ year + '.nc'
met_nc = Dataset(metFile, 'r')
met_lat = met_nc.variables['latitude'][:]
met_lon = met_nc.variables['longitude'][:]
met_time= met_nc.variables['time'][:]
sp_slice = met_nc.variables['sp'][180,:,:]
met_nc.close()

if sanityCheck:
	m = Basemap(projection='robin',lon_0=0,resolution='c')
	lons, lats = np.meshgrid(met_lon, met_lat)
	x,y=m(lons,lats)
	m.drawcoastlines()
	#m.fillcontinents(color='coral',lake_color='aqua')
	# draw parallels and meridians.
	#plt.pcolormesh(x,y, sp_slice)
	c = plt.contour(x, y, sp_slice )
	#m.drawmapboundary(fill_color='aqua')
	plt.title("ERA-Interim Grid")
	plt.show(block=False)
	#plt.close()

# Get emissions file (use yearly nc files)
if inventory == 'GFED4s':
	fire_file = dataDir + 'GFED4s/GFED4.1s_' + species + '_'+ year + '.nc'
elif inventory == 'FINN':
	fire_file = dataDir + 'FINN/FINN_' + species + '_'+ year + '.nc'
elif inventory == 'HMS':
	fire_file = dataDir + 'HMS/HYSPLITPoints_' + species  + '.nc'
else:
	print 'Unknown emissions type'

# Pull the time first, get an array of the years, find the indicies, then pull only
# these from the emissions variable. For all but HMS this should load of time indicies
# indo the workspace.
fire_nc   = Dataset(fire_file, 'r')
fire_time = fire_nc.variables['time']
fire_lat  = fire_nc.variables['latitude'][:]
fire_lon  = fire_nc.variables['longitude'][:]

# preserve for writing nc data
fire_time_units = fire_time.units


print '----------------------------------------------------------------------'
print 'Working on loading the really big emissions file for year ' + year
print '----------------------------------------------------------------------'
timeStart = timer.time()

# TODO: Burn area

# Handle different time format of the monthly data
if 'monthly' in species:

	fire_time = fire_time[:]

	si = species.index('_') + 1
	sf = len(species)
	grid_area = fire_nc.variables['grid_area'][:]


	if getBurnArea:
		emissions = fire_nc.variables['burn_area_fraction'][:]
		emissions_units = 'm**2'
		species = 'monthly_burned_area'

	else:
		emissions = fire_nc.variables[species[si:sf]][:]
		emissions_units = fire_nc.variables[species[si:sf]].units
		# The for loop has gotten rid of the area component so we need to remove
		# from the units
		emissions_units = re.sub('m\*\*-2 ','',emissions_units)

	for n in range(len(fire_time)):
		emissions[n,:,:] = emissions[n,:,:] * grid_area


else:

	# When daily data subset the arrays?
	t_temp, month_temp, year_temp = cnm.get_era_interim_time(fire_time)
	timeIndex = np.where(year_temp == int(year))[0]

	# Now subset the time array to the year of interest only
	fire_time = fire_time[timeIndex]

	# Remove variables that are no longer needed.
	del t_temp, month_temp, year_temp

	# Get daily emissions data
	emissions_units = fire_nc.variables[species].units
	emissions = fire_nc.variables[species][timeIndex,:,:]



dt = (timer.time() - timeStart) / 60.
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to load emissions array into workspace'
print '----------------------------------------------------------------------'


print '----------------------------------------------------------------------'
print 'Working on shifting emissions from -180:180 grid to 0:360 grid.'
print '----------------------------------------------------------------------'
regridStart = timer.time()

emissions_new_grid = np.zeros(emissions.shape) # same size overall
fire_lon_new = np.zeros(fire_lon.shape)

eastOfGreenwich = np.where(fire_lon >= 0)[0]
westOfGreenwich = np.where(fire_lon < 0)[0]

# Assign the eastern hemisphere
fire_lon_new[range(len(eastOfGreenwich))] = fire_lon[eastOfGreenwich]
emissions_new_grid[:,:,range(len(eastOfGreenwich))] = emissions[:,:,eastOfGreenwich]
# assign the western hemisphere
assignment = np.where(fire_lon_new==0)[0]
fire_lon_new[assignment] = fire_lon[westOfGreenwich] + 360.
emissions_new_grid[:,:, assignment] = emissions[:,:,westOfGreenwich]

# Replace old grids and longitude with new
emissions = emissions_new_grid
fire_lon = fire_lon_new


dt = (timer.time() - regridStart) / 60.
print '-----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to regrid (lons) of the emissions array'
print '-----------------------------------------------------------------------'

######################################################################################
# Subset this grid, by the maximum and minimum bounds of fires, for HMS inventory. No
# need to loop all over the world when there are only recorded fires in North America.
# This dramatically improves performance when HMS fires are those being regridded.
######################################################################################
if inventory == 'HMS':
	minLat = 10.
	maxLat = 85.
	minLon = 190.
	maxLon = 325.

	emissions, fire_lat, fire_lon = cnm.mask2dims(emissions, fire_lon, fire_lat,\
												  0, minLon, maxLon, minLat, maxLat)


######################################################################################
# Sanity plot for interactive analysis mode
######################################################################################
if sanityCheck:

	emissions_new_sum = np.sum(emissions, axis=0)
	emissions_ma = np.ma.masked_where(emissions_new_sum==0,emissions_new_sum, copy=True)
	m = Basemap(projection='robin',lon_0=0,resolution='c')
	lons, lats = np.meshgrid(fire_lon, fire_lat)
	x,y=m(lons,lats)

	m = Basemap(projection='robin', lon_0=0,resolution='c')
	m.drawcoastlines()
	#m.fillcontinents(color='coral',lake_color='aqua')
	# draw parallels and meridians.
	m.drawparallels(np.arange(-90.,120.,30.))
	m.drawmeridians(np.arange(0.,360.,60.))
	c = m.pcolor(x, y, emissions_ma)
	#m.drawmapboundary(fill_color='aqua')
	plt.title("REGRID GRID")
	plt.show()

# How make an array that represets all dates in this year to index the yearly data
# This is needed because HMS is missing dates
nTStep = len(fire_time)

# NOTE: float32 specified otherwise addition leads to small errors and mass
# NOTE: is not conserved in the regridding process.
fire_new_grid = np.zeros((nTStep, len(met_lat), len(met_lon)), dtype='float')
deltaMass = np.zeros(nTStep, dtype='float')

# For each time step, find the best match in terms of min lon lat
# distance
# TODO: Only loop over non zero locations? Over subregion?
timeStart = timer.time()

for t in range(nTStep): # nTStep
	#print str(t) + ' time step'


	# Where does this date for this year fall in fire_time?
	if 'monthly' in species:
		tMask = t
	else:
		# Make sure the date exists in both MET and the fire inventory
		tMask = met_time[t] == fire_time[:]

	# See if the date exists, if it does not that means there is no fire emission
	# data for this date and we can move forward with the time loop.
	if (np.sum(tMask) > 0) or (nTStep == 12):

		print str(np.round(float(t)/nTStep*100.)) + " % complete"

		# Loop over fire lon values
		for x in range(len(fire_lon)):
			# to make it a point, loop over each lat for each lon
			for y in range(len(fire_lat)):

				# y and x represent the emissions data location we are
				# trying to assign!

				# Difference in regular grid arrahmsys, units of deg
				dx = np.abs(fire_lon[x] - met_lon)
				dy = np.abs(fire_lat[y] - met_lat)

				# Index of best match
				xi = np.argmin(dx)
				yi = np.argmin(dy)

				# assign the data, new grid is t referenced, emissions tMask
				fire_new_grid[t, yi, xi] = (fire_new_grid[t, yi, xi] + emissions[tMask, y, x])

		# How much mass is lost or gained?
		deltaMass[t] = np.sum(fire_new_grid[t, :, :]) - np.sum(emissions[tMask, :, :])
	else:
		print str(met_time[t]) + '  - missing in fire_emissions'

dt = (timer.time() - timeStart) / 60.
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to run the re-grid loop'
print '----------------------------------------------------------------------'

######################################################################################
# Write the NETCDF data. Make sure to include all relevant units for a given species!
######################################################################################
print 'Working on writing the output as netCDF data'
nLat = len(met_lat)
nLon = len(met_lon)
nTime = nTStep # should match number of days in calendar year of chosen year

# TODO: Make this save name (and whole section) dynamic based on the passed met grid
# TODO: to use

# Handle emissions specific labels
if inventory == 'GFED4s':
	outputFile = dataDir + 'GFED4s/' + 'GFED4.1s_' + METGrid + '_' + species + '_' +\
		         str(year) + '.nc'

elif inventory == 'FINN':
	outputFile = dataDir + 'FINN/' + 'FINN_' + METGrid + '_' + species + '_' +\
		         str(year) + '.nc'

elif inventory == 'HMS':
	outputFile = dataDir + 'HMS/' + 'HYSPLITPoints_' + METGrid + '_' + species +'_' +\
	             str(year) + '.nc'


ncFile = Dataset(outputFile, 'w', format='NETCDF4')
ncFile.description = 'Daily ' + inventory + ' regridded to era-interim grid.'
ncFile.location = 'Global'
ncFile.createDimension('time',  nTime )
ncFile.createDimension('latitude', nLat )
ncFile.createDimension('longitude', nLon )

VAR_ = ncFile.createVariable(species, 'f4',('time','latitude','longitude'))

if inventory != 'HMS':
	VAR_.units = emissions_units
else:
	VAR_.units = 'smoke production duration hours'

# Create time variable
time_ = ncFile.createVariable('time', 'i4', ('time',))
time_.units = fire_time_units

# Keep track of the change in mass
deltaMass_ = ncFile.createVariable('deltaMass', 'f4', ('time',))
deltaMass_.units = 'change in grams whole grid each timestep due to regridding'

# create lat variable
latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
latitude_.units = 'degrees north'

# create longitude variable
longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
longitude_.units = 'degrees east'

# Write the actual data to these dimensions
VAR_[:]       = fire_new_grid[:]
latitude_[:]  = met_lat
longitude_[:] = met_lon
time_[:]      = fire_time[:]
deltaMass_[:] = deltaMass[:]

ncFile.close()

print 'Script executed without error!'

