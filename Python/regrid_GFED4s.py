#!/usr/bin/env python2


###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# The purpose of this script is to regrid GFED data to ecmwf era-interim grid
# (or any lat lon passed) and save as a single daily nc file. 

# Example using scipy 
# http://christopherbull.com.au/python/scipy-interpolate-griddata/


import sys # for reading command line arguments
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
import numpy as np
import scipy.interpolate
import time as timer


print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	year      =  str(sys.argv[1]) 
	species   =  str(sys.argv[2])

else: 
	# Development environment. Set variables manually here. 
	year     = str(2003)
	species  = 'C' # 'C' , 'DM', 'small_fire_fraction' (These have daily fraction est.)


# Consider making this directory a command line argument
dataDir  = '/barnes-scratch/sbrey/'
sanityCheck = False

# Get the era-interim (or other MET grid). They ALL live on the same grid, except
# time
# TODO: Handle time matching in this script as well!
metFile = dataDir + 'era_interim_nc_daily/sp_'+ year + '.nc' 
met_nc = Dataset(metFile, 'r')
met_lat = met_nc.variables['latitude'][:]
met_lon = met_nc.variables['longitude'][:]
met_time= met_nc.variables['time']
sp_slice = met_nc.variables['sp'][180,:,:]


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
fire_file = dataDir + 'GFED4s/GFED4.1s_' + species + '_'+ year + '.nc'
fire_nc   = Dataset(fire_file, 'r')
fire_lat  = fire_nc.variables['latitude'][:] 
fire_lon  = fire_nc.variables['longitude'][:] 
fire_time = fire_nc.variables['time']

print '----------------------------------------------------------------------'
print 'Working on loading the really big emissions file' 
print '----------------------------------------------------------------------'
timeStart = timer.time()
emissions = fire_nc.variables[species][:]
dt = (timer.time() - timeStart) / 60.
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to load emissions array into workspace'
print '----------------------------------------------------------------------'


print '----------------------------------------------------------------------'
print 'Working on shifting emissions from -180:180 grid to 0:360 grid.' 
print '----------------------------------------------------------------------'
regridStart = timer.time()

emissions_new_grid = np.zeros(emissions.shape)
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
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to regrid the emissions array'
print '----------------------------------------------------------------------'

##############
# Sanity plot 
##############
if sanityCheck:

	emissions_new_sum = np.sum(emissions_new_grid, axis=0)
	emissions_ma = np.ma.masked_where(emissions_new_sum==0,emissions_new_sum, copy=True)
	m = Basemap(projection='robin',lon_0=0,resolution='c')
	lons, lats = np.meshgrid(fire_lon_new, fire_lat)
	x,y=m(lons,lats)

	m = Basemap(projection='robin', lon_0=0,resolution='c')
	m.drawcoastlines()
	#m.fillcontinents(color='coral',lake_color='aqua')
	# draw parallels and meridians.
	m.drawparallels(np.arange(-90.,120.,30.))
	m.drawmeridians(np.arange(0.,360.,60.))
	c = m.pcolor(x, y, emissions_ma)
	#m.drawmapboundary(fill_color='aqua')
	plt.title("GFED GRID")
	plt.show()



# Get the total number of time steps
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
	#print str(float(t)/nTStep*100.) + " % complete"
	# Loop over fire lon values
	for x in range(len(fire_lon)):
		# to make it a point, loop over each lat for each lon
		for y in range(len(fire_lat)):
	
			# y and x represent the emissions data location we are
			# trying to assign!
	
			# Difference in regular grid arrays, units of deg
			dx = np.abs(fire_lon[x] - met_lon)
			dy = np.abs(fire_lat[y] - met_lat)
	
			# Index of best match
			xi = np.argmin(dx)
			yi = np.argmin(dy)
			#xi = np.where(dx == dx.min())[0][0]
			#yi = np.where(dy == dy.min())[0][0]

			# assign the data
			fire_new_grid[t, yi, xi] = (fire_new_grid[t, yi, xi] + emissions[t, y, x])

	# How much mass is lost or gained?
	deltaMass[t] = np.sum(fire_new_grid[t, :, :]) - np.sum(emissions[t, :, :])

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
nTime = len(fire_time)

outputFile = dataDir + 'GFED4s/' + 'GFED4.1s_METGrid_' + species + '_' +\
		str(year) + '.nc'


ncFile = Dataset(outputFile, 'w', format='NETCDF4')
ncFile.description = 'Daily GFED regridded to era-interim grid.'
ncFile.location = 'Global'
ncFile.createDimension('time',  nTime )
ncFile.createDimension('latitude', nLat )
ncFile.createDimension('longitude', nLon )

VAR_ = ncFile.createVariable(species,\
	            'f4',('time','latitude','longitude'))
VAR_.units = 'g ' + species + ' per grid cell per day'

# Create time variable	
time_ = ncFile.createVariable('time', 'i4', ('time',))
time_.units = fire_time.units

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



