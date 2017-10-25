#!/usr/bin/env python2

###############################################################################
# ------------------------- Description ---------------------------------------
###############################################################################
# This script is designed to create a mask of stagnation days as defined by
# Wang and Angell [1999]. The three variables read in and assessed for
# stagnation condition are 500mb , 1000mb (later taken as SLP) geostrophic
# winds, and precipitation.
# TODO: UPDATE RH calculation to better estimate
# TODO: Make a dynamic argument to choose all domain or "_NA_" domain.
# TODO: the current working version ony uses merged years NA domain.
# TODO: Make block identification not terrible.


# Set the directory where the data structure starts
# data creation pipeline:
#    get_era_interim_data.py
#       average6hourlyData.py,
#           merge_nc_data.py



import os
import os.path
import numpy as np
import sys
from netCDF4 import Dataset
import time as timer
from datetime import timedelta
import datetime
import matplotlib.pyplot as plt
plt.ioff() # makes it so figures only show if you say .show()
import cesm_nc_manager as cnm
from mpl_toolkits.basemap import Basemap, cm
from skimage import measure
from skimage import filters


# Figure out what machine this code is running on. Set file paths.
drive = cnm.getDrive()
dataDirBase = drive + "era_interim_nc_daily_merged/"

################################################################################
# ------------------------- Functions  -----------------------------------------
################################################################################


def find_blocking_days(sdFactor=0.5, startDate="2003-01-01", endDate="2016-12-31",\
					      minDays=3, plotBlocks=False, minBlobSize=15.):
	"""
	This function finds blocking days based on a very simple definition.
	Blocking days are defined as days when the 500 mb geopotential height
	is one standard deviation above the jDay mean (1979-2016) for at least
	five days. This function takes one argument.
		Arguments:
			sdFactor:  Number multiplied by monthly std when setting
					    the threshold.
			startDate: The first date to create blocking event mask for.
			endDate:   The last date to create blocking event mask for.
			minDays:   The minimum number of consecutive days required for high z
					    values to be considered a block.
			plotBlocks: True or False. If true the z climatology, daily value, and
			            identified area of blocking are plotted and saved. SLOW.
			minBlobSize: The minimum size of a blob in terms of degrees. Value
						   is squared.


		return:
			blocking days: An array of 1 and 0 indicating if blocking
					exists (where equal to 1).
	"""

	# get the start and end times
	startDate = datetime.datetime.strptime(startDate, '%Y-%m-%d')
	endDate   = datetime.datetime.strptime(endDate, '%Y-%m-%d')

	all_Z500_dir = drive + 'era_interim_nc_daily/'
	z_nc    = Dataset(all_Z500_dir + 'z_all_daily.nc', 'r')
	z       = z_nc.variables['z500']
	lat = z_nc.variables['latitude'][:]
	lon = z_nc.variables['longitude'][:]
	time = z_nc.variables['time']


	dx = np.abs(np.mean(np.diff(lon)))
	dy = np.abs(np.mean(np.diff(lat)))

	# Translates minBlobSpan degree argument to indecies needed to make this many
	# degrees in our latitude x longitude gridded data.
	blobSpan = np.round(minBlobSize / dx)

#	# For test plotting make bounds
#	minLat     = lat.min()
#	maxLat     = lat.max()
#	minLon     = lon[241]
#	maxLon     = lon[479]

#	map = Basemap(projection='robin',llcrnrlat=minLat, urcrnrlat=maxLat,\
#                 llcrnrlon=minLon, urcrnrlon=maxLon,resolution='c',\
#	    	       lon_0=0, lat_0=90)

	map = Basemap(projection='ortho',lon_0=-105,lat_0=60,resolution='l')


	# grid coords for mesh plotting of values.
	lons, lats = np.meshgrid(lon, lat)
	x, y = map(lons, lats)

	# Make a nice month and time array for masking
	t, month, year = cnm.get_era_interim_time(time)
	nTime = len(t)

	day = []
	for i in range(len(t)):
		day.append(t[i].day)
	day = np.array(day)

#	# Now sinces they are annoying, lets ignore Feb 29 all the time. Meaning
#	# we are getting rid of it in the time and z arrays.
#	notLeapDayMask = (month != 2) & (day != 29)
#
#	t = t[notLeapDayMask]
#	month = month[notLeapDayMask]
#	year  = year[notLeapDayMask]
#	day   = day[notLeapDayMask]
#	nTime = len(t)
#
#	if np.sum((month == 2) & (day == 29)) > 0:
#		raise ValueError('There is still a February 29th in the time series.')

	# Create Julian day specific threshold values based on that JDay
	# mean and sd for the ~39 years of reanalysis I am working with
	# the z_thresh will be equal spatial same shape as z but with
	# Julian day time axis.
	jDays = np.arange(1, 367)
	nJDays = len(jDays)

	# Create an array of the Julian dates associated with the time
	# axis of the Z data.
	t_jDay = []
	for i in range(len(t)):
		thisJDay = t[i].timetuple().tm_yday
		t_jDay.append(thisJDay)
	t_jDay = np.array(t_jDay)

	# Create the threshold mask.
	# NOTE: If we save the spatial_mean and spatial_std, masks could be created
	# NOTE: later on much easier.
	# NOTE: Some years June 1 has a different julain day because of leap year.
	# NOTE: so these methods will smooth things out a bit.
	spatial_mean = np.zeros((nJDays, z.shape[1], z.shape[2]))
	spatial_std  = np.zeros((nJDays, z.shape[1], z.shape[2]))
	z_thresh     = np.zeros((nJDays, z.shape[1], z.shape[2]))

	# Sample statistics summary.
	# To see if JDay 366 has way bigger std. Turns out it is smaller. Probably
	# not getting the true variability properly represented.
	julianDaySD = np.zeros((nJDays))
	nSamples    = np.zeros(nJDays)

	for i in range(nJDays):
		# Find mask for this jDay days in the record, should be ~40
		dayMask             = jDays[i] == t_jDay
		nSamples[i]         = np.sum(dayMask)

		spatial_mean[i,:,:] = np.mean(z[dayMask,:,:], axis=0)
		spatial_std[i,:,:]  = np.std(z[dayMask,:,:], axis=0)
		julianDaySD[i]      = np.mean(spatial_std[i,:,:])
		z_thresh[i,:,:]     = spatial_mean[i,:,:] + (spatial_std[i,:,:] * sdFactor)

	# Only create a blocking even mask for dates between the date arguments
	analysisDateIndex = np.where((t >= startDate) & (t <= endDate))[0]
	nAnalaysisDays = len(analysisDateIndex)

	# Create an array where each days blocking mask summary can be saved.
	blocking_mask = np.zeros(shape=(nAnalaysisDays, z.shape[1], z.shape[2]), dtype=bool)

	for i in range(nAnalaysisDays):

		# for this analysis date (i), where does that put us on jDays?
		print 'working on: ' + str(t[analysisDateIndex[i]])
		jDayIndex = np.where(t_jDay[analysisDateIndex[i]] == jDays)[0][0]

		# Here, the numbers are in reference to days in the past from
		# day we are tying to make a mask for
		high_z_masks = np.zeros((minDays, z.shape[1], z.shape[2]))
		for d in range(minDays):
			high_z_masks[d,:,:] = z[analysisDateIndex[i]-d, :, :] >= z_thresh[jDayIndex-d,:,:]

		# Figure out where these 2D arrays are all true, sum over days dimension
		block_count = np.sum(high_z_masks, axis=0)

		# Turn the boolean into a numeric 1 or 0 array
		ridgeMask = np.array(block_count == minDays, dtype=int)
		blocking_mask[i, :, :] = ridgeMask

		# For plotting, mask out where there are no blocks, easy to plot blocks
		ridgeMask_ma = np.ma.masked_where(ridgeMask==0, ridgeMask)
		ridge_z_values_ma = np.ma.masked_where(ridgeMask==0, z[analysisDateIndex[i],:,:])

		# Show this dates Z for plotting, divide by 100 m to show decameters
		# the way they do at: http://weather.rap.ucar.edu/upper/upaCNTR_500.gif
		todays_z = z[analysisDateIndex[i],:,:] / 100.
		todays_climo_z = spatial_mean[jDayIndex,:,:] / 100.

		#########################################################################
		# Set a minimum size requirement for block 'features'.
		# This of course means features have to be identified.
		# http://www.scipy-lectures.org/packages/scikit-image/auto_examples/plot_labels.html
		#########################################################################
		# TODO: Try 15 x 15 deg min size to dfine a block. This means we have to
		# TODO: check the blobs. I think we should do this from the centriod.
		im = ridgeMask
		blobs = im == 1 # Bool condition of blobs is where ridgemask == 1 by def

		all_labels = measure.label(blobs)
		blobs_labels = measure.label(blobs, background=0, connectivity=2)
		uniqueBobIDs = np.unique(blobs_labels)

		# Check each blob for minimum size requirement
		for b in uniqueBobIDs:
			if b !=0: # 0 is background so skip
				blobMask = b == blobs_labels
				blobArea = np.sum(blobMask)
				if blobArea < blobSpan**2:
					# I do not want this to remain a blob
					blobs_labels[blobMask] = 0


		# Mask non-blobs for plotting
		blobs_labels_ma = np.ma.masked_where(blobs_labels==0, blobs_labels)

		# TODO: !!!!!!!!!
		# Go through unique blob labels. The sum of the blog label met is the
		# total area. Use that as a cuttoff and get rid of blobs that are too
		# small.


		if plotBlocks:

			# Plot the mean field for this day also.

#			# plotting subset indicies
#			lon_g, lat_g = np.meshgrid(lon, lat)
			lat_i = (lat > 20)
#			m = (lon_g > 180.) & (lon_g < 360.) & (lat_g > 8.) & (lat_g < 80.)

			fig = plt.figure(figsize=(12,12))
			map.drawcoastlines()
			map.drawstates()
			map.drawcountries()

			# Plot the julain day climatology in pcolor shaded.
			c_z_climotology =  map.pcolor(x[lat_i,:], y[lat_i,:], todays_climo_z[lat_i,:])
			bar = plt.colorbar(c_z_climotology)

			# Show the identified features
			c_peaks = map.pcolor(x[lat_i,:], y[lat_i,:], blobs_labels_ma[lat_i,:], cmap="spectral")


			# Plot the daily heights as contours using the same colorbar
			c_height = map.contour(x[lat_i,:], y[lat_i,:], todays_z[lat_i,:], linewidths=4)
			plt.clabel(c_height, inline=1, fontsize=10)


			# shade out the area where we define a block using semi-transparent
			# shading.
			c_ridge = map.pcolor(x[lat_i,:], y[lat_i,:], ridgeMask_ma[lat_i,:],
						hatch="/.", alpha=0.)

			dateString = str(t[analysisDateIndex[i]])[0:10]
			plt.title('Date: ' + dateString + \
			          ' Julain day = ' + str(jDays[jDayIndex]))
			plt.savefig('../Figures/block_test/z_show_' + dateString\
			            + '_sd='+str(sdFactor)+\
						'_days='+str(minDays)+'_minBlobSize='+str(minBlobSize)+\
						'.png')
			plt.close(fig)





	# Finally close the very large nc file connection.
	z_nc.close()

	return blocking_mask

# TODO: Make years of analysis arguments? Rather than assume NA and 2003-2016?
# TODO: 'NA' needs to be changed to 'west' as it only covers 30. - 49.5 N and
# TODO: 234. - 258.75 E.
def make_era_interim_met_masks(windSfcLim=8., wind500Lim=13., precLim=0.01,
							   TThresh=297.039, RHThresh=25., windThresh=6.7056,
							   writeNC=True, region = "_"):
	"""
	This function takes limits and creates masks (1 condition is true, 0 condition
	is not true) for different meteorology event or threshold types.
		Argument Units:
			windSfcLim  = m/s
			wind500Lim  = m/s
			precLim     = inches/day

			# For single variable thresholds
			TThresh    = K
			RHThresh   = %
			windThresh = m/s
			writeNC    = True masks written to nc file. If false
					       they are not.
			region     = regions nc data are presliced into formatted as
			             '_regionName_'. '_' is no region and correct filename
						    formatting.

		These are the defualt definitions of stagnation defined:
		http://www.arl.noaa.gov/documents/reports/atlas.pdf

	"""

	startTime = timer.time()
	#############################################################################
	# Load surface winds
	# NOTE: x = store_x/scale + offset.
	#############################################################################


	u10_nc = Dataset(dataDirBase + 'u10' + region + '2003_2016.nc', 'r')
	u10    = u10_nc.variables['u10'][:]
	#u10_   = u10[:] #/ u10.scale_factor + u10.scale_factor????
	u10_nc.close()

	v10_nc = Dataset(dataDirBase + 'v10' + region + '2003_2016.nc', 'r')
	v10    = v10_nc.variables['v10'][:]
	v10_nc.close()

	sfc_wind = np.sqrt(v10**2 + u10**2)

	###############################################################################
	# Load 500 mb winds
	###############################################################################
	v_nc    = Dataset(dataDirBase + 'v'+ region +'2003_2016.nc', 'r')
	level   = v_nc.variables['level']
	level_i = np.where(level[:] == 500)[0][0]
	v       = v_nc.variables['v'][:,level_i,:,:]
	v_nc.close()

	u_nc    = Dataset(dataDirBase + 'u'+ region +'2003_2016.nc', 'r')
	u       = u_nc.variables['u'][:,level_i,:,:]
	u_nc.close()

	upper_wind = np.sqrt(v**2 + u**2)

	###############################################################################
	# Get precipitation
	###############################################################################
	tp_nc = Dataset(dataDirBase + 'tp'+ region +'2003_2016.nc', 'r')
	tp_meters = tp_nc.variables['tp'] # meters per calendar date
	inchPerM = 39.3701      # [inch/m]
	tp = tp_meters[:] * inchPerM
	latitude = tp_nc.variables['latitude']
	longitude = tp_nc.variables['longitude']
	time = tp_nc.variables['time']

	# build the individual masks, first tp (total precipitation)
	mask_sfc  = np.array(sfc_wind < windSfcLim, dtype=bool)
	mask_500  = np.array(upper_wind < wind500Lim, dtype=bool)
	mask_tp   = np.array(tp < precLim, dtype=bool)

	# Combined stagnation mask
	stagnation_mask = np.array(mask_sfc & mask_500 & mask_tp, dtype=int)


	###############################################################################
	# Sanity check the output before writing the mask to an nc file
	###############################################################################
	if tp[mask_tp].max() >= precLim:
		print 'The maximum value of precip on stangation days exceeds threshold!'
		raise ValueError("This means creating the mask has failed!!!!!")


	###############################################################################
	# Now make individual masks for high wind, T, and low prec and RH days
	# http://w1.weather.gov/glossary/index.php?word=Red%20Flag%20Warning
	# For red flat warning:
	# 	T > 75 F = 297.039 K
	#	RH% <= 25%
	#	surface wind >= 15 mph = 6.7056 m/s
	###############################################################################
	t2m_nc = Dataset(dataDirBase + 't2m'+ region +'2003_2016.nc', 'r')
	t2m = t2m_nc.variables['t2m'][:]
	t2m_nc.close()

	d2m_nc = Dataset(dataDirBase + 'd2m'+ region +'2003_2016.nc', 'r')
	d2m = d2m_nc.variables['d2m'][:]
	d2m_nc.close()

	# Calculate (and save) relative humidity
	T_C = t2m - 273.15
	Td_C = d2m - 273.15

	# http://andrew.rsmas.miami.edu/bmcnoldy/Humidity.html
	# TODO: Confirm this formula and estimate for RH
	RH = 100.*(np.exp((17.625*Td_C) / (243.04+Td_C))/np.exp((17.625*T_C)/(243.04+T_C)))

	# Accepted approximation
	# http://journals.ametsoc.org/doi/pdf/10.1175/BAMS-86-2-225
	#RH = 100. - 5. * (T_C - Td_C)

	# Make the masks
	high_wind_mask  = np.array(sfc_wind >= windThresh, dtype=int)
	low_RH_mask     = np.array(RH < RHThresh, dtype=int)
	high_T_mask     = np.array(t2m >= TThresh, dtype=int)
	blocking_mask   = np.array(find_blocking_days(), dtype=int) # Calls comlex function
	low_precip_mask = np.array(tp < precLim, dtype=int)

	writingComplete = timer.time()
	dt = (writingComplete - startTime) / 60.
	print '----------------------------------------------------------------------'
	print 'It took ' + str(dt) + ' minutes to create the requested masks.'
	print '----------------------------------------------------------------------'

	###############################################################################
	# Write the stagnation mask as daily netCDF data
	###############################################################################
	if writeNC:

		saveName = dataDirBase + 'met_event_masks' + region + '2003_2016.nc'

		ncFile = Dataset(saveName, 'w', format='NETCDF4')
		ncFile.description = 'Masks indicating threshold conditions'
		ncFile.location = 'Global'
		ncFile.createDimension('time', len(time[:]) )
		ncFile.createDimension('latitude', len(latitude[:]) )
		ncFile.createDimension('longitude', len(longitude[:]) )

		# Create variables on the dimension they live on

		# Stagnation
		stagnation_mask_ = ncFile.createVariable('stagnation_mask', 'i', ('time','latitude','longitude'))
		stagnation_mask_.units = 'limts = surface wind >= ' + str(windSfcLim) +\
				         ' 500 mb wind lim < ' +str(wind500Lim) + 'precip < ' + str(precLim)
		stagnation_mask_[:] = stagnation_mask

		# Precip
		low_precip_mask
		low_precip_mask_ = ncFile.createVariable('low_precip_mask', 'i', ('time','latitude','longitude'))
		low_precip_mask_.units = 'days precip < ' + str(precLim) + ' inches/day'
		low_precip_mask_[:] = low_precip_mask[:]

		# wind
		high_wind_mask_ = ncFile.createVariable('high_wind_mask', 'i', ('time','latitude','longitude'))
		high_wind_mask_.units = 'days wind > ' + str(windThresh) + ' m/s'
		high_wind_mask_[:] = high_wind_mask[:]

		# RH
		low_RH_mask_ = ncFile.createVariable('low_RH_mask', 'i', ('time','latitude','longitude'))
		low_RH_mask_.units = 'RH% less than ' + str(RHThresh)
		low_RH_mask_[:] = low_RH_mask[:]

		# Temperature
		high_T_mask_ = ncFile.createVariable('high_T_mask', 'i', ('time','latitude','longitude'))
		high_T_mask_.units = 'T >= ' + str(TThresh)
		high_T_mask_[:] = high_T_mask[:]

		# Blocking days
		blocking_mask_ = ncFile.createVariable('blocking_mask', 'i', ('time','latitude','longitude'))
		blocking_mask_.units = 'Location date has z > monthly mean + std for today and yesterday.'
		blocking_mask_[:] = blocking_mask[:]

		# dimension values assignments
		time_ = ncFile.createVariable('time', 'i4', ('time',))
		time_.units = time.units
		time_[:] = time[:]

		latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
		latitude_.units = latitude.units
		latitude_[:]     = latitude[:]

		longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
		longitude_.units = longitude.units
		longitude_[:]    = longitude[:]

		ncFile.close()

		#######################################################################
		# Save RH as its own met variable!!!
		#######################################################################
		saveName = dataDirBase + 'RH_NA_2003_2016.nc'

		ncFile = Dataset(saveName, 'w', format='NETCDF4')
		ncFile.description = 'Mask indicating stangation days'
		ncFile.location = 'Global'
		ncFile.createDimension('time', len(time[:]) )
		ncFile.createDimension('latitude', len(latitude[:]) )
		ncFile.createDimension('longitude', len(longitude[:]) )

		RH_ = ncFile.createVariable('RH', 'i', ('time','latitude','longitude'))
		RH_.units = '%'
		RH_[:] = RH

		# dimension values assignments
		time_ = ncFile.createVariable('time', 'i4', ('time',))
		time_.units = time.units
		time_[:] = time[:]

		latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
		latitude_.units = latitude.units
		latitude_[:] = latitude[:]

		longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
		longitude_.units = longitude.units
		longitude_[:]    = longitude[:]

		ncFile.close()

		dt = (timer.time() - writingComplete) / 60.
		print '----------------------------------------------------------------------'
		print 'It took ' + str(dt) + ' minutes to write the data as nc files.'
		print '----------------------------------------------------------------------'



