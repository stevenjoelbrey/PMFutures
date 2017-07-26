#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script is designed to create a mask of stagnation days as defined by 
# Wang and Angell [1999]. The three variables read in and assessed for 
# stagnation condition are 500mb , 1000mb (later taken as SLP) geostrophic 
# winds, and precipitation. 
# TODO: UPDATE RH calculation to better estimate

# Set the directory where the data structure starts 
# data creation pipeline get_era_interim_data.py -> average6hourlyData.py, 
# merge_nc_data.py
dataDirBase= "/barnes-scratch/sbrey/era_interim_nc_daily_merged/"

import os
import os.path
import numpy as np
import sys
from netCDF4 import Dataset 
import time as timer
from datetime import timedelta
import datetime
import matplotlib.pyplot as plt

def find_blocking_days(sdFactor=1.):
	"""
	This function finds blocking days based on a very simple definition.
	Blocking days are defined as days when the 500 mb geopotential height
	is one standard deviation above the monthly mean for at least two days.
	This function takes one argument.
		Arguments:
			sdFactor: Number multiplied by monthly std when setting
					the threshold. 
		return:
			blocking days: An array of 1 and 0 indicating if blocking
					exists.
	"""

	z_nc    = Dataset(dataDirBase + 'z_NA_2003_2016.nc', 'r')
	level   = z_nc.variables['level']
	level_i = np.where(level[:] == 500)[0][0]
	z       = z_nc.variables['z'][:,level_i,:,:]
	#lat = z_nc.variables['latitude'][:]
	#lon = z_nc.variables['longitude'][:]
	#x,y=np.meshgrid(lon,lat)
	time = z_nc.variables['time'][:]
	z_nc.close()

	# Make hourly time array into a nice datetime object 
	t0 = datetime.datetime(year=1900, month=1, day=1,\
                               hour=0, minute=0, second=0)

	# Make a nice month and time array for masking 
	nTime = len(time)
	t = []
	month = []
	for i in range(nTime):
		dt = timedelta(hours=int(time[i]))
		t_new  = t0 + dt 
		t.append(t_new)
		month.append(t_new.month)
	t = np.array(t)
	month = np.array(month)
	
	# Create month specific monthly threshold values 
	z_thresh = np.zeros((12, z.shape[1], z.shape[2]))
	for m in np.unique(month):
		monthMask       = m == month
		spatial_mean    = np.mean(z[monthMask,:,:], axis=0)
		spatial_std     = np.std(z[monthMask,:,:], axis=0)
		z_thresh[m-1,:,:] = spatial_mean + spatial_std * sdFactor

	
	# Because criteria needs two days of information, have to start at
	# the second day to look back into the past one day
	index = np.arange(1, nTime)
	blocking_mask = np.zeros(z.shape, dtype=bool)

	for i in index:

		# Need to know what month this day belongs to to choose
		# the correct climatological threshold.
		monthIndex = month[i] - 1

		high_z_yesterday       = z[i-1, :, :] >= z_thresh[monthIndex,:,:]
		high_z_today           = z[i, :, :] >= z_thresh[monthIndex,:,:]
		blocking_mask[i, :, :] = high_z_yesterday & high_z_today
		
		#fig = plt.figure()
		#c=plt.pcolor(x,y, blocking_mask[i, :, :])
		#bar = plt.colorbar(c)
		#plt.savefig('../Figures/block_test/block_test_'+str(i)+'.png')
		#plt.close()

		#fig = plt.figure()
		#c=plt.contour(x,y, z[i, :, :], vmin=z.min(), vmax=z.max())
		#bar = plt.colorbar(c)
		#plt.savefig('../Figures/block_test/z_'+str(i)+'.png')
		#plt.close()

	return blocking_mask


def make_era_interim_met_masks(windSfcLim=8., wind500Lim=13., precLim=0.01,
                               TThresh=297.039, RHThresh=25., windThresh=6.7056,
			       writeNC=True):
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

		These are the defualt definitions of stagnation defined:
		http://www.arl.noaa.gov/documents/reports/atlas.pdf

	"""

	startTime = timer.time()
	###############################################################################
	# Load surface winds
	# NOTE: x = store_x/scale + offset.
	###############################################################################
	u10_nc = Dataset(dataDirBase + 'u10_NA_2003_2016.nc', 'r')
	u10    = u10_nc.variables['u10'][:]
	#u10_   = u10[:] #/ u10.scale_factor + u10.scale_factor????
	u10_nc.close()

	v10_nc = Dataset(dataDirBase + 'v10_NA_2003_2016.nc', 'r')
	v10    = v10_nc.variables['v10'][:]
	v10_nc.close()

	sfc_wind = np.sqrt(v10**2 + u10**2)

	###############################################################################
	# Load 500 mb winds
	###############################################################################
	v_nc    = Dataset(dataDirBase + 'v_NA_2003_2016.nc', 'r')
	level   = v_nc.variables['level']
	level_i = np.where(level[:] == 500)[0][0]
	v       = v_nc.variables['v'][:,level_i,:,:]
	v_nc.close()

	u_nc    = Dataset(dataDirBase + 'u_NA_2003_2016.nc', 'r')
	u       = u_nc.variables['u'][:,level_i,:,:]
	u_nc.close()

	upper_wind = np.sqrt(v**2 + u**2)

	###############################################################################
	# Get precipitation
	###############################################################################
	tp_nc = Dataset(dataDirBase + 'tp_NA_2003_2016.nc', 'r')
	tp_meters = tp_nc.variables['tp'] # meters per calendar date
	inchPerM = 39.3701      # [inch/m]
	tp = tp_meters[:] * inchPerM
	latitude = tp_nc.variables['latitude']
	longitude = tp_nc.variables['longitude']
	time = tp_nc.variables['time']

	# Create numpy arrays to store mask of stagnation events
	stagnationMask = np.zeros(tp.shape, dtype=int)

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
	t2m_nc = Dataset(dataDirBase + 't2m_NA_2003_2016.nc', 'r')
	t2m = t2m_nc.variables['t2m'][:]
	t2m_nc.close()

	d2m_nc = Dataset(dataDirBase + 'd2m_NA_2003_2016.nc', 'r')
	d2m = d2m_nc.variables['d2m'][:]
	d2m_nc.close()

	# Calculate (and save) relative humidity
	T_C = t2m - 273.15
	Td_C = d2m - 273.15

	# http://andrew.rsmas.miami.edu/bmcnoldy/Humidity.html
	# TODO: Confirm this formula and estimate for RH
	#RH = 100.*(np.exp((17.625*Td_C) / (243.04+Td_C))/np.exp((17.625*T_C)/(243.04+T_C)))

	# Accepted approximation
	# http://journals.ametsoc.org/doi/pdf/10.1175/BAMS-86-2-225
	RH = 100. - 5. * (T_C - Td_C)

	# Make the masks 
	high_wind_mask  = np.array(sfc_wind >= windThresh, dtype=int)
	low_RH_mask     = np.array(RH < RHThresh, dtype=int)
	high_T_mask     = np.array(t2m >= TThresh, dtype=int)
	blocking_mask   = np.array(find_blocking_days(), dtype=int)
	
	writingComplete = timer.time()
	dt = (writingComplete - startTime) / 60. 
	print '----------------------------------------------------------------------'
	print 'It took ' + str(dt) + ' minutes to create the requested masks.'
	print '----------------------------------------------------------------------'

	###############################################################################
	# Write the stagnation mask as daily netCDF data
	###############################################################################
	if writeNC:

		saveName = dataDirBase + 'met_event_masks_NA_2003_2016.nc'

		ncFile = Dataset(saveName, 'w', format='NETCDF4')
		ncFile.description = 'Masks indicating threshold conditions'
		ncFile.location = 'Global'
		ncFile.createDimension('time', len(time[:]) )
		ncFile.createDimension('latitude', len(latitude[:]) )
		ncFile.createDimension('longitude', len(longitude[:]) )

		# Create variables on the dimension they live on 
		stagnation_mask_ = ncFile.createVariable('stagnation_mask', 'i', ('time','latitude','longitude'))
		stagnation_mask_.units = 'limts = surface wind >= ' + str(windSfcLim) +\
				         ' 500 mb wind lim < ' +str(wind500Lim) + 'precip < ' + str(precLim)  	

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
		time_[:]         = time[:]

		latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
		latitude_.units = latitude.units
		latitude_[:]     = latitude[:]	


		longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
		longitude_.units = longitude.units
		longitude_[:]    = longitude[:]

		ncFile.close()

		dt = (timer.time() - writingComplete) / 60. 
		print '----------------------------------------------------------------------'
		print 'It took ' + str(dt) + ' minutes to write the data as nc files.'
		print '----------------------------------------------------------------------'

