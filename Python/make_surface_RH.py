#!/usr/bin/env python2
# -*- coding: utf-8 -*-


###############################################################################
# ------------------------- Description ---------------------------------------
###############################################################################
# This script will calculate global ground level RH% using ecmwf ERA interim
# grids. It will create yearly files that will require merging by
# merge_yearly_nc.py
#
# Follows ----------------------------------------
# 	- get_ERA_Interim_data.py -> average6HourlyData.py
# Precedes
#	- merge_yearly_nc.py (combine the yearly RH files made by this script)
#	- make_era_interim_met_masks.py
#	- any script that uses RH% in calculation or analysis
#
################################################################################
# This is the link with the description of how to best calculate surface RH%
# when using ECMWF reanalysis data.
# https://software.ecmwf.int/wiki/display/CKB/Do+ERA+datasets+contain+parameters+for+near-surface+humidity
# Page 94 of this pdf for tetons equation version used:
# https://www.ecmwf.int/sites/default/files/elibrary/2016/16648-part-iv-physical-processes.pdf
################################################################################

import cesm_nc_manager as cnm
import os
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import matplotlib.ticker as tkr
import datetime
import time as timer
import os.path

drive = cnm.getDrive()
dataDir   = os.path.join(drive,"era_interim_nc_6_hourly/")
years = np.arange(2002, 2003)

# Create surface RH for each year of chosen range
for year in years:

	year = str(year)

	print 'Calculating 2-meter RH% for year: ' + year

	# Load surface pressure, temperature, and dew point temperature
	d2m_nc = Dataset(os.path.join(dataDir, 'd2m_' + year + '.nc'),'r')
	d2m = d2m_nc.variables['d2m'][:]
	time = d2m_nc.variables['time']
	latitude = d2m_nc.variables['latitude']
	longitude = d2m_nc.variables['longitude']

	t2m_nc = Dataset(os.path.join(dataDir, 't2m_' + year + '.nc'),'r')
	t2m = t2m_nc.variables['t2m'][:]
	t2m_nc.close()

	sp_nc = Dataset(os.path.join(dataDir, 'sp_' + year + '.nc'),'r')
	sp = sp_nc.variables['sp'][:]
	sp_nc.close()


	if (sp.shape == t2m.shape == sp.shape) != True:
		raise ValueError('The dimensions of the loaded variables do not match.')

	# Relative humidity is the ratio of the water vapor pressure to the maximum
	# water vapor pressure at a given temperature.
	# RH = 100 * es(Td)/es(T)

	def e_sat(T, surface='water'):
		"""Calculates saturation vapor pressure for a given temperature. Calculation
		can be made with respected to 'water' and 'ice'"""

		if surface == 'water':
			a1 = 611.21 # Pa
			a3 = 17.502 # dimensionless
			a4 = 32.19  # K
		elif surface == 'ice':
			a1 = 611.21 # Pa
			a3 = 22.587
			a4 = -0.7   # K

		# Does not change with respect to surface state
		To = 273.16 # K

		e = a1 * np.exp(a3 * ( (T - To) / (T - a4) ) )

		return e

	# Calculate relative humiditty everywhere!
	RH = 100. * e_sat(d2m, 'water') / e_sat(t2m, 'water')

	################################################################################
	# Write the RH calculation!
	################################################################################

	saveName = dataDir + 'rh2m_' + year +'.nc'

	ncFile = Dataset(saveName, 'w', format='NETCDF4')
	ncFile.description = 'Relative humidity calculated using es(Td)/es(T) where es istetens eqn'
	ncFile.location = 'Global'
	ncFile.createDimension('time', len(time[:]) )
	ncFile.createDimension('latitude', len(latitude[:]) )
	ncFile.createDimension('longitude', len(longitude[:]) )

	# Create variables on the dimension they live on
	time_ = ncFile.createVariable('time', 'i4', ('time',))
	time_.units = time.units
	time_.calendar = time.calendar
	time_[:] = time[:]

	# create lat variable
	latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
	latitude_.units = latitude.units
	latitude_[:] = latitude

	# create longitude variable
	longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
	longitude_.units = longitude.units
	longitude_[:] = longitude

	# Create the RH variable as a floating point value
	RH_ = ncFile.createVariable('rh2m', 'f4', ('time','latitude','longitude'))
	RH_.units = '%'

	RH_[:] = RH

	ncFile.close()
