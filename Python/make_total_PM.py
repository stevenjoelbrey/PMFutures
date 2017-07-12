#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to calculate total PM concentration. 
# From Lamarque et al 2012
# Please refer to the Make_PM25_README in /FIRE_EMissions and the recomended 
# publications for a better understanding of the formula in this script. 

# TODO: Save a separate file that saves out fire emitted PM25 seperately?

import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:

	scenario    = sys.argv[1]
	print 'Using arguments passed via command line.'
	print sys.argv

else:
	scenario    = '2000Firev1'
	print 'Using default: ' + scenario

import os
import numpy as np
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import cesm_nc_manager as cnm
import time as timer
import os.path

dataDirBase= os.path.join("/pierce-scratch","mariavm")

startTime = timer.time()

# species| Formula                       | emitted_by_fire_module
#--------------------------------------------------------------------
# CB1    | C, hydrophpbic black carbon   | y, as BC
# CB2    | C, hydrophilic black carbon   | y, as BC
# DST01  | AlSiO5                        | n
# DST02  | AlSiO5                        | n
# DST03  | AlSiO5						 | n
# DST04  | AlSiO5                        | n
# NH4    |         						 | n	
# NH4NO3 |                               | n
# OC1    | C, hydrophobic organic carbon | y, as OC
# OC2    | C, hydrophilic organic carbon | y, as OC
# SOA    | C12                           | n
# SO4    |                               | n
# SSLT01 | NaCl                          | n
# SSLT02 | NaCl                          | n
# SSLT03 | NaCl                          | n
# SSLT04 | NaCl                          | n

# To calculate PM2.5 = 
# 	  SO4 + 
#         CB1 + CB2 + 
#	  1.8*(OC1 + OC2) + 
#	  (SOAT + SOAM + SOAI + SOAB + SOAX) +
#         SSLT01 + SSLT02
#         NH4NO3 + 
#         DST01 + DST02 + 


def getSelf(dataDirBase, scenario, NCVariable):
	"""Quick function for returning nc file connection."""
	ncFile = cnm.makeAQNCFile(dataDirBase, NCVariable, scenario, tStep="daily")
	nc     = Dataset(ncFile, 'r')
	ncVar  = nc.variables[NCVariable]
	return ncVar 

# Get sulfate
SO4_SRF = getSelf(dataDirBase, scenario, 'SO4_SRF')

# Get black carbon
# BC = CB1 + CB2
CB1_SRF = getSelf(dataDirBase, scenario, 'CB1_SRF')
CB2_SRF = getSelf(dataDirBase, scenario, 'CB2_SRF')

# Get secondary organic aerosol (SOA)
# SOA = SOAB + SOAI + SOAM + SOAX + SOAT
SOAB_SRF = getSelf(dataDirBase, scenario, 'SOAB_SRF')
SOAI_SRF = getSelf(dataDirBase, scenario, 'SOAI_SRF')
SOAM_SRF = getSelf(dataDirBase, scenario, 'SOAM_SRF')
SOAX_SRF = getSelf(dataDirBase, scenario, 'SOAX_SRF')
SOAT_SRF = getSelf(dataDirBase, scenario, 'SOAT_SRF')

# Get sea salt aerosol
# FineSSLT = SSLT01 + SSLT02
SSLT01_SRF = getSelf(dataDirBase, scenario, 'SSLT01_SRF')
SSLT02_SRF = getSelf(dataDirBase, scenario, 'SSLT02_SRF')

# Get dust
# FineDST = DST01 + DST02
DST01_SRF = getSelf(dataDirBase, scenario, 'DST01_SRF')
DST02_SRF = getSelf(dataDirBase, scenario, 'DST02_SRF')

# Organic Aerosol (OA)
# Organic aerosol (OA) = OC1+OC2
OC1_SRF = getSelf(dataDirBase, scenario, 'OC1_SRF')
OC2_SRF = getSelf(dataDirBase, scenario, 'OC2_SRF')

# Get NH4NO3 (Ammonium nitrate)
NH4NO3_SRF = getSelf(dataDirBase, scenario, 'NH4NO3_SRF')

# All needed nc connections are made. Now do the math for total PM2.5
t   = getSelf(dataDirBase, scenario, 'date')
n_t = len(t)

# Make an array to store the PM2.5
PM25     = np.zeros(SO4_SRF.shape)
FirePM25 = np.zeros(SO4_SRF.shape)

# Start putting the massive data arrays into the work space
SO4 = SO4_SRF[:]

BC  = CB1_SRF[:] + CB2_SRF[:]

OC = OC1_SRF[:] + OC2_SRF[:] 

SOA = SOAT_SRF[:] + SOAM_SRF[:] +  SOAI_SRF[:] + SOAB_SRF[:] + SOAX_SRF[:] 

NH4NO3 = NH4NO3_SRF[:]

FineDST = DST01_SRF[:] + DST02_SRF[:]

FineSSLT = SSLT01_SRF[:] + SSLT02_SRF[:]

# Calculate PM2.5 in units of kg/kg (all of these arrays are kg/kg)
PM25 = SO4 + BC + 1.8*OC + SOA + NH4NO3 + FineDST + FineSSLT

# Calculate PM2.5 emitted directly by fires
FirePM25 = BC + OC


additionComplete = timer.time()
dt = (additionComplete - startTime) / 60. 
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to calculate [PM25] from very big arrays'
print '----------------------------------------------------------------------'


# TODO: estimate the concentrations
#print 'Working on converting PM25 kg/kg to ug/m3'
# Now I want PM25 in units of ug/m3 using the idea gas law for dry air
T  = getSelf(dataDirBase, scenario, 'T')
PS = getSelf(dataDirBase, scenario, 'PS')

R = 287.04 # [J/K/kg ]
Rho = PS[:] / (R * T[:,25,:,:])

# Mass conversion variables
gramsPerKg  = 1000.
ugPerGram   = 1e6

# Calculate mass concentration from kg/kg
PM25_ugPerm3      = PM25 * Rho * gramsPerKg * ugPerGram
fire_PM25_ugPerm3 = FirePM25 * Rho * gramsPerKg * ugPerGram

# Use PS top get spatial dimensions of data 
ncFile = cnm.makeAQNCFile(dataDirBase, "SO4_SRF", scenario, tStep="daily")
nc     = Dataset(ncFile, 'r')
lat    = nc.variables['lat'][:]
lon    = nc.variables['lon'][:]
time   = nc.variables['time'][:]
nc.close()

###############################################################################
# Write total fire PM25 as daily netCDF data
###############################################################################	
# TODO: write cesm function that is designed to save a single netcdf variable
# TODO: this can be used in a lot of different places already.
print 'Working on writing the total fire PM NETCDF output.'

saveName = cnm.makeAQNCFile(dataDirBase, 'FirePM25', scenario, 'daily')
ncFile = Dataset(saveName, 'w', format='NETCDF4')
ncFile.description = 'Total surface layer fire emitted PM2.5 as defined in README_CESM'
ncFile.location = 'Global'
ncFile.createDimension('time', n_t )
ncFile.createDimension('lat', len(lat) )
ncFile.createDimension('lon', len(lon) )

# Create variables on the dimension they live on 
PMVar = ncFile.createVariable('FirePM25', 'f4', ('time','lat','lon'))
PMVar.units = 'kg/kg'
PMVar.description = 'kg of emitted species per kg of air in surface layer. P=Rho*R*T to get dry air concentration estimate.'

time_var = ncFile.createVariable('time', 'i4', ('time',))
time_var.units = 'days from origin'

latitude = ncFile.createVariable('lat', 'f4', ('lat',))
latitude.units = 'degrees north'
 
longitude = ncFile.createVariable('lon', 'f4', ('lon',))
longitude.units = 'degrees east'

# Write the actual data to these dimensions
PMVar[:]         = PM25
latitude[:]      = lat
longitude[:]     = lon
time_var[:]      = time

ncFile.close()

###############################################################################
# Write total PM25 as daily netCDF data
###############################################################################	
# TODO: write cesm function that is designed to save a single netcdf variable
# TODO: this can be used in a lot of different places already.
print 'Working on writing the fire emitted PM NETCDF output.'

saveName = cnm.makeAQNCFile(dataDirBase, 'PM25', scenario, 'daily')
ncFile = Dataset(saveName, 'w', format='NETCDF4')
ncFile.description = 'PM25 in surface layer'
ncFile.location = 'Global'
ncFile.createDimension('time', n_t )
ncFile.createDimension('lat', len(lat) )
ncFile.createDimension('lon', len(lon) )

# Create variables on the dimension they live on 
PMVar = ncFile.createVariable('PM25', 'f4', ('time','lat','lon'))
PMVar.units = 'kg/kg'
PMVar.description = 'kg PM25 per kg of air in surface layer. P=Rho*R*T to get dry air concentration estimate.'

time_var = ncFile.createVariable('time', 'i4', ('time',))
time_var.units = 'days from origin'

latitude = ncFile.createVariable('lat', 'f4', ('lat',))
latitude.units = 'degrees north'
 
longitude = ncFile.createVariable('lon', 'f4', ('lon',))
longitude.units = 'degrees east'

# Write the actual data to these dimensions
PMVar[:]         = FirePM25
latitude[:]      = lat
longitude[:]     = lon
time_var[:]      = time

ncFile.close()


###############################################################################
# Write total PM25_conc as daily netCDF data
###############################################################################	
# TODO: write cesm function that is designed to save a single netcdf variable
# TODO: this can be used in a lot of different places already.
print 'Working on writing the total PM25_conc NETCDF output.'

saveName = cnm.makeAQNCFile(dataDirBase, 'PM25_conc', scenario, 'daily')
ncFile = Dataset(saveName, 'w', format='NETCDF4')
ncFile.description = 'PM25 concentration in surface layer'
ncFile.location = 'Global'
ncFile.createDimension('time', n_t )
ncFile.createDimension('lat', len(lat) )
ncFile.createDimension('lon', len(lon) )

# Create variables on the dimension they live on 
PMVar = ncFile.createVariable('PM25_conc', 'f4', ('time','lat','lon'))
PMVar.units = 'ug/m3'
PMVar.description = 'ug/m3 in surface layer air volume. P=Rho*R*T of surface used to get from kg/kg to mass concentration.'

firePMVar = ncFile.createVariable('fire_PM25_conc', 'f4', ('time','lat','lon'))
firePMVar.units = 'ug/m3'
firePMVar.description = 'fire emitted PM ug/m3 in surface layer air volume. P=Rho*R*T of surface used to get from kg/kg to mass concentration.'

time_var = ncFile.createVariable('time', 'i4', ('time',))
time_var.units = 'days from origin'

latitude = ncFile.createVariable('lat', 'f4', ('lat',))
latitude.units = 'degrees north'
 
longitude = ncFile.createVariable('lon', 'f4', ('lon',))
longitude.units = 'degrees east'

# Write the actual data to these dimensions
PMVar[:]         = PM25_ugPerm3
firePMVar[:]     = fire_PM25_ugPerm3   
latitude[:]      = lat
longitude[:]     = lon
time_var[:]      = time

ncFile.close()

writingComplete = timer.time()
dt = (writingComplete - startTime) / 60. 
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to run the entire script.'
print '----------------------------------------------------------------------'




