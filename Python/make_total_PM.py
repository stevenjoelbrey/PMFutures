#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to calculate total PM concentration. 
# From Lamarque et al 2012
# Please refer to the README_CESM and the recomended publications for a better
# understanding of the formula in this script. 

import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	scenario    = sys.argv[1]

else:
	print 'Using default '
	scenario    = '2000Base'

import os
import numpy as np
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import cesm_nc_manager as cnm
import time as timer

startTime = timer.time()

# species| Formula                       | emitted_by_fire_module
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

# To calculate PM2.5,
#		BC = CB1 + CB2
#		SOA = SOAB + SOAI + SOAM + SOAX + SOAT
#		FineSSLT = SSLT01 + SSLT02
#		FineDST = DST01 + DST02
#		Organic aerosol (OA) = OC1 + OC2 + SOA
	
#		PM25 = BC + 1.2*OA + SO4 + NH3NO4 + FineSSLT + FineDST


def getSelf(scenario, species):
	"""Quick function for returning nc file connection."""
	ncFile = cnm.makeAQNCFile(species, scenario, 'daily')
	nc     = Dataset(ncFile, 'r')
	ncVar  = nc.variables[species]
	return ncVar 

# Get black carbon
# BC = CB1 + CB2
CB1_SRF = getSelf(scenario, 'CB1_SRF')
CB2_SRF = getSelf(scenario, 'CB2_SRF')

# Get secondary organic aerosol (SOA)
# SOA = SOAB + SOAI + SOAM + SOAX + SOAT
SOAB_SRF = getSelf(scenario, 'SOAB_SRF')
SOAI_SRF = getSelf(scenario, 'SOAI_SRF')
SOAM_SRF = getSelf(scenario, 'SOAM_SRF')
SOAX_SRF = getSelf(scenario, 'SOAX_SRF')
SOAT_SRF = getSelf(scenario, 'SOAT_SRF')

# Get sea salt aerosol
# FineSSLT = SSLT01 + SSLT02
SSLT01_SRF = getSelf(scenario, 'SSLT01_SRF')
SSLT02_SRF = getSelf(scenario, 'SSLT02_SRF')

# Get dust
# FineDST = DST01 + DST02
DST01_SRF = getSelf(scenario, 'DST01_SRF')
DST02_SRF = getSelf(scenario, 'DST02_SRF')

# Organic Aerosol (OA)
# Organic aerosol (OA) = OC1+OC2+SOA
OC1_SRF = getSelf(scenario, 'OC1_SRF')
OC2_SRF = getSelf(scenario, 'OC2_SRF')

# Get sulfate
SO4_SRF = getSelf(scenario, 'SO4_SRF')

# Get NH4NO3 (Ammonium nitrate)
# NOTE: There is a typo in the README_CESM. NH4NO3 is written NH3NO4
NH4NO3_SRF = getSelf(scenario, 'NH4NO3_SRF')

# All needed nc connections are made. Now do the math for total PM2.5
t   = getSelf(scenario, 'date')
n_t = len(t)

# Make an array to store the PM2.5
PM25     = np.zeros(SO4_SRF.shape)
FirePM25 = np.zeros(SO4_SRF.shape)

for i in range(n_t): #n_t

	#print i

	BC  = CB1_SRF[i,:,:] + CB2_SRF[i,:,:]

	FineSSLT = SSLT01_SRF[i,:,:] + SSLT02_SRF[i,:,:]

	FineDST = DST01_SRF[i,:,:] + DST02_SRF[i,:,:]
	
	SOA = SOAB_SRF[i,:,:] + SOAI_SRF[i,:,:] + SOAM_SRF[i,:,:]+\
          SOAX_SRF[i,:,:] + SOAT_SRF[i,:,:]

	OA = OC1_SRF[i,:,:] + OC2_SRF[i,:,:] + SOA

	SO4 = SO4_SRF[i,:,:]

	NH4NO3 = NH4NO3_SRF[i,:,:]

	# Calculate PM2.5 in units of kg/kg (all of these arrays are kg/kg)
	PM25[i,:,:] = BC + 1.2 * OA + SO4 + NH4NO3 + FineSSLT + FineDST

	# Calculate PM2.5 emitted directly by fires
	FirePM25[i,:,:] = BC + OC1_SRF[i,:,:] + OC2_SRF[i,:,:]
	

#print 'Working on converting PM25 kg/kg to ug/m3'
# Now I want PM25 in units of ug/m3 using the idea gas law for dry air
#T  = getSelf(scenario, 'T')
#PS = getSelf(scenario, 'PS')

#R = 287.04 # [J/K/kg ]
#Rho = PS[:] / (R * T[:,25,:,:])

# Mass conversion variables
#gramsPerKg  = 1000.
#ugPerGram   = 1e6

#PM25_ugPerm3 = (PM25 * Rho) * (gramsPerKg * ugPerGram)


# Use PS top get spatial dimensions of data 
ncFile = cnm.makeAQNCFile("PS", scenario, 'daily')
nc     = Dataset(ncFile, 'r')
lat    = nc.variables['lat'][:]
lon    = nc.variables['lon'][:]
time   = nc.variables['time'][:]
nc.close()

###############################################################################
# Write total PM25 as daily netCDF data
###############################################################################	
# TODO: write cesm function that is designed to save a single netcdf variable
# TODO: this can be used in a lot of different places already.
print 'Working on writing the total PM NETCDF output.'

saveName = cnm.makeAQNCFile('PM25', scenario, 'daily')
ncFile = Dataset(saveName, 'w', format='NETCDF4')
ncFile.description = 'Total surface layer PM2.5 as defined in README_CESM'
ncFile.location = 'Global'
ncFile.createDimension('time', n_t )
ncFile.createDimension('lat', len(lat) )
ncFile.createDimension('lon', len(lon) )

# Create variables on the dimension they live on 
PMVar = ncFile.createVariable('PM25', 'f4', ('time','lat','lon'))
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

saveName = cnm.makeAQNCFile('FirePM25', scenario, 'daily')
ncFile = Dataset(saveName, 'w', format='NETCDF4')
ncFile.description = 'Fire emitted PM2.5 in surface layer'
ncFile.location = 'Global'
ncFile.createDimension('time', n_t )
ncFile.createDimension('lat', len(lat) )
ncFile.createDimension('lon', len(lon) )

# Create variables on the dimension they live on 
PMVar = ncFile.createVariable('FirePM25', 'f4', ('time','lat','lon'))
PMVar.units = 'kg/kg'
PMVar.description = 'kg of fire emitted PM per kg of air in surface layer. P=Rho*R*T to get dry air concentration estimate.'

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


writingComplete = timer.time()
dt = (writingComplete - startTime) / 60. 
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to run the entire script.'
print '----------------------------------------------------------------------'


