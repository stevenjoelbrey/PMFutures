# PMFutures Data Pipeline

This document is designed to show the data structure of this project. The idea is that 
this document can serve as an easy guide to reproducing all of the data used by R and 
Python analysis scripts in the event that 1) the data is lost and I (steven brey) need to 
re-create everything fast, 2) someone wants to reproduce this work, and is lost in 
figuring out where the data come from.

~ indicates "PMFutures/" top directory. 
/scratch = a directory not under version control where I keep very large datasets on 
either my local machine or local cluster. These are never under version control. 

Each major type of data (e.g. emissions, meteorology, fire records) will contain the 
following information

	- Datasource (url or other resource description)
	- How the data was downloaded, script of manual. If script, what script and where 
	  does it live?
	- Location of this raw, not yet touched or changed by this product data, and 
	  version control information. 
	- Various scripts used to transform/subset/comb data and where that output lives and
	  what scripts use that final version of production data
	  
	
### GFED4s 

Datasource: http://www.geo.vu.nl/~gwerf/GFED/GFED4/*hdf5

Datasource README: http://www.geo.vu.nl/~gwerf/GFED/GFED4/Readme.pdf

	- These data were downloaded by manually clicking each hdf5 file
	- location: /scratch/GFED4s/
		- This directory contains hdf5 files that contain daily emissions 2003-2016
		
	- ~Python/GFED4s_to_nc.py
		- This script reads yearly hdf5 files and saves out yearly and multiyear .nc
		  files of a single emissions species in the hdf5 file, e.g. C or DM. This is
		  where daily emission estimate arrays are applied to the GFED4s emissions, based
		  on  
		  
	- ~Python/regrid_fire_emissions.py
		- This script takes the GFED4s nc data and regrids to ecmwf 0.75 deg grid based
		  on nearest neighbor method. 'MEG_GRID' string is inserted to these .nc files,
		  which also live in /scratch/GFED4s/
		

### ECMWF Reanalysis (era-interim)

Datasource: http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=sfc/

and:		http://apps.ecmwf.int/datasets/data/interim-full-daily/levtype=pl/
	
	- ~Python/get_ERA_Interim_data.py 
	    - these data are downloaded using the script mentioned above. 
		- this script uses the ecmwfapi which allows much easier downloading than
		  the web interface. 
		- requests are made for one variable at a time. The variable name requested 
	      matches the name of the downloaded file name, the downloaded file name 
	      contains the year of the data. e.g. z_2003.nc contains 6-hourly data of 'z'
		- these 6-hourly data are saved in scratch/era_interim_nc_6_hourly/. These
		  files are strait from the web and have not been altered by my project.

	- ~Python/average6HourlyData.py 
		- this averages the 6hourly data saved in scratch/era_interim_nc_6_hourly/ and
		  saves them to scratch/era_interim_nc_daily/
		  
	- ~Python/merge_year_ny.py 
		- This script is used to merge yearly daily era-interim nc files or GFED4s
          fire emissions files. These merged files are saved in 
          scratch/era_interim_nc_daily_merged/
          
    - ~Python/get_all_era_interim_z.py 
    	- unrelated to the rest of the data pipeline structure. This is used for getting
    	  as much z data as is possible in order to create a julain day 500 mb height 
    	  climatology for the US.
  
### FINN (Fire INventory from NCAR)

The gridded version of FINN used for this work was created by the GEOS-Chem Support Team. 
I used the HEMCO/FINN/v2015-02 version of these data. 

18 Feb 2015
GEOS-Chem Support Team
geos-chem-support@as.harvard.edu
 
"The Fire INventory from NCAR version 1.0 (FINNv1) provides daily, 1-km
resolution, global estimates of the trace gas and particle emissions from
open burning of biomass, which includes wildfire, agricultural fires, and
prescribed burning and does not include biofuel use and trash burning"
(Wiedinmyer et al., 2011, citation below).	  

##### GEOS-Chem wiki:

  http://wiki.geos-chem.org/FINN_biomass_burning_emissions

##### FINN reference paper:

  Wiedinmyer, C., S. K. Akagi, R .J. Yokelson, L. K. Emmons, J. A. Al-Saadi,
  J. J. Orlando, and A. J. Soja. "The Fire INventory from NCAR (FINN): a
  high resolution global model to estimate the emissions from open burning:",
  Geosci. Model Dev., 4, 625641, 2011.
	 
	
TODO: 
	- ~Python/process_FINN.py
		- This script processes FINN data organized for GEOS-Chem and makes it match
		  the GFED emissions format already used by this work. 
	- ~Python/regrid_fire_emissions.py
		- This script takes fire emissions and places them onto the meteorology grid of 
		  my choosing. 	 
	 