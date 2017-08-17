# PMFutures Data Pipeline

This document is designed to show the data structure of this project. The idea is that 
this document can serve as an easy guide to reproducing all of the data used by R and 
Python analysis scripts in the event that 1) the data is lost and I (steven brey) need to 
re-create everything fast, 2) someone wants to reproduce this work, and is lost in 
figuring out where the data come from.

Each major type of data (e.g. emissions, meteorology, fire records) will contain the 
following information

	- Datasource (url or other resource description)
	- How the data was downloaded, script of manual. If script, what script and where 
	  does it live?
	- Location of this raw, not yet touched or changed by this product data 
	- Various scripts used to transform/subset/comb data and where that output lives and
	  what scripts use that final version of production data
	  
	
# GFED4s 

## Datasource: http://www.geo.vu.nl/~gwerf/GFED/GFED4/*hdf5
	-	These data were downloaded by manually clicking each hdf5 file


# ECMWF Reanalysis (era-interim)

 
	  
	 