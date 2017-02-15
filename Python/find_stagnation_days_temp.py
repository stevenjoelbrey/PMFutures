def findStagnationDays(scenario):
	
	# For development mode only
	

	# Load daily surface pressure to en-route to make sea surface 
	# geostrophic winds 
	#PSLFile  = makeAQNCFile('PSL', '2000Base', 'hourly')
	#PSL_nc   = Dataset(PSLFile, 'r')
	#PSL      = PSL_nc.variables['PSL'][0:24,:,:]
	#PS.close()

	# We need ground level geopotential, which we will use to calc-
	# ulate surface elevation. 
	#srf_geo = getUSGSTopo()
	#zg = srf_geo['topo']  # meters
	#H = 8000.             # meters

	# This equal only works for small Zg, so it will not work over the 
	# mountinous regions of the U.S. 
	# TODO: How do I get Tv for a layer where the model only has T
	# TODO: above ground?? 
	#P0 = P_ground[1000,:,:] * np.exp(zg / H) 

	# Definition of geostrophic winds
	# Phi = geopotential [kg/m**2/s**2]
	# dy  = grid spacing difference in latitude in units of m  
	# dx  = grid spacing difference in longitude units of m 
	# f = 2 * omega * np.sin(lat), array that matches length of lat dim


	# Load geopotential height on pressure surfaces for calculating 
	# geostrophic winds
	Z3File = makeAQNCFile('Z3_P', '2050RCP45', 'daily')
	Z3_nc   = Dataset(Z3File, 'r')
	plevel  = Z3_nc.variables['plevel'][:]  
	lat     = Z3_nc.variables['lat'][:]
	lon     = Z3_nc.variables['lon'][:]
	# We are going to use 1000 hPa as a proxy for sea level pressure.
	slpIndex = np.where(plevel == 1000)[0][0]

	# Get the data 
	Z1000 = Z3_nc.variables['Z3'][10:20, slpIndex, :, :]

	
	# Create the coriolis parameter and spherical coords constants
	omega = 7.292e-5
	mPerDegLat = 111. * 1000 
	g = 9.81 

	###############################
	# Create a grid of dx and dy 
	###############################

	# Express the location where the geowind will be calcuated, between
	# grid points. 
	lonMid = lon[:-1] + (lon[1:] - lon[:-1]) / 2
	latMid = lat[:-1] + (lat[1:] - lat[:-1]) / 2
	nrows = len(latMid)
	ncols = len(lonMid)

	# Coriolis at the mid points where geo winds calculated 
	fMid     = 2. * omega * np.sin( (np.pi/180.) * latMid )

	# This is where grid distances in meters will be stored
    dx = np.zeros( (nrows, ncols) )
	dy = dx

	for i in range(nrows):	
		dx[i,:] =  np.diff(lon) * (mPerDegLat * np.cos(latMid[i]*np.pi/180.))
	
	for j in range(ncols):
		dy[:,j] = np.abs(np.diff(lat)) * mPerDegLat	

	dx    = np.diff(lon) * np.cos( (np.pi/180.) * lat) * mPerDegLat
	dy    = np.diff(lat) * mPerDefLat


	u_g  = np.zeros( (nrows, ncols) )
	v_g  = u_g

	for i in range(ncols):
		u_g[:,i] = (-1) * (g / fMid[i]) * np.diff(Z1000[8, :, i]) / dy[:, i]


	
	v_g = (1./fMid)  * np.diff(Z1000) / dx 

	v_g = u_g**2 


	# Place back on the original latitude longitude grid 

	# TODO: Compare 1000 mb geostrophic winds to U and V 1000mb winds as a sanity check


	# Find dates that sea-level geo winds are less than 8 m/s. Flag.

	# Load 500 hPa geostrophic winds 

	# Find locations where 500 hPa geostrophic winds are less then 13 m/s 

	# Find locations where there is less then 0.01 inches of precipitation

	# Return an equal size array equal to one where all of these conditions are
    # met and 0 where they are not. Consider returning all individual masks if
	# useful elsewhere. 
